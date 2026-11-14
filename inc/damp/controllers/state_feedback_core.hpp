// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file state_feedback_core.hpp
 * @brief The runtime the LQ laws share: detail::StateFeedbackCore, one tick with every feature
 *        (integral action, a filter, a steady-state target, input registers) optional
 *
 * Each runtime set is its own header: lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp (design function,
 * result, runtime). The types a design hands this runtime are in state_feedback_types.hpp; the design
 * itself is detail/state_feedback_design.hpp and Law<kRes, T> detail/lq_law.hpp.
 */

#include "damp/backend.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"

namespace damp {
namespace detail {

/**
 * @brief The integral action of a law: the integrators, the integrated outputs Ci (on the
 *        measurement) and their feedthrough Ci·D, the per-sample step, the error held for the next
 *        predict, and with an estimator the measurement held from the last update (NYH = NY; 0: none).
 */
template<size_t NI, size_t NC, size_t NU, size_t NYH, typename T, typename CiMatrix = Matrix<NI, NC, T>, typename CiuMatrix = Matrix<NI, NU, T>>
struct IntegralAction {
    using HeldY = maybe<(NYH > 0), ColVec<NYH, T>, OffSlot::held_y>;

    ColVec<NI, T>                   xi{};    ///< Integral of the tracking error
    [[no_unique_address]] CiMatrix  Ci{};    ///< Integrated outputs: e = r − Ci·(measurement)
    [[no_unique_address]] CiuMatrix Ci_u{};  ///< Ci·D: the integrated outputs' feedthrough of the input
    T                               step{1}; ///< Per-sample integrator gain: the design's Ts (xi += Ts·e ≈ ∫e dt)
    ColVec<NI, T>                   e{};     ///< Error held from control for the next predict, before the input's feedthrough
    [[no_unique_address]] HeldY     y{};     ///< Measurement from the last update (with an estimator)
};

/** @brief A gain design result other than StateFeedbackResult: K, success and Ts (PlaceResult, LQRResult). */
template<typename R>
concept gain_result = requires(const R& r) {
    r.success;
    r.Ts;
    r.K;
    requires MatrixLike<damp::remove_cvref_t<decltype(r.K)>>;
} && !is_state_feedback_result<R>::value;

/** @brief A gain design is discrete: its @c discrete flag (LQRResult), else a sample time Ts > 0 (PlaceResult, StateFeedbackResult). */
template<typename R>
[[nodiscard]] constexpr bool discrete_gain(const R& r) {
    if constexpr (requires { r.discrete; }) {
        return r.discrete;
    } else {
        return r.Ts != decltype(r.Ts){0};
    }
}

/** @brief Where a state-feedback law keeps its matrices: in the object, any value at run time. */
struct HeldGains {
    static constexpr bool constant = false; ///< The matrices are members

    template<size_t R, size_t C, typename T>
    using K = Matrix<R, C, T>; ///< The gain
    template<size_t R, size_t C, typename T>
    using Ci = Matrix<R, C, T>; ///< The integrated outputs
    template<size_t R, size_t C, typename T>
    using Ci_u = Matrix<R, C, T>; ///< Their feedthrough
    template<size_t R, size_t C, typename T>
    using TargetD = Matrix<R, C, T>; ///< The target per unit d̂
    template<size_t R, size_t C, typename T>
    using TargetR = Matrix<R, C, T>; ///< The target per unit r
};

/**
 * @brief The implementation the state-feedback runtimes share: u = u* − K·(x̂ − x*), every feature
 *        optional and free when off.
 *
 * LQR, LQI, LQG, LQGI and OFLQG wrap it, each exposing one update, one control and one step.
 *
 * | Feature | Parameter | Off |
 * | --- | --- | --- |
 * | integral action on r − Ci·y | NI, NY | 0 |
 * | estimator (FixedGainObserver with design::kalman's gain, disturbance model) | Est | NoEstimator: the state is measured |
 * | steady-state target for r, under the disturbance model if there is one (SteadyStateTarget) | NZ | 0 |
 * | input registers: the inputs applied over the last NRU periods, held by the law (a measured state) | NRU | 0 |
 *
 * One tick is @ref step(): @ref predict (the time update with the input @ref back_calculate
 * recorded, else the last command), then @ref update with the measurement at sample k (the
 * state itself without an estimator), then @ref control, which returns u_k. Discrete-time
 * index: predict uses the input held over [k−1, k); update consumes m_k; control emits u_k.
 * The low-latency ISR is update(m_k) → control → write → back_calculate(u) → predict(): the
 * trailing predict is the next sample's time update. For a measured-state LQR, predict is
 * empty, so step(x) equals update then control.
 *
 * With input registers (NRU > 0, a design on a plant whose input arrives late) the law's model is
 * [x; v(k−1); …; v(k−NRU)]. On a measured state update takes the plant state x, and the law holds
 * the inputs it applied: @ref predict shifts the input back_calculate recorded (else the last
 * command) into the first register, so it runs every tick, as step does. A filter carries them as
 * states of its model, shifted by its own predict. The reference point of a register is u*: an
 * input held there is the one the law would send at the reference, so x_ref is the plant's.
 *
 * The law:
 * - about a reference point (x*, u*): u = u* − K·[x̂ − x*; xi]; x* = 0, u* = 0 is the regulator;
 * - with a target and no disturbance model (TrackedOutputs): the target's steady state for r,
 *   x_t = Nₓ·r and u_t = Nᵤ·r, adds to (x*, u*);
 * - with a disturbance model (ND > 0): the target is the steady state that holds the controlled
 *   outputs at r under d̂, taken about (x*, u*); the estimator runs in deviations from it
 *   (SteadyStateTarget::anchor), and u = u_a − K·(x̂ − x_a) at the deviation the command meets;
 * - with integral action: xi integrates r − Ci·y, xi += Ts·(r − Ci·y) per tick (xi ≈ ∫e dt),
 *   y = C·x + D·u the output over the period, u the input applied over it (the design's model):
 *   control holds the part it knows, r − Ci·C·x (from the measured state, or from the
 *   measurement less the feedthrough D·u of the input it was taken under), and @ref predict
 *   subtracts Ci·D·u for the input applied, then advances xi unless the applied input shows the
 *   actuator clamping against it (conditional integration).
 *
 * @see "Optimal Control" (Anderson & Moore, 1990), §4.1 and §8
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §3.5: back-calculation anti-windup
 *
 * @tparam NX  Plant states
 * @tparam NU  Inputs
 * @tparam T   Scalar
 * @tparam NI  Integrators (0: none)
 * @tparam NY  Measured outputs (with no estimator, the outputs the integrators read through C)
 * @tparam Est Estimator: FixedGainObserver (the fixed Kalman gain of design::kalman), or NoEstimator for a measured state
 * @tparam NZ  Controlled outputs of the steady-state target (0: none)
 * @tparam Store Where the matrices live: HeldGains (members), or constants of the type (Law of a `constexpr` design)
 * @tparam NRU Input registers the law holds, NRU·NU states after the plant's (0: none)
 */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Est, size_t NZ, typename Store = HeldGains, size_t NRU = 0>
struct StateFeedbackCore {
    using Traits = detail::estimator_traits<Est>;
    static constexpr bool   kEstimated = Traits::present;         ///< A filter estimates the state
    static constexpr size_t NXL = NX + (NRU * NU);                ///< The law's model: the plant, then its input registers
    static constexpr size_t ND = Traits::ND;                      ///< Disturbance states of the estimator
    static constexpr size_t NXE = kEstimated ? Traits::NXE : NXL; ///< The filter's model before the disturbance (the law's without one)
    static constexpr size_t NIn = kEstimated ? NY : NX;           ///< What update takes: y, or the plant state
    static constexpr size_t NM = kEstimated ? NY : NXL;           ///< What the integrators read: y, or the law's state
    static constexpr size_t NR = (NI > 0) ? NI : NZ;              ///< The reference r: integrated or tracked outputs (0: none)
    static_assert(!kEstimated || Traits::NXE == NXL, "The estimator is on the law's model.");
    static_assert(NI == 0 || NZ == 0, "Integral action and a disturbance model are two offset-free choices: lqgi or oflqg, not both.");
    static_assert(ND == 0 || NZ > 0, "A disturbance model tracks through its steady-state target.");
    static constexpr bool kRecord = (NI > 0) || kEstimated || (NRU > 0); ///< The law keeps its command and the applied input
    using Record = damp::conditional_t<kRecord, detail::AppliedRecord<NU, T>, detail::NoAppliedRecord>;

    using Target = detail::maybe<(NZ > 0), SteadyStateTarget<NXE, NU, ND, NZ, T, typename Store::template TargetD<NXE + NU, ND, T>, typename Store::template TargetR<NXE + NU, NZ, T>>, detail::OffSlot::target>;
    using Integral = detail::maybe<(NI > 0), detail::IntegralAction<NI, NM, NU, (kEstimated ? NY : 0), T, typename Store::template Ci<NI, NM, T>, typename Store::template Ci_u<NI, NU, T>>, detail::OffSlot::integral>;
    using HeldX = detail::maybe<!kEstimated, ColVec<NXL, T>, detail::OffSlot::held_x>;

    [[no_unique_address]] typename Store::template K<NU, NXL + NI, T> K{};         ///< Gain on [x; xi]
    [[no_unique_address]] Est                                         estimator{}; ///< The filter (none for a measured state)
    [[no_unique_address]] Target                                      target{};    ///< The steady-state target (NZ > 0)
    [[no_unique_address]] Integral                                    integral{};  ///< The integrators (NI > 0)

    constexpr StateFeedbackCore() = default;

    /** @brief From a gain on the law's state (no estimator): K = @p K_. */
    constexpr explicit StateFeedbackCore(const Matrix<NU, NXL + NI, T>& K_)
        requires(!kEstimated && !Store::constant)
        : K(K_) {}

    /**
     * @brief From a gain design (design::place's PlaceResult, LQRResult; any scalar): converts K.
     *
     * A failed or continuous (not discrete) result does not compile when constant-evaluated; at run
     * time either leaves K = 0.
     */
    template<typename R>
        requires(detail::gain_result<R> && !kEstimated && NI == 0 && NZ == 0 && NRU == 0 && !Store::constant)
    constexpr StateFeedbackCore(const R& result) // NOLINT
    {
        if (usable(result, detail::discrete_gain(result))) {
            K = result.K.template as<T>();
        }
    }

    /**
     * @brief From an LQ design (any scalar): the gain, and the filter, integrators and target the
     *        result carries.
     *
     * A failed or continuous result (Ts = 0) does not compile when constant-evaluated; at run
     * time either leaves K = 0.
     */
    template<size_t NXR, size_t NUR, typename U, size_t NIR, size_t NYR, typename Kalman, size_t NZR, size_t NRR>
        requires(!Store::constant)
    constexpr StateFeedbackCore(const design::StateFeedbackResult<NXR, NUR, U, NIR, NYR, Kalman, NZR, NRR>& result) { // NOLINT
        constexpr bool same = NXR == NX && NUR == NU && NIR == NI && NYR == NY && NZR == NZ && NRR == NRU && damp::is_same_v<Est, detail::runtime_estimator_t<Kalman, T>>;
        static_assert(same, "This law's sizes are not the design's: let them deduce, LQG law{res}, or Law<res, T> for a member.");
        if constexpr (same) {
            if (usable(result, result.Ts != U{0})) {
                K = result.K.template as<T>();
            }
            if constexpr (kEstimated) {
                estimator = Est(result.kalman);
            }
            if constexpr (NI > 0 && kEstimated) {
                integral.Ci = result.integral.Ci.template as<T>(); // on the measurement y
            } else if constexpr (NI > 0) {
                integral.Ci = result.integral.Ci_x.template as<T>(); // on the measured state: Ci·C
            }
            if constexpr (NI > 0) {
                integral.Ci_u = result.integral.Ci_u.template as<T>();
                integral.step = static_cast<T>(result.integral.step);
            }
            if constexpr (NZ > 0) {
                target = Target(result.target);
            }
        }
    }

    // ---- the tick ------------------------------------------------------------------------

    /**
     * @brief Record the input the plant received, @p u_sat (after the actuator's clamp), and the
     *        feedforward @p u_ff added outside the law, for the next @ref predict.
     *
     * The clamp's excess is u_sat − (u_cmd + u_ff), u_cmd the law's last command, so a
     * feedforward is not taken for saturation. The estimator moves on @p u_sat. A law that keeps
     * no memory (no integrator or estimator) ignores it.
     */
    constexpr void back_calculate(const ColVec<NU, T>& u_sat, const ColVec<NU, T>& u_ff = ColVec<NU, T>{}) {
        if constexpr (kRecord) {
            applied_.record(u_sat, u_ff);
        } else {
            static_cast<void>(u_sat);
            static_cast<void>(u_ff);
        }
    }

    /**
     * @brief Time update with the input the plant received over the last period (back_calculate's,
     *        else the last command): the integrators advance unless that would wind them up
     *        against a clamp, the estimator predicts.
     */
    constexpr void predict() {
        if constexpr (kRecord) {
            applied_.take();
        }
        if constexpr (NI > 0) {
            integrate();
        }
        if constexpr (kEstimated) {
            estimator.predict(applied_.applied);
        }
        if constexpr (NRU > 0 && !kEstimated) {
            for (size_t i = NXL; i-- > NX + NU;) { // oldest first: each register takes the one before it
                x_(i) = x_(i - NU);
            }
            for (size_t j = 0; j < NU; ++j) {
                x_(NX + j) = applied_.applied(j);
            }
        }
    }

    /**
     * @brief Measurement update: the filter's with the outputs @p m (the integrators read them at
     *        the next control), or, with no filter, the plant state @p m held for the next control
     *        (the input registers after it are the law's own).
     */
    constexpr void update(const ColVec<NIn, T>& m) {
        if constexpr (kEstimated) {
            if constexpr (NI > 0) {
                integral.y = m;
            }
            estimator.update(m);
        } else if constexpr (NRU == 0) {
            x_ = m;
        } else {
            for (size_t i = 0; i < NX; ++i) {
                x_(i) = m(i);
            }
        }
    }

    /**
     * @brief The law at the last update's measurement, about the reference point (@p x_ref, @p u_ref):
     *        integral action on @p r (NI > 0), or the steady-state target for @p r (NZ > 0, under d̂ when ND > 0).
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref) {
        const ColVec<NXL, T> ref = on_law_model(x_ref, u_ref);
        if constexpr (ND > 0) {
            return anchored_law(r, ref, u_ref);
        } else if constexpr (NZ > 0) {
            const auto p = target.equilibrium(r); // on the law's model: each register at u_t
            return reference_law(law_state(), ColVec<NXL, T>(ref + p.x), ColVec<NU, T>(u_ref + p.u));
        } else if constexpr (NI > 0 && kEstimated) {
            // y was measured under the input the last predict applied: take its feedthrough out.
            return integral_law(r, law_state(), ref, u_ref, integral.y, ColVec<NI, T>(integral.Ci_u * applied_.applied));
        } else if constexpr (NI > 0) {
            return integral_law(r, x_, ref, u_ref, x_, ColVec<NI, T>{});
        } else {
            static_cast<void>(r);
            return reference_law(law_state(), ref, u_ref);
        }
    }

    /** @brief One tick: predict(), update(@p m), then control(@p r, @p x_ref, @p u_ref). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NIn, T>& m) {
        predict();
        update(m);
        return control(r, x_ref, u_ref);
    }

    // ---- state ---------------------------------------------------------------------------

    /** @brief The command the last control returned: the input a predict without back_calculate uses. */
    [[nodiscard]] constexpr const ColVec<NU, T>& last_command() const
        requires(kRecord)
    {
        return applied_.command;
    }

    /** @brief The input the last predict used: back_calculate's, else the command before it. */
    [[nodiscard]] constexpr const ColVec<NU, T>& applied() const
        requires(kRecord)
    {
        return applied_.applied;
    }

    /** @brief The estimated plant state x̂ (absolute). */
    [[nodiscard]] constexpr ColVec<NX, T> state_estimate() const
        requires(kEstimated)
    {
        return plant_estimate();
    }

    /** @brief The estimated disturbance d̂. */
    [[nodiscard]] constexpr ColVec<ND, T> disturbance_estimate() const
        requires(ND > 0)
    {
        return estimator.disturbance_estimate();
    }

    /** @brief The gain K on [x; xi]. */
    [[nodiscard]] constexpr const Matrix<NU, NXL + NI, T>& getK() const { return K; }

    /** @brief Clear the integrators and the input memory (an estimate and its anchor stay). */
    constexpr void reset() {
        if constexpr (kRecord) {
            applied_ = Record{};
        }
        if constexpr (kEstimated) {
            estimator.reset();
        }
        if constexpr (NI > 0) {
            integral.xi = ColVec<NI, T>{};
            integral.e = ColVec<NI, T>{};
            if constexpr (kEstimated) {
                integral.y = ColVec<NY, T>{};
            }
        }
        if constexpr (NRU > 0 && !kEstimated) {
            x_ = ColVec<NXL, T>{};
        }
    }

private:
    [[no_unique_address]] Record applied_{}; ///< The last command and the input back_calculate recorded
    [[no_unique_address]] HeldX  x_{};       ///< The measured state from the last update (no estimator)

    /** @brief Whether a design's gain runs on this tick: it succeeded and is discrete. Otherwise it does not compile when constant-evaluated (the error names its status) and the law keeps K = 0. */
    template<typename R>
    [[nodiscard]] static constexpr bool usable(const R& result, bool discrete) {
        const bool ok = result.success && discrete;
        if (!ok && damp::is_constant_evaluated()) {
            if constexpr (requires { result.status; }) {
                detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_sample_time : result.status);
            } else {
                detail::refuse_failed_design(); // a gain result with no status
            }
        }
        return ok;
    }

    [[nodiscard]] constexpr ColVec<NX, T> plant_estimate() const
        requires(kEstimated)
    {
        return detail::leading<NX>(estimator.state_estimate());
    }

    [[nodiscard]] constexpr ColVec<NXL, T> law_state() const {
        if constexpr (kEstimated) {
            return detail::leading<NXL>(estimator.state_estimate());
        } else {
            return x_;
        }
    }

    /** @brief The reference point on the law's model: x*, then u* in each input register. */
    [[nodiscard]] static constexpr ColVec<NXL, T> on_law_model(const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref) {
        if constexpr (NRU == 0) {
            static_cast<void>(u_ref);
            return x_ref;
        } else {
            ColVec<NXL, T> out{};
            for (size_t i = 0; i < NX; ++i) {
                out(i) = x_ref(i);
            }
            for (size_t i = NX; i < NXL; ++i) {
                out(i) = u_ref((i - NX) % NU);
            }
            return out;
        }
    }

    [[nodiscard]] constexpr auto integrators() const {
        if constexpr (NI > 0) {
            return integral.xi;
        } else {
            return ColVec<0, T>{};
        }
    }

    [[nodiscard]] constexpr ColVec<NU, T> record(const ColVec<NU, T>& u) {
        if constexpr (kRecord) {
            applied_.command = u;
        }
        return u;
    }

    /** @brief u = u* − K·z, z = [x − x*; xi], x and x* on the law's model. */
    [[nodiscard]] constexpr ColVec<NU, T> reference_law(const ColVec<NXL, T>& x, const ColVec<NXL, T>& x_ref, const ColVec<NU, T>& u_ref) {
        const ColVec<NXL + NI, T> z = detail::stacked(ColVec<NXL, T>(x - x_ref), integrators());
        return record(ColVec<NU, T>(u_ref - (K * z)));
    }

    /** @brief The law, then the integrators' error from @p r and the measurement, with the feedthrough @p fed (Ci·D·u) it carries added back. */
    template<size_t NMeas>
    [[nodiscard]] constexpr ColVec<NU, T> integral_law(
        const ColVec<NI, T>&    r,
        const ColVec<NXL, T>&   x,
        const ColVec<NXL, T>&   x_ref,
        const ColVec<NU, T>&    u_ref,
        const ColVec<NMeas, T>& meas,
        const ColVec<NI, T>&    fed
    ) {
        const ColVec<NU, T> u = reference_law(x, x_ref, u_ref);
        integral.e = ColVec<NI, T>((r - (integral.Ci * meas)) + fed);
        return u;
    }

    /** @brief The estimator anchored on the target for @p r about (x*, u*), then u = u_a − K·(the deviation the command meets). */
    [[nodiscard]] constexpr ColVec<NU, T> anchored_law(const ColVec<NZ, T>& r, const ColVec<NXL, T>& x_ref, const ColVec<NU, T>& u_ref)
        requires(ND > 0)
    {
        target.anchor(estimator, r, x_ref, u_ref);
        return record(ColVec<NU, T>(estimator.anchor_input() - (K * detail::leading<NXL>(estimator.at_command()))));
    }

    /** @brief Integrate the held error less the applied input's feedthrough Ci·D·u, unless the clamp would wind it up: integrator i holds when −K(j, i)·eᵢ points into a clamped input j. */
    constexpr void integrate()
        requires(NI > 0)
    {
        const ColVec<NI, T> e = integral.e - (integral.Ci_u * applied_.applied);
        ColVec<NU, T>       excess{};
        for (size_t j = 0; j < NU; ++j) {
            const T u_req = applied_.command(j) + applied_.feedforward(j);
            const T diff = applied_.applied(j) - u_req;
            const T tol = T{4} * eps<T>() * (damp::abs(applied_.applied(j)) + damp::abs(u_req));
            excess(j) = (damp::abs(diff) > tol) ? diff : T{0};
        }
        for (size_t i = 0; i < NI; ++i) {
            bool winds_up = false;
            for (size_t j = 0; j < NU; ++j) {
                winds_up = winds_up || (K(j, NXL + i) * e(i) * excess(j) > T{0});
            }
            if (!winds_up) {
                integral.xi(i) += integral.step * e(i);
            }
        }
        integral.e = ColVec<NI, T>{};
    }
};

} // namespace detail
} // namespace damp
