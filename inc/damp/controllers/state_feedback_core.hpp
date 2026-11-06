// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file state_feedback_core.hpp
 * @brief What the LQ runtimes share: detail::StateFeedbackCore, the design result StateFeedbackResult,
 *        the weights, IntegralAction, and the design machinery behind lqr, lqi, lqg, lqgi and oflqg
 *
 * Each runtime set is its own header: lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp (design function,
 * result, runtime). This header is what they include.
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/options.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/luenberger.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {

/** @brief The LQ design result the runtimes run (defined below, with the design machinery). */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ>
struct StateFeedbackResult;

/** @brief The filter half of an estimated design (defined in kalman.hpp). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND>
struct KalmanResult;

} // namespace design

/**
 * @brief The estimator of a law that measures its state: none, an empty type, so a
 *        `[[no_unique_address]]` member of it takes no storage.
 */
struct NoEstimator {
    /** @brief Any scalar: still nothing. */
    template<typename U>
    [[nodiscard]] constexpr NoEstimator as() const {
        return {};
    }
};

namespace detail {

/** @brief What an estimator type carries: none (NoEstimator). */
template<typename E>
struct estimator_traits {
    static constexpr bool   present = false; ///< A filter runs
    static constexpr size_t NXE = 0;         ///< Filter states before the disturbance
    static constexpr size_t ND = 0;          ///< Disturbance states
};

/** @brief A Luenberger observer's sizes. */
template<size_t NXK, size_t NU, size_t NY, typename T, size_t ND_>
struct estimator_traits<Luenberger<NXK, NU, NY, T, ND_>> {
    static constexpr bool   present = true;  ///< A filter runs
    static constexpr size_t NXE = NXK - ND_; ///< Filter states before the disturbance
    static constexpr size_t ND = ND_;        ///< Disturbance states
};

/**
 * @brief The integral action of a law: the integrators, the integrated outputs Ci (on the
 *        measurement) and their feedthrough Ci·D, the per-sample step, the error held for the next
 *        predict, and with an estimator the measurement held from the last update (NYH = NY; 0: none).
 */
template<size_t NI, size_t NC, size_t NU, size_t NYH, typename T>
struct IntegralAction {
    using HeldY = maybe<(NYH > 0), ColVec<NYH, T>, OffSlot::held_y>;

    ColVec<NI, T>               xi{};    ///< Integral of the tracking error
    Matrix<NI, NC, T>           Ci{};    ///< Integrated outputs: e = r − Ci·(measurement)
    Matrix<NI, NU, T>           Ci_u{};  ///< Ci·D: the integrated outputs' feedthrough of the input
    T                           step{1}; ///< Per-sample integrator gain: the design's Ts (xi += Ts·e ≈ ∫e dt)
    ColVec<NI, T>               e{};     ///< Error held from control for the next predict, before the input's feedthrough
    [[no_unique_address]] HeldY y{};     ///< Measurement from the last update (with an estimator)
};

/** @brief Not a StateFeedbackResult. */
template<typename R>
struct is_state_feedback_result : std::false_type {};

/** @brief A StateFeedbackResult. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ>
struct is_state_feedback_result<design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ>> : std::true_type {};

/** @brief A gain design result other than StateFeedbackResult: K, success and Ts (PlaceResult, LQRResult). */
template<typename R>
concept gain_result = requires(const R& r) {
    r.success;
    r.Ts;
    r.K;
    requires MatrixLike<std::remove_cvref_t<decltype(r.K)>>;
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

/** @brief The runtime of a design's filter half: none for a measured state. */
template<typename Kalman, typename T>
struct runtime_estimator {
    using type = NoEstimator; ///< NoEstimator
};

/** @brief A KalmanResult runs as the Luenberger observer of its sizes (its fixed gain M), in the runtime's scalar @p T. */
template<size_t NXK, size_t NU, size_t NY, typename U, size_t NW, size_t NV, size_t ND, typename T>
struct runtime_estimator<design::KalmanResult<NXK, NU, NY, U, NW, NV, ND>, T> {
    using type = Luenberger<NXK, NU, NY, T, ND>; ///< The observer and the disturbance
};

/** @brief The runtime of a design's filter half (@p Kalman: KalmanResult or NoEstimator), in the scalar @p T. */
template<typename Kalman, typename T>
using runtime_estimator_t = runtime_estimator<Kalman, T>::type;

} // namespace detail

namespace detail {

/**
 * @brief The implementation the state-feedback runtimes share: u = u* − K·(x̂ − x*), every feature
 *        optional and free when off.
 *
 * LQR, LQI, LQG, LQGI and OFLQG wrap it, each exposing one update, one control and one step.
 *
 * | Feature | Parameter | Off |
 * | --- | --- | --- |
 * | integral action on r − Ci·y | NI, NY | 0 |
 * | estimator (Luenberger with design::kalman's gain, disturbance model) | Est | NoEstimator: the state is measured |
 * | steady-state target for r, under the disturbance model if there is one (SteadyStateTarget) | NZ | 0 |
 *
 * One tick is @ref step(): @ref predict (the time update with the input @ref back_calculate
 * recorded, else the last command), then @ref update with the measurement at sample k (the
 * state itself without an estimator), then @ref control, which returns u_k. Discrete-time
 * index: predict uses the input held over [k−1, k); update consumes m_k; control emits u_k.
 * The low-latency ISR is update(m_k) → control → write → back_calculate(u) → predict(): the
 * trailing predict is the next sample's time update. For a measured-state LQR, predict is
 * empty, so step(x) equals update then control.
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
 * @tparam Est Estimator: Luenberger (the fixed Kalman gain of design::kalman), or NoEstimator for a measured state
 * @tparam NZ  Controlled outputs of the steady-state target (0: none)
 */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Est, size_t NZ>
struct StateFeedbackCore {
    using Traits = detail::estimator_traits<Est>;
    static constexpr bool   kEstimated = Traits::present;         ///< A filter estimates the state
    static constexpr size_t NXL = NX;                             ///< The law's model: the plant
    static constexpr size_t ND = Traits::ND;                      ///< Disturbance states of the estimator
    static constexpr size_t NXE = kEstimated ? Traits::NXE : NXL; ///< The filter's model before the disturbance (the law's without one)
    static constexpr size_t NM = kEstimated ? NY : NX;            ///< The measurement: y, or the state itself
    static constexpr size_t NR = (NI > 0) ? NI : NZ;              ///< The reference r: integrated or tracked outputs (0: none)
    static_assert(!kEstimated || Traits::NXE == NXL, "The estimator is on the law's model.");
    static_assert(NI == 0 || NZ == 0, "Integral action and a disturbance model are two offset-free choices: lqgi or oflqg, not both.");
    static_assert(ND == 0 || NZ > 0, "A disturbance model tracks through its steady-state target.");
    static constexpr bool kRecord = (NI > 0) || kEstimated; ///< The law keeps its command and the applied input
    using Record = std::conditional_t<kRecord, detail::AppliedRecord<NU, T>, detail::NoAppliedRecord>;

    using Target = detail::maybe<(NZ > 0), SteadyStateTarget<NXE, NU, ND, NZ, T>, detail::OffSlot::target>;
    using Integral = detail::maybe<(NI > 0), detail::IntegralAction<NI, NM, NU, (kEstimated ? NY : 0), T>, detail::OffSlot::integral>;
    using HeldX = detail::maybe<!kEstimated, ColVec<NX, T>, detail::OffSlot::held_x>;

    Matrix<NU, NXL + NI, T>        K{};         ///< Gain on [x; xi]
    [[no_unique_address]] Est      estimator{}; ///< The filter (none for a measured state)
    [[no_unique_address]] Target   target{};    ///< The steady-state target (NZ > 0)
    [[no_unique_address]] Integral integral{};  ///< The integrators (NI > 0)

    constexpr StateFeedbackCore() = default;

    /** @brief From a gain on the law's state (no estimator): K = @p K_. */
    constexpr explicit StateFeedbackCore(const Matrix<NU, NXL + NI, T>& K_)
        requires(!kEstimated)
        : K(K_) {}

    /**
     * @brief From a gain design (design::place's PlaceResult, LQRResult; any scalar): converts K.
     *
     * A failed or continuous (not discrete) result does not compile when constant-evaluated; at run
     * time either leaves K = 0.
     */
    template<typename R>
        requires(detail::gain_result<R> && !kEstimated && NI == 0 && NZ == 0)
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
    template<size_t NXR, size_t NUR, typename U, size_t NIR, size_t NYR, typename Kalman, size_t NZR>
    constexpr StateFeedbackCore(const design::StateFeedbackResult<NXR, NUR, U, NIR, NYR, Kalman, NZR>& result) { // NOLINT
        constexpr bool same = NXR == NX && NUR == NU && NIR == NI && NYR == NY && NZR == NZ && std::is_same_v<Est, detail::runtime_estimator_t<Kalman, T>>;
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
    }

    /**
     * @brief Measurement update: the filter's with the outputs @p m (the integrators read them at
     *        the next control), or, with no filter, the state @p m held for the next control.
     */
    constexpr void update(const ColVec<NM, T>& m) {
        if constexpr (kEstimated) {
            if constexpr (NI > 0) {
                integral.y = m;
            }
            estimator.update(m);
        } else {
            x_ = m;
        }
    }

    /**
     * @brief The law at the last update's measurement, about the reference point (@p x_ref, @p u_ref):
     *        integral action on @p r (NI > 0), or the steady-state target for @p r (NZ > 0, under d̂ when ND > 0).
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref) {
        if constexpr (ND > 0) {
            return anchored_law(r, x_ref, u_ref);
        } else if constexpr (NZ > 0) {
            const auto p = target.equilibrium(r);
            return reference_law(law_state(), ColVec<NX, T>(x_ref + p.x), ColVec<NU, T>(u_ref + p.u));
        } else if constexpr (NI > 0 && kEstimated) {
            // y was measured under the input the last predict applied: take its feedthrough out.
            return integral_law(r, plant_estimate(), x_ref, u_ref, integral.y, ColVec<NI, T>(integral.Ci_u * applied_.applied));
        } else if constexpr (NI > 0) {
            return integral_law(r, x_, x_ref, u_ref, x_, ColVec<NI, T>{});
        } else {
            static_cast<void>(r);
            return reference_law(law_state(), x_ref, u_ref);
        }
    }

    /** @brief One tick: predict(), update(@p m), then control(@p r, @p x_ref, @p u_ref). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& m) {
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
    }

private:
    [[no_unique_address]] Record applied_{}; ///< The last command and the input back_calculate recorded
    [[no_unique_address]] HeldX  x_{};       ///< The measured state from the last update (no estimator)

    /** @brief Whether a design's gain runs on this tick: it succeeded and is discrete. Otherwise it does not compile when constant-evaluated (the error names its status) and the law keeps K = 0. */
    template<typename R>
    [[nodiscard]] static constexpr bool usable(const R& result, bool discrete) {
        const bool ok = result.success && discrete;
        if (!ok && std::is_constant_evaluated()) {
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

    [[nodiscard]] constexpr ColVec<NX, T> law_state() const {
        if constexpr (kEstimated) {
            return plant_estimate();
        } else {
            return x_;
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

    /** @brief u = u* − K·z, z = [x − x*; xi]. */
    [[nodiscard]] constexpr ColVec<NU, T> reference_law(const ColVec<NX, T>& x, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref) {
        const ColVec<NXL + NI, T> z = detail::stacked(ColVec<NX, T>(x - x_ref), integrators());
        return record(ColVec<NU, T>(u_ref - (K * z)));
    }

    /** @brief The law, then the integrators' error from @p r and the measurement, with the feedthrough @p fed (Ci·D·u) it carries added back. */
    template<size_t NMeas>
    [[nodiscard]] constexpr ColVec<NU, T> integral_law(
        const ColVec<NI, T>&    r,
        const ColVec<NX, T>&    x,
        const ColVec<NX, T>&    x_ref,
        const ColVec<NU, T>&    u_ref,
        const ColVec<NMeas, T>& meas,
        const ColVec<NI, T>&    fed
    ) {
        const ColVec<NU, T> u = reference_law(x, x_ref, u_ref);
        integral.e = ColVec<NI, T>((r - (integral.Ci * meas)) + fed);
        return u;
    }

    /** @brief The estimator anchored on the target for @p r about (x*, u*), then u = u_a − K·(the deviation the command meets). */
    [[nodiscard]] constexpr ColVec<NU, T> anchored_law(const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref)
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

namespace detail {

/** @brief The runtime an LQ design runs as, in the scalar @p T (specialized by lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp). */
template<typename R, typename T>
struct lq_runtime {};

} // namespace detail

/** @brief Law<kRes, T> of an LQ or gain design: the runtime it runs as (LQR, LQI, LQG, LQGI, OFLQG), default-constructed from the design. */
template<const auto& Res, typename T>
    requires requires {
        typename detail::lq_runtime<std::remove_cvref_t<decltype(Res)>, T>::type;
        Res.success;
    }
struct Law<Res, T> : detail::lq_runtime<std::remove_cvref_t<decltype(Res)>, T>::type {
    using Runtime = detail::lq_runtime<std::remove_cvref_t<decltype(Res)>, T>::type; ///< The runtime the design runs as

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res.template as<T>()} {}

    using Runtime::operator=;
};


} // namespace damp

namespace damp {

namespace design {

namespace detail {

/**
 * @brief An option left out of StateFeedbackOptions: that feature is off.
 *
 * Bare braces for the options, `{.noise = …}`, give every field this type, and the
 * conversion says to name the struct: `design::StateFeedbackOptions{.noise = …}`.
 */
struct Unset {
    constexpr Unset() = default;

    /** @brief Bare braces reach this: it does not compile, and says to name the options struct. */
    template<typename U>
    constexpr Unset(const U& value) { // NOLINT(google-explicit-constructor): implicit, so bare braces reach the message
        static_cast<void>(value);
        static_assert(sizeof(U) == 0, "Name the options struct: design::StateFeedbackOptions{.noise = ...}. Bare braces cannot carry an option's type.");
    }
};

/**
 * @brief The features of a state-feedback law, one named field each; leave a field out and the
 *        feature is off.
 *
 * Write it by name, with the fields you use, in this order:
 *
 * | Field | Value | The law gets |
 * | --- | --- | --- |
 * | `.noise` | design::KalmanNoise{.process, .measurement} | a Kalman filter: the law measures y, not x (LQG) |
 * | `.offset_free` | design::IntegralAction{}, or a disturbance model: InputDisturbance, OutputDisturbance, PlantDisturbance or DisturbanceModel | no steady-state error on the outputs r tracks (LQI, LQGI, OFLQG): by integrators on r − y, or by a disturbance the filter estimates (needs `.noise`) |
 * | `.track` | design::TrackedOutputs{} | the outputs r holds through the model's steady-state target, with no offset-free action (LQR, LQG) |
 *
 * @code
 * design::StateFeedbackOptions{.noise = design::KalmanNoise{.process = Qn, .measurement = Rn}, .offset_free = design::InputDisturbance{}}
 * @endcode
 *
 * Each field's type is a template argument deduced from what you put in it, so the sizes an
 * option carries (the outputs listed, a latency) size the result.
 */
template<typename NoiseT = Unset, typename OffsetFreeT = Unset, typename TrackT = Unset>
struct StateFeedbackOptions {
    NoiseT      noise{};       ///< design::KalmanNoise{.process, .measurement}: estimate the state with a Kalman filter
    OffsetFreeT offset_free{}; ///< design::IntegralAction{}, or a disturbance model the filter estimates: no steady-state error
    TrackT      track{};       ///< design::TrackedOutputs{}: the outputs r holds through the steady-state target
};

} // namespace detail

/**
 * @brief Offset-free by integral action: integrators on r − y, every measured output
 *        (`IntegralAction{}`) or the listed ones (`IntegralAction<2>{{0, 3}}`, in that order:
 *        the references r follow it).
 * @tparam NI Outputs integrated, listed in @ref outputs (0: every measured output)
 */
template<size_t NI = 0>
struct IntegralAction {
    damp::array<size_t, NI> outputs{}; ///< Indices of the integrated outputs, each once (none listed: all)
};

/**
 * @struct LQRWeights
 * @brief The cost of a state-feedback design: what to penalize and how much
 *
 * Weights in the plant's own time domain: a cost integral on a continuous plant (with a Ts
 * argument it is discretized exactly), per step on a discrete one. Braces size it from the
 * call: `{.Q = Q, .R = R}`, and `.Q_integral` with integral action.
 *
 * | Field | Penalizes | Size |
 * | --- | --- | --- |
 * | Q | the plant states x | NX × NX, positive semidefinite |
 * | R | the inputs u | NU × NU, positive definite |
 * | N | the cross term 2xᵀNu (default zero) | NX × NU |
 * | Q_integral | the integrals ∫e of the integrated outputs (lqi, lqgi) | NI × NI, positive definite |
 *
 * The integrator weight is block-diagonal with Q: the closed loop still couples them, so
 * raising Q_integral speeds the integral action and moves the regulator's poles too. An integrator Q_integral does not weight is a mode at z = 1
 * the cost does not see, with no stabilizing gain, so a design with integral action and a
 * Q_integral that is not positive definite (left out of the braces, for one) is refused with
 * DesignStatus::bad_weights, as is an R that is not positive definite, on every route.
 */
template<size_t NX, size_t NU, typename T = double, size_t NI = 0>
struct LQRWeights {
    Matrix<NX, NX, T> Q{};          ///< On the plant states x (positive semidefinite)
    Matrix<NU, NU, T> R{};          ///< On the inputs u (positive definite)
    Matrix<NX, NU, T> N{};          ///< Cross term between x and u (default zero)
    Matrix<NI, NI, T> Q_integral{}; ///< On the integrals ∫e of the integrated outputs (positive definite)
};

/** @brief Sizes from the matrices: `LQRWeights{Q, R}`. */
template<size_t NX, size_t NU, typename T>
LQRWeights(Matrix<NX, NX, T>, Matrix<NU, NU, T>) -> LQRWeights<NX, NU, T>;

/** @brief Sizes from the matrices: `LQRWeights{Q, R, N}`. */
template<size_t NX, size_t NU, typename T>
LQRWeights(Matrix<NX, NX, T>, Matrix<NU, NU, T>, Matrix<NX, NU, T>) -> LQRWeights<NX, NU, T>;

/** @brief Sizes from the matrices: `LQRWeights{Q, R, N, Q_integral}`. */
template<size_t NX, size_t NU, typename T, size_t NI>
LQRWeights(Matrix<NX, NX, T>, Matrix<NU, NU, T>, Matrix<NX, NU, T>, Matrix<NI, NI, T>) -> LQRWeights<NX, NU, T, NI>;

/** @brief The integral action of a design: the integrated outputs, their maps of the state and the input, and the per-sample step. */
template<size_t NI, size_t NX, size_t NU, size_t NY, typename T = double>
struct IntegralDesign {
    Matrix<NI, NY, T> Ci{};    ///< The integrated outputs: e = r − Ci·y, y = C·x + D·u
    Matrix<NI, NX, T> Ci_x{};  ///< Ci·C: the integrated outputs as a map of the plant state (what a law that measures its state integrates)
    Matrix<NI, NU, T> Ci_u{};  ///< Ci·D: the integrated outputs' feedthrough of the input applied over the period
    T                 step{1}; ///< Per-sample integrator gain: the design's Ts (xi += Ts·e ≈ ∫e dt)

    /** @brief The same integral action in the scalar U. */
    template<typename U>
    [[nodiscard]] constexpr IntegralDesign<NI, NX, NU, NY, U> as() const {
        return {Ci.template as<U>(), Ci_x.template as<U>(), Ci_u.template as<U>(), static_cast<U>(step)};
    }
};

/**
 * @struct StateFeedbackResult
 * @brief A state-feedback design: the gain, and the parts the law's options add
 *
 * Its template arguments are the runtime's: LQR, LQI, LQG, LQGI or OFLQG, the one its design names, deduces them from it. The
 * gain is on the law's state [x; xi]: the plant, then the integrators (NI > 0).
 *
 * | Member | What it is | Present |
 * | --- | --- | --- |
 * | K, S, e | the gain, the Riccati solution, the closed-loop poles | always |
 * | success, status, Ts | ok or why not; the law's sample time | always |
 * | kalman | the Kalman filter (design::kalman's result) | lqg, lqgi, oflqg |
 * | integral | the integrated outputs Ci, Ci·C, Ci·D and the step | NI > 0 |
 * | target | the steady-state target (design::steady_state_target's result) | oflqg; lqr and lqg with TrackedOutputs (NZ > 0) |
 *
 * @tparam NX     Plant states
 * @tparam NI     Integrators (0: none)
 * @tparam NY     Measured outputs (NX for a law without integrators or estimator)
 * @tparam Kalman The filter's design: KalmanResult, or NoEstimator for a measured state
 * @tparam NZ     Controlled outputs of the target (0: none)
 *
 * @see "Optimal Control" (Anderson & Moore, 1990), §4 and §8
 */
template<size_t NX, size_t NU, typename T = double, size_t NI = 0, size_t NY = NX, typename Kalman = NoEstimator, size_t NZ = 0>
struct [[nodiscard]] StateFeedbackResult {
    using Filter = ::damp::detail::estimator_traits<::damp::detail::runtime_estimator_t<Kalman, T>>;
    static constexpr bool   kEstimated = Filter::present; ///< A Kalman filter estimates the state
    static constexpr size_t NXL = NX;                     ///< The law's model: the plant
    static constexpr size_t NK = NXL + NI;                ///< The gain's state: the plant and the integrators
    static constexpr size_t ND = Filter::ND;              ///< Disturbance states of the filter
    using Integral = ::damp::detail::maybe<(NI > 0), IntegralDesign<NI, NX, NU, NY, T>, ::damp::detail::OffSlot::design_integral>;
    using Target = ::damp::detail::maybe<(NZ > 0), SteadyStateTargetResult<NXL, NU, ND, NZ, T>, ::damp::detail::OffSlot::design_target>;

    Matrix<NU, NK, T>                 K{};                                 ///< Optimal gain: u = −K·[x; xi]
    Matrix<NK, NK, T>                 S{};                                 ///< Riccati solution (positive semidefinite)
    damp::array<damp::complex<T>, NK> e{};                                 ///< Closed-loop poles of the law's model
    bool                              success{false};                      ///< status == DesignStatus::ok
    DesignStatus                      status{DesignStatus::not_converged}; ///< ok, or the first failure: the gain's, the filter's, the target's
    T                                 Ts{0};                               ///< Sample time of the law: 0 continuous (e in the s-plane), > 0 discrete
    [[no_unique_address]] Kalman      kalman{};                            ///< The Kalman filter (lqg, lqgi, oflqg)
    [[no_unique_address]] Integral    integral{};                          ///< The integral action (NI > 0)
    [[no_unique_address]] Target      target{};                            ///< The steady-state target (NZ > 0)

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same design in the scalar @p U: `res.as<float>()` for the embedded runtime. */
    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        using Out = StateFeedbackResult<NX, NU, U, NI, NY, std::remove_cvref_t<decltype(kalman.template as<U>())>, NZ>;
        Out out{K.template as<U>(), S.template as<U>(), array_as<damp::complex<U>>(e), success, status, static_cast<U>(Ts), kalman.template as<U>()};
        if constexpr (NI > 0) {
            out.integral = integral.template as<U>();
        }
        if constexpr (NZ > 0) {
            out.target = target.template as<U>();
        }
        return out;
    }

    /**
     * @brief Check if the closed loop is stable: every pole strictly inside the unit circle
     *        (Ts > 0), or strictly in the left half plane (Ts = 0). False for a failed design.
     */
    [[nodiscard]] constexpr bool is_stable() const {
        if (!success) {
            return false; // a failed design has no closed loop (e is zero)
        }
        for (size_t i = 0; i < NK; ++i) {
            if (Ts > T{0} ? (e[i].abs() >= T{1}) : (e[i].real() >= T{0})) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief The law as a discrete state-space block, for closing the loop in analysis
     *
     * The realization of the runtime's tick, by the law's features:
     *
     * | Law | State | Inputs | Output |
     * | --- | --- | --- | --- |
     * | measured state | none | x | u = −Kx |
     * | measured state, integrators | ξ | [r; x] | u |
     * | Kalman filter | x̂(k−1\|k−1) | y | u |
     * | Kalman filter, integrators | [x̂(k\|k−1); ξ; u_prev] | [r; y] | u |
     *
     * With integrators, @f$ u = -K_x x - K_i \xi @f$ reads the current ξ, then
     * @f$ \xi \leftarrow \xi + h(r - C_i (C x + D u)) @f$ (h = integral.step), the measurement
     * less the feedthrough of the input it was taken under standing in for C x with a filter. With a filter it is the
     * current estimator: predict with the last command, the measurement update with the
     * innovation gain M (KalmanResult::M), then the law,
     * @f[
     *   A_c = (I-MC)(A-BK) + MDK,\quad B_c = M,\quad C_c = -K A_c,\quad D_c = -KM.
     * @f]
     * A law with a disturbance model holds more than these
     * states, so it has no such block. The steady-state target of TrackedOutputs is a
     * feedforward on r outside this block: the block is the regulator (r = 0).
     */
    [[nodiscard]] constexpr auto to_ss() const
        requires(ND == 0)
    {
        if constexpr (!kEstimated && NI == 0) {
            return StateSpace<0, NX, NU, T>{.A = Matrix<0, 0, T>{}, .B = Matrix<0, NX, T>{}, .C = Matrix<NU, 0, T>{}, .D = -K, .Ts = Ts};
        } else if constexpr (!kEstimated) {
            return integral_ss();
        } else if constexpr (NI == 0) {
            return estimator_ss();
        } else {
            return integral_estimator_ss();
        }
    }

private:
    /**
     * @brief The law with integrators on a measured state: u = C_c ξ + D_c [r; x], C_c = −K_i,
     *        D_c = [0, −K_x], and ξ ← ξ + h·(r − Ci·C·x − Ci·D·u): A_c = I − h·Ci·D·C_c,
     *        B_c = [hI, −h·Ci·C] − h·Ci·D·D_c.
     */
    [[nodiscard]] constexpr StateSpace<NI, NI + NX, NU, T> integral_ss() const {
        const T                 h = integral.step;
        const Matrix<NU, NI, T> Cc = -K.template block<NU, NI>(0, NX);
        Matrix<NU, NI + NX, T>  Dc{};
        Dc.template block<NU, NX>(0, NI) = -K.template block<NU, NX>(0, 0);
        Matrix<NI, NI + NX, T> Bc{};
        Bc.template block<NI, NI>(0, 0) = Matrix<NI, NI, T>::identity() * h;
        Bc.template block<NI, NX>(0, NI) = -integral.Ci_x * h;
        const Matrix<NI, NU, T> hCiD = integral.Ci_u * h;
        const Matrix<NI, NI, T> Ac = Matrix<NI, NI, T>::identity() - (hCiD * Cc);
        return {.A = Ac, .B = Matrix<NI, NI + NX, T>(Bc - (hCiD * Dc)), .C = Cc, .D = Dc, .Ts = Ts};
    }

    /** @brief The current-estimator law: the filter's estimate x̂(k−1|k−1) is the state. */
    [[nodiscard]] constexpr StateSpace<NX, NY, NU, T> estimator_ss() const {
        const auto& sys = kalman.sys;
        const auto  ImMC = Matrix<NX, NX, T>::identity() - (kalman.M * sys.C);
        const auto  Ac = (ImMC * (sys.A - (sys.B * K))) + ((kalman.M * sys.D) * K);
        return {.A = Ac, .B = kalman.M, .C = -(K * Ac), .D = -(K * kalman.M), .Ts = Ts};
    }

    /** @brief The current-estimator law with integrators: the state is [x̂(k|k−1); ξ; u_prev], the inputs [r; y] (to_ss()). */
    [[nodiscard]] constexpr StateSpace<NX + NI + NU, NI + NY, NU, T> integral_estimator_ss() const {
        constexpr size_t        NS = NX + NI + NU;
        const auto&             sys = kalman.sys;
        const auto&             M = kalman.M;
        const Matrix<NU, NX, T> Kx = K.template block<NU, NX>(0, 0);
        const Matrix<NU, NI, T> Ki = K.template block<NU, NI>(0, NX);
        const auto              ImMC = Matrix<NX, NX, T>::identity() - (M * sys.C);
        const auto              MD = M * sys.D;
        const auto              KxImMC = Kx * ImMC;
        const auto              AmBKx = sys.A - (sys.B * Kx);
        const auto              KxM = Kx * M;
        const T                 h = integral.step;

        Matrix<NS, NS, T> Ac{};
        Ac.template block<NX, NX>(0, 0) = AmBKx * ImMC;
        Ac.template block<NX, NI>(0, NX) = -(sys.B * Ki);
        Ac.template block<NX, NU>(0, NX + NI) = -(AmBKx * MD);
        Ac.template block<NI, NI>(NX, NX) = Matrix<NI, NI, T>::identity();
        Ac.template block<NU, NX>(NX + NI, 0) = -KxImMC;
        Ac.template block<NU, NI>(NX + NI, NX) = -Ki;
        Ac.template block<NU, NU>(NX + NI, NX + NI) = Kx * MD;
        Matrix<NS, NI + NY, T> Bc{};
        Bc.template block<NI, NI>(NX, 0) = Matrix<NI, NI, T>::identity() * h;
        Bc.template block<NX, NY>(0, NI) = AmBKx * M;
        Bc.template block<NI, NY>(NX, NI) = -integral.Ci * h;
        Bc.template block<NU, NY>(NX + NI, NI) = -KxM;
        Matrix<NU, NS, T> Cc{};
        Cc.template block<NU, NX>(0, 0) = -KxImMC;
        Cc.template block<NU, NI>(0, NX) = -Ki;
        Cc.template block<NU, NU>(0, NX + NI) = Kx * MD;
        Matrix<NU, NI + NY, T> Dc{};
        Dc.template block<NU, NY>(0, NI) = -KxM;
        //! ξ ← ξ + h·(r − Ci·y + Ci·D·u_prev − Ci·D·u): the measurement less the feedthrough of the
        //! input it was taken under, then the feedthrough of the input the law applies, u = C_c s + D_c [r; y].
        const Matrix<NI, NU, T> hCiD = integral.Ci_u * h;
        Ac.template block<NI, NU>(NX, NX + NI) = hCiD;
        const Matrix<NI, NS, T>      fed_a = hCiD * Cc;
        const Matrix<NI, NI + NY, T> fed_b = hCiD * Dc;
        for (size_t i = 0; i < NI; ++i) {
            for (size_t j = 0; j < NS; ++j) {
                Ac(NX + i, j) -= fed_a(i, j);
            }
            for (size_t j = 0; j < NI + NY; ++j) {
                Bc(NX + i, j) -= fed_b(i, j);
            }
        }
        return {.A = Ac, .B = Bc, .C = Cc, .D = Dc, .Ts = Ts};
    }
};

namespace detail {

/** @brief Whether an option was left out. */
template<typename O>
inline constexpr bool is_unset = std::is_same_v<O, Unset>;

/** @brief Not integral action. */
template<typename>
inline constexpr bool is_integral_action = false;

/** @brief IntegralAction. */
template<size_t NI>
inline constexpr bool is_integral_action<IntegralAction<NI>> = true;

/** @brief Ci, an NI × NY matrix (matlab::lqi's form). */
template<size_t R, size_t C, typename U>
inline constexpr bool is_integral_action<Matrix<R, C, U>> = true;

/** @brief Not TrackedOutputs. */
template<typename>
inline constexpr bool is_tracked_outputs = false;

/** @brief TrackedOutputs. */
template<size_t NZ>
inline constexpr bool is_tracked_outputs<TrackedOutputs<NZ>> = true;

/** @brief A disturbance model with states: every disturbance spec but NoDisturbance. */
template<typename Spec>
inline constexpr bool is_disturbance_model = is_disturbance<Spec>::value && !std::is_same_v<Spec, NoDisturbance>;

/** @brief The fields of an options type: none for anything but StateFeedbackOptions. */
template<typename Options>
struct options_traits {
    static constexpr bool named = false; ///< Not a StateFeedbackOptions
    using noise = Unset;                 ///< .noise
    using offset_free = Unset;           ///< .offset_free
    using disturbance = Unset;           ///< .offset_free when it is a disturbance model
    using integrate = Unset;             ///< .offset_free when it is integral action
    using track = Unset;                 ///< .track
};

/** @brief The field types of a StateFeedbackOptions. */
template<typename NoiseT, typename OffsetFreeT, typename TrackT>
struct options_traits<StateFeedbackOptions<NoiseT, OffsetFreeT, TrackT>> {
    static constexpr bool named = true;                                                            ///< A StateFeedbackOptions
    using noise = NoiseT;                                                                          ///< .noise
    using offset_free = OffsetFreeT;                                                               ///< .offset_free
    using disturbance = std::conditional_t<is_disturbance_model<OffsetFreeT>, OffsetFreeT, Unset>; ///< .offset_free when it is a disturbance model
    using integrate = std::conditional_t<is_integral_action<OffsetFreeT>, OffsetFreeT, Unset>;     ///< .offset_free when it is integral action
    using track = TrackT;                                                                          ///< .track
};

/** @brief How many outputs a selection integrates on a plant with NY outputs: none for anything but a selection. */
template<typename Sel, size_t NY>
struct integrated_count {
    static constexpr bool   selection = false; ///< Not a selection
    static constexpr size_t value = 0;         ///< NI
};

/** @brief Every output (none listed), or the listed ones. */
template<size_t NI, size_t NY>
struct integrated_count<IntegralAction<NI>, NY> {
    static constexpr bool   selection = true;            ///< A selection
    static constexpr size_t value = (NI == 0) ? NY : NI; ///< NI
};

/** @brief The rows of Ci, an NI × NY matrix. */
template<size_t R, size_t C, typename U, size_t NY>
struct integrated_count<Matrix<R, C, U>, NY> {
    static constexpr bool   selection = (C == NY); ///< A selection when it has a column per output
    static constexpr size_t value = R;             ///< NI
};

/** @brief NI of a design's options. */
template<size_t NY, typename Options>
inline constexpr size_t integrator_count = integrated_count<typename options_traits<Options>::integrate, NY>::value;

/** @brief The weights of a design with @p Options; not deduced, so braces take their sizes from the plant. */
template<size_t NX, size_t NU, size_t NY, typename T, typename Options>
using weights_for = std::type_identity_t<LQRWeights<NX, NU, T, integrator_count<NY, Options>>>;

/** @brief Whether a disturbance spec fits the plant: a DisturbanceModel on its states and outputs, in its scalar. */
template<typename Spec, size_t NX, size_t NY, typename T>
inline constexpr bool disturbance_fits = true;

/** @brief A DisturbanceModel fits a plant of its sizes and scalar. */
template<size_t NXD, size_t NYD, size_t ND, size_t NZ, typename U, size_t NX, size_t NY, typename T>
inline constexpr bool disturbance_fits<DisturbanceModel<NXD, NYD, ND, NZ, U>, NX, NY, T> = (NXD == NX) && (NYD == NY) && std::is_same_v<U, T>;

/** @brief Each field of the options is its own kind of option: a message per mistake. @return Every field is */
template<size_t NX, size_t NY, typename T, typename Options>
consteval bool fields_valid() {
    using O = options_traits<Options>;
    using F = O::offset_free;
    constexpr bool noise = is_unset<typename O::noise> || is_kalman_noise<typename O::noise>::value;
    constexpr bool choice = is_unset<F> || is_integral_action<F> || is_disturbance_model<F>;
    constexpr bool fits = disturbance_fits<typename O::disturbance, NX, NY, T>;
    constexpr bool integrate = is_unset<typename O::integrate> || integrated_count<typename O::integrate, NY>::selection;
    constexpr bool track = is_unset<typename O::track> || is_tracked_outputs<typename O::track>;
    constexpr bool one_reference = is_unset<typename O::track> || is_unset<F>;
    static_assert(noise, "The noise is design::KalmanNoise{.process = ..., .measurement = ...}.");
    static_assert(
        choice,
        "options.offset_free is design::IntegralAction{} (integrators on r - y), or a disturbance model the Kalman filter estimates: "
        "design::InputDisturbance{}, OutputDisturbance{}, PlantDisturbance{} or DisturbanceModel{.Bd, .Cd}. Leave it out for neither."
    );
    static_assert(fits, "The DisturbanceModel is on the plant's states and outputs, in the plant's scalar.");
    static_assert(integrate, "The integrated outputs Ci·y: Ci is an NI x NY matrix, a column per measured output.");
    static_assert(track, "options.track is design::TrackedOutputs{} (every output) or TrackedOutputs<NZ>{{i, ...}}.");
    static_assert(one_reference, "Integral action and a disturbance model already hold r: .track is for a law with neither (lqr, lqg).");
    return noise && choice && fits && integrate && track && one_reference;
}

/** @brief The features go together and the noise is sized to the plant: a message per mistake. @return They do and it is */
template<size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Options>
consteval bool features_valid() {
    using O = options_traits<Options>;
    constexpr bool   kNoise = !is_unset<typename O::noise>;
    constexpr size_t ND = [] {
        if constexpr (is_disturbance<typename O::disturbance>::value) {
            return disturbance_traits<typename O::disturbance, NU, NY, NW>::states;
        } else {
            return size_t{0};
        }
    }();
    constexpr bool estimated = kNoise || is_unset<typename O::disturbance>;
    constexpr bool sized = !kNoise || std::is_same_v<typename O::noise, KalmanNoise<NW, NV, T, ND>>;
    static_assert(estimated, "A disturbance model is estimated by the Kalman filter: give its noise, with .disturbance the disturbance's random walk.");
    static_assert(
        sized || !estimated,
        "The noise is sized to the plant: .process NW x NW (the columns of G), .measurement NV x NV (the columns of H), in the plant's scalar; with a disturbance model, also .disturbance ND x ND (the random walk of d), else none."
    );
    return estimated && sized;
}

/** @brief The spec of the options' steady-state target: the disturbance model, else the tracked outputs; NoDisturbance without either. */
template<typename Options>
[[nodiscard]] constexpr auto disturbance_of(const Options& opts) {
    if constexpr (!is_unset<typename options_traits<Options>::disturbance>) {
        return opts.offset_free;
    } else if constexpr (!is_unset<typename options_traits<Options>::track>) {
        return opts.track;
    } else {
        static_cast<void>(opts);
        return NoDisturbance{};
    }
}

/** @brief The integral action of the options; Unset without it. */
template<typename Options>
[[nodiscard]] constexpr auto integral_action_of(const Options& opts) {
    if constexpr (is_unset<typename options_traits<Options>::integrate>) {
        static_cast<void>(opts);
        return Unset{};
    } else {
        return opts.offset_free;
    }
}

/** @brief No integral action: Ci has no rows. */
template<size_t NI, size_t NY, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NI, NY, T>> integrated_selector(const Unset& none) {
    static_cast<void>(none);
    return Matrix<NI, NY, T>{};
}

/** @brief Ci of integral action: the identity for every output, else the listed rows; nullopt for an index out of range or listed twice. */
template<size_t NI, size_t NY, typename T, size_t NS>
[[nodiscard]] constexpr damp::optional<Matrix<NI, NY, T>> integrated_selector(const IntegralAction<NS>& sel) {
    if constexpr (NS == 0) {
        static_cast<void>(sel);
        return Matrix<NY, NY, T>::identity();
    } else {
        const auto rows = selected_outputs<NY>(sel.outputs);
        if (!rows) {
            return damp::nullopt;
        }
        Matrix<NI, NY, T> Ci{};
        for (size_t i = 0; i < NI; ++i) {
            Ci(i, (*rows)[i]) = T{1};
        }
        return Ci;
    }
}

/** @brief Ci as given, in the plant's scalar. */
template<size_t NI, size_t NY, typename T, typename U>
[[nodiscard]] constexpr damp::optional<Matrix<NI, NY, T>> integrated_selector(const Matrix<NI, NY, U>& Ci) {
    return Ci.template as<T>();
}

/** @brief The discrete model of [x; xi] under the regulator (r = 0): xi[k+1] = xi[k] − h·Ci·(C·x[k] + Dy·u[k]). */
template<size_t NX, size_t NU, size_t NI, typename T>
struct IntegralModel {
    Matrix<NX + NI, NX + NI, T> A{}; ///< [A 0; −h·Ci·C I]
    Matrix<NX + NI, NU, T>      B{}; ///< [B; −h·Ci·Dy]
};

/** @brief Augment a discrete model with the integrators of the selected outputs (see IntegralModel). */
template<size_t NX, size_t NU, size_t NY, size_t NI, typename T>
[[nodiscard]] constexpr IntegralModel<NX, NU, NI, T> integral_model(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NY, NX, T>& C,
    const Matrix<NY, NU, T>& Dy,
    const Matrix<NI, NY, T>& Ci,
    T                        h
) {
    IntegralModel<NX, NU, NI, T> m{};
    m.A.template block<NX, NX>(0, 0) = A;
    m.A.template block<NI, NX>(NX, 0) = -(Ci * C) * h;
    m.A.template block<NI, NI>(NX, NX) = Matrix<NI, NI, T>::identity(); // integrator memory
    m.B.template block<NX, NU>(0, 0) = B;
    m.B.template block<NI, NU>(NX, 0) = -(Ci * Dy) * h;
    return m;
}

/**
 * @brief The continuous cost over one sample on [x; xi], for the runtime's integrator
 *
 * Between samples the plant runs continuously on the held input and the integrator ramps on the
 * error sampled at the start, ẋi = −Ci·(Cm·x_k + Dm·u). The Van Loan model is [x; xi; x_k], x_k
 * held over the period; folding x_k onto x gives the exact discrete weights of
 * ∫ ([x; xi]ᵀQ[x; xi] + 2[x; xi]ᵀNu + uᵀRu) dt on [x; xi] and u.
 */
template<size_t NX, size_t NU, size_t NY, size_t NI, typename T>
[[nodiscard]] constexpr damp::optional<LQRCost<NX + NI, NU, T>> sampled_integral_cost(
    const Matrix<NX, NX, T>&       A,
    const Matrix<NX, NU, T>&       B,
    const Matrix<NY, NX, T>&       Cm,
    const Matrix<NY, NU, T>&       Dm,
    const Matrix<NI, NY, T>&       Ci,
    const LQRCost<NX + NI, NU, T>& weights,
    T                              Ts
) {
    constexpr size_t NZ = NX + NI; // [x; xi], the weighted part
    constexpr size_t NS = NZ + NX; // [x; xi; x_k]
    const auto       CiC = Ci * Cm;

    Matrix<NS, NS, T> As{};
    Matrix<NS, NU, T> Bs{};
    As.template block<NX, NX>(0, 0) = A;
    for (size_t i = 0; i < NI; ++i) { // the integrator sees the sampled output
        for (size_t j = 0; j < NX; ++j) {
            As(NX + i, NZ + j) = -CiC(i, j);
        }
    }
    Bs.template block<NX, NU>(0, 0) = B;
    Bs.template block<NI, NU>(NX, 0) = -(Ci * Dm);
    Matrix<NS, NS, T> Qs{};
    Matrix<NS, NU, T> Ns{};
    Qs.template block<NZ, NZ>(0, 0) = weights.Q;
    Ns.template block<NZ, NU>(0, 0) = weights.N;
    const auto cost_opt = discretize_lqr_cost(As, Bs, Qs, weights.R, Ts, Ns);
    if (!cost_opt) {
        return damp::nullopt;
    }
    const auto& cost = *cost_opt;

    //! [x; xi; x_k] = P·[x; xi]: x_k is x at the sample
    Matrix<NS, NZ, T> P{};
    P.template block<NZ, NZ>(0, 0) = Matrix<NZ, NZ, T>::identity();
    P.template block<NX, NX>(NZ, 0) = Matrix<NX, NX, T>::identity();
    const Matrix<NZ, NS, T> Pt = P.t();
    const Matrix<NZ, NZ, T> Qf = Pt * cost.Q * P;
    const Matrix<NZ, NZ, T> Qft = Qf.t();
    return LQRCost<NZ, NU, T>{
        .Q = Matrix<NZ, NZ, T>((Qf + Qft) * static_cast<T>(0.5)),
        .R = cost.R,
        .N = Matrix<NZ, NU, T>(Pt * cost.N),
    };
}

/** @brief The weights on [x; xi]: Q and Q_integral on the diagonal, N on the plant's rows. */
template<size_t NX, size_t NU, typename T, size_t NI>
[[nodiscard]] constexpr LQRCost<NX + NI, NU, T> stacked_weights(const LQRWeights<NX, NU, T, NI>& w) {
    LQRCost<NX + NI, NU, T> out{.R = w.R};
    out.Q.template block<NX, NX>(0, 0) = w.Q;
    out.Q.template block<NI, NI>(NX, NX) = w.Q_integral;
    out.N.template block<NX, NU>(0, 0) = w.N;
    return out;
}

/** @brief The gain on the law's state. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, size_t NK>
[[nodiscard]] constexpr StateFeedbackResult<NX, NU, T, NI, NY> with_gain(StateFeedbackResult<NX, NU, T, NI, NY> out, const LQRResult<NK, NU, T>& lq) {
    out.K = lq.K;
    out.S = lq.S;
    out.e = lq.e;
    return ::damp::detail::set_status(out, lq.status);
}

/**
 * @brief The law on a measured state: the gain on [x; xi] for the cost of @p w, and the
 *        integrators of @p Ci (NI > 0). A continuous plant with no Ts and no integrators gets the
 *        continuous law (CARE).
 */
template<bool Discrete, size_t NI, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto law_design(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const LQRWeights<NX, NU, T, NI>& w, const Matrix<NI, NY, T>& Ci) {
    using Result = StateFeedbackResult<NX, NU, T, NI, ((NI > 0) ? NY : NX)>;
    if (!mat::cholesky(w.R).has_value()) {
        return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // R not positive definite, on every route
    }
    if constexpr (Discrete && NI == 0) {
        if (sys.Ts == T{0}) {
            const auto lq = continuous_lqr(sys.A, sys.B, w.Q, w.R, w.N);
            return with_gain(Result{}, lq);
        }
    }
    const auto m = sampled_model(sys, Ts, Discrete);
    if (m.status != DesignStatus::ok) {
        return ::damp::detail::refused<Result>(m.status);
    }
    const T h = Discrete ? sys.Ts : Ts;
    if constexpr (NI == 0) {
        LQRCost<NX, NU, T> cd{w.Q, w.R, w.N};
        if constexpr (!Discrete) {
            const auto c = discretize_lqr_cost(sys.A, sys.B, w.Q, w.R, Ts, w.N);
            if (!c) {
                return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // the plant was sampled: a non-finite weight
            }
            cd = *c;
        }
        return with_gain(Result{.Ts = h}, discrete_lqr(m.law.A, m.law.B, cd));
    } else {
        if (!mat::cholesky(w.Q_integral).has_value()) {
            return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // an unweighted integrator: a mode at z = 1 the cost does not see
        }
        // Forward-Euler integrator xi += h·e (MATLAB®'s discrete lqi): xi ≈ ∫e dt, so Q_integral
        // means the same on either route. The runtime integrates e = r − Ci·(C·x + D·u), u the
        // input applied over the period, and so does this model on both routes.
        const auto              im = integral_model(m.law.A, m.law.B, m.law.C, m.law.D, Ci, h);
        LQRCost<NX + NI, NU, T> cd = stacked_weights(w);
        if constexpr (!Discrete) {
            const auto c = sampled_integral_cost(sys.A, sys.B, sys.C, sys.D, Ci, cd, Ts);
            if (!c) {
                return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // the plant was sampled: a non-finite weight
            }
            cd = *c;
        }
        Result out{.Ts = h};
        out.integral = IntegralDesign<NI, NX, NU, NY, T>{Ci, Ci * sys.C, Ci * sys.D, h};
        return with_gain(out, discrete_lqr(im.A, im.B, cd));
    }
}

/** @brief The design on checked options: the law, then the filter and the target they ask for; status the first failure. */
template<bool Discrete, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Options>
[[nodiscard]] constexpr auto designed(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const weights_for<NX, NU, NY, T, Options>& w, const Options& opts) {
    using O = options_traits<Options>;
    constexpr size_t NI = integrator_count<NY, Options>;
    constexpr bool   kEstimated = !is_unset<typename O::noise>;
    constexpr bool   kTarget = !is_unset<typename O::disturbance> || !is_unset<typename O::track>;
    const auto       spec = disturbance_of(opts);

    const auto Ci = integrated_selector<NI, NY, T>(integral_action_of(opts));
    auto       law = law_design<Discrete, NI>(sys, Ts, w, Ci ? *Ci : Matrix<NI, NY, T>{});
    if (!Ci) {
        law = ::damp::detail::refused<decltype(law)>(DesignStatus::bad_weights); // an integrated output out of range or listed twice
    }
    if constexpr (!kEstimated && !kTarget) {
        return law;
    } else {
        constexpr size_t NDS = disturbance_traits<std::remove_cvref_t<decltype(spec)>, NU, NY, NW>::listed;
        constexpr size_t NZS = (NDS == 0) ? NY : NDS; // the spec's outputs, or every output
        constexpr size_t NZ = kTarget ? NZS : 0;
        const auto       filter = [&] {
            if constexpr (kEstimated) {
                return kalman_design<Discrete>(sys, Ts, opts.noise, spec);
            } else {
                return NoEstimator{};
            }
        }();
        using Result = StateFeedbackResult<NX, NU, T, NI, (kEstimated ? NY : NX), std::remove_cvref_t<decltype(filter)>, NZ>;
        Result out{law.K, law.S, law.e, law.success, law.status, law.Ts, filter};
        if constexpr (NI > 0) {
            out.integral = law.integral;
        }
        DesignStatus status = law.status;
        if constexpr (kEstimated) {
            status = (status == DesignStatus::ok) ? filter.status : status;
        }
        if constexpr (kTarget) {
            out.target = target_design<Discrete>(sys, Ts, spec);
            status = (status == DesignStatus::ok) ? out.target.status : status;
        }
        return ::damp::detail::set_status(out, status);
    }
}

/**
 * @brief The state-feedback design on @p sys in its own time domain (@p Discrete) or sampled at
 *        @p Ts, once the options pass their checks: a mistake stops here with its message alone.
 */
template<bool Discrete, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Options>
[[nodiscard]] constexpr auto state_feedback_design(
    const StateSpace<NX, NU, NY, T, NW, NV>&   sys,
    T                                          Ts,
    const weights_for<NX, NU, NY, T, Options>& w,
    const Options&                             opts
) {
    if constexpr (!options_traits<Options>::named) {
        static_assert(options_traits<Options>::named, "The options are design::StateFeedbackOptions{.noise = ..., .offset_free = ...}: one named field per feature.");
        return StateFeedbackResult<NX, NU, T>{};
    } else if constexpr (!fields_valid<NX, NY, T, Options>()) {
        return StateFeedbackResult<NX, NU, T>{}; // fields_valid says which field
    } else if constexpr (!features_valid<NU, NY, T, NW, NV, Options>()) {
        return StateFeedbackResult<NX, NU, T>{}; // features_valid says why
    } else {
        return designed<Discrete>(sys, Ts, w, opts);
    }
}

} // namespace detail

namespace detail {

/**
 * @brief The linear-quadratic state-feedback law for a plant, in the plant's own time domain
 *
 * The law u = −K·x minimizes the weighted cost of @p weights; @p options adds the features of the
 * law, one named field each (StateFeedbackOptions). A discrete plant (@c sys.Ts > 0) gets the
 * discrete law (DARE) at its own Ts; a continuous plant with no options gets the continuous law
 * (CARE, Ts = 0), for analysis: only a discrete law deploys. For a digital law from a continuous
 * plant, pass the sample time: state_feedback(sys, Ts, weights, options).
 *
 * @f[
 *   J = \sum_k \bigl( x_k^\top Q x_k + u_k^\top R u_k + 2 x_k^\top N u_k + \xi_k^\top Q_\xi \xi_k \bigr)
 * @f]
 *
 * With `.offset_free = IntegralAction{}`, the integrators ξ advance by xi += Ts·(r − Ci·y)
 * (MATLAB®'s discrete lqi, so ξ ≈ ∫e dt and Q_integral means what it does on the sampled route).
 * An input delay of up to one sample is part of the plant: design on the model
 * discretize(plant, Ts, InputDelay{θ}) returns, its register state a state like any other.
 * The other options are state_feedback(sys, Ts, weights, options)'s.
 *
 * Refused with DesignStatus::invalid_sample_time for a negative or non-finite sys.Ts, or a
 * continuous plant with any option; with DesignStatus::bad_weights for a listed output out
 * of range or listed twice.
 *
 * @note Compare with MATLAB®'s lqr(sys, Q, R, N), lqi, lqg and lqgtrack, which also follow the
 *       plant's domain.
 *
 * @param sys      Plant (A, B, C, D, G, H and Ts)
 * @param weights  {.Q, .R, .N, .Q_integral}, in the plant's time domain
 * @param options  design::StateFeedbackOptions{…}: the law's features (default none)
 * @return StateFeedbackResult: the runtime its design names runs it
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Options>
[[nodiscard]] constexpr auto state_feedback(
    const StateSpace<NX, NU, NY, T, NW, NV>&   sys,
    const weights_for<NX, NU, NY, T, Options>& weights,
    const Options&                             options
) {
    return state_feedback_design<true>(sys, T{0}, weights, options);
}

/**
 * @brief The sampled-data state-feedback law: a continuous plant and its cost, the digital law at
 *        @p Ts
 *
 * Discretizes the dynamics (ZOH) and the cost integral exactly (Van Loan, discretize_lqr_cost),
 * then solves the discrete problem, so the weights mean the same at any rate. @p options adds
 * the law's features, one named field each:
 *
 * | Option | The law |
 * | --- | --- |
 * | `.noise = KalmanNoise{.process, .measurement}` | estimates x with a Kalman filter from y (design::kalman(sys, Ts, noise)): the process noise a density, the measurement a variance per reading |
 * | `.offset_free = IntegralAction{}` or `IntegralAction<NI>{{i, …}}` | integrates r − Ci·y: ξ̇ = r − Ci·y, so ξ = ∫e dt and Q_integral holds at any rate |
 * | `.offset_free = InputDisturbance{}`, … | estimates a disturbance d and holds r through the steady-state target under it (OFLQG, below; needs `.noise`) |
 *
 * @par Integral action
 * The integrator advances by xi += Ts·(r − Ci·y) on each sampled error, and the cost integral
 * is taken over each sample with the plant running on the held input and ξ ramping on the held
 * error (Van Loan), so the closed-loop poles on the result are those of the deployed loop.
 * Constant references and disturbances leave no steady-state error on the integrated outputs.
 *
 * @par Input delay
 * A delay of up to one sample is a model of the plant, not an option of the law: design on
 * discretize(plant, Ts, InputDelay{θ}), whose register state v(k−1) the law measures or the
 * filter estimates like any other state.
 *
 * @par Offset-free control by a disturbance model
 * A disturbance model in `.offset_free` (InputDisturbance, OutputDisturbance, PlantDisturbance,
 * DisturbanceModel) augments the filter's model with integrating disturbances,
 * @f[
 *   x_{k+1} = A x_k + B u_k + B_d d_k,\qquad d_{k+1} = d_k,\qquad y_k = C x_k + C_d d_k + D u_k ,
 * @f]
 * the filter estimates @f$ [\hat x;\, \hat d] @f$ (the random walk of d is noise.disturbance, positive definite),
 * and the steady-state target holds the controlled outputs z (the rows the spec lists, all
 * without a list) at r under d̂:
 * @f[
 *   \begin{bmatrix} A-I & B \\ C_z & D_z \end{bmatrix}
 *   \begin{bmatrix} x_t \\ u_t \end{bmatrix}
 *   =
 *   \begin{bmatrix} -B_d\hat d \\ r - C_{d,z}\hat d \end{bmatrix},
 * @f]
 * so the law is @f$ u = u_t - K(\hat x - x_t) @f$ with the plain gain. Where d enters: an input
 * disturbance adds to u (B_d = B), an output disturbance is a bias on each reading (C_d = I),
 * PlantDisturbance is the plant's own G and H, DisturbanceModel is explicit. The augment must be
 * detectable (ND ≤ NY), else DesignStatus::not_detectable; a target needs a unique steady state
 * for NZ = NU, else DesignStatus::singular.
 *
 * Refused with DesignStatus::invalid_sample_time when @p Ts is not finite and positive, or when
 * @p sys is already discrete, even at the same Ts (a second discretization is the bug this
 * signature rules out); with DesignStatus::bad_weights for a listed output out of range or listed twice. The status is
 * the first failure: the gain's, the filter's, then the target's.
 *
 * @code
 * const auto res = design::state_feedback(plant, Ts, {.Q = Q, .R = R},
 *     design::detail::StateFeedbackOptions{.noise = design::KalmanNoise{.process = Qn, .measurement = Rn, .disturbance = Qd},
 *                                          .offset_free = design::InputDisturbance{}});
 * static constinit OFLQG law{res.as<float>()}; // the public form: design::oflqg
 * // each tick: auto u = law.step(r, y);
 * @endcode
 *
 * @note Compare with MATLAB®'s lqrd(sys, Q, R, Ts); lqi, lqg and lqgtrack on a continuous plant,
 *       followed by c2d; and for the offset-free law, its mpc with a disturbance model, without
 *       a QP.
 * @see "Optimal Control" (Anderson & Moore, 1990), §4 and §8
 * @see "Computing Integrals Involving the Matrix Exponential" (Van Loan, 1978)
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings,
 *      2003), AIChE Journal 49(2)
 * @see "Offset-free reference tracking with model predictive control" (Maeder, Borrelli &
 *      Morari, 2009), Automatica 45(9)
 *
 * @param sys      Continuous plant (sys.Ts = 0; A, B, C, D, G, H)
 * @param Ts       Sampling time of the law [s]
 * @param weights  Continuous {.Q, .R, .N, .Q_integral}: a cost integral
 * @param options  design::StateFeedbackOptions{…}: the law's features (default none)
 * @return StateFeedbackResult on the discretized plant: the runtime its design names runs it
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Options>
[[nodiscard]] constexpr auto state_feedback(
    const StateSpace<NX, NU, NY, T, NW, NV>&   sys,
    std::type_identity_t<T>                    Ts,
    const weights_for<NX, NU, NY, T, Options>& weights,
    const Options&                             options
) {
    return state_feedback_design<false>(sys, Ts, weights, options);
}

} // namespace detail

} // namespace design

} // namespace damp
