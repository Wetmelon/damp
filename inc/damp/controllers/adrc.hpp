// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file adrc.hpp
 * @brief Active Disturbance Rejection Control design and runtime
 *
 * Linear ADRC: an extended-state observer (ESO) estimates the plant output, its derivatives
 * and the total disturbance f on the model y⁽ⁿ⁾ = f + b₀u, and a PD law on the estimate cancels
 * f̂. The observer is @ref design::eso / @ref ESO; this header is that observer plus the PD law.
 * The runtime's ESO is the discrete one of Miklosovic, Radke and Gao: the ESO's linear
 * dynamics discretized by ZOH, in current-estimator form, with its poles at the image
 * @f$e^{-\omega_o T_s}@f$ of the continuous design's poles at @f$-\omega_o@f$.
 */

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/eso.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"

namespace damp {

namespace design {

/**
 * @struct ADRCResult
 * @brief Active Disturbance Rejection Control design result
 *
 * Contains the computed observer and controller gains for ADRC: the continuous ESO gains β
 * (for analysis, to_ss / to_tf) and, for Ts > 0, the discrete ESO the runtime ticks,
 * z̄_k = Φ ẑ_{k−1} + Γ u_{k−1}, ẑ_k = z̄_k + L (y_k − z̄_k,₀).
 * Use .as\<U\>() to convert for type conversion (e.g., double to float).
 *
 * @see "From PID to Active Disturbance Rejection Control" (Han, 1998)
 * @see "Discrete implementation and generalization of the extended state observer"
 *      (Miklosovic, Radke & Gao, 2006), ACC
 */
template<size_t NX, typename T = double>
struct ADRCResult {
    static_assert(NX == 1 || NX == 2, "ADRC supports 1st- and 2nd-order plants only");

    T wc{}; ///< Controller bandwidth
    T wo{}; ///< FixedGainObserver bandwidth
    T b0{}; ///< Plant gain

    ColVec<NX + 1, T> beta{}; ///< ESO gains

    T Kp{}; ///< Unscaled proportional gain (wc or wc²); tick divides by b0
    T Kd{}; ///< Unscaled derivative gain (0 or 2·wc); tick divides by b0

    Matrix<NX + 1, NX + 1, T> Phi{};   ///< ZOH state map Φ of the ESO (zero when Ts = 0)
    ColVec<NX + 1, T>         Gamma{}; ///< ZOH input map Γ on u_{k−1}
    ColVec<NX + 1, T>         L{};     ///< Current-estimator gain on the innovation y_k − z̄_k,₀

    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the design failed
    T            Ts{0};                               ///< [s] Sample time of the law: 0 continuous (analysis), > 0 the runtime's tick

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return ADRCResult<NX, U>{
            static_cast<U>(wc),
            static_cast<U>(wo),
            static_cast<U>(b0),
            beta.template as<U>(),
            static_cast<U>(Kp),
            static_cast<U>(Kd),
            Phi.template as<U>(),
            ColVec<NX + 1, U>(Gamma.template as<U>()),
            ColVec<NX + 1, U>(L.template as<U>()),
            success,
            status,
            static_cast<U>(Ts),
        };
    }

    /** @brief The observer this law ticks: @ref design::eso of the same ω_o, b₀ and T_s. */
    [[nodiscard]] constexpr ESOResult<NX, T> as_eso() const {
        return {
            .wo = wo,
            .b0 = b0,
            .beta = beta,
            .Phi = Phi,
            .Gamma = Gamma,
            .L = L,
            .success = success,
            .status = status,
            .Ts = Ts,
        };
    }

    /**
     * @brief Linear ADRC as a 2-DOF compensator: u from [r; y]
     *
     * Continuous ESO + PD (unsaturated). State is the ESO [ẑ₁, (ẑ₂,) f̂].
     * Failed designs return a zero system.
     */
    [[nodiscard]] constexpr StateSpace<NX + 1, 2, 1, T> to_ss_2dof() const {
        StateSpace<NX + 1, 2, 1, T> sys{};
        if (!success || !(b0 > T{0})) {
            return sys;
        }
        if constexpr (NX == 1) {
            sys.A(0, 0) = -(beta[0] + Kp);
            sys.A(1, 0) = -beta[1];
            sys.B(0, 0) = Kp;
            sys.B(0, 1) = beta[0];
            sys.B(1, 1) = beta[1];
            sys.C(0, 0) = -Kp / b0;
            sys.C(0, 1) = -T{1} / b0;
            sys.D(0, 0) = Kp / b0;
        } else {
            sys.A(0, 0) = -beta[0];
            sys.A(0, 1) = T{1};
            sys.A(1, 0) = -(beta[1] + Kp);
            sys.A(1, 1) = -Kd;
            sys.A(2, 0) = -beta[2];
            sys.B(1, 0) = Kp;
            sys.B(0, 1) = beta[0];
            sys.B(1, 1) = beta[1];
            sys.B(2, 1) = beta[2];
            sys.C(0, 0) = -Kp / b0;
            sys.C(0, 1) = -Kd / b0;
            sys.C(0, 2) = -T{1} / b0;
            sys.D(0, 0) = Kp / b0;
        }
        return sys;
    }

    /**
     * @brief SISO loop map C(s) = −U/Y at r = 0
     *
     * Unity-negative-feedback convention (same as PID / lead / PR `to_ss`):
     * u = C(s)(r − y). Bode / `series` / `feedback` consume this map.
     */
    [[nodiscard]] constexpr StateSpace<NX + 1, 1, 1, T> to_ss() const {
        const auto                  full = to_ss_2dof();
        StateSpace<NX + 1, 1, 1, T> sys{};
        sys.A = full.A;
        for (size_t i = 0; i < NX + 1; ++i) {
            sys.B(i, 0) = full.B(i, 1);
        }
        for (size_t j = 0; j < NX + 1; ++j) {
            sys.C(0, j) = -full.C(0, j);
        }
        sys.D(0, 0) = -full.D(0, 1);
        return sys;
    }

    /**
     * @brief Transfer function of @ref to_ss (Leverrier)
     *
     * NX = 1: C(s) = (K_p β₂ + (K_p β₁ + β₂)s) / (b₀ s (s + β₁ + K_p)).
     */
    [[nodiscard]] constexpr TransferFunction<NX + 2, NX + 2, T> to_tf() const {
        return ::damp::to_transfer_function(to_ss());
    }
};

/**
 * @brief Active Disturbance Rejection Control design
 *
 * The observer is @ref design::eso (poles at −ω_o, plant gain b₀, the discrete ESO at Ts > 0).
 * Controller gains are unscaled (b0 is applied only in the runtime law u = (u0 − f̂)/b0):
 * - 1st order: Kp = wc, Kd = 0
 * - 2nd order: Kp = wc², Kd = 2·wc
 *
 * Only plant orders NX = 1 and NX = 2 are supported (same surface as typical
 * linear ADRC toolboxes / MATLAB® discrete ADRC). Higher-order ADRC is out of
 * scope.
 *
 * Ts = 0 is the continuous design, for analysis (to_ss, to_tf); Φ, Γ and L stay zero.
 *
 * Requires wc > 0, wo > 0, b0 > 0 (all finite), else refused with invalid_parameter and
 * zero gains so the runtime cannot divide by an invalid b0; a negative or non-finite Ts
 * is refused with invalid_sample_time.
 *
 * @param wc Controller bandwidth
 * @param wo Observer bandwidth (same as @ref design::eso)
 * @param b0 Plant gain
 * @param Ts [s] Sample time of the runtime (0 = continuous)
 *
 * @return ADRCResult with computed gains
 *
 * @see design::eso — the observer this law ticks
 * @see "From PID to Active Disturbance Rejection Control" (Han, 1998)
 * @see "Discrete implementation and generalization of the extended state observer"
 *      (Miklosovic, Radke & Gao, 2006), ACC: the current discrete ESO and its gains
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr ADRCResult<NX, T> adrc(T wc, T wo, T b0, T Ts = T{0}) {
    static_assert(NX == 1 || NX == 2, "ADRC supports 1st- and 2nd-order plants only");

    const bool ok = damp::isfinite(wc) && damp::isfinite(wo) && damp::isfinite(b0) && (wc > T{0})
                 && (wo > T{0}) && (b0 > T{0});
    if (!ok) {
        return ::damp::detail::refused<ADRCResult<NX, T>>(DesignStatus::invalid_parameter);
    }

    const auto e = eso<NX>(wo, b0, Ts);
    if (!e.success) {
        return ::damp::detail::refused<ADRCResult<NX, T>>(e.status);
    }

    return ADRCResult<NX, T>{
        .wc = wc,
        .wo = e.wo,
        .b0 = e.b0,
        .beta = e.beta,
        .Kp = (NX == 1) ? wc : (wc * wc),
        .Kd = (NX == 1) ? T{0} : (T{2} * wc),
        .Phi = e.Phi,
        .Gamma = e.Gamma,
        .L = e.L,
        .success = true,
        .status = DesignStatus::ok,
        .Ts = e.Ts,
    };
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Active Disturbance Rejection Control (ADRC)
 *
 * Discrete ADRC with an extended-state observer (ESO) for SISO systems.
 * Supports both 1st-order and 2nd-order plant models.
 * Based on the work of Jingqing Han and Zhiqiang Gao. The observer is @ref ESO
 * (design::eso): one mat-vec per tick; this runtime adds the PD law that cancels f̂.
 *
 * References:
 * - Han, J. (1998). From PID to Active Disturbance Rejection Control. IEEE Transactions on Industrial Electronics, 45(5), 900-906.
 * - Gao, Z. (2006). Active Disturbance Rejection Control: A Paradigm Shift in Feedback Control System Design. American Control Conference.
 *
 * @see "Discrete implementation and generalization of the extended state observer"
 *      (Miklosovic, Radke & Gao, 2006), ACC
 * @tparam NX Plant order: 1 (1st-order) or 2 (2nd-order) only
 * @tparam T Scalar type (default: float)
 */
template<size_t NX, typename T = float, typename EsoStore = detail::HeldEso<NX, T>>
class ADRCController {
    static_assert(NX == 1 || NX == 2, "ADRC supports 1st- and 2nd-order plants only");

    template<size_t, typename, typename>
    friend class ADRCController; // cross-precision converting ctor

    T                    b0{1.0f}; ///< input gain
    T                    Kp{1.0f}; ///< proportional gain
    T                    Kd{1.0f}; ///< derivative gain
    ESO<NX, T, EsoStore> eso_{};   ///< Discrete ESO (Miklosovic–Radke–Gao current estimator)

public:
    constexpr ADRCController() = default;

    /** @brief The gains, its ESO's maps being constants of the type (Law of a `constexpr` design). */
    constexpr ADRCController(T b0_, T Kp_, T Kd_)
        requires(EsoStore::constant)
        : b0(b0_), Kp(Kp_), Kd(Kd_) {}

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T> && !EsoStore::constant)
    constexpr ADRCController(const design::ADRCResult<NX, U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : ADRCController(result.template as<T>()) {}

    /**
     * @brief From the design result, which carries the sample time.
     *
     * A failed or continuous (Ts = 0) design is ill-formed under constant evaluation; at run
     * time it leaves the controller invalid (valid() false, output 0).
     */
    constexpr ADRCController(const design::ADRCResult<NX, T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        requires(!EsoStore::constant)
        : b0(result.b0),
          Kp(result.Kp),
          Kd(result.Kd),
          eso_(result.as_eso()) {
        if ((!result.success || !(result.Ts > T{0})) && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_sample_time : result.status); // constant-evaluated (constinit / constexpr): a failed or continuous design does not compile
        }
    }

    template<typename U>
        requires(!EsoStore::constant)
    constexpr ADRCController(const ADRCController<NX, U>& other)
        : b0(static_cast<T>(other.b0)),
          Kp(static_cast<T>(other.Kp)),
          Kd(static_cast<T>(other.Kd)),
          eso_(other.eso_) {}

    /**
     * @brief Reference-tracking control step (Gao's linear ADRC).
     *
     * Runs the discrete extended-state observer one step (predict with the previously-applied
     * command, correct with y_k), then forms the control law
     *
     *     u = (Kp·(r − ẑ₁) − Kd·ẑ₂ − f̂) / b0
     *
     * where `f̂ = z_{NX+1}` is the ESO's total-disturbance estimate. Cancelling
     * `f̂` is what gives ADRC its disturbance rejection and (integral-free) zero
     * steady-state error. The two design knobs are the controller and observer
     * bandwidths `wc`/`wo` (see @ref design::adrc); there is no separate
     * integrator gain.
     *
     * The ESO uses the previously-applied command. After computing a new command this
     * overload optimistically stores the unsaturated value; if a downstream stage
     * clamps the command, follow up with `back_calculate(u_sat)` so the next ESO tick
     * uses what was actually applied.
     *
     * Failed or continuous designs (success = false, b0 ≤ 0 or Ts = 0) emit zero and hold the
     * ESO state, so the control law never divides by an invalid b0 or runs without the observer.
     *
     * @param r  Reference (setpoint).
     * @param y  Measurement (plant output).
     * @return Control command `u`.
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        if (!valid()) {
            return T{0};
        }

        const auto& z = eso_.step(y); // current estimator on the last applied command

        // --- Control law: PD on the estimated state, minus the disturbance.
        const T u0 = (Kp * (r - z[0])) - (Kd * z[1]); // Kd == 0 for NX == 1
        const T u = (u0 - z[NX]) / b0;                // cancel f̂ = z_{NX+1}

        eso_.back_calculate(u);
        return u;
    }

    /**
     * @brief Anti-windup hook for cascade-level saturation propagation.
     *
     * Overrides the stored command with the actually-realized value so
     * the next ESO tick reflects what the plant received, not what the
     * controller wanted. This is the natural ADRC anti-windup behavior:
     * the ESO's disturbance estimate self-corrects once it sees the true
     * input, without an explicit integrator unwind step.
     *
     * Same signature as the other laws' back_calculate. The ESO needs only the applied
     * input, and a feedforward added downstream is already part of @p u_sat, so @p u_ff
     * does not enter.
     *
     * @param u_sat Command actually applied after downstream clamping.
     * @param u_ff  Feedforward added to this law's command before the clamp (unused here).
     */
    constexpr void back_calculate(T u_sat, [[maybe_unused]] T u_ff = T{0}) { eso_.back_calculate(u_sat, u_ff); }

    /** @brief A successful discrete (Ts > 0) design with b0 > 0; gates control(). */
    [[nodiscard]] constexpr bool valid() const { return eso_.valid() && b0 > T{0}; }

    /** @brief ESO state [ŷ, (ŷ̇,) f̂] — plot these when commissioning. */
    [[nodiscard]] constexpr const ColVec<NX + 1, T>& observer() const { return eso_.state(); }

    /** @brief Command the ESO used last tick (post-`back_calculate` if you called it). */
    [[nodiscard]] constexpr T last_command() const { return eso_.last_command(); }

    /** @brief Zero the ESO state and the held input. */
    constexpr void reset() { eso_.reset(); }

    /** @brief Start the ESO from @p z0 with no held input: e.g. [y₀, 0, 0] avoids the first-command kick from ẑ = 0. */
    constexpr void reset(const ColVec<NX + 1, T>& z0) { eso_.reset(z0); }
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<size_t NX, typename T>
ADRCController(const design::ADRCResult<NX, T>&) -> ADRCController<NX, T>;

namespace detail {

/** @brief The runtime an ADRC design runs as: none for anything else. */
template<const auto& Res, typename T, typename R = damp::remove_cvref_t<decltype(Res)>>
struct adrc_runtime {};

/** @brief design::adrc: with members, or its ESO on the design's constants. */
template<const auto& Res, typename T, size_t NX, typename U>
struct adrc_runtime<Res, T, design::ADRCResult<NX, U>> {
    using type = ADRCController<NX, T>;                              ///< With members
    using on_constants = ADRCController<NX, T, ConstantEso<Res, T>>; ///< Its ESO on the design's constants
};

/** @brief An ADRC design. */
template<const auto& Res, typename T>
concept adrc_design = requires { typename adrc_runtime<Res, T>::type; };

} // namespace detail

/** @brief Law<kRes, T> of a design::adrc initialized at start-up: the ADRCController it runs as. A runtime of the same shape assigned replaces it. */
template<const auto& Res, typename T>
    requires detail::adrc_design<Res, T>
struct Law<Res, T> : detail::adrc_runtime<Res, T>::type {
    using Runtime = detail::adrc_runtime<Res, T>::type; ///< The runtime the design runs as

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res} {}

    using Runtime::operator=;
};

/** @brief Law<kRes, T> of a `constexpr` design::adrc: the same law, its ESO on the design's matrices as constants of the type. */
template<const auto& Res, typename T>
    requires(detail::adrc_design<Res, T> && detail::constant_design<Res>)
struct Law<Res, T> : detail::adrc_runtime<Res, T>::on_constants {
    using Runtime = detail::adrc_runtime<Res, T>::type; ///< The runtime with members: the type of a member a redesign assigns

    static_assert(detail::constant_design_runs(Res.success, Res.status, Res.Ts));

    /** @brief The design's law. */
    constexpr Law() : detail::adrc_runtime<Res, T>::on_constants(static_cast<T>(Res.b0), static_cast<T>(Res.Kp), static_cast<T>(Res.Kd)) {}
};
} // namespace damp
