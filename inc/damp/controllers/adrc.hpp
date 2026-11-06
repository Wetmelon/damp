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
 * f̂. The runtime's ESO is the discrete one of Miklosovic, Radke and Gao: the ESO's linear
 * dynamics discretized by ZOH, in current-estimator form, with its poles at the image
 * @f$e^{-\omega_o T_s}@f$ of the continuous design's poles at @f$-\omega_o@f$.
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
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
 * ẑ_k = F ẑ_{k−1} + G_u u_{k−1} + G_y y_k.
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
    T wo{}; ///< Luenberger bandwidth
    T b0{}; ///< Plant gain

    ColVec<NX + 1, T> beta{}; ///< ESO gains

    T Kp{}; ///< Unscaled proportional gain (wc or wc²); tick divides by b0
    T Kd{}; ///< Unscaled derivative gain (0 or 2·wc); tick divides by b0

    Matrix<NX + 1, NX + 1, T> F{};  ///< Discrete ESO state map (I − L C)Φ; eigenvalues e^{−ω_o T_s} (zero when Ts = 0)
    ColVec<NX + 1, T>         Gu{}; ///< Discrete ESO input map (I − L C)Γ on u_{k−1}
    ColVec<NX + 1, T>         Gy{}; ///< Discrete ESO measurement gain L on y_k (current estimator)

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
            F.template as<U>(),
            ColVec<NX + 1, U>(Gu.template as<U>()),
            ColVec<NX + 1, U>(Gy.template as<U>()),
            success,
            status,
            static_cast<U>(Ts),
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

namespace detail {

/**
 * @brief Fill the discrete ESO of @p r (F, Gu, Gy) from its b0, wo and Ts > 0: the ZOH model of
 *        the chain of integrators and the current-estimator gain with every pole at e^{−ω_o T_s}.
 * @see "Discrete implementation and generalization of the extended state observer"
 *      (Miklosovic, Radke & Gao, 2006), ACC
 */
template<size_t NX, typename T>
constexpr void discrete_eso(ADRCResult<NX, T>& r) {
    constexpr size_t       n = NX + 1;
    StateSpace<n, 1, 1, T> eso{};
    for (size_t i = 0; i < NX; ++i) {
        eso.A(i, i + 1) = T{1}; // ż_i = ẑ_{i+1}; f̂ enters the highest derivative
    }
    eso.B(NX - 1, 0) = r.b0;
    eso.C(0, 0) = T{1};
    const auto d_opt = ::damp::detail::discretize_zoh_impl(eso, r.Ts); // exact: A is nilpotent
    if (!d_opt) {
        ::damp::detail::set_status(r, DesignStatus::invalid_parameter); // b0 and Ts are checked finite first
        return;
    }
    const auto& d = *d_opt;

    const T      h = r.Ts;
    const T      lambda = damp::exp(-r.wo * h);
    const T      m = T{1} - lambda;
    ColVec<n, T> L{};
    if constexpr (NX == 1) {
        L(0) = m * (T{1} + lambda);
        L(1) = (m * m) / h;
    } else {
        L(0) = m * (T{1} + lambda + (lambda * lambda));
        L(1) = (T{3} * m * m * (T{1} + lambda)) / (T{2} * h);
        L(2) = (m * m * m) / (h * h);
    }

    // (I − L C) M = M − L (C M): C picks row 0
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            r.F(i, j) = d.A(i, j) - (L(i) * d.A(0, j));
        }
        r.Gu(i) = d.B(i, 0) - (L(i) * d.B(0, 0));
    }
    r.Gy = L;
}

} // namespace detail

/**
 * @brief Active Disturbance Rejection Control design
 *
 * Designs ESO gains using pole placement for observer poles at -wo.
 * Places all poles at -wo for the extended system.
 * Controller gains are unscaled (b0 is applied only in the runtime law
 * u = (u0 − f̂)/b0):
 * - 1st order: Kp = wc, Kd = 0
 * - 2nd order: Kp = wc², Kd = 2·wc
 *
 * Only plant orders NX = 1 and NX = 2 are supported (same surface as typical
 * linear ADRC toolboxes / MATLAB® discrete ADRC). Higher-order ADRC is out of
 * scope.
 *
 * For Ts > 0 the design also builds the discrete ESO the runtime ticks (Miklosovic, Radke &
 * Gao, 2006). The ESO's linear dynamics @f$\dot z = A z + B u + \beta (y - C z)@f$, with
 * @f$A@f$ the chain of integrators, are discretized by ZOH (@ref discretize), exact for this
 * model, to @f$(\Phi, \Gamma)@f$. The current estimator
 * @f[
 *   \bar z_k = \Phi \hat z_{k-1} + \Gamma u_{k-1},\qquad
 *   \hat z_k = \bar z_k + L\,(y_k - C \bar z_k)
 * @f]
 * has error dynamics @f$(I - L C)\Phi@f$, and @f$L@f$ places all its poles at
 * @f$\lambda = e^{-\omega_o T_s}@f$, the image of the continuous poles at @f$-\omega_o@f$:
 * @f[
 *   n = 2:\ L = \begin{bmatrix} 1-\lambda^2 \\ (1-\lambda)^2/T_s \end{bmatrix},\qquad
 *   n = 3:\ L = \begin{bmatrix} 1-\lambda^3 \\ \tfrac{3}{2}(1-\lambda)^2(1+\lambda)/T_s \\
 *              (1-\lambda)^3/T_s^2 \end{bmatrix}.
 * @f]
 * The result carries @f$F = (I - LC)\Phi@f$, @f$G_u = (I - LC)\Gamma@f$ and @f$G_y = L@f$, so
 * the tick is @f$\hat z_k = F \hat z_{k-1} + G_u u_{k-1} + G_y y_k@f$, mat-vec only. Unlike
 * a forward-Euler step of the continuous ESO, whose poles 1 − ω_o T_s leave the unit circle
 * for ω_o T_s > 2, these poles stay at e^{−ω_o T_s} for any ω_o T_s.
 * Ts = 0 is the continuous design, for analysis (to_ss, to_tf); F, G_u and G_y stay zero.
 *
 * Requires wc > 0, wo > 0, b0 > 0 (all finite), else refused with invalid_parameter and
 * zero gains so the runtime cannot divide by an invalid b0; a negative or non-finite Ts
 * is refused with invalid_sample_time.
 *
 * @param wc Controller bandwidth
 * @param wo Luenberger bandwidth
 * @param b0 Plant gain
 * @param Ts [s] Sample time of the runtime (0 = continuous)
 *
 * @return ADRCResult with computed gains
 *
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
    if (!(Ts >= T{0}) || !damp::isfinite(Ts)) {
        return ::damp::detail::refused<ADRCResult<NX, T>>(DesignStatus::invalid_sample_time);
    }

    ColVec<NX + 1, T> beta{};
    constexpr size_t  n = NX + 1;
    for (size_t i = 1; i <= n; ++i) {
        T binom = T{1};
        for (size_t k = 1; k <= i; ++k) {
            binom *= static_cast<T>(n - k + 1);
            binom /= static_cast<T>(k);
        }
        beta[i - 1] = binom * damp::pow(wo, static_cast<int>(i));
    }

    ADRCResult<NX, T> out{
        .wc = wc,
        .wo = wo,
        .b0 = b0,
        .beta = beta,
        .Kp = (NX == 1) ? wc : (wc * wc),
        .Kd = (NX == 1) ? T{0} : (T{2} * wc),
        .success = true,
        .status = DesignStatus::ok,
        .Ts = Ts,
    };
    if (Ts > T{0}) {
        detail::discrete_eso(out);
    }
    return out;
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Active Disturbance Rejection Control (ADRC)
 *
 * Discrete ADRC with an extended-state observer (ESO) for SISO systems.
 * Supports both 1st-order and 2nd-order plant models.
 * Based on the work of Jingqing Han and Zhiqiang Gao. The ESO is the ZOH-discretized current
 * estimator the design builds (design::adrc): one mat-vec per tick.
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
template<size_t NX, typename T = float>
class ADRCController {
    static_assert(NX == 1 || NX == 2, "ADRC supports 1st- and 2nd-order plants only");

    template<size_t, typename>
    friend class ADRCController; // cross-precision converting ctor

    T b0{1.0f}; ///< input gain
    T Kp{1.0f}; ///< proportional gain
    T Kd{1.0f}; ///< derivative gain

    Matrix<NX + 1, NX + 1, T> F{};           ///< Discrete ESO state map
    ColVec<NX + 1, T>         Gu{};          ///< Discrete ESO input map on u_{k−1}
    ColVec<NX + 1, T>         Gy{};          ///< Discrete ESO measurement gain on y_k
    ColVec<NX + 1, T>         z{};           ///< ESO state: [ŷ, ŷ̇, ..., f̂]
    T                         u_prev_{T{0}}; ///< Last applied command (post-saturation)
    bool                      valid_{false}; ///< A successful discrete (Ts > 0) design with b0 > 0; gates control()

public:
    constexpr ADRCController() = default;
    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!std::is_same_v<U, T>)
    constexpr ADRCController(const design::ADRCResult<NX, U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : ADRCController(result.template as<T>()) {}

    /**
     * @brief From the design result, which carries the sample time.
     *
     * A failed or continuous (Ts = 0) design is ill-formed under constant evaluation; at run
     * time it leaves the controller invalid (valid() false, output 0).
     */
    constexpr ADRCController(const design::ADRCResult<NX, T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : b0(result.b0),
          Kp(result.Kp),
          Kd(result.Kd),
          F(result.F),
          Gu(result.Gu),
          Gy(result.Gy),
          valid_(result.success && result.b0 > T{0} && result.Ts > T{0}) {
        if ((!result.success || !(result.Ts > T{0})) && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_sample_time : result.status); // constant-evaluated (constinit / constexpr): a failed or continuous design does not compile
        }
    }

    template<typename U>
    constexpr ADRCController(const ADRCController<NX, U>& other)
        : b0(static_cast<T>(other.b0)),
          Kp(static_cast<T>(other.Kp)),
          Kd(static_cast<T>(other.Kd)),
          F(other.F.template as<T>()),
          Gu(other.Gu.template as<T>()),
          Gy(other.Gy.template as<T>()),
          z(other.z.template as<T>()),
          u_prev_(static_cast<T>(other.u_prev_)),
          valid_(other.valid_) {}

    /**
     * @brief Reference-tracking control step (Gao's linear ADRC).
     *
     * Runs the discrete extended-state observer one step, ẑ_k = F ẑ_{k−1} + G_u u_{k−1} +
     * G_y y_k (predict with the previously-applied command, correct with y_k), then forms the
     * control law
     *
     *     u = (Kp·(r − ẑ₁) − Kd·ẑ₂ − f̂) / b0
     *
     * where `f̂ = z_{NX+1}` is the ESO's total-disturbance estimate. Cancelling
     * `f̂` is what gives ADRC its disturbance rejection and (integral-free) zero
     * steady-state error. The two design knobs are the controller and observer
     * bandwidths `wc`/`wo` (see @ref design::adrc); there is no separate
     * integrator gain.
     *
     * The ESO uses the previously-applied command (stored in `u_prev_`). After
     * computing a new command this overload optimistically stores the
     * unsaturated value in `u_prev_`; if a downstream stage clamps the command,
     * follow up with `back_calculate(u_sat)` so the next ESO tick uses
     * what was actually applied.
     *
     * Failed or continuous designs (success = false, b0 ≤ 0 or Ts = 0) emit zero and hold the
     * ESO state, so the control law never divides by an invalid b0 or runs without the observer.
     *
     * @param r  Reference (setpoint).
     * @param y  Measurement (plant output).
     * @return Control command `u`.
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        if (!valid_) {
            return T{0};
        }

        // --- Discrete ESO (ZOH, current estimator) on the last applied command.
        // Plant model: y^(NX) = f + b0·u, with z = [ŷ, ŷ̇, …, ŷ^(NX-1), f̂].
        z = ColVec<NX + 1, T>((F * z) + (Gu * u_prev_) + (Gy * y));

        // --- Control law: PD on the estimated state, minus the disturbance.
        const T u0 = (Kp * (r - z[0])) - (Kd * z[1]); // Kd == 0 for NX == 1
        const T u = (u0 - z[NX]) / b0;                // cancel f̂ = z_{NX+1}

        u_prev_ = u;
        return u;
    }

    /**
     * @brief Anti-windup hook for cascade-level saturation propagation.
     *
     * Overrides the stored `u_prev_` with the actually-realized command so
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
    constexpr void back_calculate(T u_sat, [[maybe_unused]] T u_ff = T{0}) { u_prev_ = u_sat; }

    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /** @brief ESO state [ŷ, (ŷ̇,) f̂] — plot these when commissioning. */
    [[nodiscard]] constexpr const ColVec<NX + 1, T>& observer() const { return z; }

    /** @brief Command the ESO used last tick (post-`back_calculate` if you called it). */
    [[nodiscard]] constexpr T last_command() const { return u_prev_; }

    constexpr void reset() {
        z = ColVec<NX + 1, T>{};
        u_prev_ = T{0};
    }
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<size_t NX, typename T>
ADRCController(const design::ADRCResult<NX, T>&) -> ADRCController<NX, T>;
} // namespace damp
