// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file stsmc.hpp
 * @brief Super-twisting (second-order sliding-mode) controller — continuous
 *        finite-time rejection of a Lipschitz matched disturbance.
 */

#include "damp/math/math.hpp"
namespace damp {

namespace design {

/**
 * @struct STSMCResult
 * @brief Super-twisting (second-order sliding-mode) controller design result.
 *
 * Holds the two super-twisting gains plus the optional shaping parameters. Use
 * .as\<U\>() to convert for embedded deployment.
 *
 * @see "Strict Lyapunov Functions for the Super-Twisting Algorithm" (Moreno &
 *      Osorio, 2012), IEEE TAC
 */
template<typename T = double>
struct STSMCResult {
    T    k1{};           ///< Continuous gain on |s|^½·sign(s)
    T    k2{};           ///< Integral gain on sign(s)
    T    lambda{};       ///< Sliding-surface slope for the (r, y) convenience form (s = λe + ė)
    T    k_lin{};        ///< Generalized-STA linear gain (0 = classic super-twisting)
    T    epsilon{};      ///< Boundary-layer thickness for the sign function (0 = true sign)
    bool success{false}; ///< true if gains are valid (k1, k2 > 0)

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return STSMCResult<U>{
            static_cast<U>(k1), static_cast<U>(k2), static_cast<U>(lambda), static_cast<U>(k_lin),
            static_cast<U>(epsilon), success
        };
    }
};

/**
 * @brief Synthesize super-twisting gains from a disturbance-derivative bound.
 *
 * The super-twisting algorithm rejects a matched disturbance d(t) acting on the
 * sliding dynamics ṡ = u + d, provided the disturbance is Lipschitz with a known
 * bound on its rate, |ḋ| ≤ L. A standard sufficient gain choice (Levant /
 * Moreno-Osorio) is
 *
 *     k₁ = 1.5·√L,   k₂ = 1.1·L
 *
 * scaled by an optional @p gain_margin ≥ 1 for robustness to model error. The
 * unique feature versus first-order SMC: the control is *continuous* (the
 * discontinuity is hidden under the integrator), so it rejects the disturbance
 * with finite-time convergence without the boundary-layer chattering of
 * `SMCController`.
 *
 * For a noisy measurement or a higher-order/under-modelled actuator (e.g. a
 * hydraulic swash-plate stage), set @p k_lin > 0 to use the generalized
 * super-twisting algorithm, which adds linear damping terms (`+ k_lin·s`) that
 * stiffen convergence and de-sensitize the loop to the `|s|^½` term amplifying
 * sensor noise — the usual cure when textbook super-twisting still chatters on
 * real hardware.
 *
 * @note Discretized with explicit Euler on the integral state. Keep gains matched
 *       to the sample rate; if `s` is differentiated from a noisy signal, filter
 *       it (or supply `s` directly) — noise on `s` is the dominant real-world
 *       chattering source.
 *
 * @tparam T Scalar type
 * @param disturbance_bound L, the bound on |ḋ| (must be > 0)
 * @param lambda            Sliding-surface slope for the (r, y) form (0 if you supply s directly)
 * @param k_lin             Generalized-STA linear gain (≥ 0; 0 = classic super-twisting)
 * @param epsilon           Boundary-layer thickness for the sign function (≥ 0; 0 = true sign)
 * @param gain_margin       Robustness multiplier on both gains (≥ 1)
 * @return STSMCResult with computed k₁, k₂ and the shaping parameters
 *
 * @see "Sliding order and sliding accuracy in sliding mode control" (Levant, 2003)
 * @see smc() for the first-order (boundary-layer) sliding-mode controller
 */
template<typename T = double>
[[nodiscard]] constexpr STSMCResult<T> stsmc(
    T disturbance_bound,
    T lambda = T{0},
    T k_lin = T{0},
    T epsilon = T{0},
    T gain_margin = T{1}
) {
    if (disturbance_bound <= T{0} || gain_margin < T{1} || k_lin < T{0} || epsilon < T{0}) {
        return STSMCResult<T>{};
    }
    const T L = disturbance_bound;
    return STSMCResult<T>{
        .k1 = gain_margin * static_cast<T>(1.5) * damp::sqrt(L),
        .k2 = gain_margin * static_cast<T>(1.1) * L,
        .lambda = lambda,
        .k_lin = k_lin,
        .epsilon = epsilon,
        .success = true,
    };
}

/**
 * @brief Super-twisting controller from gains you specify directly.
 *
 * Bypasses the disturbance-bound formula for users who tune k₁/k₂ themselves.
 * Named @c stsmc_gains (not an overload of @ref stsmc) so a two-arg call is not
 * ambiguous with @ref stsmc's `(L, λ)` form under default arguments.
 *
 * @param k1 Continuous gain on |s|^½·sign(s) (> 0)
 * @param k2 Integral gain on sign(s) (> 0)
 * @param lambda  Sliding-surface slope for the (r, y) form
 * @param k_lin      Generalized-STA linear gain (≥ 0)
 * @param epsilon Boundary-layer thickness (≥ 0)
 */
template<typename T = double>
[[nodiscard]] constexpr STSMCResult<T>
stsmc_gains(T k1, T k2, T lambda = T{0}, T k_lin = T{0}, T epsilon = T{0}) {
    if (!(k1 > T{0} && k2 > T{0} && k_lin >= T{0} && epsilon >= T{0})) {
        return STSMCResult<T>{};
    }
    return STSMCResult<T>{
        .k1 = k1,
        .k2 = k2,
        .lambda = lambda,
        .k_lin = k_lin,
        .epsilon = epsilon,
        .success = true,
    };
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Super-twisting controller (second-order sliding mode).
 *
 * Continuous-control sliding mode for systems of relative degree 1 in the
 * sliding variable (ṡ = u + d). The control law is
 *
 *     u   = −k₁·(|s|^½·sign(s) + k_lin·s) + v
 *     v̇   = −k₂·(sign(s) + 3·k_lin·|s|^½·sign(s) + 2·k_lin²·s)
 *
 * With k_lin = 0 this is the classic super-twisting algorithm (Levant):
 * `u = −k₁·|s|^½·sign(s) + v`, `v̇ = −k₂·sign(s)`. The integral state v is what
 * makes u continuous — finite-time rejection of a Lipschitz matched disturbance
 * with none of the first-order boundary-layer chattering. With k_lin > 0 the
 * generalized form adds linear damping for noisy / higher-order plants.
 *
 * Two entry points per step:
 *  - `control(s, Ts)`    — canonical: you form the sliding variable s yourself.
 *  - `control(r, y, Ts)` — convenience: builds s = λ·e + ė (ė by backward
 *                          difference), matching `SMCController`'s surface.
 *
 * The integral state is updated by explicit Euler. The output is the
 * sliding-mode control contribution; add any equivalent/feedforward control
 * separately.
 *
 * Example: reject a sinusoidal load disturbance on an integrator plant
 * @code
 * using namespace damp;
 * // |ḋ| ≤ 2, 1 kHz loop; supply s directly.
 * constexpr auto art = design::stsmc(2.0);
 * static_assert(art.success);
 * STSMCController u_st(art.as<float>());
 * // In the loop: float u = u_st.control_sliding(s, 1e-3f);  // s = your sliding variable
 * @endcode
 *
 * @tparam T Scalar type (default: float)
 * @see "Strict Lyapunov Functions for the Super-Twisting Algorithm" (Moreno &
 *      Osorio, 2012), IEEE TAC, https://doi.org/10.1109/TAC.2012.2186179
 */
template<typename T = float>
class STSMCController {
public:
    constexpr STSMCController() = default;

    constexpr STSMCController(const design::STSMCResult<T>& result, T Ts = T{0})
        : k1_(result.k1), k2_(result.k2), lambda_(result.lambda), k_lin_(result.k_lin), epsilon_(result.epsilon), Ts_(Ts), valid_(result.success) {}

    template<typename U>
    constexpr STSMCController(const STSMCController<U>& other)
        : k1_(static_cast<T>(other.k1_)),
          k2_(static_cast<T>(other.k2_)),
          lambda_(static_cast<T>(other.lambda_)),
          k_lin_(static_cast<T>(other.k_lin_)),
          epsilon_(static_cast<T>(other.epsilon_)),
          Ts_(static_cast<T>(other.Ts_)),
          v_(static_cast<T>(other.v_)),
          e_prev_(static_cast<T>(other.e_prev_)),
          first_(other.first_),
          valid_(other.valid_) {}

    /// Super-twisting control from the sliding variable @p s (canonical form).
    /// Named distinctly from control(r, y) because both would otherwise be
    /// two-scalar overloads.
    /// @param s  Sliding variable
    /// @param Ts Sample time [s] for the explicit-Euler integral state
    [[nodiscard]] constexpr T control_sliding(T s, T Ts) {
        if (!valid_) {
            return T{0};
        }
        const T abs_s = damp::abs(s);
        // sign(s), optionally softened by a boundary layer to trade a little
        // accuracy for less numerical chatter near s = 0.
        const T sign_s = (epsilon_ > T{0}) ? (s / (abs_s + epsilon_)) : static_cast<T>(damp::sgn(s));
        const T sqrt_s = damp::sqrt(abs_s);

        const T phi1 = (sqrt_s * sign_s) + (k_lin_ * s);
        const T u = (-k1_ * phi1) + v_;

        // φ₂ = sign(s) + 3·k_lin·|s|^½·sign(s) + 2·k_lin²·s  (= classic sign(s) when k_lin = 0,
        // keeping the k₁ = 1.5√L / k₂ = 1.1L formulas valid). Explicit-Euler integral.
        const T phi2 = sign_s + (T{3} * k_lin_ * sqrt_s * sign_s) + (T{2} * k_lin_ * k_lin_ * s);
        if (Ts > T{0}) {
            v_ += -k2_ * phi2 * Ts; // explicit-Euler integral; skip on a stalled/non-positive dt
        }

        return u;
    }

    /// Sliding-variable form using the stored sample time (set at construction).
    [[nodiscard]] constexpr T control_sliding(T s) { return control_sliding(s, Ts_); }

    /// Convenience: build s = λ·e + ė (ė via backward difference) then apply STA.
    /// @note Assumes the control enters the surface dynamics with the canonical
    ///       sign (ṡ = u + d). If your plant gives ṡ = −u + d, negate the output
    ///       (or your λ/s), or use the `control_sliding(s, Ts)` form with a surface
    ///       you sign yourself.
    [[nodiscard]] constexpr T control(T r, T y, T Ts) {
        // Guard the backward-difference divide; a non-positive Ts holds and emits nothing.
        if (!valid_ || Ts <= T{0}) {
            return T{0};
        }
        const T e = r - y;
        if (first_) {
            e_prev_ = e; // seed: no first-tick derivative kick into the surface
            first_ = false;
        }
        const T de = (e - e_prev_) / Ts;
        e_prev_ = e;
        return control_sliding((lambda_ * e) + de, Ts);
    }

    /**
     * @brief 2-arg fixed-rate form using the stored sample time.
     *
     * The SisoController concept shape: construct with Ts and call
     * control(r, y) each tick. With no stored Ts it emits zero and holds
     * state, exactly like the 3-arg form with Ts = 0.
     */
    [[nodiscard]] constexpr T control(T r, T y) { return control(r, y, Ts_); }

    /// Current value of the integral state v (the disturbance estimate, −d̂).
    [[nodiscard]] constexpr T integral_state() const { return v_; }

    [[nodiscard]] constexpr bool valid() const { return valid_; }

    constexpr void reset() {
        v_ = T{0};
        e_prev_ = T{0};
        first_ = true;
    }

private:
    template<typename>
    friend class STSMCController; // cross-precision converting ctor

    T    k1_{T{0}};
    T    k2_{T{0}};
    T    lambda_{T{0}};
    T    k_lin_{T{0}};
    T    epsilon_{T{0}};
    T    Ts_{T{0}}; ///< Stored sample time for the 2-arg control(r, y) / 1-arg control_sliding(s) forms [s]
    T    v_{T{0}};
    T    e_prev_{T{0}};
    bool first_{true};
    bool valid_{false};
};

} // namespace damp
