// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqi.hpp
 * @brief Linear-Quadratic-Integral design and runtime controller
 *
 * @code
 * using namespace damp;
 * // Discrete plant with Ts already set on sys
 * constexpr auto res = design::discrete_lqi(
 *     sys,
 *     Matrix<NX + NY, NX + NY>::identity(), // Q on [x; xi]
 *     Matrix<NU, NU>::identity());          // R
 * static_assert(res.success);
 * LQI ctrl(res);
 * auto u = ctrl.control(r, y, x); // advances internal xi
 *
 * // Integrate only sensors 0 and 3 of six: Q is on [x; xi] with two integrators.
 * constexpr auto res2 = design::discrete_lqi(sys6, Q_aug2, R, {0, 3});
 * LQI ctrl2(res2);
 * const auto u_req = ctrl2.feedback(r2, y6, x); // r2: references for sensors 0 and 3
 * const auto u_app = clamp_each(u_req);          // the actuator's limits
 * ctrl2.commit(u_app);                           // integrators hold when the clamp would wind them up
 * @endcode
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/stability.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace detail {

/// NI×NY selector of the first NI outputs (the identity when NI = NY).
template<size_t NI, size_t NY, typename T>
[[nodiscard]] constexpr Matrix<NI, NY, T> leading_output_selector() {
    Matrix<NI, NY, T> Ci{};
    for (size_t i = 0; i < NI && i < NY; ++i) {
        Ci(i, i) = T{1};
    }
    return Ci;
}

/// NI×NY selector picking the listed outputs; nullopt when an index is outside [0, NY).
template<size_t NI, size_t NY, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NI, NY, T>> output_selector(const size_t (&outputs)[NI]) {
    Matrix<NI, NY, T> Ci{};
    for (size_t i = 0; i < NI; ++i) {
        if (outputs[i] >= NY) {
            return damp::nullopt;
        }
        Ci(i, outputs[i]) = T{1};
    }
    return Ci;
}

} // namespace detail

namespace design {

/**
 * @struct LQIResult
 * @brief LQI design result
 *
 * Optimal gain on the augmented state [x; xi] for discrete LQI tracking,
 * with Riccati solution S and closed-loop poles of the augmented system.
 * Use .as<float>() for deploy.
 *
 * @see discrete_lqi()
 * @see "Optimal Control" (Anderson & Moore, 1990), §4 (LQR on the augmented plant)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NI = NY>
struct LQIResult {
    Matrix<NU, NX + NI, T>            K{};                                                    ///< Optimal gain: u = -K*[x; xi]
    Matrix<NX + NI, NX + NI, T>       S{};                                                    ///< Riccati equation solution
    ColVec<NX + NI, damp::complex<T>> e{};                                                    ///< Closed-loop poles
    bool                              success{false};                                         ///< true if DARE converged
    Matrix<NI, NY, T>                 Ci{damp::detail::leading_output_selector<NI, NY, T>()}; ///< Integrated outputs: xi += r − Ci·y

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return LQIResult<NX, NU, NY, U, NI>{
            K.template as<U>(),
            S.template as<U>(),
            e.template as<damp::complex<U>>(),
            success,
            Ci.template as<U>(),
        };
    }
};

/**
 * @brief Discrete Linear-Quadratic-Integral (LQI) design for output tracking
 *
 * Augments the discrete plant with an output-error integrator and solves DARE
 * LQR on the augmented pair (A_aug, B_aug):
 *
 *     A_aug = |  A      0 |      B_aug = | B |
 *             | −Ci·C   I |              | 0 |
 *
 *     xi[k+1] = xi[k] + (r[k] − Ci·y[k])
 *     u[k]    = −K [x[k]; xi[k]]
 *
 * Ci (NI×NY) picks the integrated outputs; it is the identity in the overload
 * without it (every output integrated). Q is the cost on [x; xi] (size
 * (NX+NI)×(NX+NI)); R is the input cost. On a model with delay states
 * (discretize with InputDelay / OutputDelay) Q may instead be sized to the
 * plant, (NX_plant+NI)×(NX_plant+NI): the delay states between the plant block
 * and the integrators get no cost. Constant references and constant
 * output disturbances yield zero steady-state error on the integrated outputs
 * when the design succeeds and the loop is stable. A selection the inputs
 * cannot hold apart (more integrators than inputs, or a rank-deficient Ci) is
 * not stabilizable and reports success = false.
 *
 * @note Compare with MATLAB®'s lqi(sys, Q, R) — exposed as the lqi() alias in matlab.hpp.
 *
 * @see dare() for the Riccati solve
 * @see discrete_lqr() for pure regulation without integral action
 * @see "Optimal Control" (Anderson & Moore, 1990), §4
 *
 * @param sys  Discrete-time state-space plant (uses A, B, C)
 * @param Q    Augmented state cost ((NX+NI)×(NX+NI), or plant-sized on a delayed model; positive semidefinite)
 * @param R    Input cost (NU×NU, positive definite)
 * @param Ci   Integrated-output selector (NI×NY): integrator i follows r_i − (Ci·y)_i
 * @return LQIResult with K, S, closed-loop poles, Ci, and success
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI, size_t NQA = NX + NI>
[[nodiscard]] constexpr LQIResult<NX, NU, NY, T, NI> discrete_lqi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NQA, NQA, T>&               Q,
    const Matrix<NU, NU, T>&                 R,
    const Matrix<NI, NY, T>&                 Ci
) {
    const auto                  Q_aug = damp::detail::pad_integral_weight<NX, NI>(Q); // delay states: no cost
    Matrix<NX + NI, NX + NI, T> A_aug{};
    Matrix<NX + NI, NU, T>      B_aug{};
    A_aug.template block<NX, NX>(0, 0) = sys.A;
    A_aug.template block<NI, NX>(NX, 0) = -(Ci * sys.C);
    A_aug.template block<NI, NI>(NX, NX) = Matrix<NI, NI, T>::identity(); // integrator memory
    B_aug.template block<NX, NU>(0, 0) = sys.B;

    LQIResult<NX, NU, NY, T, NI> out{};
    const auto                   dare_opt = dare(A_aug, B_aug, Q_aug, R);
    if (!dare_opt) {
        return out;
    }
    const auto K_opt = lqr_gain(A_aug, B_aug, dare_opt.value(), R);
    if (!K_opt) {
        return out;
    }
    out.K = K_opt.value();
    out.S = dare_opt.value();
    out.e = stability::closed_loop_poles(A_aug, B_aug, out.K);
    out.Ci = Ci;
    out.success = true;
    return out;
}

/// @copydoc discrete_lqi(const StateSpace<NX, NU, NY, T, NW, NV>&, const Matrix<NQA, NQA, T>&, const Matrix<NU, NU, T>&, const Matrix<NI, NY, T>&)
/// Every output is integrated (Ci = I).
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0, size_t NQA = NX + NY>
[[nodiscard]] constexpr LQIResult<NX, NU, NY, T> discrete_lqi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NQA, NQA, T>&               Q,
    const Matrix<NU, NU, T>&                 R
) {
    return discrete_lqi(sys, Q, R, Matrix<NY, NY, T>::identity());
}

/**
 * @brief Discrete LQI integrating the listed outputs only
 *
 * `discrete_lqi(sys, Q, R, {0, 3})` integrates outputs 0 and 3 of @p sys; the
 * runtime then takes references for those two outputs and the full measurement
 * vector. An index outside [0, NY) reports success = false.
 *
 * @param sys      Discrete-time state-space plant
 * @param Q        Augmented state cost on [x; xi] ((NX+NI)×(NX+NI), or plant-sized on a delayed model)
 * @param R        Input cost (NU×NU)
 * @param outputs  Indices of the integrated outputs, in integrator order
 * @return LQIResult with Ci built from @p outputs
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI, size_t NQA = NX + NI>
[[nodiscard]] constexpr LQIResult<NX, NU, NY, T, NI> discrete_lqi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NQA, NQA, T>&               Q,
    const Matrix<NU, NU, T>&                 R,
    const size_t (&outputs)[NI]
) {
    const auto Ci = damp::detail::output_selector<NI, NY, T>(outputs);
    if (!Ci) {
        return LQIResult<NX, NU, NY, T, NI>{};
    }
    return discrete_lqi(sys, Q, R, Ci.value());
}
} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-Quadratic-Integral (LQI) controller
 *
 * Output tracking controller with integral action: u = -K * [x; xi]
 * where xi integrates the output error (r − Ci·y); Ci selects the integrated
 * outputs (the identity when every output is integrated).
 * Provides zero steady-state error for constant references and disturbances.
 *
 * control(r, y, x) is the self-contained tick. When the actuator can clamp, split
 * it: feedback(r, y, x) returns the raw command, and commit(u_applied) advances
 * the integrators except those that would drive a clamped input further into its
 * limit (conditional integration), so the integral does not wind up.
 *
 * @tparam NX Number of states
 * @tparam NU Number of control inputs
 * @tparam NY Number of outputs
 * @tparam T  Scalar type (default: float — embedded runtime)
 * @tparam NI Number of integrated outputs (default: NY)
 */
template<size_t NX, size_t NU, size_t NY, typename T = float, size_t NI = NY>
struct LQI {
    Matrix<NU, NX + NI, T> K{};                                                    ///< Optimal gain: u = -K*[x; xi]
    ColVec<NI, T>          xi{};                                                   ///< Integral of tracking error, owned by the controller
    Matrix<NI, NY, T>      Ci{damp::detail::leading_output_selector<NI, NY, T>()}; ///< Integrated outputs: xi += r − Ci·y

    constexpr LQI() = default;
    constexpr LQI(const Matrix<NU, NX + NI, T>& K_) : K(K_) {} // NOLINT

    /// From design result (any scalar); converts K and Ci via @c .as\<T\>().
    template<typename U>
    constexpr LQI(const design::LQIResult<NX, NU, NY, U, NI>& result) // NOLINT
        : K(result.K.template as<T>()), Ci(result.Ci.template as<T>()) {}

    template<typename U>
    constexpr LQI(const LQI<NX, NU, NY, U, NI>& other) // NOLINT
        : K(other.K.template as<T>()), xi(other.xi.template as<T>()), Ci(other.Ci.template as<T>()) {}

    /**
     * @brief Compute control from a caller-supplied augmented state
     *
     * Computes u = -K * x_aug where x_aug = [x; xi]. Use this when you maintain
     * the integral state yourself; it does not touch the internal integrator.
     *
     * @param x_aug Augmented state vector [x; xi]
     * @return Control input vector u
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NX + NI, T>& x_aug) const {
        return ColVec<NU, T>(-K * x_aug);
    }

    /**
     * @brief Compute control with the controller's own integral action
     *
     * Uses the internal integrator state: u = -[Kx Ki]·[x; xi] with the current
     * xi, then advances xi[k+1] = xi[k] + (r − Ci·y) (matching the −C augmentation
     * sign convention). Gives zero steady-state error to constant references
     * without the caller tracking the integral.
     *
     * Reference-first argument order (reference, measurement, state), matching
     * the library-wide control(r, feedback...) convention.
     *
     * @param r Reference for the integrated outputs
     * @param y Measured output
     * @param x Current (estimated or measured) plant state
     * @return Control input vector u
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NI, T>& r, const ColVec<NY, T>& y, const ColVec<NX, T>& x) {
        return control(r, y, x, ColVec<NX, T>{});
    }

    /**
     * @brief Compute control against a full state reference (trajectory tracking)
     *
     * u = -Kx·(x − x_ref) − Ki·xi, then xi[k+1] = xi[k] + (r − Ci·y). The state-feedback term is
     * the feedback half of a 2-DOF law (see LQR::control(const ColVec<NX,T>&, const
     * ColVec<NX,T>&) for the u_ff input feedforward it omits). Without @p x_ref the feedback
     * regulates the whole state toward the origin, so on a moving reference it opposes the
     * commanded motion and the integrator alone drags the plant along — a velocity-proportional
     * lag plus the overshoot that windup dumps at arrival. Supplying a feasible reference state
     * (e.g. the profile's position and velocity, consistent with C·x_ref = r) makes the feedback
     * act on trajectory error; the integrator then supplies the quasi-static feedforward, so
     * unlike the pure regulator there is no droop — but a fast-changing u_ff (the acceleration
     * phase of a move) still lags unless fed forward on the input.
     *
     * Reference-first argument order, matching the library-wide control(r, feedback...)
     * convention.
     *
     * @param r     Reference for the integrated outputs (drives the integrator)
     * @param y     Measured output
     * @param x     Current (estimated or measured) plant state
     * @param x_ref State reference the feedback regulates toward (must be feasible)
     * @return Control input vector u
     */
    [[nodiscard]] constexpr ColVec<NU, T>
    control(const ColVec<NI, T>& r, const ColVec<NY, T>& y, const ColVec<NX, T>& x, const ColVec<NX, T>& x_ref) {
        const ColVec<NU, T> u = feedback(r, y, x, x_ref);
        commit(u);
        return u;
    }

    /**
     * @brief Tracking law on a model with delay states: @p x_ref covers the NR plant states only.
     *
     * @p x is still the whole state. The past-input states' reference is the past u*, zero for
     * this law (the integrator supplies the feedforward). Past-output states (OutputDelay) also
     * regulate toward zero here, since this runtime carries no model to rebuild the reference's
     * past outputs; LQGI, which does, fills them in.
     */
    template<size_t NR>
        requires(NR < NX)
    [[nodiscard]] constexpr ColVec<NU, T>
    control(const ColVec<NI, T>& r, const ColVec<NY, T>& y, const ColVec<NX, T>& x, const ColVec<NR, T>& x_ref) {
        return control(r, y, x, detail::zero_extend<NX>(x_ref));
    }

    /**
     * @brief First half of a saturation-aware tick: the raw (pre-saturation) command
     *
     * u = -Kx·(x − x_ref) − Ki·xi at the current integral. The tracking error is held
     * for commit, which decides whether to integrate it.
     * @code
     * const auto u_req = ctrl.feedback(r, y, x);
     * const auto u     = clamp_each(u_req + u_ff); // feedforward, then the actuator's limits
     * ctrl.commit(u, u_ff);
     * @endcode
     *
     * @param r     Reference for the integrated outputs
     * @param y     Measured output
     * @param x     Current (estimated or measured) plant state
     * @param x_ref State reference the feedback regulates toward (default: origin)
     * @return Raw control input u (pre-saturation)
     */
    [[nodiscard]] constexpr ColVec<NU, T> feedback(
        const ColVec<NI, T>& r,
        const ColVec<NY, T>& y,
        const ColVec<NX, T>& x,
        const ColVec<NX, T>& x_ref = ColVec<NX, T>{}
    ) {
        const auto Kx = K.template block<NU, NX>(0, 0);
        const auto Ki = K.template block<NU, NI>(0, NX);
        u_cmd_ = ColVec<NU, T>(-((Kx * (x - x_ref)) + (Ki * xi)));
        e_ = ColVec<NI, T>(r - (Ci * y));
        return u_cmd_;
    }

    /// feedback on a model with delay states: @p x_ref covers the NR plant states only (see control).
    template<size_t NR>
        requires(NR < NX)
    [[nodiscard]] constexpr ColVec<NU, T>
    feedback(const ColVec<NI, T>& r, const ColVec<NY, T>& y, const ColVec<NX, T>& x, const ColVec<NR, T>& x_ref) {
        return feedback(r, y, x, detail::zero_extend<NX>(x_ref));
    }

    /**
     * @brief Second half of the tick: integrate the held error unless the clamp would wind it up
     *
     * The clamp excess is s = u_applied − (u_req + u_ff), with u_req from feedback. Integrator
     * i advances by its error eᵢ unless, on some clamped input j (sⱼ ≠ 0), its next-tick effect
     * −Kᵢ(j, i)·eᵢ points away from the feasible side (−Kᵢ(j, i)·eᵢ·sⱼ < 0) — deeper into the
     * limit. Without a clamp (s = 0) every integrator advances, as in control(r, y, x).
     *
     * @param u_applied Input the plant actually received (after the clamp, including u_ff)
     * @param u_ff      Feedforward added to the command before the clamp (default: none)
     */
    constexpr void commit(const ColVec<NU, T>& u_applied, const ColVec<NU, T>& u_ff = ColVec<NU, T>{}) {
        ColVec<NU, T> excess{};
        for (size_t j = 0; j < NU; ++j) {
            const T u_req = u_cmd_(j) + u_ff(j);
            const T diff = u_applied(j) - u_req;
            const T tol = T{4} * eps<T>() * (damp::abs(u_applied(j)) + damp::abs(u_req));
            excess(j) = (damp::abs(diff) > tol) ? diff : T{0};
        }
        for (size_t i = 0; i < NI; ++i) {
            bool winds_up = false;
            for (size_t j = 0; j < NU; ++j) {
                winds_up = winds_up || (K(j, NX + i) * e_(i) * excess(j) > T{0});
            }
            if (!winds_up) {
                xi(i) += e_(i);
            }
        }
        e_ = ColVec<NI, T>{};
        u_cmd_ = ColVec<NU, T>{};
    }

    /// Clear the integral state and any held tick.
    constexpr void reset() {
        xi = ColVec<NI, T>{};
        e_ = ColVec<NI, T>{};
        u_cmd_ = ColVec<NU, T>{};
    }

    [[nodiscard]] constexpr const Matrix<NU, NX + NI, T>& getK() const { return K; }

private:
    ColVec<NU, T> u_cmd_{}; ///< Raw command from the last feedback
    ColVec<NI, T> e_{};     ///< Tracking error held from feedback for commit
};

} // namespace damp