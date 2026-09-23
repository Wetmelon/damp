// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file foc.hpp
 * @brief Machine-agnostic dq current regulator and related closed-form helpers.
 *
 * Stage-2 electrical pipeline primitive: PI + decoupling FF → clamped Vdq.
 * Frame transforms and modulation belong to the machine adapter. Plant torque /
 * capability maps live in kit headers (spm.hpp / ipm.hpp), not here.
 */

#include <limits>

#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/pid_design.hpp" // ::damp::design::pi_pole_placement_first_order
#include "damp/math/transforms.hpp"

namespace damp::motor {

/**
 * @defgroup foc_design Field-Oriented Control Design
 * @brief Closed-form maps for current loops and SVPWM geometry.
 */

/**
 * @brief Constant-power torque ceiling Tₘₐₓ = P_rated / |ω|.
 * @ingroup foc_design
 * @param P_rated   [W] rated / base-speed power
 * @param omega_rad [rad/s] mechanical speed
 * @param T_peak    [Nm] peak torque in the constant-torque region
 * @param omega_min [rad/s] speed floor (avoid /0)
 * @return Tₘₐₓ [Nm]
 */
template<typename T = double>
[[nodiscard]] constexpr T constant_power_torque_limit(T P_rated, T omega_rad, T T_peak, T omega_min = static_cast<T>(1e-3)) {
    const T w = damp::abs(omega_rad);
    if (w <= omega_min) {
        return T_peak;
    }
    return damp::min(T_peak, P_rated / w);
}

/**
 * @brief SVPWM voltage-circle radius Vₘₐₓ = m · Vdc / √3.
 * @ingroup foc_design
 * @param Vdc            [V] DC-bus voltage
 * @param max_modulation [-] fraction of the linear circle (default 1)
 * @return Vₘₐₓ [V]
 */
template<typename T = double>
[[nodiscard]] constexpr T voltage_circle_radius(T Vdc, T max_modulation = T{1}) {
    return max_modulation * Vdc * damp::numbers::inv_sqrt3_v<T>;
}

/**
 * @brief Decoupling + back-EMF feedforward at R = 0 (reference currents).
 * @ingroup foc_design
 *
 * @f$ v_d = -\omega L_q i_q,\quad v_q = \omega(L_d i_d + \lambda) @f$.
 * @p lambda is the flux used for FF (PM linkage or ACIM-referred rotor flux).
 */
template<typename T = double>
[[nodiscard]] constexpr DirectQuadrature<T> decoupling_feedforward(
    T omega_e, const DirectQuadrature<T>& Ldq, T lambda, const DirectQuadrature<T>& Idq_ref
) {
    return DirectQuadrature<T>{
        -(omega_e * Ldq.q * Idq_ref.q),
        omega_e * ((Ldq.d * Idq_ref.d) + lambda),
    };
}

/**
 * @brief Clamped dq voltage command from @ref FOController::current_controller.
 * @tparam T Scalar type
 */
template<typename T = float>
struct DqCommand {
    DirectQuadrature<T> Vdq = {};             ///< [V] dq voltage target (|Vdq| ≤ Vmax)
    bool                is_saturated = false; ///< |Vdq| hit the voltage circle
    T                   v_excess = T{0};      ///< |Vdq|/Vmax before limit; 0 if unlimited
};

/**
 * @brief dq current regulator: PI + decoupling FF → clamped Vdq.
 *
 * Park / inv-Park and modulation stay on the machine adapter. Tune with
 * @ref ::damp::design::pi_pole_placement_first_order via @ref tune.
 *
 * @tparam T Scalar type (embedded default float)
 */
template<typename T = float>
struct FOController {
    using DQ = DirectQuadrature<T>;

    DQ Ldq = {};
    T  R = {};
    T  lambda = {}; ///< [Wb] flux for back-EMF FF (PM λ or ACIM-referred λr)

    /// Fraction of Vdc/√3 available to the regulator (1 = full linear SVPWM circle).
    T max_modulation = T{1};

    /// Optional R·i + L·Δi/Ts feedforward; retune PI if enabled (off by default).
    bool plant_inversion_ff = false;

    ContinuousPID<T> dctrl = {};
    ContinuousPID<T> qctrl = {};

    constexpr FOController() = default;

    /**
     * @brief Construct from R–L–λ and seed current-loop bandwidth.
     * @param Ldq_in     [H] dq inductances
     * @param R_in       [ohm] phase resistance
     * @param lambda_in  [Wb] FF flux
     * @param current_bw [rad/s] closed-loop current bandwidth for @ref tune
     */
    constexpr FOController(DQ Ldq_in, T R_in, T lambda_in, T current_bw = T{1000}, T b = T{1})
        : Ldq(Ldq_in), R(R_in), lambda(lambda_in) {
        tune(current_bw, T{1}, b);
    }

    /**
     * @brief Place both dq PI loops on the R–L plant.
     * @param omega_bw [rad/s] closed-loop bandwidth
     * @param zeta     [-] damping (default 1)
     * @param b        [-] setpoint weight: 1 = PI, 0 = I-P
     */
    constexpr void tune(T omega_bw, T zeta = T{1}, T b = T{1}) {
        const auto d = ::damp::design::pi_pole_placement_first_order(Ldq.d, R, omega_bw, zeta, b);
        const auto q = ::damp::design::pi_pole_placement_first_order(Ldq.q, R, omega_bw, zeta, b);

        dctrl.Kp = d.Kp;
        dctrl.Ki = d.Ki;
        dctrl.Kbc = d.Kbc;
        dctrl.b = d.b;

        qctrl.Kp = q.Kp;
        qctrl.Ki = q.Ki;
        qctrl.Kbc = q.Kbc;
        qctrl.b = q.b;
    }

    /// Bumpless disable: integrators track @p Vdq_track.
    constexpr void disable(const DQ& Vdq_track = {}) {
        dctrl.disable(Vdq_track.d);
        qctrl.disable(Vdq_track.q);
    }

    /// Resume Auto on both dq loops.
    constexpr void enable() {
        dctrl.enable();
        qctrl.enable();
    }

    /// Clear both dq integrators.
    constexpr void reset() {
        dctrl.reset();
        qctrl.reset();
    }

    /**
     * @brief One current-loop step → clamped dq voltage.
     *
     * Decoupling FF → optional plant-inversion FF → PI → circular |Vdq|≤Vmax + AW.
     *
     * @param Idq_ref [A] dq current reference
     * @param Idq     [A] measured dq current
     * @param omega_e [rad/s] electrical speed
     * @param Ts      [s] sample time
     * @param Vmax    [V] voltage-circle radius (default: no limit). Non-positive
     *                forces zero voltage and unwinds the PI (collapsed bus).
     * @return Clamped Vdq plus saturation flags
     */
    [[nodiscard]] DqCommand<T> current_controller(
        const DQ& Idq_ref,
        const DQ& Idq,
        const T   omega_e,
        const T   Ts,
        const T   Vmax = std::numeric_limits<T>::max()
    ) {
        DQ Vff = decoupling_feedforward(omega_e, Ldq, lambda, Idq_ref);

        if (plant_inversion_ff) {
            Vff.d += R * Idq.d;
            Vff.q += R * Idq.q;
            if (Ts > T{0}) {
                Vff.d += Ldq.d * (Idq_ref.d - Idq.d) / Ts;
                Vff.q += Ldq.q * (Idq_ref.q - Idq.q) / Ts;
            }
        }

        const T ud = dctrl.control(Idq_ref.d, Idq.d, Ts);
        const T uq = qctrl.control(Idq_ref.q, Idq.q, Ts);
        DQ      Vdq{Vff.d + ud, Vff.q + uq};

        DqCommand<T> cmd;
        // Non-positive ceiling: the bus cannot support a voltage. Unwind the PI
        // toward zero; do not treat it as "no limit".
        if (!(Vmax > T{0})) {
            dctrl.back_calculate(ud, T{0}, Ts);
            qctrl.back_calculate(uq, T{0}, Ts);
            cmd.is_saturated = true;
            cmd.Vdq = {};
            return cmd;
        }

        const T Vmag = Vdq.abs();
        const T Vcap = std::numeric_limits<T>::max();
        if (Vmag > Vmax) {
            const T scale = Vmax / Vmag;
            // Scale the summed command onto the circle, but charge the integrator
            // only for the PI's own scaled output. Feedforward (back-EMF) saturation
            // belongs to field weakening, not the current integrator.
            dctrl.back_calculate(ud, scale * ud, Ts);
            qctrl.back_calculate(uq, scale * uq, Ts);
            Vdq = Vdq * scale;
            cmd.is_saturated = true;
            cmd.v_excess = Vmag / Vmax;
        } else if (Vmax < Vcap && Vmag > T{0}) {
            cmd.v_excess = Vmag / Vmax;
        }

        cmd.Vdq = Vdq;
        return cmd;
    }
};

} // namespace damp::motor
