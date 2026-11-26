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
 * capability maps live in kit headers (spm.hpp; the IPM one is not yet public), not here.
 */


#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/math/transforms.hpp"

namespace damp::motor {

namespace detail {

/**
 * @brief The PI of an R–L current loop 1/(L s + R) at the bandwidth ω: Kp = ω·L, Kᵢ = ω·R
 *
 * The PI's zero at −Kᵢ/Kp = −R/L cancels the plant's pole, so the loop is ω/s and the closed loop
 * from the reference, with b = 1, is ω/(s + ω): one pole at −ω, whatever ω is against R/L. It is
 * the internal-model current controller of Harnefors and Nee, and the formula ODrive tunes its
 * current loop with. The cancelled pole stays in the response to a voltage disturbance, which
 * decays at R/L; a law that places it with the pair is design::pid(plant, Ts,
 * PolePlacement{.omega, .zeta}). With b = 0 (I-P) the reference sees the pole too: ω·(R/L)/((s +
 * R/L)(s + ω)).
 *
 * @param L     [H] inductance (finite, > 0)
 * @param R     [Ω] resistance (finite, ≥ 0; 0 gives Kᵢ = 0, the loop ω/s on a pure inductance)
 * @param omega [rad/s] closed-loop bandwidth ω (finite, > 0)
 * @param b     [-] proportional setpoint weight (1 = PI, 0 = I-P)
 * @return The gains (tracking time 0: the integral holds against a clamp), or nullopt for a number
 *         outside its domain
 *
 * @see "Model-based current control of AC machines using the internal model control method"
 *      (Harnefors & Nee, 1998), IEEE Trans. Ind. Appl. 34(1)
 */
template<typename T>
[[nodiscard]] constexpr damp::optional<design::PIDGains<T>> current_loop_pi(T L, T R, T omega, T b) {
    if (!damp::finite_positive(L) || !damp::finite_non_negative(R) || !damp::finite_positive(omega) || !damp::isfinite(b)) {
        return damp::nullopt;
    }
    return design::PIDGains<T>{.Kp = omega * L, .Ki = omega * R, .b = b};
}

} // namespace detail

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
 * Park / inv-Park and modulation stay on the machine adapter. @ref tune sets both PIs to the
 * bandwidth ω: Kp = ω·L, Kᵢ = ω·R, the closed loop ω/(s + ω) (detail::current_loop_pi).
 *
 * @tparam T Scalar type (embedded default float)
 */
template<typename T = float>
struct FOController {
    using DQ = DirectQuadrature<T>;

    DQ Ldq = {};
    T  R = {};
    T  lambda = {}; ///< [Wb] flux for back-EMF FF (PM λ or ACIM-referred λr)

    /** @brief Optional R·i + L·Δi/Ts feedforward; retune PI if enabled (off by default). */
    bool plant_inversion_ff = false;

    PIController<T> dctrl = {}; ///< d-axis current PI, ticked at the measured Ts (gains from @ref tune; set_gains / set_limits to change them)
    PIController<T> qctrl = {}; ///< q-axis current PI, the same

    constexpr FOController() noexcept = default;

    /**
     * @brief Construct from R–L–λ and seed current-loop bandwidth.
     * @param Ldq_in     [H] dq inductances
     * @param R_in       [ohm] phase resistance
     * @param lambda_in  [Wb] FF flux
     * @param current_bw [rad/s] closed-loop current bandwidth for @ref tune
     * @param b          [-] PI setpoint weight (1 = PI, 0 = I-P)
     */
    constexpr FOController(DQ Ldq_in, T R_in, T lambda_in, T current_bw = T{1000}, T b = T{1}) noexcept
        : Ldq(Ldq_in), R(R_in), lambda(lambda_in) {
        tune(current_bw, b);
    }

    /**
     * @brief Set both dq PI loops to the bandwidth ω on the R–L plant: Kp = ω·L, Kᵢ = ω·R
     *
     * Each PI's zero cancels its axis's R/L pole, so each closed loop is ω/(s + ω) at any ω
     * (detail::current_loop_pi). Each loop takes the gains with set_gains, its integral and limits
     * kept. A refused input (L ≤ 0, R < 0, ω ≤ 0 or a number not finite) keeps the previous gains
     * and returns false; under constant evaluation it does not compile.
     *
     * @param omega_bw [rad/s] closed-loop bandwidth ω
     * @param b        [-] setpoint weight: 1 = PI, 0 = I-P
     * @return true if both loops were set
     */
    constexpr bool tune(T omega_bw, T b = T{1}) noexcept {
        const auto d = detail::current_loop_pi(Ldq.d, R, omega_bw, b);
        const auto q = detail::current_loop_pi(Ldq.q, R, omega_bw, b);
        if (!d || !q) {
            if (damp::is_constant_evaluated()) {
                ::damp::detail::refuse_failed_design(::damp::design::DesignStatus::invalid_parameter); // constant-evaluated (constinit / constexpr): a refused input does not compile
            }
            return false;
        }

        PIController<T> dn = dctrl; // both or neither
        PIController<T> qn = qctrl;
        if (!dn.set_gains(*d) || !qn.set_gains(*q)) {
            if (damp::is_constant_evaluated()) {
                ::damp::detail::refuse_failed_design(::damp::design::DesignStatus::invalid_parameter); // a coefficient overflows
            }
            return false;
        }
        dctrl = damp::move(dn);
        qctrl = damp::move(qn);
        return true;
    }

    /**
     * @brief In place of @ref current_controller while another source drives the bridge
     *
     * Each PI tracks its share of the applied @p Vdq, the voltage less the feedforward the
     * current controller would add, so the next current_controller at these currents, speed and
     * Ts starts from @p Vdq (PIController::track).
     */
    constexpr void track(const DQ& Vdq, const DQ& Idq_ref, const DQ& Idq, const T omega_e, const T Ts) {
        const DQ Vff = feedforward(Idq_ref, Idq, omega_e, Ts);
        dctrl.track(Vdq.d - Vff.d, Idq_ref.d, Idq.d, Ts);
        qctrl.track(Vdq.q - Vff.q, Idq_ref.q, Idq.q, Ts);
    }

    /** @brief Clear both dq integrators. */
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
    [[nodiscard]] constexpr DqCommand<T> current_controller(
        const DQ& Idq_ref,
        const DQ& Idq,
        const T   omega_e,
        const T   Ts,
        const T   Vmax = damp::numeric_limits<T>::max()
    ) {
        const DQ Vff = feedforward(Idq_ref, Idq, omega_e, Ts);

        const T ud = dctrl.step(Idq_ref.d, Idq.d, Ts);
        const T uq = qctrl.step(Idq_ref.q, Idq.q, Ts);
        DQ      Vdq{Vff.d + ud, Vff.q + uq};

        DqCommand<T> cmd;
        // Non-positive ceiling: the bus cannot support a voltage. Unwind the PI
        // toward zero; do not treat it as "no limit".
        if (!damp::finite_positive(Vmax)) {
            dctrl.back_calculate(T{0});
            qctrl.back_calculate(T{0});
            cmd.is_saturated = true;
            cmd.Vdq = {};
            return cmd;
        }

        const T Vmag = Vdq.abs();
        const T Vcap = damp::numeric_limits<T>::max();
        if (Vmag > Vmax) {
            const T scale = Vmax / Vmag;
            // Scale the summed command onto the circle, but charge the integrator
            // only for the PI's own scaled output. Feedforward (back-EMF) saturation
            // belongs to field weakening, not the current integrator.
            dctrl.back_calculate(scale * ud);
            qctrl.back_calculate(scale * uq);
            Vdq = Vdq * scale;
            cmd.is_saturated = true;
            cmd.v_excess = Vmag / Vmax;
        } else if (Vmax < Vcap && Vmag > T{0}) {
            cmd.v_excess = Vmag / Vmax;
        }

        cmd.Vdq = Vdq;
        return cmd;
    }

private:
    /** @brief The decoupling feedforward, and with plant_inversion_ff R·i + L·Δi/Ts. */
    [[nodiscard]] constexpr DQ feedforward(const DQ& Idq_ref, const DQ& Idq, const T omega_e, const T Ts) const {
        DQ Vff = decoupling_feedforward(omega_e, Ldq, lambda, Idq_ref);
        if (plant_inversion_ff) {
            Vff.d += R * Idq.d;
            Vff.q += R * Idq.q;
            if (Ts > T{0}) {
                Vff.d += Ldq.d * (Idq_ref.d - Idq.d) / Ts;
                Vff.q += Ldq.q * (Idq_ref.q - Idq.q) / Ts;
            }
        }
        return Vff;
    }
};

} // namespace damp::motor
