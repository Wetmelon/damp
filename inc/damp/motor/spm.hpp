// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file damp/motor/spm.hpp
 * @brief SPM design helpers (surface PM, id = 0 envelope + datasheet maps).
 *
 * Closed-form Kt/λ/Kv conversions, id=0 torque, and voltage-limited iq ceiling
 * for surface machines. Used with @ref foc.hpp for FOC teaching and host sizing.
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/math/transforms.hpp"

namespace damp::motor {
namespace spm {

/**
 * @brief Kₜ = 1½ p λ [Nm/A] (amplitude-invariant).
 */
template<typename T = double>
[[nodiscard]] constexpr T torque_constant_from_flux(T pole_pairs, T lambda) {
    return static_cast<T>(1.5) * pole_pairs * lambda;
}

/**
 * @brief λ from datasheet Kₜ (amplitude / peak-per-phase convention).
 */
template<typename T = double>
[[nodiscard]] constexpr T flux_from_torque_constant(T pole_pairs, T Kt) {
    return Kt / (static_cast<T>(1.5) * pole_pairs);
}

/**
 * @brief Kₜ from hobby Kᵥ [RPM/V] (peak line-to-line / bus-volt sense).
 * @see Krishnan, PM Synchronous and BLDC Motor Drives, Ch. 9
 */
template<typename T = double>
[[nodiscard]] constexpr T torque_constant_from_Kv(T Kv) {
    constexpr T c = (T{30} / damp::numbers::pi_v<T>)*(damp::numbers::sqrt3_v<T> / T{2});
    return c / Kv;
}

/**
 * @brief λ from Kᵥ via @ref torque_constant_from_Kv.
 */
template<typename T = double>
[[nodiscard]] constexpr T flux_from_Kv(T pole_pairs, T Kv) {
    return flux_from_torque_constant(pole_pairs, torque_constant_from_Kv(Kv));
}

/**
 * @brief Motor constant Kₘ = Kₜ / √R [Nm/√W].
 */
template<typename T = double>
[[nodiscard]] constexpr T motor_constant(T Kt, T R) {
    return Kt / damp::sqrt(R);
}

/**
 * @brief i_q = T_e / Kₜ for id = 0.
 */
template<typename T = double>
[[nodiscard]] constexpr T iq_from_torque(T Te, T pole_pairs, T lambda) {
    return Te / torque_constant_from_flux(pole_pairs, lambda);
}

/**
 * @brief Steady-state Vdq (rotor frame).
 */
template<typename T = double>
[[nodiscard]] constexpr DirectQuadrature<T> steady_state_vdq(
    T omega_e, const DirectQuadrature<T>& Ldq, T R, T lambda, const DirectQuadrature<T>& Idq
) {
    return DirectQuadrature<T>{
        (R * Idq.d) - (omega_e * Ldq.q * Idq.q),
        (R * Idq.q) + (omega_e * ((Ldq.d * Idq.d) + lambda)),
    };
}

/**
 * @brief |Vdq| at a steady-state operating point.
 */
template<typename T = double>
[[nodiscard]] constexpr T steady_state_voltage_magnitude(
    T omega_e, const DirectQuadrature<T>& Ldq, T R, T lambda, const DirectQuadrature<T>& Idq
) {
    const auto v = steady_state_vdq(omega_e, Ldq, R, lambda, Idq);
    return damp::hypot(v.d, v.q);
}

/**
 * @brief Electromagnetic torque of an SPM (i_d = 0): T_e = Kₜ i_q.
 */
template<typename T = double>
[[nodiscard]] constexpr T torque(T pole_pairs, T lambda, T iq) {
    return torque_constant_from_flux(pole_pairs, lambda) * iq;
}

/**
 * @brief Max positive i_q for an SPM (i_d = 0) at electrical speed.
 */
template<typename T = double>
[[nodiscard]] constexpr T max_iq(T omega_e, T Vmax, T Imax, T Lq, T R, T lambda) {
    if (!(Imax > T{0}) || !(Vmax > T{0}) || !(Lq > T{0}) || !(lambda > T{0})) {
        return T{0};
    }
    if (!(damp::abs(omega_e) > T{0})) {
        if (R > T{0}) {
            return damp::min(Imax, Vmax / R);
        }
        return Imax;
    }

    const T we = damp::abs(omega_e);
    const T a = (we * Lq) * (we * Lq) + (R * R);
    const T b = T{2} * R * we * lambda;
    const T c = (we * lambda) * (we * lambda) - (Vmax * Vmax);
    if (!(a > T{0})) {
        return T{0};
    }
    const T disc = (b * b) - (T{4} * a * c);
    if (disc < T{0}) {
        return T{0};
    }
    const T iq_pos = (-b + damp::sqrt(disc)) / (T{2} * a);
    if (!(iq_pos > T{0})) {
        return T{0};
    }
    return damp::min(Imax, iq_pos);
}

/**
 * @brief One sample of the SPM torque–speed envelope (MTPA i_d=0 branch).
 */
template<typename T = double>
struct TorqueSpeedPoint {
    T    omega_mech{};
    T    omega_elec{};
    T    iq{};
    T    Te{};
    T    Vmag{};
    bool voltage_limited{false};
    bool current_limited{false};

    template<typename U>
    [[nodiscard]] constexpr TorqueSpeedPoint<U> as() const {
        return TorqueSpeedPoint<U>{
            static_cast<U>(omega_mech),
            static_cast<U>(omega_elec),
            static_cast<U>(iq),
            static_cast<U>(Te),
            static_cast<U>(Vmag),
            voltage_limited,
            current_limited,
        };
    }
};

/**
 * @brief SPM max-torque point at a mechanical speed (i_d = 0 only).
 */
template<typename T = double>
[[nodiscard]] constexpr TorqueSpeedPoint<T>
max_torque_at_speed(T omega_mech, T Vmax, T Imax, T Lq, T R, T lambda, T pole_pairs) {
    if (!(pole_pairs > T{0})) {
        return TorqueSpeedPoint<T>{};
    }
    const T omega_elec = damp::abs(omega_mech) * pole_pairs;
    const T iq = max_iq(omega_elec, Vmax, Imax, Lq, R, lambda);
    const T Te = torque(pole_pairs, lambda, iq);
    const T Vmag = steady_state_voltage_magnitude(
        omega_elec, DirectQuadrature<T>{Lq, Lq}, R, lambda, DirectQuadrature<T>{T{0}, iq}
    );
    const bool current_limited = (iq >= (Imax * (T{1} - static_cast<T>(1e-9))));
    bool       voltage_limited = !current_limited && (iq > T{0});
    if (!(iq > T{0}) && (omega_elec > T{0})) {
        voltage_limited = true;
    }
    return TorqueSpeedPoint<T>{
        .omega_mech = omega_mech,
        .omega_elec = omega_elec,
        .iq = iq,
        .Te = Te,
        .Vmag = Vmag,
        .voltage_limited = voltage_limited,
        .current_limited = current_limited,
    };
}

/**
 * @brief Fixed table of SPM torque–speed envelope samples (host / LUT bake).
 */
template<std::size_t N, typename T = double>
[[nodiscard]] constexpr damp::array<TorqueSpeedPoint<T>, N>
torque_speed_envelope(T omega_mech_max, T Vmax, T Imax, T Lq, T R, T lambda, T pole_pairs) {
    damp::array<TorqueSpeedPoint<T>, N> table{};
    if (N == 0) {
        return table;
    }
    if (!(omega_mech_max > T{0})) {
        for (std::size_t i = 0; i < N; ++i) {
            table[i] = max_torque_at_speed(T{0}, Vmax, Imax, Lq, R, lambda, pole_pairs);
        }
        return table;
    }
    for (std::size_t i = 0; i < N; ++i) {
        const T frac = (N == 1) ? T{0} : static_cast<T>(i) / static_cast<T>(N - 1);
        table[i] = max_torque_at_speed(frac * omega_mech_max, Vmax, Imax, Lq, R, lambda, pole_pairs);
    }
    return table;
}

} // namespace spm
} // namespace damp::motor
