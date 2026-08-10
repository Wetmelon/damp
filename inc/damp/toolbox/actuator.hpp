// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file toolbox/actuator.hpp
 * @brief Bridge from motion-profile output to drive-native servoactuator commands.
 *
 * The trajectory generators (@ref damp::TrajectoryState, the bank in
 * @c damp/trajectory/polynomial.hpp, and the task-space planners in
 * @c damp/trajectory/cartesian_move.hpp / @c damp/trajectory/topp.hpp) all emit
 * motion in SI joint units — position in rad (or m for a prismatic/linear axis),
 * velocity in rad/s, acceleration in rad/s². A field servo drive instead wants a
 * per-axis setpoint triple in its own mechanical units: `input_pos` [turns],
 * `input_vel` / vel_ff [turns/s], and `input_torque` / torque_ff [Nm at the motor].
 *
 * This header is the thin, allocation-free glue between the two:
 *
 * - ServoCommand — the drive-native `{position, velocity, torque}` triple.
 * - ServoAxis — a single transmission: one signed ratio that maps the joint
 *   unit to motor turns, plus a homing offset. The torque reflection is *derived*
 *   from that same ratio by power conservation (see below), so an axis is fully
 *   described by one number.
 * - ConstantInertiaFeedforward — the per-axis torque-feedforward policy that
 *   matches a decoupled drivetrain (single actuators, Cartesian / CoreXY / gantry /
 *   polar): `τ = J·a + b·v + τ_c·sign(v) + g`. This is the standard inertia
 *   feedforward model, extended with viscous + Coulomb friction
 *   and a settable gravity term. Coupled rigid-body torque for a serial arm (the
 *   configuration-dependent mass matrix `M(q)`, Coriolis, gravity) is a separate
 *   future dynamics policy — it is *not* a per-axis constant and is out of scope here.
 * - ServoBank — maps a whole `array<TrajectoryState, NAxes>` (one
 *   `TrajectoryBank::eval(t)` / `step(dt)` call) straight to an
 *   `array<ServoCommand, NAxes>` you write to the drives.
 *
 * ### Power-consistent transmission
 *
 * For an ideal (lossless, backlash-free) transmission, one signed ratio
 * @p turns_per_unit `r` describes the whole axis. Kinematics scale forward:
 * @f[
 *   p_\text{drive} = r\,q + p_0, \qquad v_\text{drive} = r\,\dot q .
 * @f]
 * Torque reflects *backward*, and the scale follows from power balance
 * @f$\tau_\text{motor}\,\omega_\text{motor} = \tau_\text{joint}\,\dot q@f$ with
 * @f$\omega_\text{motor} = 2\pi r\,\dot q@f$:
 * @f[
 *   \tau_\text{motor} = \frac{\tau_\text{joint}}{2\pi\,r}.
 * @f]
 * The sign of @p r (drive direction) therefore reflects into the torque sign
 * automatically — no separate direction flag is needed. For a rotary @f$n{:}1@f$
 * gearbox on a joint in radians, @f$r = n/2\pi@f$ turns per rad; for a leadscrew of
 * lead @f$L@f$ [m/rev] on a linear axis in metres, @f$r = 1/L@f$ turns per metre.
 *
 * @note Ideal transmission. Reflection assumes lossless, rigid, one-piece
 *       transmission (no efficiency, compliance, or backlash). Add motor rotor
 *       inertia by reflecting it to the joint side, @f$J_\text{joint} =
 *       J_\text{load} + (2\pi r)^2 J_\text{motor}@f$, when populating the
 *       feedforward inertia.
 *
 * @see damp/trajectory/trajectory_types.hpp for the source @ref damp::TrajectoryState.
 * @see damp/toolbox/scaling.hpp for the affine-calibration primitives this builds on.
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/trajectory/trajectory_types.hpp"

namespace damp {

/**
 * @brief A drive-native servoactuator setpoint: position, velocity, torque.
 *
 * Units are the drive's own — typically {turns, turns/s, Nm}. Stream
 * these to the three feedforward inputs of a closed-loop position drive
 * (`input_pos`, `vel_ff`, `torque_ff`).
 *
 * @tparam T Scalar type
 */
template<typename T = double>
struct ServoCommand {
    T position{T{0}}; ///< Commanded position [drive units, e.g. turns]
    T velocity{T{0}}; ///< Velocity feedforward [drive units/s, e.g. turns/s]
    T torque{T{0}};   ///< Torque feedforward at the motor [Nm]

    /// Precision rebind.
    template<typename U>
    [[nodiscard]] constexpr ServoCommand<U> as() const {
        return ServoCommand<U>{static_cast<U>(position), static_cast<U>(velocity), static_cast<U>(torque)};
    }
};

/**
 * @brief One servoactuator transmission: SI joint unit ⟷ drive (motor) units.
 *
 * Holds a single signed kinematic ratio turns_per_unit and a homing offset.
 * Forward kinematics scale joint position/velocity to motor turns; torque reflects
 * backward through the same ratio by power conservation (see the file header). Build
 * one with @ref rotary_gearbox / @ref linear_screw rather than setting fields by hand.
 *
 * @tparam T Scalar type
 */
template<typename T = double>
struct ServoAxis {
    T turns_per_unit{T{1}};  ///< Signed ratio r: motor turns per joint unit (rad or m)
    T position_offset{T{0}}; ///< Added to the position command [turns] (home/zero offset)

    /// Joint position [rad or m] → drive position command [turns].
    [[nodiscard]] constexpr T position_command(T q) const { return (turns_per_unit * q) + position_offset; }

    /// Joint velocity [rad/s or m/s] → velocity feedforward [turns/s].
    [[nodiscard]] constexpr T velocity_command(T v) const { return turns_per_unit * v; }

    /// Joint-side effort (torque [Nm] or force [N]) → motor torque feedforward [Nm].
    /// τ_motor = τ_joint / (2π·r); sign of r reflects into the torque sign.
    [[nodiscard]] constexpr T reflect_torque(T joint_effort) const {
        return joint_effort / (T{2} * damp::numbers::pi_v<T> * turns_per_unit);
    }

    /// Convert a full joint TrajectoryState (position + velocity) plus an
    /// already-computed joint-side effort into a drive command. Torque comes from
    /// the feedforward policy, not the profile, so it is passed in separately.
    [[nodiscard]] constexpr ServoCommand<T> to_command(const TrajectoryState<T>& s, T joint_effort = T{0}) const {
        return ServoCommand<T>{position_command(s.position), velocity_command(s.velocity), reflect_torque(joint_effort)};
    }

    /// Precision rebind.
    template<typename U>
    [[nodiscard]] constexpr ServoAxis<U> as() const {
        return ServoAxis<U>{static_cast<U>(turns_per_unit), static_cast<U>(position_offset)};
    }
};

/**
 * @brief Build a ServoAxis for a rotary joint behind a gearbox.
 *
 * Joint coordinate is in radians. @p gear_ratio is the reduction
 * @f$n@f$ (motor revs per joint rev); the motor turns per joint radian is
 * @f$r = n/2\pi@f$. @p direction (±1) flips the motor's sense relative to the joint.
 *
 * @param gear_ratio     Reduction n:1 (motor revs per output rev), > 0
 * @param direction      +1 or −1 (motor rotation sense vs. joint)
 * @param offset_turns   Homing offset added to the position command [turns]
 */
template<typename T = double>
[[nodiscard]] constexpr ServoAxis<T> rotary_gearbox(T gear_ratio, T direction = T{1}, T offset_turns = T{0}) {
    return ServoAxis<T>{direction * gear_ratio / (T{2} * damp::numbers::pi_v<T>), offset_turns};
}

/**
 * @brief Build a ServoAxis for a linear axis driven by a leadscrew/belt.
 *
 * Joint coordinate is in metres. @p lead is the linear travel per motor
 * revolution [m/rev], so the motor turns per metre is @f$r = 1/L@f$. With this ratio
 * the torque reflection @f$\tau_\text{motor} = F\,L/2\pi@f$ falls out of
 * ServoAxis::reflect_torque automatically (F is the joint-side force [N]).
 *
 * @param lead           Linear travel per motor revolution [m/rev], ≠ 0
 * @param direction      +1 or −1 (motor rotation sense vs. travel)
 * @param offset_turns   Homing offset added to the position command [turns]
 */
template<typename T = double>
[[nodiscard]] constexpr ServoAxis<T> linear_screw(T lead, T direction = T{1}, T offset_turns = T{0}) {
    return ServoAxis<T>{direction / lead, offset_turns};
}

/**
 * @brief Per-axis decoupled torque feedforward: `τ = J·a + b·v + τ_c·sign(v) + g`.
 *
 * The right torque model for a decoupled drivetrain — a single actuator, or an
 * orthogonal machine (Cartesian / CoreXY / gantry / polar) where each axis sees a
 * roughly constant effective inertia. Per axis:
 * @f[
 *   \tau_i = J_i\,\ddot q_i + b_i\,\dot q_i + \tau_{c,i}\,\operatorname{sign}(\dot q_i) + g_i .
 * @f]
 * @ref inertia (@f$J_i@f$, the joint-side reflected inertia incl. reflected rotor
 * inertia) is the standard inertia feedforward term; @ref viscous and
 * @ref coulomb add the friction terms an inertia-only model misses; @ref gravity is
 * a per-axis hold-torque you update each tick (default 0). Effort is returned in
 * joint-side units (Nm for a rotary joint, N for a linear one) — the owning
 * ServoAxis reflects it to motor Nm.
 *
 * @note A serial arm is not decoupled: its effective joint inertia is the
 *       configuration-dependent mass matrix `M(q)` with Coriolis/centrifugal
 *       coupling. That needs a rigid-body (RNEA) dynamics policy fed by per-link
 *       inertial parameters; it is intentionally not modelled here. For an arm,
 *       leave torque at zero (pos/vel feedforward only) or drive @ref gravity from
 *       an external gravity-compensation term.
 *
 * @tparam NAxes Number of axes
 * @tparam T     Scalar type
 */
template<size_t NAxes, typename T = double>
struct ConstantInertiaFeedforward {
    damp::array<T, NAxes> inertia{}; ///< Joint-side effective inertia Jᵢ [kg·m² or kg]
    damp::array<T, NAxes> viscous{}; ///< Viscous friction bᵢ [Nm·s/rad or N·s/m]
    damp::array<T, NAxes> coulomb{}; ///< Coulomb friction τ_c,ᵢ [Nm or N]
    damp::array<T, NAxes> gravity{}; ///< Gravity / hold torque gᵢ [Nm or N], settable per tick

    /// Joint-side feedforward effort for axis @p i from its motion state.
    [[nodiscard]] constexpr T operator()(size_t i, const TrajectoryState<T>& s) const {
        return (inertia[i] * s.acceleration) + (viscous[i] * s.velocity) + (coulomb[i] * damp::sgn(s.velocity)) + gravity[i];
    }
};

/**
 * @brief A bank of ServoAxis transmissions: maps a synchronized multi-axis
 *        TrajectoryState array straight to drive commands in one call.
 *
 * Pairs one-to-one with TrajectoryBank (or any planner producing an
 * `array<TrajectoryState, NAxes>`). The end-to-end loop body is:
 * @code
 *   auto states = bank.step(dt);            // array<TrajectoryState, N> in SI joint units
 *   auto cmds   = servos.command(states, ff); // array<ServoCommand, N> in drive units
 *   for (size_t i = 0; i < N; ++i) {
 *       drive[i].input_pos    = cmds[i].position;   // turns
 *       drive[i].vel_ff       = cmds[i].velocity;   // turns/s
 *       drive[i].torque_ff    = cmds[i].torque;     // Nm
 *   }
 * @endcode
 *
 * @tparam NAxes Number of axes
 * @tparam T     Scalar type
 */
template<size_t NAxes, typename T = double>
class ServoBank {
public:
    using State = TrajectoryState<T>;
    using StateArray = damp::array<State, NAxes>;
    using CommandArray = damp::array<ServoCommand<T>, NAxes>;

    constexpr ServoBank() = default;
    constexpr explicit ServoBank(const damp::array<ServoAxis<T>, NAxes>& axes) : axes_(axes) {}

    [[nodiscard]] constexpr const ServoAxis<T>& axis(size_t i) const { return axes_[i]; }
    constexpr void                              set_axis(size_t i, const ServoAxis<T>& a) { axes_[i] = a; }

    /// Map all axes with a torque-feedforward policy `T effort = policy(i, state)`.
    template<typename TorquePolicy>
    [[nodiscard]] constexpr CommandArray command(const StateArray& states, const TorquePolicy& policy) const {
        CommandArray out{};
        for (size_t i = 0; i < NAxes; ++i) {
            out[i] = axes_[i].to_command(states[i], policy(i, states[i]));
        }
        return out;
    }

    /// Map all axes with no torque feedforward (position + velocity only).
    [[nodiscard]] constexpr CommandArray command(const StateArray& states) const {
        CommandArray out{};
        for (size_t i = 0; i < NAxes; ++i) {
            out[i] = axes_[i].to_command(states[i]);
        }
        return out;
    }

    /// Precision rebind.
    template<typename U>
    [[nodiscard]] constexpr ServoBank<NAxes, U> rebind() const {
        damp::array<ServoAxis<U>, NAxes> a{};
        for (size_t i = 0; i < NAxes; ++i) {
            a[i] = axes_[i].template as<U>();
        }
        return ServoBank<NAxes, U>{a};
    }

private:
    damp::array<ServoAxis<T>, NAxes> axes_{};
};

} // namespace damp
