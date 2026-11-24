// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file steady_state_target.hpp
 * @brief The steady-state target a law tracks through: design::steady_state_target → SteadyStateTargetResult → SteadyStateTarget
 *
 * The target holds the controlled outputs @f$ z = C_z x + D_z u + C_{d,z} \hat d @f$ (rows of
 * the measured outputs) at r under the estimated disturbance:
 * @f[
 *   \begin{bmatrix} A-I & B \\ C_z & D_z \end{bmatrix}
 *   \begin{bmatrix} x_t \\ u_t \end{bmatrix}
 *   =
 *   \begin{bmatrix} -B_d\hat d \\ r - C_{d,z}\hat d \end{bmatrix},
 * @f]
 * solved once at design time for unit d and r, so the runtime target is two matrix–vector
 * products (with NZ < NU the minimum-norm solution). More held outputs than inputs does not
 * compile: no steady state holds them all. Shared by the laws that track through a target
 * (OFLQG, OFMPC); a law with a disturbance model re-anchors its estimator (a FixedGainObserver
 * with design::kalman's fixed gain) on it. A MeasuredDisturbance d, read each tick, enters
 * beside d̂ through its own B_d (none on the readings): the right side is then under [d̂; d].
 *
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings,
 *      2003), AIChE Journal 49(2)
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/law_detail.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/options.hpp"
#include "damp/design/status.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/matrix/svd.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @brief The outputs a law without integrators or a disturbance model holds at r through its
 *        steady-state target (design::state_feedback and design::mpc, with or without a KalmanNoise,
 *        and design::steady_state_target, the target alone): `TrackedOutputs{{0}}`, in that order
 *        (the references r follow it). Left out, the law is a regulator and takes no r.
 *
 * The target is the model's: x* = Nₓ·r, u* = Nᵤ·r. A model error or a constant disturbance
 * leaves an offset; integral action (lqi, lqgi) or a disturbance model (oflqg, an offset-free
 * mpc, whose spec's outputs are the held ones) removes it. An MPC's QP still weighs every output,
 * those it does not list toward C·x*. Beside a disturbance spec it does not compile: the spec's
 * .outputs are the held ones. Nor does a list longer than the plant's inputs, or TrackedOutputs{}
 * on a plant with more outputs than inputs: a steady state holds at most one output per input.
 *
 * @see "Feedback Control of Dynamic Systems" (Franklin, Powell & Emami-Naeini), the reference
 *      input of state-space design (N_x, N_u)
 *
 * @ref outputs is an array, so `TrackedOutputs{{0, 2}}` takes NZ from its list; `TrackedOutputs{}` is
 * every measured output.
 *
 * @tparam NZ Outputs held at r, listed in @ref outputs (0: every measured output)
 */
template<size_t NZ = 0>
struct TrackedOutputs {
    using option_kind = kind::TrackedOutputs; ///< Its kind of design option
    size_t outputs[NZ]{};                     ///< Indices of the measured outputs held at r, each once
};

/** @brief Every measured output held at r: no list. */
template<>
struct TrackedOutputs<0> {
    using option_kind = kind::TrackedOutputs; ///< Its kind of design option
};

/**
 * @struct SteadyStateTargetResult
 * @brief The steady-state target maps: [x_t; u_t] per unit disturbance and per unit reference
 *
 * With a MeasuredDisturbance, target_d is per unit [d̂; d], the estimated disturbance's columns
 * then the measured one's NM, as every law's target is: SteadyStateTarget::equilibrium(r, [d̂; d])
 * gives the (x_t, u_t) the law holds under that d.
 *
 * @tparam NX States of the filter's model before the disturbance (the plant)
 * @tparam ND Disturbance columns: the spec's states, then a MeasuredDisturbance's inputs (0: the target from r alone)
 * @tparam NZ Controlled outputs (the listed measured outputs, or all of them)
 */
template<size_t NX, size_t NU, size_t ND, size_t NZ, typename T = double>
struct [[nodiscard]] SteadyStateTargetResult {
    /** @brief The maps in another scalar: as<U>()'s type. */
    template<typename U>
    using as_t = SteadyStateTargetResult<NX, NU, ND, NZ, U>;

    Matrix<NX + NU, ND, T> target_d{};                          ///< [x_t; u_t] per unit estimated disturbance, then per unit measured disturbance
    Matrix<NX + NU, NZ, T> target_r{};                          ///< [x_t; u_t] per unit controlled-output reference
    bool                   success{false};                      ///< status == DesignStatus::ok
    DesignStatus           status{DesignStatus::not_converged}; ///< ok, or why the target is refused

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr SteadyStateTargetResult<NX, NU, ND, NZ, U> as() const {
        return SteadyStateTargetResult<NX, NU, ND, NZ, U>{target_d.template as<U>(), target_r.template as<U>(), success, status};
    }
};

namespace detail {

/** @brief A spec that lists outputs: a disturbance spec or TrackedOutputs with an .outputs list (none: every output). */
template<typename S>
concept ListsOutputs = requires(const S& spec) { spec.outputs; };

/** @brief The outputs a disturbance spec or TrackedOutputs lists, as an array: none when it lists none (every output). */
template<typename Spec>
[[nodiscard]] constexpr auto listed_outputs(const Spec& spec) {
    if constexpr (ListsOutputs<Spec>) {
        constexpr size_t        NZ = sizeof(spec.outputs) / sizeof(size_t);
        damp::array<size_t, NZ> out{};
        for (size_t i = 0; i < NZ; ++i) {
            out[i] = spec.outputs[i];
        }
        return out;
    } else {
        static_cast<void>(spec);
        return damp::array<size_t, 0>{};
    }
}

/** @brief The disturbance channels of NoDisturbance: none. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXE, size_t NWE, size_t NVE>
[[nodiscard]] constexpr DisturbanceChannels<NX, NXE, NY, 0, T> disturbance_channels(
    const NoDisturbance&,
    const StateSpace<NX, NU, NY, T, NW, NV>&,
    const StateSpace<NXE, NU, NY, T, NWE, NVE>&,
    T
) {
    return {};
}

/** @brief The disturbance channels of TrackedOutputs: none. */
template<size_t NZ, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXE, size_t NWE, size_t NVE>
[[nodiscard]] constexpr DisturbanceChannels<NX, NXE, NY, 0, T> disturbance_channels(
    const TrackedOutputs<NZ>&,
    const StateSpace<NX, NU, NY, T, NW, NV>&,
    const StateSpace<NXE, NU, NY, T, NWE, NVE>&,
    T
) {
    return {};
}

/** @brief Tracked outputs: no disturbance states. */
template<size_t NZ, size_t NU, size_t NY, size_t NW>
struct disturbance_traits<TrackedOutputs<NZ>, NU, NY, NW> {
    static constexpr size_t states = 0;  ///< ND
    static constexpr size_t listed = NZ; ///< NZ listed
};

/** @brief Not TrackedOutputs. */
template<typename>
inline constexpr bool is_tracked_outputs = false;

/** @brief TrackedOutputs. */
template<size_t NZ>
inline constexpr bool is_tracked_outputs<TrackedOutputs<NZ>> = true;

/** @brief Whether a design's options name its held outputs in one place: no TrackedOutputs beside a disturbance spec, whose .outputs are the held ones. */
template<typename... Opts>
inline constexpr bool held_in_one_place = (option_count<kind::TrackedOutputs, Opts...> == 0) || !(DisturbanceSpec<Opts> || ...);

/** @brief The controlled outputs of a target under @p Spec: those it lists, else every output. */
template<typename Spec, size_t NU, size_t NY, size_t NW>
inline constexpr size_t target_outputs = (disturbance_traits<Spec, NU, NY, NW>::listed == 0) ? NY : disturbance_traits<Spec, NU, NY, NW>::listed;

/**
 * @brief The outputs a target under @p Aim (target_spec_of) holds at r: those a disturbance spec or
 *        a TrackedOutputs lists, every output when it lists none; none for a regulator
 *        (NoDisturbance), which takes no r.
 */
template<typename Aim, size_t NY>
inline constexpr size_t reference_outputs = damp::is_same_v<Aim, NoDisturbance> ? 0 : target_outputs<Aim, 0, NY, 0>;

/**
 * @brief How many outputs the target of a design with @p Opts holds, at r or at 0: those a
 *        disturbance spec lists, else those a TrackedOutputs lists, every one of the NY when its
 *        list is empty; with neither the design is a regulator, whose target under a
 *        MeasuredDisturbance holds every output at 0, or an LQI's @p NI integrated ones, and
 *        which has none without one.
 */
template<size_t NU, size_t NY, size_t NI, typename... Opts>
consteval size_t held_count() {
    using Spec = option_t<kind::Disturbance, NoDisturbance, Opts...>;
    using Held = option_t<kind::TrackedOutputs, NoDisturbance, Opts...>;
    if constexpr (DisturbanceSpec<Spec>) {
        return target_outputs<Spec, NU, NY, 0>;
    } else if constexpr (is_tracked_outputs<Held>) {
        return target_outputs<Held, NU, NY, 0>;
    } else if constexpr (option_count<kind::MeasuredDisturbance, Opts...> > 0) {
        return (NI > 0) ? NI : NY;
    } else {
        return 0;
    }
}

/**
 * @brief The compile-time check of the outputs a design's target holds, the same in every design
 *        with one (design::state_feedback, design::mpc, design::steady_state_target): no
 *        TrackedOutputs beside a disturbance spec, and at most NU outputs held. One TrackedOutputs
 *        at most is the design's option set's check.
 *
 * A steady state holds at most one output per input: with more, none holds them all, and the
 * nearest one in least squares holds none of them, so a design that asks for more does not
 * compile and names the outputs it holds. A regulator (no spec, no TrackedOutputs) has a target
 * only under a MeasuredDisturbance, which holds every output at 0, or, on an LQI, its @p NI
 * integrated outputs (0: no IntegralAction).
 */
template<size_t NU, size_t NY, size_t NI, typename... Opts>
constexpr void check_held_outputs() {
    static_assert(
        held_in_one_place<Opts...>,
        "A disturbance spec holds the outputs it lists at r: give them as its .outputs, InputDisturbance{{i, ...}} or "
        "DisturbanceModel{.Bd = ..., .outputs = {i, ...}}, not as TrackedOutputs."
    );
    static_assert(
        held_count<NU, NY, NI, Opts...>() <= NU,
        "No steady state holds more outputs than the plant has inputs: list the ones the target holds, at most one per input, "
        "TrackedOutputs{{i, ...}} (r = 0 for a regulator under a MeasuredDisturbance) or a disturbance spec's .outputs = {i, ...}; "
        "an LQI under a MeasuredDisturbance holds its IntegralAction's .outputs."
    );
}

/**
 * @brief The spec a target is designed on, from a design's disturbance spec and TrackedOutputs, the
 *        one pick of every design: the disturbance model, which holds the outputs it lists, else
 *        the outputs a TrackedOutputs @p held lists (every output when it lists none), else
 *        NoDisturbance, a regulator, which holds nothing at r. Anything else in either place (a
 *        NoDisturbance, an IntegralAction, an option left out) is not a spec.
 */
template<typename Spec, typename Held>
[[nodiscard]] constexpr auto target_spec_of(const Spec& spec, const Held& held) {
    static_cast<void>(spec);
    static_cast<void>(held);
    if constexpr (DisturbanceSpec<Spec>) {
        return spec;
    } else if constexpr (is_tracked_outputs<Held>) {
        return held;
    } else {
        return NoDisturbance{};
    }
}

/** @brief The options steady_state_target takes: the model features of the design (a KalmanNoise is ignored). */
using target_options = option_set<kind::Disturbance, kind::TrackedOutputs, kind::MeasuredDisturbance, kind::KalmanNoise, kind::InputDelay, kind::OutputDelay>;

/**
 * @brief The target maps of a model: [A − I, B; C_z, D_z]·[x_t; u_t] = [−B_d·d; r − C_{d,z}·d]
 *        solved for unit d and r, @p M the left side (see steady_state_target): square with
 *        NZ = NU, the minimum-norm solution with fewer. More (NZ > NU) does not reach it from a
 *        design: check_held_outputs refuses it at compile time.
 */
template<size_t NXE, size_t NU, size_t NZ, size_t ND, typename T>
[[nodiscard]] constexpr SteadyStateTargetResult<NXE, NU, ND, NZ, T> target_solve(const Matrix<NXE + NZ, NXE + NU, T>& M, const Matrix<NXE, ND, T>& Bd, const Matrix<NZ, ND, T>& Cdz) {
    using Result = SteadyStateTargetResult<NXE, NU, ND, NZ, T>;
    //! d enters the state update through B_d and the controlled outputs through their C_d rows:
    //! z = C_z x + D_z u + C_{d,z} d is held at r.
    Matrix<NXE + NZ, ND, T> rhs_d{};
    rhs_d.template block<NXE, ND>(0, 0) = -Bd;
    for (size_t i = 0; i < NZ; ++i) {
        for (size_t j = 0; j < ND; ++j) {
            rhs_d(NXE + i, j) = -Cdz(i, j);
        }
    }
    Matrix<NXE + NZ, NZ, T> rhs_r{};
    rhs_r.template block<NZ, NZ>(NXE, 0) = Matrix<NZ, NZ, T>::identity();
    Result out{};
    if constexpr (NZ == NU) {
        const auto td = mat::solve(M, rhs_d);
        const auto tr = mat::solve(M, rhs_r);
        if (!td || !tr) {
            return refused<Result>(DesignStatus::singular); // no unique steady state for these outputs
        }
        out.target_d = td.value();
        out.target_r = tr.value();
    } else {
        const auto Mp = mat::pseudo_inverse(M);
        if (!Mp) {
            return refused<Result>(DesignStatus::not_converged); // its sweeps ran out: M is finite here
        }
        out.target_d = *Mp * rhs_d;
        out.target_r = *Mp * rhs_r;
        {
            //! NZ < NU: a steady state holds r only if [A − I, B; C_z, D_z] has full row rank. The
            //! minimum-norm solution must then solve the system to roundoff (relative to ‖M‖·‖T‖,
            //! which carries the conditioning); a rank-deficient M leaves an O(1) residual.
            const T tol = static_cast<T>(NXE + NU) * T{64} * eps<T>();
            const T res_d = Matrix<NXE + NZ, ND, T>((M * out.target_d) - rhs_d).norm();
            const T res_r = Matrix<NXE + NZ, NZ, T>((M * out.target_r) - rhs_r).norm();
            if (!(res_d <= tol * ((M.norm() * out.target_d.norm()) + rhs_d.norm())) || !(res_r <= tol * ((M.norm() * out.target_r.norm()) + rhs_r.norm()))) {
                return refused<Result>(DesignStatus::singular); // no steady state holds every r
            }
        }
    }
    return set_status(out, DesignStatus::ok);
}

/** @brief [A − I, B; C_z, D_z] of a model, C_z and D_z its outputs @p Cz, @p Dz held at r. */
template<size_t NXE, size_t NU, size_t NZ, typename T>
[[nodiscard]] constexpr Matrix<NXE + NZ, NXE + NU, T> target_system(
    const Matrix<NXE, NXE, T>& A,
    const Matrix<NXE, NU, T>&  B,
    const Matrix<NZ, NXE, T>&  Cz,
    const Matrix<NZ, NU, T>&   Dz
) {
    Matrix<NXE + NZ, NXE + NU, T> M{};
    M.template block<NXE, NXE>(0, 0) = A - Matrix<NXE, NXE, T>::identity();
    M.template block<NXE, NU>(0, NXE) = B;
    M.template block<NZ, NXE>(NXE, 0) = Cz;
    M.template block<NZ, NU>(NXE, NXE) = Dz;
    return M;
}

/**
 * @brief The target maps on the filter's model of @p sys under the options (see
 *        steady_state_target), in the plant's time domain (@p D plant) or sampled at @p Ts,
 *        with the delay states of @p delay. @p measured, a measured disturbance's channel on
 *        the same model (measured_channel), adds its NM columns after the spec's: target_d is
 *        then [x_t; u_t] per unit [d̂; d], and the measured d reaches no reading directly.
 */
template<Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Spec, typename Delay = NoDelay, size_t NM = 0>
[[nodiscard]] constexpr auto target_design(
    const StateSpace<NX, NU, NY, T, NW, NV>&               sys,
    T                                                      Ts,
    const Spec&                                            spec,
    const Delay&                                           delay = Delay{},
    const Matrix<NX + delay_states<Delay, NU, NY>, NM, T>& measured = Matrix<NX + delay_states<Delay, NU, NY>, NM, T>{}
) {
    constexpr size_t NXE = NX + delay_states<Delay, NU, NY>;
    constexpr size_t ND = disturbance_traits<Spec, NU, NY, NW>::states;
    constexpr size_t NZ = target_outputs<Spec, NU, NY, NW>;
    using Result = SteadyStateTargetResult<NXE, NU, ND + NM, NZ, T>;
    const auto m = sampled_model(sys, Ts, D, delay);
    if (m.status != DesignStatus::ok) {
        return refused<Result>(m.status);
    }
    const auto rows = selected_outputs<NY>(listed_outputs(spec));
    if (!rows) {
        return refused<Result>(DesignStatus::bad_weights); // a controlled output out of range or listed twice
    }
    const auto ch = on_late_readings(spec, disturbance_channels(spec, sys, m.law, Ts), sys, Ts, delay); // the filter's: d on a late reading as it was when taken
    if (!mat::all_of(ch.Bd, [](T v) { return damp::isfinite(v); }) || !mat::all_of(ch.Cd, [](T v) { return damp::isfinite(v); })) {
        return refused<Result>(DesignStatus::bad_weights); // a non-finite B_d or C_d
    }

    //! [A − I, B; C_z, D_z]: z the listed rows of the law's model, C_{d,z} the same rows of C_d
    Matrix<NZ, NXE, T>     Cz{};
    Matrix<NZ, NU, T>      Dz{};
    Matrix<NZ, ND + NM, T> Cdz{};
    for (size_t i = 0; i < NZ; ++i) {
        for (size_t j = 0; j < NXE; ++j) {
            Cz(i, j) = m.law.C((*rows)[i], j);
        }
        for (size_t j = 0; j < NU; ++j) {
            Dz(i, j) = m.law.D((*rows)[i], j);
        }
        for (size_t j = 0; j < ND; ++j) {
            Cdz(i, j) = ch.Cd((*rows)[i], j);
        }
    }
    return target_solve<NXE, NU>(target_system(m.law.A, m.law.B, Cz, Dz), side_by_side(ch.Bd, measured), Cdz);
}

/**
 * @brief The target of a law whose r is not the target's (an LQ regulator, an LQI) under a
 *        measured disturbance alone: the steady state that holds @p S·y at 0 under d on the
 *        model of @p sys with @p delay, @p measured the disturbance's channel on it. Only its
 *        map of d is kept (target_r has no columns): an LQ regulator holds every output at 0
 *        (S = I), an LQI its integrated outputs (S = Ci), whose integrators then hold r.
 */
template<Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NS, typename Delay, size_t NM>
[[nodiscard]] constexpr SteadyStateTargetResult<NX + delay_states<Delay, NU, NY>, NU, NM, 0, T> held_at_zero(
    const StateSpace<NX, NU, NY, T, NW, NV>&               sys,
    T                                                      Ts,
    const Matrix<NS, NY, T>&                               S,
    const Delay&                                           delay,
    const Matrix<NX + delay_states<Delay, NU, NY>, NM, T>& measured
) {
    constexpr size_t NXE = NX + delay_states<Delay, NU, NY>;
    using Result = SteadyStateTargetResult<NXE, NU, NM, 0, T>;
    const auto m = sampled_model(sys, Ts, D, delay);
    if (m.status != DesignStatus::ok) {
        return refused<Result>(m.status);
    }
    const auto full = target_solve<NXE, NU>(target_system(m.law.A, m.law.B, Matrix<NS, NXE, T>(S * m.law.C), Matrix<NS, NU, T>(S * m.law.D)), measured, Matrix<NS, NM, T>{});
    Result     out{};
    out.target_d = full.target_d;
    return set_status(out, full.status);
}

/**
 * @brief The target held_target returns for the pick @p Aim on the model of a plant of NX states,
 *        NU inputs, NY outputs and NW noise inputs with @p Delay, under NM measured disturbances:
 *        [x_t; u_t] per unit [d̂; d] (the spec's states, then NM), per unit r on the outputs it holds.
 */
template<typename Aim, size_t NX, size_t NU, size_t NY, size_t NW, typename Delay, size_t NM, typename T>
using held_target_t = SteadyStateTargetResult<NX + delay_states<Delay, NU, NY>, NU, disturbance_traits<Aim, NU, NY, NW>::states + NM, reference_outputs<Aim, NY>, T>;

/**
 * @brief The target a design holds on its pick @p aim (target_spec_of), @p measured a measured
 *        disturbance's channel on the model of @p sys with @p delay: target_design on a
 *        disturbance spec or a TrackedOutputs; a regulator's (NoDisturbance) the steady state
 *        under d that holds every output at 0 (held_at_zero), or none without a measured d, the
 *        target x = 0, u = 0, with no map at all. design::mpc's and design::steady_state_target's;
 *        design::state_feedback's is the same, an LQI holding its integrated outputs at 0.
 */
template<Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Aim, typename Delay, size_t NM>
[[nodiscard]] constexpr held_target_t<Aim, NX, NU, NY, NW, Delay, NM, T> held_target(
    const StateSpace<NX, NU, NY, T, NW, NV>&               sys,
    T                                                      Ts,
    const Aim&                                             aim,
    const Delay&                                           delay,
    const Matrix<NX + delay_states<Delay, NU, NY>, NM, T>& measured
) {
    if constexpr (!damp::is_same_v<Aim, NoDisturbance>) {
        return target_design<D>(sys, Ts, aim, delay, measured);
    } else if constexpr (NM > 0) {
        return held_at_zero<D>(sys, Ts, Matrix<NY, NY, T>::identity(), delay, measured);
    } else {
        static_cast<void>(sys);
        static_cast<void>(Ts);
        static_cast<void>(delay);
        static_cast<void>(measured);
        SteadyStateTargetResult<NX + delay_states<Delay, NU, NY>, NU, 0, 0, T> none{};
        return set_status(none, DesignStatus::ok);
    }
}

/** @brief target_with once its options fit: held_target under the measured disturbance's channel on the model of @p sys with @p delay. */
template<Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Spec, typename Delay, size_t NM>
[[nodiscard]] constexpr auto held_target_under(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const Spec& spec, const Delay& delay, const MeasuredDisturbance<NX, NM, T>& measured) {
    const auto channel = measured_channel<D>(sys, Ts, measured, delay);
    using Result = held_target_t<Spec, NX, NU, NY, NW, Delay, NM, T>;
    if (!channel) {
        const auto m = sampled_model(sys, Ts, D, delay);
        return refused<Result>((m.status != DesignStatus::ok) ? m.status : DesignStatus::bad_weights); // the model's refusal first, as the designs'; else a non-finite B_d
    }
    return held_target<D>(sys, Ts, spec, delay, *channel);
}

/**
 * @brief steady_state_target on its options, on the model of @p sys with @p delay: held_target on
 *        the outputs the laws hold (target_spec_of: the disturbance spec's, else the
 *        TrackedOutputs', else a regulator's under a MeasuredDisturbance) and, after the spec's
 *        columns, the MeasuredDisturbance's channel on that model (measured_channel), the target
 *        every design with the same options holds. Without one (NM = 0) the channel has no
 *        columns: the spec's target alone. With neither a spec nor a TrackedOutputs nor a
 *        MeasuredDisturbance there is nothing to hold, and it does not compile.
 * @return refused as the model is (a sample time, a delay), else with DesignStatus::bad_weights
 *         for a non-finite B_d
 */
template<Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Delay, typename... Opts>
[[nodiscard]] constexpr auto target_with(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const Delay& delay, const Opts&... opts) {
    target_options::check<Opts...>();
    check_held_outputs<NU, NY, 0, Opts...>();
    constexpr bool   kModelFits = check_delays<D, Opts...>() && check_disturbance<NX, NY, T, Opts...>();
    const auto       spec = target_spec_of(option<kind::Disturbance>(NoDisturbance{}, opts...), option<kind::TrackedOutputs>(NoDisturbance{}, opts...));
    const auto       measured = measured_of<NX, T>(opts...);
    constexpr size_t NM = measured_inputs<damp::remove_cvref_t<decltype(measured)>>;
    static_assert(
        !damp::is_same_v<damp::remove_cvref_t<decltype(spec)>, NoDisturbance> || NM > 0,
        "steady_state_target holds the outputs a TrackedOutputs or a disturbance spec lists: TrackedOutputs{} for every output, "
        "TrackedOutputs{{i, ...}} for some. With neither it is a regulator's, which holds nothing at r: give it a MeasuredDisturbance "
        "for the steady state under d."
    );
    if constexpr (!kModelFits) {
        static_cast<void>(sys);
        static_cast<void>(Ts);
        static_cast<void>(delay);
        return SteadyStateTargetResult<NX, NU, 0, 0, T>{}; // the message says why: nothing further
    } else {
        return held_target_under<D>(sys, Ts, spec, delay, measured);
    }
}

} // namespace detail

/**
 * @brief The steady-state target of a law on a discrete plant: the maps from the estimated
 *        disturbance and the reference to [x_t; u_t]
 *
 * Solves, for unit d and r,
 * [A − I, B; C_z, D_z]·[x_t; u_t] = [−B_d d; r − C_{d,z} d] on the filter's model of
 * design::kalman with the same options: C_z and D_z are the rows of the measured outputs the
 * laws hold at r, those the disturbance spec lists, else those a TrackedOutputs lists, all of
 * them when the one given lists none; B_d and C_d the spec's channels, C_{d,z} the same rows of
 * C_d (none without a spec: the target from r alone). With neither a spec nor a TrackedOutputs
 * it is a regulator's target, as design::state_feedback and design::mpc then are: no map of r,
 * and under a MeasuredDisturbance the steady state that holds every output at 0. Without a
 * MeasuredDisturbance too it would hold nothing, and does not compile: TrackedOutputs{} holds
 * every output at r. With NZ = NU the system is square and
 * solved; with NZ < NU, the minimum-norm solution. With NZ > NU no steady state holds every r,
 * and the design does not compile: a plant with more outputs than inputs lists the ones r
 * holds, at most NU, in a TrackedOutputs or the spec's .outputs.
 * Offset-free control needs the steady state to exist and the augment detectable
 * (design::kalman checks it).
 *
 * A MeasuredDisturbance{.Bd} among @p opts, a disturbance read each tick, adds its own B_d to
 * the state rows' right side (none to the outputs'), its NM columns after the spec's: target_d
 * is then [x_t; u_t] per unit [d̂; d], the target design::state_feedback and design::mpc hold
 * with the same options, and a design without one is the type it was.
 *
 * Refused with DesignStatus::singular when NZ = NU and the steady state is not unique, or
 * NZ < NU and no steady state holds every r (the system is rank deficient); with
 * DesignStatus::bad_weights for a controlled output out of range or listed twice, or a
 * non-finite entry of the plant, B_d or C_d, or of a MeasuredDisturbance's B_d. An option of
 * a kind given twice, a TrackedOutputs beside a disturbance spec, or more held outputs than
 * inputs does not compile, as in design::state_feedback and design::mpc.
 *
 * @see "Offset-free reference tracking with model predictive control" (Maeder, Borrelli &
 *      Morari, 2009), Automatica 45(9)
 * @see "Disturbance modeling for offset-free linear model predictive control" (Muske &
 *      Badgwell, 2002), Journal of Process Control 12(5)
 *
 * @param sys   Discrete-time plant (Ts > 0)
 * @param opts  The model features of the design, in any order: a disturbance spec or a
 *              TrackedOutputs, a MeasuredDisturbance (one of the three at least; a KalmanNoise is
 *              accepted and ignored)
 * @return SteadyStateTargetResult with the maps and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto steady_state_target(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Opts&... opts) {
    return detail::target_with<detail::Domain::plant>(sys, T{0}, detail::NoDelay{}, opts...);
}

/**
 * @brief The steady-state target of a law on a continuous plant sampled at @p Ts (see
 *        steady_state_target(sys, opts…)): on the discrete filter's model of
 *        design::kalman(sys, Ts, noise, opts…).
 *
 * With InputDelay<N>{θ} among @p opts the model is the plant whose input arrives θ late, the
 * plant then its delay states v(k−1) … v(k−N); at a steady state each holds u_t. With
 * OutputDelay<M>{θ}, the model of design::kalman(sys, Ts, noise, …, OutputDelay<M>{θ}), its
 * readings late: the target is on that model, so each of the readings' delay states holds the
 * steady output, the state a law with that filter anchors on.
 *
 * A MeasuredDisturbance's B_d is continuous, ZOH'd with the plant (d held over the period). It
 * enters the plant's states, none of an InputDelay's delay states, and an OutputDelay's as the
 * plant moved under it until each reading was taken: the channel design::kalman's filter takes.
 *
 * @param sys   Continuous plant (sys.Ts = 0)
 * @param Ts    Sampling time of the law [s]
 * @param opts  The model features of the design, in any order: a disturbance spec or a
 *              TrackedOutputs, a MeasuredDisturbance, an InputDelay, an OutputDelay
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto steady_state_target(const StateSpace<NX, NU, NY, T, NW, NV>& sys, damp::type_identity_t<T> Ts, const Opts&... opts) {
    return detail::target_with<detail::Domain::sampled>(sys, Ts, detail::delays_of(opts...), opts...);
}

} // namespace design

/** @brief A state and the input at it: a steady state (SteadyStateTarget), or a point of a reference trajectory (OutputInverse). */
template<size_t NX, size_t NU, typename T>
struct TargetPoint {
    ColVec<NX, T> x{}; ///< x_t
    ColVec<NU, T> u{}; ///< u_t
};

/**
 * @brief The steady-state target at run time: (x_t, u_t) from the reference and the estimated
 *        disturbance, two matrix–vector products, and the re-anchoring of a FixedGainObserver on it.
 *
 * Shared by the laws that track through a target (OFLQG, MPC, OFMPC): equilibrium() for a law
 * on the absolute state, anchor() for one on the estimator's deviation (ND > 0).
 *
 * @tparam NX States of the filter's model before the disturbance; ND Disturbance states;
 *            NZ Controlled outputs (0: a regulator's target, TargetR then any empty type)
 */
template<size_t NX, size_t NU, size_t ND, size_t NZ, typename T = float, typename TargetD = Matrix<NX + NU, ND, T>, typename TargetR = Matrix<NX + NU, NZ, T>>
struct SteadyStateTarget {
    [[no_unique_address]] TargetD target_d{}; ///< [x_t; u_t] per unit d̂
    [[no_unique_address]] TargetR target_r{}; ///< [x_t; u_t] per unit r

    constexpr SteadyStateTarget() = default;

    /** @brief From design::steady_state_target (any scalar). A failed design does not compile when constant-evaluated. */
    template<typename U>
        requires(damp::is_same_v<TargetD, Matrix<NX + NU, ND, T>> && damp::is_same_v<TargetR, Matrix<NX + NU, NZ, T>>)
    constexpr SteadyStateTarget(const design::SteadyStateTargetResult<NX, NU, ND, NZ, U>& result) // NOLINT
        : target_d(result.target_d.template as<T>()), target_r(result.target_r.template as<T>()) {
        if (!result.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /** @brief The steady state (x_t, u_t) that holds z = @p r under the disturbance @p d (a regulator's, NZ = 0: under d alone). */
    [[nodiscard]] constexpr TargetPoint<NX, NU, T> equilibrium(const ColVec<NZ, T>& r, const ColVec<ND, T>& d = ColVec<ND, T>{}) const {
        ColVec<NX + NU, T> s{};
        if constexpr (ND > 0 && NZ > 0) {
            s = (target_d * d) + (target_r * r);
        } else if constexpr (ND > 0) {
            static_cast<void>(r);
            s = target_d * d; // no map of r: TargetR may hold nothing
        } else if constexpr (NZ > 0) {
            static_cast<void>(d);
            s = target_r * r;
        } else {
            static_cast<void>(r);
            static_cast<void>(d);
        }
        TargetPoint<NX, NU, T> p{};
        for (size_t i = 0; i < NX; ++i) {
            p.x(i) = s(i);
        }
        for (size_t i = 0; i < NU; ++i) {
            p.u(i) = s(NX + i);
        }
        return p;
    }

    /**
     * @brief Move @p estimator's anchor to the target for @p r about the reference point
     *        (@p x_ref, @p u_ref) under its d̂: x_a = x_t + x*, u_a = u_t + u*.
     */
    template<typename Estimator>
        requires(ND > 0)
    constexpr void anchor(
        Estimator&           estimator,
        const ColVec<NZ, T>& r,
        const ColVec<NX, T>& x_ref = ColVec<NX, T>{},
        const ColVec<NU, T>& u_ref = ColVec<NU, T>{}
    ) const {
        const auto p = equilibrium(r, estimator.disturbance_estimate());
        estimator.anchor(ColVec<NX, T>(p.x + x_ref), ColVec<NU, T>(p.u + u_ref), x_ref, u_ref);
    }

    /**
     * @brief anchor(estimator, r, x_ref, u_ref) under the measured disturbance @p d too: the target
     *        for @p r under [d̂; d], its last NM columns d's, and the estimator's anchor input [u_a; d],
     *        its measured inputs after u (design::kalman with a MeasuredDisturbance), so its next
     *        predict moves the deviation by B_d·(d⁺ − d) as d changes.
     */
    template<typename Estimator, size_t NM>
        requires(NM > 0 && NM < ND)
    constexpr void anchor(Estimator& estimator, const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d) const {
        const auto               d_hat = estimator.disturbance_estimate();
        const auto               p = equilibrium(r, ::damp::detail::stacked(d_hat, d));
        const ColVec<NU + NM, T> u_a = ::damp::detail::stacked(ColVec<NU, T>(p.u + u_ref), d);
        estimator.anchor(ColVec<NX, T>(p.x + x_ref), u_a, x_ref, ::damp::detail::stacked(u_ref, ColVec<NM, T>{}));
    }
};

/** @brief Deduce the target from its design, keeping the result's scalar: SteadyStateTarget target{res.as<float>()}. */
template<size_t NX, size_t NU, size_t ND, size_t NZ, typename T>
SteadyStateTarget(const design::SteadyStateTargetResult<NX, NU, ND, NZ, T>&) -> SteadyStateTarget<NX, NU, ND, NZ, T>;

} // namespace damp
