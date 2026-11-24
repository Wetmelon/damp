// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file disturbance.hpp
 * @brief The disturbance a design models: NoDisturbance, InputDisturbance, OutputDisturbance, PlantDisturbance or DisturbanceModel, and the channels they enter through; and a MeasuredDisturbance, one it is given
 *
 * A disturbance spec is a model feature, recognized by type among design::steady_state_target's
 * options, and the disturbance option of design::kalman, design::state_feedback and design::mpc with a KalmanNoise: ND integrating
 * disturbances d_{k+1} = d_k enter the filter's model through B_d and the readings through C_d.
 * A spec's .outputs are the ones the target holds at r, every measured output when it lists none,
 * and at most one per input: a longer list, or none on a plant with more outputs than inputs,
 * does not compile. On readings that arrive late (OutputDelay) d is a state of the plant held
 * over the period, so each reading holds d as it was when the reading was taken (held_channel);
 * an InputDisturbance's adds to u, which the model carries there. A MeasuredDisturbance is the other kind, an input the
 * law reads each tick: the option of design::state_feedback, design::mpc, design::kalman and
 * design::steady_state_target, alone or beside a disturbance spec, whose d̂ then estimates only
 * what it does not measure.
 */


#include "damp/backend.hpp"
#include "damp/controllers/detail/law_detail.hpp"
#include "damp/design/options.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @brief No disturbance model (the default): the plain filter on x.
 */
struct NoDisturbance {};

/**
 * @brief An input disturbance: d adds to u, B_d = B and C_d = D, one disturbance per input.
 *
 * The usual offset-free model for a load or an actuator bias (Pannocchia & Rawlings, 2003).
 *
 * @ref outputs is an array, so `InputDisturbance{{0, 2}}` takes NZ from its list; `InputDisturbance{}` holds
 * every measured output.
 * @tparam NZ Outputs held offset-free, listed in @ref outputs (0: every measured output)
 */
template<size_t NZ = 0>
struct InputDisturbance {
    size_t outputs[NZ]{}; ///< Indices of the measured outputs held offset-free, each once (an array, so bare braces give NZ)
};

/** @brief An input disturbance, every measured output held offset-free: no list. */
template<>
struct InputDisturbance<0> {};

/**
 * @brief An output disturbance: a bias on each reading, B_d = 0 and C_d = I, one per output.
 *
 * The model of MATLAB®'s default MPC estimator (a step on each measured output).
 *
 * @ref outputs is an array, so `OutputDisturbance{{0, 2}}` takes NZ from its list; `OutputDisturbance{}` holds
 * every measured output.
 * @tparam NZ Outputs held offset-free, listed in @ref outputs (0: every measured output)
 */
template<size_t NZ = 0>
struct OutputDisturbance {
    size_t outputs[NZ]{}; ///< Indices of the measured outputs held offset-free, each once (an array, so bare braces give NZ)
};

/** @brief An output disturbance, every measured output held offset-free: no list. */
template<>
struct OutputDisturbance<0> {};

/**
 * @brief The plant's own disturbance channel: B_d = G and C_d = H, one disturbance per column of
 *        G (H needs as many columns).
 *
 * For a plant whose G and H describe how an unmeasured force enters the dynamics and the
 * sensors (a force on the carriage that an accelerometer also reads). The process and
 * measurement noise still enter through G and H.
 *
 * @ref outputs is an array, so `PlantDisturbance{{0, 2}}` takes NZ from its list; `PlantDisturbance{}` holds
 * every measured output.
 * @tparam NZ Outputs held offset-free, listed in @ref outputs (0: every measured output)
 */
template<size_t NZ = 0>
struct PlantDisturbance {
    size_t outputs[NZ]{}; ///< Indices of the measured outputs held offset-free, each once (an array, so bare braces give NZ)
};

/** @brief The plant's own disturbance channel, every measured output held offset-free: no list. */
template<>
struct PlantDisturbance<0> {};

/**
 * @brief An explicit disturbance model: d enters the plant's state through B_d and the readings
 *        through C_d.
 *
 * B_d is in the plant's time domain (continuous on a continuous plant, discretized with it by
 * ZOH; per step on a discrete one), one row per plant state. On readings that arrive late
 * (OutputDelay) each holds what the plant moved under d until it was taken, and C_d·d as it was
 * then: the filter's model is the delayed discretization of [A, B_d; 0, 0] read through [C, C_d].
 *
 * @tparam NX Plant states; NY Measured outputs (0: d reaches no reading, C_d = 0); ND Disturbance states
 *
 * @ref outputs is an array, so `DisturbanceModel{.Bd = Bd, .Cd = Cd, .outputs = {0}}` takes NZ from
 * its list, and NX, NY and ND from the matrices; left out, every measured output is held. With
 * @ref Cd left out, d does not reach the readings: `DisturbanceModel{.Bd = sys.B}` is
 * DisturbanceModel<NX, 0, ND>, and the design reads C_d as zero on the plant's outputs.
 * @tparam NZ Outputs held offset-free, listed in @ref outputs (0: every measured output)
 */
template<size_t NX, size_t NY = 0, size_t ND = 1, size_t NZ = 0, typename T = double>
struct DisturbanceModel {
    Matrix<NX, ND, T> Bd{};          ///< How d enters the state update
    Matrix<NY, ND, T> Cd{};          ///< How d enters the readings directly
    size_t            outputs[NZ]{}; ///< Indices of the measured outputs held offset-free, each once (an array, so bare braces give NZ)
};

/** @brief An explicit disturbance model holding every measured output offset-free: no list. */
template<size_t NX, size_t NY, size_t ND, typename T>
struct DisturbanceModel<NX, NY, ND, 0, T> {
    Matrix<NX, ND, T> Bd{}; ///< How d enters the state update
    Matrix<NY, ND, T> Cd{}; ///< How d enters the readings directly
};

/**
 * @brief A measured disturbance: a known input the controller does not command, entering the
 *        plant's state through @c Bd (a load torque, a feed rate, the ambient temperature), read
 *        each tick.
 *
 * A trailing option of design::state_feedback, design::mpc, design::kalman and
 * design::steady_state_target. The law's target is the steady state under it, steady_state_target's
 * with the same options (target_d per unit [d̂; d]): u = u_t(r, d) − K·(x − x_t(r, d)) on an LQ law
 * and the QP posed about that target on an MPC, which also holds d over its horizon; the runtime
 * takes it each tick, after the reference and its feedforward: control(r, d),
 * control(r, x_ref, u_ref, d), step(r, d, m), step(r, x_ref, u_ref, d, m). A filter (a
 * KalmanNoise) takes it as a known input after u, its model's B being [B, B_d] and its noise
 * unchanged, so beside a disturbance spec its d̂ estimates only the part of a load that is not
 * measured.
 *
 * A regulator (no r: no disturbance spec, TrackedOutputs or IntegralAction) holds every output at
 * 0 under d, where a steady state can. With as many outputs as inputs and no unique steady state
 * it is refused (DesignStatus::singular). With more outputs than inputs no steady state holds
 * them all (a motor law that reads both speed and current under a load torque), and the design
 * does not compile: list the ones it holds, at most one per input, `TrackedOutputs{{i, …}}`, and
 * pass r = 0.
 *
 * B_d is in the plant's time domain: continuous on a continuous plant, discretized with it by ZOH
 * (d held over the period), per step on a discrete one. On a plant whose input arrives late
 * (InputDelay) d enters the plant, not the delay states: it is not a command. On readings that
 * arrive late (OutputDelay) it reaches them as the plant moved under it until each was taken.
 * Unlike a disturbance spec it is not rescaled and does not reach the readings directly.
 *
 * `MeasuredDisturbance{Bd}` and `MeasuredDisturbance{.Bd = Bd}` take NX, NM and the scalar from
 * the matrix.
 *
 * @note Compare with MATLAB®'s measured-disturbance channels (setmpcsignals 'MD').
 * @tparam NX Plant states; NM Measured-disturbance inputs
 */
template<size_t NX, size_t NM, typename T = double>
struct MeasuredDisturbance {
    Matrix<NX, NM, T> Bd{}; ///< How the measured disturbance enters the plant's state update
};

/** @brief Sizes from the matrix: `MeasuredDisturbance{Bd}`. */
template<size_t NX, size_t NM, typename T>
MeasuredDisturbance(Matrix<NX, NM, T>) -> MeasuredDisturbance<NX, NM, T>;

namespace detail {

/** @brief Not a measured disturbance. */
template<typename>
struct is_measured_disturbance : damp::false_type {};

/** @brief MeasuredDisturbance, as a kind of trailing option (options.hpp). */
template<size_t NX, size_t NM, typename T>
struct is_measured_disturbance<MeasuredDisturbance<NX, NM, T>> : damp::true_type {};

/** @brief The inputs a measured disturbance carries: NM, none for anything else (Unset, a call without one). */
template<typename M>
inline constexpr size_t measured_inputs = 0;

/** @brief MeasuredDisturbance<NX, NM>: NM inputs. */
template<size_t NX, size_t NM, typename T>
inline constexpr size_t measured_inputs<MeasuredDisturbance<NX, NM, T>> = NM;

/** @brief The measured disturbance among a call's options, on a plant of NX states in T: none (NM = 0) without one. */
template<size_t NX, typename T, typename... Opts>
[[nodiscard]] constexpr auto measured_of(const Opts&... opts) {
    static_assert(
        ::damp::detail::option_count<is_measured_disturbance, Opts...> <= 1,
        "One MeasuredDisturbance carries every measured input: MeasuredDisturbance{.Bd = Bd}, a column of Bd per input."
    );
    using M = ::damp::detail::option_t<is_measured_disturbance, MeasuredDisturbance<NX, 0, T>, Opts...>;
    static_assert(
        damp::is_same_v<M, MeasuredDisturbance<NX, measured_inputs<M>, T>>,
        "The measured disturbance enters the plant's states: MeasuredDisturbance{.Bd = Matrix<NX, NM>{...}}, a row per plant state, in the plant's scalar."
    );
    return ::damp::detail::option<is_measured_disturbance>(MeasuredDisturbance<NX, 0, T>{}, opts...);
}

/**
 * @brief Where a disturbance enters: @c Bd on the filter's model before the disturbance (NXE
 *        states), @c Cd on the readings, and @c Bd_c on the continuous plant (the Van Loan noise
 *        of a continuous design).
 */
template<size_t NX, size_t NXE, size_t NY, size_t ND, typename T>
struct DisturbanceChannels {
    Matrix<NXE, ND, T> Bd{};   ///< Into the discrete state update of the filter's model
    Matrix<NY, ND, T>  Cd{};   ///< Into the readings
    Matrix<NX, ND, T>  Bd_c{}; ///< Into the continuous plant (Van Loan of the random walk)
};

/**
 * @brief Where a disturbance d held over the period enters the model a design is on, d entering
 *        the plant's states through @p Bd and the readings through @p Cd: B_d ZOH'd with a
 *        continuous plant at @p Ts, per step on a discrete one (@p Discrete), on the plant's
 *        states; none on an input's delay states; and on a late reading's delay states what the
 *        plant moved under d and what C_d read of it when the reading was taken.
 *
 * The model is the delayed discretization of the plant with d as a state of no dynamics,
 * [A, B_d; 0, 0] read through [C, C_d]; its d columns are those of
 * discretize(StateSpace{A, B_d, C, C_d}, Ts, OutputDelay) on the delay states: a reading taken s
 * into the period loads (C·∫₀ˢ e^{Aη} dη·B_d + C_d)·d, whatever the input's delay, so a model
 * with both delays reads exactly the same d, and @c Cd keeps C_d on an output read at once, none
 * on a late one. Without an OutputDelay @c Cd is @p Cd.
 * @return The channels, @c Bd on the model's NX + delay_states rows and @c Bd_c = @p Bd; nullopt
 *         on a non-finite entry
 */
template<bool Discrete, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND, typename Delay>
[[nodiscard]] constexpr damp::optional<DisturbanceChannels<NX, NX + delay_states<Delay, NU, NY>, NY, ND, T>> held_channel(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const Matrix<NX, ND, T>&                 Bd,
    const Matrix<NY, ND, T>&                 Cd,
    const Delay&                             delay
) {
    constexpr size_t                                                     NRY = output_delay_states<Delay>;
    DisturbanceChannels<NX, NX + delay_states<Delay, NU, NY>, NY, ND, T> out{.Cd = Cd, .Bd_c = Bd};
    const auto                                                           finite = [](T v) { return damp::isfinite(v); };
    if constexpr (ND == 0) {
        static_cast<void>(sys);
        static_cast<void>(Ts);
        static_cast<void>(delay);
        static_cast<void>(finite);
        return out; // none: no exponential to take
    } else if (!mat::all_of(Bd, finite) || !mat::all_of(Cd, finite)) {
        return damp::nullopt; // the exponential would carry it; refused before
    } else if constexpr (Discrete) {
        static_cast<void>(Ts);
        static_cast<void>(delay);
        out.Bd.template block<NX, ND>(0, 0) = Bd;
    } else if constexpr (NRY == 0) {
        static_cast<void>(delay);
        const auto d = ::damp::discretize(StateSpace<NX, ND, 1, T>{.A = sys.A, .B = Bd}, Ts);
        if (!d) {
            return damp::nullopt;
        }
        out.Bd.template block<NX, ND>(0, 0) = d->B;
    } else {
        const auto d = ::damp::discretize(StateSpace<NX, ND, NY, T>{.A = sys.A, .B = Bd, .C = sys.C, .D = Cd}, Ts, delay.output);
        if (!d) {
            return damp::nullopt;
        }
        constexpr size_t first = NX + (input_delay_states<Delay> * NU); // the late readings' delay states, after the input's
        out.Bd.template block<NX, ND>(0, 0) = d->B.template block<NX, ND>(0, 0);
        out.Bd.template block<NRY * NY, ND>(first, 0) = d->B.template block<NRY * NY, ND>(NX, 0);
        out.Cd = d->D;
    }
    return out;
}

/**
 * @brief Where a measured disturbance enters the model a design is on (held_channel, C_d = 0): its
 *        B_d on the plant's states, ZOH'd with a continuous plant at @p Ts (d held over the
 *        period), per step on a discrete one (@p Discrete); none on an input's delay states; and
 *        on a late reading's delay states what the plant moved under it until the reading was
 *        taken, C·∫₀ˢ e^{Aη} dη·B_d·d.
 * @return The channel on the model's NX + delay_states rows; nullopt on a non-finite entry
 */
template<bool Discrete, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NM, typename Delay>
[[nodiscard]] constexpr damp::optional<Matrix<NX + delay_states<Delay, NU, NY>, NM, T>> measured_channel(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const MeasuredDisturbance<NX, NM, T>&    measured,
    const Delay&                             delay
) {
    const auto channel = held_channel<Discrete>(sys, Ts, measured.Bd, Matrix<NY, NM, T>{}, delay);
    if (!channel) {
        return damp::nullopt;
    }
    return channel->Bd;
}

/** @brief @p left then @p right, side by side: the channels of an estimated disturbance, then a measured one's. */
template<size_t R, size_t CL, size_t CR, typename T>
[[nodiscard]] constexpr Matrix<R, CL + CR, T> side_by_side(const Matrix<R, CL, T>& left, const Matrix<R, CR, T>& right) {
    Matrix<R, CL + CR, T> out{};
    if constexpr (CL > 0) {
        out.template block<R, CL>(0, 0) = left;
    }
    if constexpr (CR > 0) {
        out.template block<R, CR>(0, CL) = right;
    }
    return out;
}

} // namespace detail

namespace detail {

/** @brief Not a disturbance spec. */
template<typename>
struct is_disturbance : damp::false_type {};

/** @brief NoDisturbance. */
template<>
struct is_disturbance<NoDisturbance> : damp::true_type {};

/** @brief InputDisturbance. */
template<size_t NZ>
struct is_disturbance<InputDisturbance<NZ>> : damp::true_type {};

/** @brief OutputDisturbance. */
template<size_t NZ>
struct is_disturbance<OutputDisturbance<NZ>> : damp::true_type {};

/** @brief PlantDisturbance. */
template<size_t NZ>
struct is_disturbance<PlantDisturbance<NZ>> : damp::true_type {};

/** @brief DisturbanceModel. */
template<size_t NX, size_t NY, size_t ND, size_t NZ, typename T>
struct is_disturbance<DisturbanceModel<NX, NY, ND, NZ, T>> : damp::true_type {};

/** @brief A disturbance model with states: every disturbance spec but NoDisturbance. */
template<typename Spec>
inline constexpr bool is_disturbance_model = is_disturbance<Spec>::value && !damp::is_same_v<Spec, NoDisturbance>;

/**
 * @brief The sizes of a disturbance spec on a plant with NU inputs, NY outputs and NW columns of
 *        G: @c states (ND) and @c listed (the NZ of its outputs list, 0 for all).
 */
template<typename Spec, size_t NU, size_t NY, size_t NW>
struct disturbance_traits;

/** @brief No disturbance: no states. */
template<size_t NU, size_t NY, size_t NW>
struct disturbance_traits<NoDisturbance, NU, NY, NW> {
    static constexpr size_t states = 0; ///< ND
    static constexpr size_t listed = 0; ///< NZ listed
};

/** @brief One disturbance per input. */
template<size_t NZ, size_t NU, size_t NY, size_t NW>
struct disturbance_traits<InputDisturbance<NZ>, NU, NY, NW> {
    static constexpr size_t states = NU; ///< ND
    static constexpr size_t listed = NZ; ///< NZ listed
};

/** @brief One disturbance per output. */
template<size_t NZ, size_t NU, size_t NY, size_t NW>
struct disturbance_traits<OutputDisturbance<NZ>, NU, NY, NW> {
    static constexpr size_t states = NY; ///< ND
    static constexpr size_t listed = NZ; ///< NZ listed
};

/** @brief One disturbance per column of G. */
template<size_t NZ, size_t NU, size_t NY, size_t NW>
struct disturbance_traits<PlantDisturbance<NZ>, NU, NY, NW> {
    static constexpr size_t states = NW; ///< ND
    static constexpr size_t listed = NZ; ///< NZ listed
};

/** @brief The model's own count. */
template<size_t NXD, size_t NYD, size_t ND, size_t NZ, typename T, size_t NU, size_t NY, size_t NW>
struct disturbance_traits<DisturbanceModel<NXD, NYD, ND, NZ, T>, NU, NY, NW> {
    static constexpr size_t states = ND; ///< ND
    static constexpr size_t listed = NZ; ///< NZ listed
};

/** @brief The disturbance spec among a call's options (NoDisturbance when none is given). */
template<typename... Opts>
using disturbance_option_t = ::damp::detail::option_t<is_disturbance, NoDisturbance, Opts...>;

/** @brief ND of the disturbance spec among a call's options. */
template<size_t NU, size_t NY, size_t NW, typename... Opts>
inline constexpr size_t disturbance_states = disturbance_traits<disturbance_option_t<Opts...>, NU, NY, NW>::states;
} // namespace detail

namespace detail {

/**
 * @brief The listed measured outputs, checked (in range, each once), or every output when none
 *        is listed.
 */
template<size_t NY, size_t NS>
[[nodiscard]] constexpr damp::optional<damp::array<size_t, (NS == 0) ? NY : NS>> selected_outputs(const damp::array<size_t, NS>& listed) {
    constexpr size_t        NZ = (NS == 0) ? NY : NS;
    damp::array<size_t, NZ> rows{};
    for (size_t i = 0; i < NZ; ++i) {
        if constexpr (NS == 0) {
            rows[i] = i;
        } else {
            rows[i] = listed[i];
        }
        if (rows[i] >= NY) {
            return damp::nullopt;
        }
        for (size_t j = 0; j < i; ++j) {
            if (rows[j] == rows[i]) {
                return damp::nullopt;
            }
        }
    }
    return rows;
}

/**
 * @brief An input disturbance: B_d = B of the filter's model (delayed with u), C_d = its D. On a
 *        model with input delay states its walk reaches the plant through them, so none of it acts
 *        within the period (B_d,c = 0).
 */
template<size_t NZ, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXE, size_t NWE, size_t NVE>
[[nodiscard]] constexpr DisturbanceChannels<NX, NXE, NY, NU, T> disturbance_channels(
    const InputDisturbance<NZ>&,
    const StateSpace<NX, NU, NY, T, NW, NV>&    plant,
    const StateSpace<NXE, NU, NY, T, NWE, NVE>& model,
    T
) {
    if constexpr (NXE == NX) {
        return DisturbanceChannels<NX, NXE, NY, NU, T>{.Bd = model.B, .Cd = model.D, .Bd_c = plant.B};
    } else {
        return DisturbanceChannels<NX, NXE, NY, NU, T>{.Bd = model.B, .Cd = model.D};
    }
}

/** @brief An output disturbance: a bias on each reading, B_d = 0, C_d = I. */
template<size_t NZ, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXE, size_t NWE, size_t NVE>
[[nodiscard]] constexpr DisturbanceChannels<NX, NXE, NY, NY, T> disturbance_channels(
    const OutputDisturbance<NZ>&,
    const StateSpace<NX, NU, NY, T, NW, NV>&,
    const StateSpace<NXE, NU, NY, T, NWE, NVE>&,
    T
) {
    return DisturbanceChannels<NX, NXE, NY, NY, T>{.Cd = Matrix<NY, NY, T>::identity()};
}

/** @brief The plant's channel: B_d = G of the filter's model (off the delay states), C_d = H. */
template<size_t NZ, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXE, size_t NWE, size_t NVE>
[[nodiscard]] constexpr DisturbanceChannels<NX, NXE, NY, NW, T> disturbance_channels(
    const PlantDisturbance<NZ>&,
    const StateSpace<NX, NU, NY, T, NW, NV>&    plant,
    const StateSpace<NXE, NU, NY, T, NWE, NVE>& model,
    T
) {
    static_assert(NV == NW && NWE == NW, "PlantDisturbance takes C_d = H: H needs one column per column of G.");
    return DisturbanceChannels<NX, NXE, NY, NW, T>{.Bd = model.G, .Cd = plant.H, .Bd_c = plant.G};
}

/** @brief An explicit model: B_d on the plant's states (ZOH'd on a continuous plant), C_d as given (zero left out). */
template<size_t NYD, size_t ND, size_t NZ, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXE, size_t NWE, size_t NVE>
[[nodiscard]] constexpr DisturbanceChannels<NX, NXE, NY, ND, T> disturbance_channels(
    const DisturbanceModel<NX, NYD, ND, NZ, T>& spec,
    const StateSpace<NX, NU, NY, T, NW, NV>&    plant,
    const StateSpace<NXE, NU, NY, T, NWE, NVE>&,
    T Ts
) {
    static_assert(NYD == 0 || NYD == NY, "The DisturbanceModel's C_d is on the plant's outputs, or left out.");
    DisturbanceChannels<NX, NXE, NY, ND, T> out{.Bd_c = spec.Bd};
    if constexpr (NYD > 0) {
        out.Cd = spec.Cd;
    }
    if (plant.Ts == T{0}) {
        // plant.A has already been sampled, so the exponential fails only on a non-finite B_d; it is
        // carried into Bd, which design::kalman refuses (not_detectable) and steady_state_target
        // refuses (invalid_parameter).
        const auto d = ::damp::detail::discretize_zoh_impl(StateSpace<NX, ND, 1, T>{.A = plant.A, .B = spec.Bd}, Ts);
        out.Bd.template block<NX, ND>(0, 0) = d ? d->B : spec.Bd;
    } else {
        out.Bd.template block<NX, ND>(0, 0) = spec.Bd;
    }
    return out;
}

/** @brief Not an InputDisturbance. */
template<typename>
inline constexpr bool is_input_disturbance = false;

/** @brief InputDisturbance: d adds to u. */
template<size_t NZ>
inline constexpr bool is_input_disturbance<InputDisturbance<NZ>> = true;

/**
 * @brief Whether a spec's d reaches a late reading as it was when the reading was taken, a state of
 *        the plant held over the period (held_channel): with an OutputDelay in @p Delay, every spec
 *        with states but InputDisturbance, whose d adds to u and so rides the model's own B and D.
 */
template<typename Spec, typename Delay>
inline constexpr bool held_on_late_readings = is_disturbance_model<Spec> && !is_input_disturbance<Spec> && (output_delay_states<Delay> > 0);

/**
 * @brief A disturbance's channels @p ch (disturbance_channels) on the filter's model with the
 *        delays @p delay: on readings that arrive late, the held_channel of its continuous B_d and
 *        its C_d, d reaching each reading as it was when the reading was taken. Without an
 *        OutputDelay, and for an InputDisturbance, @p ch itself.
 *
 * A non-finite B_d or C_d keeps @p ch, which carries it to the design's refusal; @p sys.A has
 * already been sampled, so the exponential fails on nothing else.
 */
template<typename Spec, size_t NX, size_t NXE, size_t NY, size_t ND, typename T, size_t NU, size_t NW, size_t NV, typename Delay>
[[nodiscard]] constexpr DisturbanceChannels<NX, NXE, NY, ND, T> on_late_readings(
    const Spec&                                    spec,
    const DisturbanceChannels<NX, NXE, NY, ND, T>& ch,
    const StateSpace<NX, NU, NY, T, NW, NV>&       sys,
    T                                              Ts,
    const Delay&                                   delay
) {
    static_cast<void>(spec);
    if constexpr (held_on_late_readings<Spec, Delay>) {
        static_assert(NXE == NX + delay_states<Delay, NU, NY>, "The channels are on the filter's model, with its delay states.");
        const auto held = held_channel<false>(sys, Ts, ch.Bd_c, ch.Cd, delay); // a late reading needs the continuous plant
        return held ? *held : ch;
    } else {
        static_cast<void>(sys);
        static_cast<void>(Ts);
        static_cast<void>(delay);
        return ch;
    }
}

/**
 * @brief Whether the disturbance augment is detectable from the readings: rank
 *        [A − I, B_d; C, C_d] = NX + ND (with (A, C) detectable, which the filter DARE checks).
 *
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings,
 *      2003), AIChE Journal 49(2)
 */
template<size_t NX, size_t NY, size_t ND, typename T>
[[nodiscard]] constexpr bool augment_detectable(const Matrix<NX, NX, T>& A, const Matrix<NY, NX, T>& C, const Matrix<NX, ND, T>& Bd, const Matrix<NY, ND, T>& Cd) {
    Matrix<NX + NY, NX + ND, T> test{};
    test.template block<NX, NX>(0, 0) = A - Matrix<NX, NX, T>::identity();
    test.template block<NX, ND>(0, NX) = Bd;
    test.template block<NY, NX>(NX, 0) = C;
    test.template block<NY, ND>(NX, NX) = Cd;
    return mat::rank(test) == NX + ND;
}

/** @brief The filter's model with the disturbance states last: [A B_d; 0 I], [B; 0], [C C_d], D. */
template<size_t NWK, size_t NXE, size_t NU, size_t NY, size_t ND, typename T, size_t NWE, size_t NVE>
[[nodiscard]] constexpr StateSpace<NXE + ND, NU, NY, T, NWK, NVE> disturbance_augment(
    const StateSpace<NXE, NU, NY, T, NWE, NVE>& model,
    const Matrix<NXE, ND, T>&                   Bd,
    const Matrix<NY, ND, T>&                    Cd,
    T                                           Ts
) {
    StateSpace<NXE + ND, NU, NY, T, NWK, NVE> aug{.D = model.D, .H = model.H, .Ts = Ts};
    aug.A.template block<NXE, NXE>(0, 0) = model.A;
    aug.A.template block<NXE, ND>(0, NXE) = Bd;
    aug.A.template block<ND, ND>(NXE, NXE) = Matrix<ND, ND, T>::identity();
    aug.B.template block<NXE, NU>(0, 0) = model.B;
    aug.C.template block<NY, NXE>(0, 0) = model.C;
    aug.C.template block<NY, ND>(0, NXE) = Cd;
    return aug;
}

} // namespace detail

} // namespace design

} // namespace damp
