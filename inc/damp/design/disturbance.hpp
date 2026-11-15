// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file disturbance.hpp
 * @brief The disturbance a design models: NoDisturbance, InputDisturbance, OutputDisturbance, PlantDisturbance or DisturbanceModel, and the channels they enter through
 *
 * A disturbance spec is a model feature, recognized by type among design::steady_state_target's
 * options, and the disturbance argument of design::kalman, design::oflqg and design::ofmpc: ND integrating
 * disturbances d_{k+1} = d_k enter the filter's model through B_d and the readings through C_d.
 */


#include "damp/backend.hpp"
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
 * @tparam NZ Outputs held offset-free, listed in @ref outputs (0: every measured output)
 */
template<size_t NZ = 0>
struct InputDisturbance {
    damp::array<size_t, NZ> outputs{}; ///< Indices of the measured outputs held offset-free, each once (none listed: all)
};

/**
 * @brief An output disturbance: a bias on each reading, B_d = 0 and C_d = I, one per output.
 *
 * The model of MATLAB®'s default MPC estimator (a step on each measured output).
 * @tparam NZ Outputs held offset-free, listed in @ref outputs (0: every measured output)
 */
template<size_t NZ = 0>
struct OutputDisturbance {
    damp::array<size_t, NZ> outputs{}; ///< Indices of the measured outputs held offset-free, each once (none listed: all)
};

/**
 * @brief The plant's own disturbance channel: B_d = G and C_d = H, one disturbance per column of
 *        G (H needs as many columns).
 *
 * For a plant whose G and H describe how an unmeasured force enters the dynamics and the
 * sensors (a force on the carriage that an accelerometer also reads). The process and
 * measurement noise still enter through G and H.
 * @tparam NZ Outputs held offset-free, listed in @ref outputs (0: every measured output)
 */
template<size_t NZ = 0>
struct PlantDisturbance {
    damp::array<size_t, NZ> outputs{}; ///< Indices of the measured outputs held offset-free, each once (none listed: all)
};

/**
 * @brief An explicit disturbance model: d enters the plant's state through B_d and the readings
 *        through C_d.
 *
 * B_d is in the plant's time domain (continuous on a continuous plant, discretized with it by
 * ZOH; per step on a discrete one), one row per plant state.
 *
 * @tparam NX Plant states; NY Measured outputs; ND Disturbance states
 * @tparam NZ Outputs held offset-free, listed in @ref outputs (0: every measured output)
 */
template<size_t NX, size_t NY, size_t ND, size_t NZ = 0, typename T = double>
struct DisturbanceModel {
    Matrix<NX, ND, T>       Bd{};      ///< How d enters the state update
    Matrix<NY, ND, T>       Cd{};      ///< How d enters the readings directly
    damp::array<size_t, NZ> outputs{}; ///< Indices of the measured outputs held offset-free, each once (none listed: all)
};

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
 * @brief An input disturbance: B_d = B of the filter's model (delayed with u), C_d = its D. On a
 *        model with input registers its walk reaches the plant through them, so none of it acts
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

/** @brief The plant's channel: B_d = G of the filter's model (off the registers), C_d = H. */
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

/** @brief An explicit model: B_d on the plant's states (ZOH'd on a continuous plant), C_d as given. */
template<size_t ND, size_t NZ, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXE, size_t NWE, size_t NVE>
[[nodiscard]] constexpr DisturbanceChannels<NX, NXE, NY, ND, T> disturbance_channels(
    const DisturbanceModel<NX, NY, ND, NZ, T>& spec,
    const StateSpace<NX, NU, NY, T, NW, NV>&   plant,
    const StateSpace<NXE, NU, NY, T, NWE, NVE>&,
    T Ts
) {
    DisturbanceChannels<NX, NXE, NY, ND, T> out{.Cd = spec.Cd, .Bd_c = spec.Bd};
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
