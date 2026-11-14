// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file state_feedback_options.hpp
 * @brief The features of a state-feedback design, one named field each (StateFeedbackOptions), and
 *        their compile-time checks: a message per mistake
 */

#include "damp/backend.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp {
namespace design {
namespace detail {

/**
 * @brief An option left out of StateFeedbackOptions: that feature is off.
 *
 * Bare braces for the options, `{.noise = …}`, give every field this type, and the
 * conversion says to name the struct: `design::StateFeedbackOptions{.noise = …}`.
 */
struct Unset {
    constexpr Unset() = default;

    /** @brief Bare braces reach this: it does not compile, and says to name the options struct. */
    template<typename U>
    constexpr Unset(const U& value) { // NOLINT(google-explicit-constructor): implicit, so bare braces reach the message
        static_cast<void>(value);
        static_assert(sizeof(U) == 0, "Name the options struct: design::StateFeedbackOptions{.noise = ...}. Bare braces cannot carry an option's type.");
    }
};

/**
 * @brief The features of a state-feedback law, one named field each; leave a field out and the
 *        feature is off.
 *
 * Write it by name, with the fields you use, in this order:
 *
 * | Field | Value | The law gets |
 * | --- | --- | --- |
 * | `.noise` | design::KalmanNoise{.process, .measurement} | a Kalman filter: the law measures y, not x (LQG) |
 * | `.offset_free` | design::IntegralAction{}, or a disturbance model: InputDisturbance, OutputDisturbance, PlantDisturbance or DisturbanceModel | no steady-state error on the outputs r tracks (LQI, LQGI, OFLQG): by integrators on r − y, or by a disturbance the filter estimates (needs `.noise`) |
 * | `.track` | design::TrackedOutputs{} | the outputs r holds through the model's steady-state target, with no offset-free action (LQR, LQG) |
 *
 * @code
 * design::StateFeedbackOptions{.noise = design::KalmanNoise{.process = Qn, .measurement = Rn}, .offset_free = design::InputDisturbance{}}
 * @endcode
 *
 * Each field's type is a template argument deduced from what you put in it, so the sizes an
 * option carries (the outputs listed, a latency) size the result.
 */
template<typename NoiseT = Unset, typename OffsetFreeT = Unset, typename TrackT = Unset>
struct StateFeedbackOptions {
    NoiseT      noise{};       ///< design::KalmanNoise{.process, .measurement}: estimate the state with a Kalman filter
    OffsetFreeT offset_free{}; ///< design::IntegralAction{}, or a disturbance model the filter estimates: no steady-state error
    TrackT      track{};       ///< design::TrackedOutputs{}: the outputs r holds through the steady-state target
};

/** @brief Whether an option was left out. */
template<typename O>
inline constexpr bool is_unset = damp::is_same_v<O, Unset>;

/** @brief Not integral action. */
template<typename>
inline constexpr bool is_integral_action = false;

/** @brief IntegralAction. */
template<size_t NI>
inline constexpr bool is_integral_action<IntegralAction<NI>> = true;

/** @brief Ci, an NI × NY matrix (matlab::lqi's form). */
template<size_t R, size_t C, typename U>
inline constexpr bool is_integral_action<Matrix<R, C, U>> = true;

/** @brief Not TrackedOutputs. */
template<typename>
inline constexpr bool is_tracked_outputs = false;

/** @brief TrackedOutputs. */
template<size_t NZ>
inline constexpr bool is_tracked_outputs<TrackedOutputs<NZ>> = true;

/** @brief The fields of an options type: none for anything but StateFeedbackOptions. */
template<typename Options>
struct options_traits {
    static constexpr bool named = false; ///< Not a StateFeedbackOptions
    using noise = Unset;                 ///< .noise
    using offset_free = Unset;           ///< .offset_free
    using disturbance = Unset;           ///< .offset_free when it is a disturbance model
    using integrate = Unset;             ///< .offset_free when it is integral action
    using track = Unset;                 ///< .track
};

/** @brief The field types of a StateFeedbackOptions. */
template<typename NoiseT, typename OffsetFreeT, typename TrackT>
struct options_traits<StateFeedbackOptions<NoiseT, OffsetFreeT, TrackT>> {
    static constexpr bool named = true;                                                             ///< A StateFeedbackOptions
    using noise = NoiseT;                                                                           ///< .noise
    using offset_free = OffsetFreeT;                                                                ///< .offset_free
    using disturbance = damp::conditional_t<is_disturbance_model<OffsetFreeT>, OffsetFreeT, Unset>; ///< .offset_free when it is a disturbance model
    using integrate = damp::conditional_t<is_integral_action<OffsetFreeT>, OffsetFreeT, Unset>;     ///< .offset_free when it is integral action
    using track = TrackT;                                                                           ///< .track
};

/** @brief How many outputs a selection integrates on a plant with NY outputs: none for anything but a selection. */
template<typename Sel, size_t NY>
struct integrated_count {
    static constexpr bool   selection = false; ///< Not a selection
    static constexpr size_t value = 0;         ///< NI
};

/** @brief Every output (none listed), or the listed ones. */
template<size_t NI, size_t NY>
struct integrated_count<IntegralAction<NI>, NY> {
    static constexpr bool   selection = true;            ///< A selection
    static constexpr size_t value = (NI == 0) ? NY : NI; ///< NI
};

/** @brief The rows of Ci, an NI × NY matrix. */
template<size_t R, size_t C, typename U, size_t NY>
struct integrated_count<Matrix<R, C, U>, NY> {
    static constexpr bool   selection = (C == NY); ///< A selection when it has a column per output
    static constexpr size_t value = R;             ///< NI
};

/** @brief NI of a design's options. */
template<size_t NY, typename Options>
inline constexpr size_t integrator_count = integrated_count<typename options_traits<Options>::integrate, NY>::value;

/** @brief The weights of a design with @p Options; not deduced, so braces take their sizes from the plant. */
template<size_t NX, size_t NU, size_t NY, typename T, typename Options>
using weights_for = damp::type_identity_t<LQRWeights<NX, NU, integrator_count<NY, Options>, T>>;

/** @brief Whether a disturbance spec fits the plant: a DisturbanceModel on its states and outputs, in its scalar. */
template<typename Spec, size_t NX, size_t NY, typename T>
inline constexpr bool disturbance_fits = true;

/** @brief A DisturbanceModel fits a plant of its sizes and scalar. */
template<size_t NXD, size_t NYD, size_t ND, size_t NZ, typename U, size_t NX, size_t NY, typename T>
inline constexpr bool disturbance_fits<DisturbanceModel<NXD, NYD, ND, NZ, U>, NX, NY, T> = (NXD == NX) && (NYD == NY) && damp::is_same_v<U, T>;

/** @brief Each field of the options is its own kind of option: a message per mistake. @return Every field is */
template<size_t NX, size_t NY, typename T, typename Options>
consteval bool fields_valid() {
    using O = options_traits<Options>;
    using F = O::offset_free;
    constexpr bool noise = is_unset<typename O::noise> || is_kalman_noise<typename O::noise>::value;
    constexpr bool choice = is_unset<F> || is_integral_action<F> || is_disturbance_model<F>;
    constexpr bool fits = disturbance_fits<typename O::disturbance, NX, NY, T>;
    constexpr bool integrate = is_unset<typename O::integrate> || integrated_count<typename O::integrate, NY>::selection;
    constexpr bool track = is_unset<typename O::track> || is_tracked_outputs<typename O::track>;
    constexpr bool one_reference = is_unset<typename O::track> || is_unset<F>;
    static_assert(noise, "The noise is design::KalmanNoise{.process = ..., .measurement = ...}.");
    static_assert(
        choice,
        "options.offset_free is design::IntegralAction{} (integrators on r - y), or a disturbance model the Kalman filter estimates: "
        "design::InputDisturbance{}, OutputDisturbance{}, PlantDisturbance{} or DisturbanceModel{.Bd, .Cd}. Leave it out for neither."
    );
    static_assert(fits, "The DisturbanceModel is on the plant's states and outputs, in the plant's scalar.");
    static_assert(integrate, "The integrated outputs Ci·y: Ci is an NI x NY matrix, a column per measured output.");
    static_assert(track, "options.track is design::TrackedOutputs{} (every output) or TrackedOutputs<NZ>{{i, ...}}.");
    static_assert(one_reference, "Integral action and a disturbance model already hold r: .track is for a law with neither (lqr, lqg).");
    return noise && choice && fits && integrate && track && one_reference;
}

/** @brief The features go together and the noise is sized to the plant: a message per mistake. @return They do and it is */
template<size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Options>
consteval bool features_valid() {
    using O = options_traits<Options>;
    constexpr bool   kNoise = !is_unset<typename O::noise>;
    constexpr size_t ND = [] {
        if constexpr (is_disturbance<typename O::disturbance>::value) {
            return disturbance_traits<typename O::disturbance, NU, NY, NW>::states;
        } else {
            return size_t{0};
        }
    }();
    constexpr bool estimated = kNoise || is_unset<typename O::disturbance>;
    constexpr bool sized = !kNoise || damp::is_same_v<typename O::noise, KalmanNoise<NW, NV, T, ND>>;
    static_assert(estimated, "A disturbance model is estimated by the Kalman filter: give its noise, with .disturbance the disturbance's random walk.");
    static_assert(
        sized || !estimated,
        "The noise is sized to the plant: .process NW x NW (the columns of G), .measurement NV x NV (the columns of H), in the plant's scalar; with a disturbance model, also .disturbance ND x ND (the random walk of d), else none."
    );
    return estimated && sized;
}

/** @brief The spec of the options' steady-state target: the disturbance model, else the tracked outputs; NoDisturbance without either. */
template<typename Options>
[[nodiscard]] constexpr auto disturbance_of(const Options& opts) {
    if constexpr (!is_unset<typename options_traits<Options>::disturbance>) {
        return opts.offset_free;
    } else if constexpr (!is_unset<typename options_traits<Options>::track>) {
        return opts.track;
    } else {
        static_cast<void>(opts);
        return NoDisturbance{};
    }
}

/** @brief The integral action of the options; Unset without it. */
template<typename Options>
[[nodiscard]] constexpr auto integral_action_of(const Options& opts) {
    if constexpr (is_unset<typename options_traits<Options>::integrate>) {
        static_cast<void>(opts);
        return Unset{};
    } else {
        return opts.offset_free;
    }
}

/** @brief No integral action: Ci has no rows. */
template<size_t NI, size_t NY, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NI, NY, T>> integrated_selector(const Unset& none) {
    static_cast<void>(none);
    return Matrix<NI, NY, T>{};
}

/** @brief Ci of integral action: the identity for every output, else the listed rows; nullopt for an index out of range or listed twice. */
template<size_t NI, size_t NY, typename T, size_t NS>
[[nodiscard]] constexpr damp::optional<Matrix<NI, NY, T>> integrated_selector(const IntegralAction<NS>& sel) {
    if constexpr (NS == 0) {
        static_cast<void>(sel);
        return Matrix<NY, NY, T>::identity();
    } else {
        const auto rows = selected_outputs<NY>(sel.outputs);
        if (!rows) {
            return damp::nullopt;
        }
        Matrix<NI, NY, T> Ci{};
        for (size_t i = 0; i < NI; ++i) {
            Ci(i, (*rows)[i]) = T{1};
        }
        return Ci;
    }
}

/** @brief Ci as given, in the plant's scalar. */
template<size_t NI, size_t NY, typename T, typename U>
[[nodiscard]] constexpr damp::optional<Matrix<NI, NY, T>> integrated_selector(const Matrix<NI, NY, U>& Ci) {
    return Ci.template as<T>();
}

} // namespace detail
} // namespace design
} // namespace damp
