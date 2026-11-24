// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private, include "damp/controllers/law.hpp"

/**
 * @file law_detail.hpp
 * @brief The parts the laws, their designs and estimators share: the applied-input record and
 *        the discrete model a design is on. Law<kRes, T> itself is controllers/law.hpp.
 */

#include "damp/backend.hpp"
#include "damp/design/options.hpp"
#include "damp/design/status.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace detail {

/**
 * @brief A law's memory of its last command and of the input back_calculate recorded since the
 *        last predict: what its time update predicts with.
 *
 * Back-calculation in the sense of anti-windup: the law learns the input that actually went
 * out (after the actuator's clamp, with any feedforward added outside it), and its own last
 * command gives the input it asked for, so the excess u_sat − (u_cmd + u_ff) is the clamp's.
 *
 * @see "Advanced PID Control" (Åström & Hägglund, 2006): back-calculation anti-windup
 * @see "Computer-Controlled Systems" (Åström & Wittenmark, 1997): the observer driven by the
 *      applied input
 */
template<size_t NU, typename T>
struct AppliedRecord {
    ColVec<NU, T> command{};       ///< The last command the law returned
    ColVec<NU, T> applied{};       ///< The input back_calculate recorded, then the one the last predict used
    ColVec<NU, T> feedforward{};   ///< The feedforward added outside the law with it
    bool          recorded{false}; ///< Whether back_calculate has recorded an input since the last predict

    /** @brief Record the input that went out and the feedforward added to the command. */
    constexpr void record(const ColVec<NU, T>& u_sat, const ColVec<NU, T>& u_ff) {
        applied = u_sat;
        feedforward = u_ff;
        recorded = true;
    }

    /** @brief The input the plant received: the recorded one, else the last command (no clamp, no feedforward); clears the record. */
    constexpr void take() {
        if (!recorded) {
            applied = command;
            feedforward = ColVec<NU, T>{};
        }
        recorded = false;
    }
};

/**
 * @brief The slot after @p i in a ring of @p last + 1 slots (i ≤ last): a compare and a wrap, no
 *        division. SmithPredictor's ring of commands advances with it.
 */
[[nodiscard]] constexpr size_t ring_next(size_t i, size_t last) {
    return (i < last) ? (i + 1) : 0;
}

/** @brief The applied-input record of a law that keeps none (no filter or integrators to feed): an empty type. */
struct NoAppliedRecord {};

/**
 * @brief Which optional member an Off stands in for. Two `[[no_unique_address]]` empty members of
 *        one object must have distinct types, so each slot is used at most once per object.
 */
enum class OffSlot : uint8_t {
    target,               ///< StateFeedbackCore::target (NZ = 0)
    integral,             ///< StateFeedbackCore::integral (NI = 0)
    design_integral,      ///< StateFeedbackResult::integral (NI = 0)
    design_target,        ///< StateFeedbackResult::target (NZ = 0)
    held_y,               ///< IntegralAction::y (no estimator: the measurement is not held)
    held_x,               ///< StateFeedbackCore::x_ (with an estimator: the state is not held)
    disturbance_lim,      ///< StateFeedbackCore::disturbance_lim (no disturbance model)
    integral_lim,         ///< StateFeedbackCore::integral_lim (NI = 0)
    design_integral_lim,  ///< StateFeedbackResult::integral_lim (NI = 0)
    tracking_time,        ///< StateFeedbackCore::tracking_time (NI = 0)
    design_tracking_time, ///< StateFeedbackResult::tracking_time (NI = 0)
    held_outputs,         ///< MPC's map from the reference to the target's outputs (every output held: r is the QP's reference)
    held_measured,        ///< The measured disturbance a filter's next predict holds (no filter, or no MeasuredDisturbance)
};

/**
 * @brief An absent optional feature: an empty type, so a `[[no_unique_address]]` member of it takes
 *        no storage. The slot is the unique tag; two Off of the same slot must not be members of
 *        one object.
 */
template<OffSlot Slot>
struct Off {};

/**
 * @brief Feature when @p present, else Off<slot> (no storage under `[[no_unique_address]]`).
 *
 * Alias the result (`using HeldY = maybe<…>`) so the member declaration is the short name.
 */
template<bool present, typename Feature, OffSlot slot>
using maybe = damp::conditional_t<present, Feature, Off<slot>>;

/** @brief The law's state [x; xi] from the state @p x and the integrators @p tail. */
template<size_t NX, size_t NI, typename T>
[[nodiscard]] constexpr ColVec<NX + NI, T> stacked(const ColVec<NX, T>& x, const ColVec<NI, T>& tail) {
    ColVec<NX + NI, T> z{};
    for_each_entry<NX, 1>([&](size_t i, size_t) { z(i) = x(i); });
    for_each_entry<NI, 1>([&](size_t i, size_t) { z(NX + i) = tail(i); });
    return z;
}

/** @brief The leading N entries of @p z: the plant part of a law's or a filter's state. */
template<size_t N, size_t M, typename T>
[[nodiscard]] constexpr ColVec<N, T> leading(const ColVec<M, T>& z) {
    static_assert(N <= M, "The leading part is at most the whole vector.");
    ColVec<N, T> x{};
    for_each_entry<N, 1>([&](size_t i, size_t) { x(i) = z(i); });
    return x;
}

/**
 * @brief A reference point (@p x_ref, @p u_ref) on a law's model: x*, then u* in each of its NRU
 *        delay states per input, the input that holds the point. Shared by the LQ and MPC runtimes.
 */
template<size_t NRU, size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr ColVec<NX + (NRU * NU), T> on_law_model(const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref) {
    if constexpr (NRU == 0) {
        static_cast<void>(u_ref);
        return x_ref;
    } else {
        ColVec<NX + (NRU * NU), T> out{};
        for (size_t i = 0; i < NX; ++i) {
            out(i) = x_ref(i);
        }
        for (size_t i = NX; i < NX + (NRU * NU); ++i) {
            out(i) = u_ref((i - NX) % NU);
        }
        return out;
    }
}

/** @brief A reference point @p ref on a filter's model of N states: the law's, then 0 on the delay states of late readings, which the law does not read. */
template<size_t N, size_t M, typename T>
[[nodiscard]] constexpr ColVec<N, T> on_filter_model(const ColVec<M, T>& ref) {
    static_assert(M <= N, "The filter's model holds the law's.");
    if constexpr (N == M) {
        return ref;
    } else {
        ColVec<N, T> out{};
        for (size_t i = 0; i < M; ++i) {
            out(i) = ref(i);
        }
        return out;
    }
}

} // namespace detail

namespace design {
namespace detail {

/** @brief The discrete model a design is on, and whether the plant was accepted. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
struct SampledModel {
    StateSpace<NX, NU, NY, T, NW, NV> law{};                    ///< The discrete model
    DesignStatus                      status{DesignStatus::ok}; ///< ok, invalid_sample_time, bad_weights or not_converged
};

/**
 * @brief Whether a design takes @p sys: in its own time domain (@p discrete) or continuous with a
 *        usable @p Ts, every entry finite.
 * @return invalid_sample_time for a plant in the other time domain or a Ts that is not finite and
 *         positive; bad_weights for a non-finite entry of A, B, C, D, G or H; else ok
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr DesignStatus plant_status(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, bool discrete) {
    if (discrete ? !damp::finite_positive(sys.Ts) : (!damp::isfinite(sys.Ts) || sys.Ts != T{0} || !damp::finite_positive(Ts))) {
        return DesignStatus::invalid_sample_time; // the plant's own time domain, one discretization
    }
    const auto finite = [](T v) { return damp::isfinite(v); };
    if (!mat::all_of(sys.A, finite) || !mat::all_of(sys.B, finite) || !mat::all_of(sys.C, finite) || !mat::all_of(sys.D, finite)
        || !mat::all_of(sys.G, finite) || !mat::all_of(sys.H, finite)) {
        return DesignStatus::bad_weights; // a non-finite plant entry, which discretize would not see on C, D or H
    }
    return DesignStatus::ok;
}

/**
 * @brief The discrete model of @p sys: itself on a discrete plant (@p discrete), else sampled at
 *        @p Ts (ZOH) by discretize.
 * @return status plant_status's; else discretize's not_converged for a model that overflows
 *         within one period
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr SampledModel<NX, NU, NY, T, NW, NV> sampled_model(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, bool discrete) {
    SampledModel<NX, NU, NY, T, NW, NV> out{.status = plant_status(sys, Ts, discrete)};
    if (out.status != DesignStatus::ok) {
        return out;
    }
    if (discrete) {
        out.law = sys;
        return out;
    }
    const auto d = ::damp::discretize(sys, Ts);
    if (!d) {
        out.status = d.status;
        return out;
    }
    out.law = *d;
    return out;
}

/** @brief No input delay: a design's model is the plant's own (sampled_model). */
struct NoDelay {};

/** @brief The input registers a delay adds to the model, per input: none without one. */
template<typename Delay>
inline constexpr size_t delay_registers = 0;

/** @brief InputDelay<N>: N registers per input. */
template<size_t N, typename S>
inline constexpr size_t delay_registers<InputDelay<N, S>> = N;

/** @brief Not an InputDelay: a kind of trailing design option (options.hpp). */
template<typename>
struct is_input_delay : damp::false_type {};

/** @brief InputDelay<N>: the delay a sampled design is on, given among its trailing options. */
template<size_t N, typename S>
struct is_input_delay<InputDelay<N, S>> : damp::true_type {};

/** @brief Not an OutputDelay: a kind of trailing design option (options.hpp). */
template<typename>
struct is_output_delay : damp::false_type {};

/** @brief OutputDelay<N>: the late reading a sampled filter design is on, given among its trailing options. */
template<size_t N, typename S>
struct is_output_delay<OutputDelay<N, S>> : damp::true_type {};

/** @brief A design's delays when its reading arrives late: the input's (NoDelay or an InputDelay), then the OutputDelay. */
template<typename In, typename Out>
struct PlantDelays {
    In  input{};  ///< NoDelay, or the InputDelay
    Out output{}; ///< The OutputDelay
};

/** @brief PlantDelays: the input's delay states, the ones a law holds. */
template<typename In, typename Out>
inline constexpr size_t delay_registers<PlantDelays<In, Out>> = delay_registers<In>;

/** @brief The delay states per output a late reading adds (NRY): none without an OutputDelay. */
template<typename Delay>
inline constexpr size_t output_delay_states = 0;

/** @brief PlantDelays with an OutputDelay<N>: N per output. */
template<typename In, size_t N, typename S>
inline constexpr size_t output_delay_states<PlantDelays<In, OutputDelay<N, S>>> = N;

/** @brief The delay states a design's delay adds after the plant's NX: NRU per input, then NRY per output. */
template<typename Delay, size_t NU, size_t NY>
inline constexpr size_t delay_states = (delay_registers<Delay> * NU) + (output_delay_states<Delay> * NY);

/** @brief The input part of a design's delay, the one its law holds: NoDelay or an InputDelay. */
template<typename Delay>
[[nodiscard]] constexpr const Delay& input_delay_of(const Delay& delay) {
    return delay;
}

/** @brief PlantDelays: its InputDelay, or NoDelay. */
template<typename In, typename Out>
[[nodiscard]] constexpr const In& input_delay_of(const PlantDelays<In, Out>& delay) {
    return delay.input;
}

/** @brief The delays among a design's options: NoDelay, the InputDelay, or PlantDelays when an OutputDelay is there. */
template<typename... Opts>
[[nodiscard]] constexpr auto delays_of(const Opts&... opts) {
    const auto input = ::damp::detail::option<is_input_delay>(NoDelay{}, opts...);
    if constexpr (::damp::detail::option_count<is_output_delay, Opts...> == 0) {
        return input;
    } else {
        const auto output = ::damp::detail::option<is_output_delay>(NoDelay{}, opts...);
        return PlantDelays<damp::remove_cvref_t<decltype(input)>, damp::remove_cvref_t<decltype(output)>>{input, output};
    }
}

/** @brief sampled_model(sys, Ts, discrete), the model of a design without a delay. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr SampledModel<NX, NU, NY, T, NW, NV> sampled_model(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, bool discrete, NoDelay) {
    return sampled_model(sys, Ts, discrete);
}

/**
 * @brief The model of a continuous @p sys sampled at @p Ts with its input @p delay late:
 *        discretize(sys, Ts, delay), the plant then the registers v(k−1) … v(k−N).
 * @return status as sampled_model's for a continuous plant (a discrete one is invalid_sample_time),
 *         then discretize's for a delay it refuses
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t N, typename S>
[[nodiscard]] constexpr SampledModel<NX + (N * NU), NU, NY, T, NW, NV> sampled_model(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    bool                                     discrete,
    const InputDelay<N, S>&                  delay
) {
    SampledModel<NX + (N * NU), NU, NY, T, NW, NV> out{.status = discrete ? DesignStatus::invalid_sample_time : plant_status(sys, Ts, false)};
    if (out.status != DesignStatus::ok) {
        return out; // a delay is sampled with the plant: there is no delayed model of a discrete one
    }
    const auto d = ::damp::discretize(sys, Ts, delay);
    if (!d) {
        out.status = d.status;
        return out;
    }
    out.law = *d;
    return out;
}

/**
 * @brief The model of a continuous @p sys sampled at @p Ts with its reading late, and its input
 *        too when @p delays holds an InputDelay: discretize(sys, Ts, [input,] output), the plant,
 *        the input's delay states, then the readings' v(k−1) … and z(k−1) … .
 * @return status as the InputDelay overload's
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename In, size_t N, typename S>
[[nodiscard]] constexpr SampledModel<NX + delay_states<PlantDelays<In, OutputDelay<N, S>>, NU, NY>, NU, NY, T, NW, NV> sampled_model(
    const StateSpace<NX, NU, NY, T, NW, NV>&  sys,
    T                                         Ts,
    bool                                      discrete,
    const PlantDelays<In, OutputDelay<N, S>>& delays
) {
    SampledModel<NX + delay_states<PlantDelays<In, OutputDelay<N, S>>, NU, NY>, NU, NY, T, NW, NV> out{
        .status = discrete ? DesignStatus::invalid_sample_time : plant_status(sys, Ts, false)
    };
    if (out.status != DesignStatus::ok) {
        return out; // a delay is sampled with the plant: there is no delayed model of a discrete one
    }
    const auto d = [&] {
        if constexpr (damp::is_same_v<In, NoDelay>) {
            return ::damp::discretize(sys, Ts, delays.output);
        } else {
            return ::damp::discretize(sys, Ts, delays.input, delays.output);
        }
    }();
    if (!d) {
        out.status = d.status;
        return out;
    }
    out.law = *d;
    return out;
}

} // namespace detail
} // namespace design

namespace detail {

/** @brief Whether a design's numbers can be part of a type: it is `constexpr`, not initialized at start-up. */
template<const auto& Res>
concept constant_design = requires { typename damp::bool_constant<(Res.success || true)>; }; // reads it: a discarded read would pass for any object

/** @brief Whether a design a Law runs on constants succeeded and is discrete. Anything else does not compile, and the error names its status. */
template<typename Ts>
[[nodiscard]] consteval bool constant_design_runs(bool success, design::DesignStatus status, Ts sample_time) {
    if (!success || !(sample_time > Ts{0})) {
        refuse_failed_design(success ? design::DesignStatus::invalid_sample_time : status);
    }
    return true;
}

/** @brief Whether a design a Law runs on constants succeeded, for a result that refuses a bad Ts itself. A failed one does not compile, and the error names its status. */
[[nodiscard]] consteval bool constant_design_runs(bool success, design::DesignStatus status) {
    if (!success) {
        refuse_failed_design(status);
    }
    return true;
}

} // namespace detail

} // namespace damp
