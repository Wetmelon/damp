// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file law_detail.hpp
 * @brief The parts the laws, their designs and estimators share: the applied-input record and
 *        the discrete model a design is on
 */

#include "damp/backend.hpp"
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
    target,          ///< StateFeedbackCore::target (NZ = 0)
    integral,        ///< StateFeedbackCore::integral (NI = 0)
    design_integral, ///< StateFeedbackResult::integral (NI = 0)
    design_target,   ///< StateFeedbackResult::target (NZ = 0)
    held_y,          ///< IntegralAction::y (no estimator: the measurement is not held)
    held_x,          ///< StateFeedbackCore::x_ (with an estimator: the state is not held)
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
    if (discrete ? !damp::finite_positive(sys.Ts) : (sys.Ts != T{0} || !damp::finite_positive(Ts))) {
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

/**
 * @brief The runtime a design result deploys as, named by the result: `Law<kRes, float> law;`
 *
 * A struct member's type cannot be deduced, so it is named by its design instead of its sizes.
 * kRes is a named result with static storage (namespace scope, or a `static` member or local).
 *
 * On a `constexpr` kRes the design's matrices are constants of the type: they take no RAM, and a
 * coefficient that is exactly 0 is no code, which removes the sparse work of a delay register, a
 * plant with an integrator or a measurement of a single state. The ticks and the numbers are the
 * runtime's own. Nothing can replace them, so a member a redesign assigns is
 * `Law<kRes, float>::Runtime law{kRes.as<float>()};`. The coefficients move from RAM to flash, as
 * constants in the code: where flash has wait states and no cache, keep the tick in RAM or tightly
 * coupled memory.
 *
 * On a kRes initialized at start-up, from values measured or identified at init, the matrices are
 * members, and assigning a runtime of the same shape replaces them.
 *
 * @code
 * inline constexpr auto kRes = design::lqg(plant, {.Q = Q, .R = R}, {.process = Qn, .measurement = Rn});
 *
 * struct Axis {
 *     Law<kRes, float> law; // the design's law, in float
 *     float            offset{0.0f};
 * };
 * @endcode
 *
 * The runtime's header defines it for the results that runtime is built from: the LQ runtimes
 * (controllers/lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp) and the observers
 * (estimation/fixed_gain_observer.hpp: FixedGainObserver, ReducedOrderObserver), design::eso (estimation/eso.hpp),
 * design::adrc (controllers/adrc.hpp) and the model inverses (design/model_inverse.hpp:
 * StateFeedforward, OutputInverse), which all run on constants when kRes is `constexpr`, and MPC
 * (controllers/mpc.hpp, ofmpc.hpp), whose QP data stays in members and whose OFMPC filter is on constants.
 * design::kalman's Law also names its recursive filter, `Law<kRes, float>::Recursive`: the KalmanFilter
 * of the same model and noise, on the same constants when kRes is `constexpr`.
 *
 * @tparam Res The design result
 * @tparam T   The runtime scalar
 */
template<const auto& Res, typename T>
struct Law {
    static_assert(sizeof(T) == 0, "Law<res, T>: res is not a design result a law is built from; include its runtime's header (controllers/lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp, mpc.hpp, ofmpc.hpp, adrc.hpp, estimation/fixed_gain_observer.hpp, estimation/eso.hpp, design/model_inverse.hpp)");
};

} // namespace damp
