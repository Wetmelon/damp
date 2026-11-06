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

#include <cstddef>
#include <cstdint>
#include <type_traits>

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
enum class OffSlot : std::uint8_t {
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
using maybe = std::conditional_t<present, Feature, Off<slot>>;

/** @brief The law's state [x; xi] from the state @p x and the integrators @p tail. */
template<size_t NX, size_t NI, typename T>
[[nodiscard]] constexpr ColVec<NX + NI, T> stacked(const ColVec<NX, T>& x, const ColVec<NI, T>& tail) {
    ColVec<NX + NI, T> z{};
    for (size_t i = 0; i < NX; ++i) {
        z(i) = x(i);
    }
    for (size_t i = 0; i < NI; ++i) {
        z(NX + i) = tail(i);
    }
    return z;
}

/** @brief The leading N entries of @p z: the plant part of a law's or a filter's state. */
template<size_t N, size_t M, typename T>
[[nodiscard]] constexpr ColVec<N, T> leading(const ColVec<M, T>& z) {
    static_assert(N <= M, "The leading part is at most the whole vector.");
    ColVec<N, T> x{};
    for (size_t i = 0; i < N; ++i) {
        x(i) = z(i);
    }
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
 * @brief The discrete model of @p sys: itself on a discrete plant (@p discrete), else sampled at
 *        @p Ts (ZOH) by discretize.
 * @return status invalid_sample_time for a plant in the other time domain or a Ts that is not
 *         finite and positive; bad_weights for a non-finite entry of A, B, C, D, G or H, on
 *         either plant; else discretize's not_converged for a model that overflows within one
 *         period.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr SampledModel<NX, NU, NY, T, NW, NV> sampled_model(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, bool discrete) {
    SampledModel<NX, NU, NY, T, NW, NV> out{};
    if (discrete ? !damp::finite_positive(sys.Ts) : (sys.Ts != T{0} || !damp::finite_positive(Ts))) {
        out.status = DesignStatus::invalid_sample_time; // the plant's own time domain, one discretization
        return out;
    }
    const auto finite = [](T v) { return damp::isfinite(v); };
    if (!mat::all_of(sys.A, finite) || !mat::all_of(sys.B, finite) || !mat::all_of(sys.C, finite) || !mat::all_of(sys.D, finite)
        || !mat::all_of(sys.G, finite) || !mat::all_of(sys.H, finite)) {
        out.status = DesignStatus::bad_weights; // a non-finite plant entry, which discretize would not see on C, D or H
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

} // namespace detail
} // namespace design

/**
 * @brief The law a design result deploys as, named by the result: `Law<kRes, float> law;`
 *
 * A struct member's type cannot be deduced, so it is named by its design instead of its sizes.
 * kRes is a named result with static storage (namespace scope, or a `static` member or local),
 * `constexpr` or initialized at start-up. A default-constructed Law is the design's law; assigning
 * a runtime of the same shape replaces it, as a redesign at run time does.
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
 * (controllers/lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp), MPC (controllers/mpc.hpp, ofmpc.hpp)
 * and the model inverses a law tracks through (design/model_inverse.hpp: StateFeedforward, OutputInverse).
 *
 * @tparam Res The design result
 * @tparam T   The runtime scalar
 */
template<const auto& Res, typename T>
struct Law {
    static_assert(sizeof(T) == 0, "Law<res, T>: res is not a design result a law is built from; include its runtime's header (controllers/lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp, mpc.hpp, ofmpc.hpp, design/model_inverse.hpp)");
};

} // namespace damp
