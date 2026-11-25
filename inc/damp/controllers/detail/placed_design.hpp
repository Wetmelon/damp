// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private: a public header that exports it provides its names

/**
 * @file placed_design.hpp
 * @brief The placed halves of a state-feedback design: the gain at Poles on the law's model
 *        (placed_law_design), and the filter of a design, the Kalman one of a KalmanNoise or the
 *        one placed at ObserverPoles (filter_design; the placed filter is
 *        estimation/detail/placed_filter.hpp's, design::observer's)
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/law_detail.hpp"
#include "damp/controllers/detail/lq_cost.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/options.hpp"
#include "damp/design/pole_placement.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/detail/placed_filter.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/observer.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {
namespace detail {

/**
 * @brief Poles lists the count its gain places: the plant's states plus the integrators. A
 *        compile error otherwise, the two counts its template arguments.
 */
template<size_t kPolesExpected, size_t kPolesGiven>
consteval bool controller_pole_count() {
    static_assert(
        kPolesExpected == kPolesGiven,
        "Poles lists one pole per plant state and per integrator of an IntegralAction (the comparison: that count == the poles listed). "
        "An InputDelay's delay states take theirs at 0, and a disturbance spec's states none."
    );
    return kPolesExpected == kPolesGiven;
}

/**
 * @brief The gain at @p placed on the law's model: on the plant, or with the integrators of @p Ci
 *        (NI > 0) and the delay states of @p delay (an InputDelay) as law_design and
 *        delayed_law_design build it, the delay states' poles at 0.
 *
 * The listed poles are the plant's and the integrators': s-plane on the sampled form, mapped to
 * z = e^{s·Ts} on the discrete model as design::place(sys, poles, Ts) maps them; z-plane on a
 * discrete plant; s-plane for the continuous gain of a continuous plant without Ts (analysis, as
 * the LQ design's CARE). The delay states hold the commands still on their way, and the gain
 * empties them (deadbeat): N·NU poles at 0 on top of the list. Placement is design::place's on
 * the law's (A, B): KNV, the Jordan construction and Ackermann's gain for one input, certified
 * against the request. With one input the gain is unique and takes any delay; with several, a
 * pole is placed at most NU times, so N ≥ 2 delay states per input (a deadbeat chain, which is
 * defective) is refused with not_placeable.
 *
 * @return The design's @p Result as law_design fills it for the LQ weights, S zero (no Riccati equation); status
 *         invalid_parameter for a pole that is not finite, sampled_model's for the plant,
 *         design::place's otherwise (not_placeable, singular)
 */
template<typename Result, Domain D, size_t NI, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NP, typename In>
[[nodiscard]] constexpr Result placed_law_design(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const PlacedPoles<NP, T>& placed, const Matrix<NI, NY, T>& Ci, const In& delay) {
    constexpr size_t NRU = input_delay_states<In>;
    constexpr size_t NXL = NX + (NRU * NU); // the law's model: the plant, then the delay states
    constexpr size_t NK = NXL + NI;
    static_assert(NP + (NRU * NU) == NK, "The poles listed are the plant's and the integrators'.");
    if (!finite_poles(placed.poles)) {
        return refused<Result>(DesignStatus::invalid_parameter); // before z = e^{s·Ts} maps −∞ to 0
    }
    if constexpr ((D == Domain::plant) && NI == 0 && NRU == 0) {
        if (sys.is_continuous()) {
            const auto gain = place_on(sys, placed.poles); // the continuous gain, for analysis
            Result     out{};
            out.K = gain.K;
            out.e = gain.e;
            return set_status(out, gain.status);
        }
    }
    const auto m = sampled_model(sys, Ts, D, delay);
    if (m.status != DesignStatus::ok) {
        return refused<Result>(m.status);
    }
    const T                  h = (D == Domain::plant) ? sys.Ts : Ts;
    Result                   out{.Ts = h};
    StateSpace<NK, NU, 1, T> law{.Ts = h};
    if constexpr (NI == 0) {
        law.A = m.law.A;
        law.B = m.law.B;
    } else {
        const auto im = integral_model(m.law.A, m.law.B, m.law.C, m.law.D, Ci, h); // as the LQ design's
        law.A = im.A;
        law.B = im.B;
        out.integral = IntegralDesign<NI, NXL, NU, NY, T>{Ci, Ci * m.law.C, Ci * m.law.D, h};
    }
    const auto gain = place_on(law, with_deadbeat<NRU * NU, D>(placed.poles, Ts));
    out.K = gain.K;
    out.e = gain.e;
    return set_status(out, gain.status);
}

/** @brief The filter spec among @p opts, for state_feedback and mpc alike: the KalmanNoise sized to the plant, or ObserverPoles (observer_spec_of). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto filter_spec_of(const Opts&... opts) {
    if constexpr (option_count<kind::ObserverPoles, Opts...> > 0) {
        return observer_spec_of<T>(the<kind::ObserverPoles>(opts...));
    } else {
        return noise_of<NX, NU, NY, T, NW, NV>(opts...);
    }
}

/**
 * @brief The filter of a design: the Kalman filter of a KalmanNoise (kalman_design), or the filter
 *        placed at ObserverPoles' PlacedPoles (placed_filter_design), the Kalman design's type.
 */
template<Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Filter, typename Spec, typename Delay>
[[nodiscard]] constexpr auto filter_design(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const Filter& filter, const Spec& spec, const Delay& delay) {
    if constexpr (OptionOf<Filter, kind::KalmanNoise>) {
        return kalman_design<D>(sys, Ts, filter, spec, delay);
    } else {
        return placed_filter_design<D>(sys, Ts, filter, spec, delay);
    }
}

} // namespace detail
} // namespace design
} // namespace damp
