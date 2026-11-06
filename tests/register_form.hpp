// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file register_form.hpp
 * @brief Test oracle: the register form of an input delay, every command in flight a state.
 *
 * The plants the delayed laws are run on, and the reference designs they are checked
 * against: a delay-augmented model is an ordinary StateSpace, and any delay-free design on it
 * plans around the delay.
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp::test {

/** @brief An input delay written as d·Ts + τ: d whole periods and a remainder 0 ≤ τ < Ts [s]. */
struct DelaySplit {
    size_t whole{0};      ///< d, the whole periods
    double fraction{0.0}; ///< [s] τ, 0 for a whole number of periods (to 16·eps)
};

/** @brief Split @p seconds into whole periods of @p Ts and the remainder; a nonpositive delay or Ts gives {0, 0}. */
[[nodiscard]] constexpr DelaySplit split_delay(double seconds, double Ts) {
    DelaySplit out{};
    if (!(seconds > 0.0) || !(Ts > 0.0)) {
        return out;
    }
    const double ratio = seconds / Ts;
    const double tol = 16.0 * eps<double>() * damp::max(1.0, ratio);
    out.whole = static_cast<size_t>(ratio + tol);
    out.fraction = (ratio - static_cast<double>(out.whole) > tol) ? (seconds - (static_cast<double>(out.whole) * Ts)) : 0.0;
    return out;
}

/**
 * @brief The plant with an input delay of @p seconds and its D past inputs as states:
 *        [x; u(k−1); …; u(k−D)], exact for a held input.
 *
 * Write the delay as d·Ts + τ, 0 ≤ τ < Ts. Within each period the plant holds u(k−d−1) for τ,
 * then u(k−d), so x(k+1) = Φx + Γ₀u(k−d) + Γ₁u(k−d−1), and the output feedthrough reads the
 * input held at the sample instant. A whole delay works on a discrete plant too. D must hold
 * the oldest input the model reads (d, or d + 1 with a fraction). The past inputs get no
 * process noise.
 */
template<size_t D, size_t NX, size_t NU, size_t NY, size_t NW, size_t NV>
[[nodiscard]] constexpr StateSpace<NX + (NU * D), NU, NY, double, NW, NV> register_form(
    const StateSpace<NX, NU, NY, double, NW, NV>& plant,
    double                                        Ts,
    double                                        seconds
) {
    const auto split = split_delay(seconds, Ts);
    const auto d = split.whole;
    const auto tau = split.fraction;

    StateSpace<NX, NU, NY, double, NW, NV> base = plant;
    Matrix<NX, NU>                         gamma0 = plant.B;
    Matrix<NX, NU>                         gamma1{};
    if (plant.Ts == 0.0) {
        base = ::damp::detail::discretize_zoh_impl(plant, Ts).value();
        gamma0 = base.B;
        if (tau > 0.0) {
            const auto head = ::damp::detail::discretize_zoh_impl(plant, tau).value();
            const auto tail = ::damp::detail::discretize_zoh_impl(plant, Ts - tau).value();
            gamma0 = tail.B;
            gamma1 = tail.A * head.B;
        }
    }
    const size_t lag_y = (tau > 0.0) ? d + 1 : d;

    StateSpace<NX + (NU * D), NU, NY, double, NW, NV> out{};
    const auto                                        place = [&](const Matrix<NX, NU>& g, size_t lag) {
        if (lag == 0) {
            out.B.template block<NX, NU>(0, 0) = g;
        } else if (lag <= D) {
            out.A.template block<NX, NU>(0, NX + (NU * (lag - 1))) = g;
        }
    };
    out.A.template block<NX, NX>(0, 0) = base.A;
    place(gamma0, d);
    place(gamma1, d + 1);
    if constexpr (D > 0) {
        out.B.template block<NU, NU>(NX, 0) = Matrix<NU, NU>::identity();
        for (size_t j = 1; j < D; ++j) {
            out.A.template block<NU, NU>(NX + (NU * j), NX + (NU * (j - 1))) = Matrix<NU, NU>::identity();
        }
    }
    out.C.template block<NY, NX>(0, 0) = base.C;
    if (lag_y == 0) {
        out.D = plant.D;
    } else if (lag_y <= D) {
        out.C.template block<NY, NU>(0, NX + (NU * (lag_y - 1))) = plant.D;
    }
    out.G.template block<NX, NW>(0, 0) = base.G;
    out.H = base.H;
    out.Ts = Ts;
    return out;
}

} // namespace damp::test
