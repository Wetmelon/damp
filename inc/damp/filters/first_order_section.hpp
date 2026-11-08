// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file first_order_section.hpp
 * @brief First-order IIR section runtime (the first-order counterpart of Biquad)
 */


#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/filters/iir_design.hpp"

namespace damp {
/**
 * @brief First-order IIR section runtime, Direct Form I.
 *
 * Runs any first-order design (design::FirstOrderCoeffs):
 *
 *   y[n] = b0·x[n] + b1·x[n−1] − a1·y[n−1]
 *
 * The design names the shape: design::highpass_1st gives the washout H(s) = s/(s+ωc)
 * (zero DC gain, "remove the DC / slow drift, keep the changes"), design::lowpass_1st
 * the first-order low-pass. A refused design (fc ≤ 0, fc ≥ Nyquist or Ts ≤ 0) does not
 * compile under constant evaluation and outputs 0 at run time.
 *
 * Use the washout to drift-compensate an integrating sensor, AC-couple a signal, or
 * extract a perturbation from its operating point (the same role the
 * extremum-seeking controller uses internally).
 *
 * @code
 * FirstOrderSection<float> hp{design::highpass_1st(2.0f, 1.0f / 1000.0f)}; // 2 Hz washout @ 1 kHz
 * float ac = hp(sample);
 * @endcode
 *
 * @see design::highpass_1st, design::lowpass_1st, Biquad
 */
template<typename T = float>
class FirstOrderSection {
public:
    constexpr FirstOrderSection() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr FirstOrderSection(const design::FirstOrderCoeffs<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : FirstOrderSection(result.template as<T>()) {}

    /** @brief Construct from first-order design coefficients (design::highpass_1st, design::lowpass_1st). */
    constexpr FirstOrderSection(const design::FirstOrderCoeffs<T>& c) // NOLINT(google-explicit-constructor): assign-from-design
        : b0_(c.b0), b1_(c.b1), a1_(c.a1) {
        if (!c.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(c.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /** @brief Process one sample (Direct Form I). */
    constexpr T operator()(T x) {
        const T y = (b0_ * x) + (b1_ * x1_) - (a1_ * y1_);
        x1_ = x;
        y1_ = y;
        return y;
    }

    constexpr void reset() {
        x1_ = T{0};
        y1_ = T{0};
    }

private:
    T b0_{T{0}}, b1_{T{0}}, a1_{T{0}};
    T x1_{T{0}}, y1_{T{0}};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
FirstOrderSection(const design::FirstOrderCoeffs<T>&) -> FirstOrderSection<T>;
} // namespace damp
