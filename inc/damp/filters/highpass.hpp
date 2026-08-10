// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file highpass.hpp
 * @brief First-order high-pass (washout) filter runtime
 */

#include "damp/filters/iir_design.hpp"

namespace damp {
/**
 * @brief First-order high-pass (washout) filter runtime.
 *
 * The everyday "remove the DC / slow drift, keep the changes" block: a
 * first-order washout with corner @p fc, designed by design::highpass_1st
 * (Tustin bilinear, same family as @ref design::lowpass_1st). Direct Form I:
 *
 *   y[n] = b0·x[n] + b1·x[n−1] − a1·y[n−1]
 *
 * with H(s) = s/(s+ωc). DC gain is zero and HF gain → 1. Invalid @p fc / @p Ts
 * (≤ 0 or fc ≥ Nyquist) yields zero coefficients (output stays 0).
 *
 * Use it to drift-compensate an integrating sensor, AC-couple a signal, or
 * extract a perturbation from its operating point (the same role the
 * extremum-seeking controller uses internally).
 *
 * @code
 * HighPass<float> hp{2.0f, 1.0f / 1000.0f}; // 2 Hz corner @ 1 kHz
 * float ac = hp(sample);
 * @endcode
 *
 * @see design::highpass_1st
 */
template<typename T = float>
class HighPass {
public:
    constexpr HighPass() = default;

    /// @param fc Corner frequency [Hz]. @param Ts Sample time [s].
    constexpr HighPass(T fc, T Ts) {
        const auto c = design::highpass_1st<T>(fc, Ts);
        b0_ = c.b0;
        b1_ = c.b1;
        a1_ = c.a1;
    }

    /// Construct from precomputed design::highpass_1st coefficients.
    constexpr explicit HighPass(const design::FirstOrderCoeffs<T>& c)
        : b0_(c.b0), b1_(c.b1), a1_(c.a1) {}

    /// Process one sample (Direct Form I).
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
} // namespace damp
