// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file moving_average.hpp
 * @brief Moving-average (boxcar / comb-notch) runtime
 */

#include <cstddef>

#include "damp/backend.hpp"

namespace damp {
/**
 * @brief Window length for a moving-average *comb* that notches @p f_notch and
 *        all its harmonics: N = round(fs / f_notch).
 *
 * An N-tap moving average has spectral nulls at every harmonic of fs/N and unity
 * gain at DC — so choosing N = fs/f_notch nulls f_notch, 2·f_notch, 3·f_notch, …
 * while preserving the DC level. Example: a HV DC bus with 2f/4f/6f… rectifier
 * ripple is cleaned with `comb_notch_window(fs, 2*f_line)` (the ripple
 * fundamental), which nulls every even line harmonic and keeps the bus average.
 *
 * @tparam T Scalar type
 * @param fs_hz      Sample rate [Hz]
 * @param f_notch_hz Lowest frequency to null [Hz] (its harmonics are nulled too)
 * @return N (≥ 1); 0 if the arguments are invalid.
 */
template<typename T>
[[nodiscard]] constexpr size_t comb_notch_window(T fs_hz, T f_notch_hz) {
    if (fs_hz <= T{0} || f_notch_hz <= T{0} || f_notch_hz >= fs_hz) {
        return 0;
    }
    return static_cast<size_t>((fs_hz / f_notch_hz) + static_cast<T>(0.5)); // round to nearest
}

/**
 * @ingroup filters
 * @brief Moving-average (boxcar) filter — also a DC-preserving harmonic-notch comb.
 *
 * Returns the mean of the most recent N samples, `y[k] = (1/N)·Σ_{i<N} x[k-i]`,
 * updated in O(1) with a running sum over a fixed-size circular buffer. Two ways
 * to think about it:
 *  - Smoother: an FIR low-pass with a flat (boxcar) impulse response.
 *  - Comb notch: spectral nulls at every harmonic of fs/N, unity gain at DC
 *    — the cheapest way to strip a periodic ripple (and *all* its harmonics)
 *    from a signal while keeping its average (see @ref comb_notch_window).
 *
 * `MaxN` bounds the buffer at compile time (allocation-free). Default active
 * window is @p MaxN; override at construction or via @ref set_window
 * (clamped to [1, MaxN]). The first N samples are a warm-up (the window fills
 * against an initial zero buffer).
 *
 * @tparam MaxN Buffer capacity (largest supported window)
 * @tparam T    Scalar type (float or double)
 */
template<size_t MaxN, typename T = float>
class MovingAverage {
public:
    static_assert(MaxN >= 1, "MovingAverage needs MaxN >= 1");

    /// Active window = MaxN.
    constexpr MovingAverage() = default;

    /// Construct with an active window of @p window samples (clamped to [1, MaxN]).
    constexpr explicit MovingAverage(size_t window) { set_window(window); }

    /// Set the active window length (clamped to [1, MaxN]) and clear state.
    constexpr void set_window(size_t window) {
        n_ = damp::clamp(window, size_t{1}, MaxN);
        reset();
    }

    /// Push one sample, return the current N-sample average.
    constexpr T operator()(T x) {
        const T old = buffer_[idx_];
        buffer_[idx_] = x;
        idx_ = (idx_ + 1) % n_;
        if (idx_ == 0) {
            // Exact resum once per full window bounds float round-off drift to a
            // single window. A plain sum (not Kahan, which -ffast-math defeats);
            // O(N) here is amortized O(1) and the WCET is O(N) either way.
            T s = T{0};
            for (size_t i = 0; i < n_; ++i) {
                s += buffer_[i];
            }
            sum_ = s;
        } else {
            sum_ += x - old;
        }
        return sum_ / static_cast<T>(n_);
    }

    constexpr void reset() {
        buffer_ = {};
        sum_ = T{0};
        idx_ = 0;
    }

    [[nodiscard]] constexpr size_t window() const { return n_; }

private:
    damp::array<T, MaxN> buffer_{};
    T                    sum_{T{0}};
    size_t               n_{MaxN}; ///< active window length
    size_t               idx_{0};
};
} // namespace damp
