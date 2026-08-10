// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file median.hpp
 * @brief Sliding-window median filter runtime
 */

#include <cstddef>

#include "damp/backend.hpp"

namespace damp {
/**
 * @brief Sliding-window median filter — nonlinear spike/outlier rejection.
 *
 * Returns the median of the last N samples. Unlike a moving average it rejects
 * impulse noise (single-sample spikes, dropouts) without smearing edges — the
 * standard despiker for noisy encoders, ADC glitches, and range finders.
 *
 * `MaxN` bounds the buffer at compile time (allocation-free). Default active
 * window is @p MaxN; override at construction or via @ref set_window
 * (clamped to [1, MaxN]). During warm-up (fewer than N samples seen) the
 * median is taken over the samples available so far.
 *
 * @tparam MaxN Buffer capacity (largest supported window)
 * @tparam T    Scalar type
 */
template<size_t MaxN, typename T = float>
class MedianFilter {
public:
    static_assert(MaxN >= 1, "MedianFilter needs MaxN >= 1");

    /// Active window = MaxN.
    constexpr MedianFilter() = default;

    /// Construct with an active window of @p window samples (clamped to [1, MaxN]).
    constexpr explicit MedianFilter(size_t window) { set_window(window); }

    /// Set the active window length (clamped to [1, MaxN]) and clear state.
    constexpr void set_window(size_t window) {
        n_ = damp::clamp(window, size_t{1}, MaxN);
        reset();
    }

    /// Push one sample, return the current windowed median.
    constexpr T operator()(T x) {
        buffer_[idx_] = x;
        idx_ = (idx_ + 1) % n_;
        if (count_ < n_) {
            ++count_;
        }

        // Insertion-sort a copy of the valid samples (n_ is small).
        damp::array<T, MaxN> s{};
        for (size_t i = 0; i < count_; ++i) {
            s[i] = buffer_[i];
        }
        for (size_t i = 1; i < count_; ++i) {
            const T key = s[i];
            size_t  j = i;
            while (j > 0 && s[j - 1] > key) {
                s[j] = s[j - 1];
                --j;
            }
            s[j] = key;
        }

        const size_t mid = count_ / 2;
        if (count_ % 2 == 0) {
            return (s[mid - 1] + s[mid]) / T{2}; // even window: mean of the two middle
        }
        return s[mid];
    }

    constexpr void reset() {
        buffer_ = {};
        idx_ = 0;
        count_ = 0;
    }

    [[nodiscard]] constexpr size_t window() const { return n_; }

private:
    damp::array<T, MaxN> buffer_{};
    size_t               n_{MaxN};  ///< active window length
    size_t               idx_{0};   ///< ring write position
    size_t               count_{0}; ///< samples seen so far (≤ n_)
};
} // namespace damp
