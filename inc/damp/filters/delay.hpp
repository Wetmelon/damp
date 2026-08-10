// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file delay.hpp
 * @brief Discrete-time delay buffer runtime
 */

#include <cstddef>

#include "damp/backend.hpp"

namespace damp {
/**
 * @brief Discrete-time delay buffer
 *
 * Circular buffer of capacity @p MaxDelay. Realizable delay is
 * @f$ 0 \le d \le \texttt{MaxDelay}-1 @f$ samples: @f$ d = \texttt{MaxDelay} @f$
 * collides read and write indices and is clamped down. Prefer sizing
 * @p MaxDelay as one more than the largest delay you need.
 */
template<size_t MaxDelay, typename T = float>
class Delay {
private:
    damp::array<T, MaxDelay> buffer_{};         ///< Circular buffer for delayed samples
    size_t                   write_idx_{0};     ///< Current write position
    size_t                   delay_samples_{1}; ///< Number of samples to delay (≤ MaxDelay−1)

public:
    /**
     * @brief Initialize delay buffer
     * @param delay_samples Samples of delay; clamped to [0, MaxDelay−1]
     */
    constexpr void init(size_t delay_samples) {
        // d == MaxDelay makes read_idx == write_idx (zero delay after write).
        constexpr size_t max_d = MaxDelay > 0 ? MaxDelay - 1 : 0;
        delay_samples_ = delay_samples > max_d ? max_d : delay_samples;
        reset();
    }

    /**
     * @brief Process input sample and return delayed output
     * @param input Current input sample
     * @return Sample from @c delay_samples ticks ago (0-delay returns previous write path)
     */
    constexpr T operator()(T input) {
        // Store current input
        buffer_[write_idx_] = input;

        // Calculate read position
        size_t read_idx = (write_idx_ + MaxDelay - delay_samples_) % MaxDelay;

        // Get delayed output
        T output = buffer_[read_idx];

        // Update write index
        write_idx_ = (write_idx_ + 1) % MaxDelay;

        return output;
    }

    /**
     * @brief Reset delay buffer to zero
     */
    constexpr void reset() {
        buffer_ = {};
        write_idx_ = 0;
    }

    /**
     * @brief Get current delay in samples
     * @return Number of delay samples
     */
    constexpr size_t get_delay_samples() const { return delay_samples_; }
};
} // namespace damp
