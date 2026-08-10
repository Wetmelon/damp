// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file impulse.hpp
 * @brief Impulse excitation — design payload and runtime generator
 */

#include <cstddef>
#include <cstdint>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp {
namespace design {

/**
 * @struct ImpulseConfig
 * @brief Configuration for a one-shot rectangular impulse (Dirac stand-in).
 *
 * Emits `amplitude` for `width_samples` sample clocks, then zero. Width 1 is
 * a unit sample pulse (scaled by A).
 *
 * @tparam T Scalar type.
 */
template<typename T = double>
struct ImpulseConfig {
    T           amplitude{T{1}};  ///< Pulse height (finite)
    std::size_t width_samples{1}; ///< Samples at amplitude (>= 1)

    [[nodiscard]] constexpr bool valid() const {
        if (!damp::isfinite(amplitude)) {
            return false;
        }
        if (width_samples == 0) {
            return false;
        }
        return true;
    }
};

/**
 * @struct ImpulseResult
 * @brief Impulse design payload.
 */
template<typename T = double>
struct ImpulseResult {
    ImpulseConfig<T> config{};
    bool             success{false};

    template<typename U>
    [[nodiscard]] constexpr ImpulseResult<U> as() const {
        return ImpulseResult<U>{
            ImpulseConfig<U>{static_cast<U>(config.amplitude), config.width_samples},
            success,
        };
    }
};

/**
 * @brief Build an impulse design payload.
 */
template<typename T = double>
[[nodiscard]] constexpr ImpulseResult<T>
impulse(const ImpulseConfig<T>& config) {
    return ImpulseResult<T>{config, config.valid()};
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Rectangular impulse (finite-width Dirac stand-in).
 *
 * Emits @c amplitude for @c width_samples ticks of step, then zero.
 * @ref done is true after the pulse has finished.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
class Impulse {
public:
    constexpr Impulse() = default;

    constexpr explicit Impulse(const design::ImpulseResult<T>& design, T /*Ts*/ = T{0})
        : amplitude_(design.config.amplitude),
          width_(design.config.width_samples),
          valid_(design.success) {}

    /**
     * @brief Emit the next impulse sample and advance the sample count.
     */
    [[nodiscard]] constexpr T step() {
        if (!valid_ || done()) {
            return T{0};
        }
        const T out = (count_ < width_) ? amplitude_ : T{0};
        ++count_;
        return out;
    }

    [[nodiscard]] constexpr bool done() const {
        return !valid_ || count_ >= width_;
    }

    constexpr void reset() { count_ = 0; }

private:
    T           amplitude_{T{0}};
    std::size_t width_{1};
    std::size_t count_{0};
    bool        valid_{false};
};

} // namespace damp
