// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file prbs.hpp
 * @brief PRBS excitation — design payload and runtime generator
 */

#include <cstddef>
#include <cstdint>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp {

namespace detail {

// Galois LFSR helpers for maximal-length PRBS (orders 2..16).

[[nodiscard]] constexpr std::uint32_t prbs_bit_mask(std::size_t order) {
    if (order == 0 || order >= 32) {
        return 0u;
    }
    return (std::uint32_t{1u} << static_cast<unsigned>(order)) - std::uint32_t{1u};
}

[[nodiscard]] constexpr std::uint32_t prbs_feedback_mask(std::size_t order) {
    switch (order) {
        case 2: return 0x3u;
        case 3: return 0x6u;
        case 4: return 0xCu;
        case 5: return 0x14u;
        case 6: return 0x30u;
        case 7: return 0x60u;
        case 8: return 0xB8u;
        case 9: return 0x110u;
        case 10: return 0x240u;
        case 11: return 0x500u;
        case 12: return 0xE08u;
        case 13: return 0x1C80u;
        case 14: return 0x3802u;
        case 15: return 0x6000u;
        case 16: return 0xD008u;
        default: return 0u;
    }
}

[[nodiscard]] constexpr std::uint32_t prbs_advance(std::uint32_t state, std::size_t order) {
    const std::uint32_t lsb = state & std::uint32_t{1u};
    state >>= 1u;
    if (lsb != 0u) {
        state ^= prbs_feedback_mask(order);
    }
    const std::uint32_t mask = prbs_bit_mask(order);
    state &= mask;
    if (state == 0u) {
        state = std::uint32_t{1u};
    }
    return state;
}

} // namespace detail

namespace design {

/**
 * @struct PRBSConfig
 * @brief Configuration for maximal-length pseudo-random binary excitation.
 *
 * The runtime generator emits `u[k] in {+amplitude, -amplitude}` from a
 * Galois LFSR. For order `n`, one period has `(2^n - 1)` chips.
 *
 * @tparam T Scalar type.
 */
template<typename T = double>
struct PRBSConfig {
    static constexpr std::size_t kMinOrder = 2;
    static constexpr std::size_t kMaxOrder = 16;

    T             amplitude{T{1}};                       ///< Output level (+/- amplitude), must be > 0
    std::size_t   lfsr_order{10};                        ///< LFSR order in [kMinOrder, kMaxOrder]
    T             clock_period_s{static_cast<T>(0.001)}; ///< Chip period in seconds (> 0)
    std::uint32_t seed{1u};                              ///< Initial LFSR state (must not map to all-zero)

    /**
     * @brief Validate PRBS configuration.
     * @return true if order/seed/clock/amplitude are valid.
     */
    [[nodiscard]] constexpr bool valid() const {
        if (!damp::finite_positive(amplitude)) {
            return false;
        }
        if (!damp::finite_positive(clock_period_s)) {
            return false;
        }
        if (lfsr_order < kMinOrder || lfsr_order > kMaxOrder) {
            return false;
        }
        const std::uint32_t mask = damp::detail::prbs_bit_mask(lfsr_order);
        if ((seed & mask) == 0u) {
            return false;
        }
        if (damp::detail::prbs_feedback_mask(lfsr_order) == 0u) {
            return false;
        }
        return true;
    }
};

/**
 * @struct PRBSResult
 * @brief PRBS design payload.
 *
 * Includes the validated configuration and sequence period (`2^n - 1` chips).
 * Use `.as<float>()` for embedded deployment.
 *
 * @tparam T Scalar type.
 */
template<typename T = double>
struct PRBSResult {
    PRBSConfig<T> config{};       ///< Validated PRBS configuration
    std::size_t   period_bits{0}; ///< Sequence period in chips: `2^order - 1`
    bool          success{false}; ///< true if `config.valid()`

    /**
     * @brief Convert design payload to another scalar type.
     * @tparam U Target scalar type.
     * @return PRBSResult\<U\> converted element-wise.
     */
    template<typename U>
    [[nodiscard]] constexpr PRBSResult<U> as() const {
        return PRBSResult<U>{
            PRBSConfig<U>{
                static_cast<U>(config.amplitude),
                config.lfsr_order,
                static_cast<U>(config.clock_period_s),
                config.seed,
            },
            period_bits,
            success,
        };
    }
};

/**
 * @brief Build a PRBS design payload from a configuration.
 *
 * @param config PRBS configuration.
 * @return PRBSResult with `success = config.valid()` and period metadata.
 */
template<typename T = double>
[[nodiscard]] constexpr PRBSResult<T>
prbs(const PRBSConfig<T>& config) {
    const bool        valid = config.valid();
    const std::size_t period_bits = valid
                                      ? ((std::size_t{1} << config.lfsr_order) - std::size_t{1})
                                      : std::size_t{0};
    return PRBSResult<T>{config, period_bits, valid};
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Maximal-length PRBS runtime generator.
 *
 * Generates a deterministic binary sequence using a Galois LFSR, mapped to
 * output levels `+amplitude` and `-amplitude`.
 *
 * The sequence period is `(2^order - 1)` chips and never enters the all-zero
 * LFSR state.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
class PRBS {
public:
    constexpr PRBS() = default;

    /**
     * @brief Construct from a PRBS design payload.
     * @param design Validated design payload.
     * @param Ts     Sample period for internal `step()` mode.
     */
    constexpr explicit PRBS(const design::PRBSResult<T>& design, T Ts)
        : config_(design.config), Ts_(Ts), valid_(design.success && Ts > T{0}) {
        reset();
    }

    /**
     * @brief Evaluate PRBS output at absolute time.
     * @param t Time in seconds.
     * @return Sequence value (`+amplitude` or `-amplitude`) at chip index `floor(t/clock_period)`.
     */
    [[nodiscard]] constexpr T step(T t) const {
        if (!config_.valid()) {
            return T{0};
        }

        const T           tc = damp::max(t, T{0});
        const std::size_t period = period_bits();
        auto              chips = static_cast<std::size_t>(damp::floor(tc / config_.clock_period_s));
        if (period > 0) {
            chips %= period;
        }

        std::uint32_t state = seeded_state();
        for (std::size_t i = 0; i < chips; ++i) {
            state = detail::prbs_advance(state, config_.lfsr_order);
        }

        return output_from_state(state);
    }

    /**
     * @brief Evaluate PRBS output at internal time and advance by `Ts`.
     * @return Current sequence value.
     */
    [[nodiscard]] constexpr T step() {
        if (!valid_) {
            return T{0};
        }
        if (done()) {
            return T{0};
        }

        const T out = output_from_state(state_);
        elapsed_chip_time_ += Ts_;

        while (elapsed_chip_time_ >= config_.clock_period_s && !done()) {
            elapsed_chip_time_ -= config_.clock_period_s;
            state_ = detail::prbs_advance(state_, config_.lfsr_order);
            ++chips_generated_;
        }

        return out;
    }

    /**
     * @brief Query whether one full PRBS period has elapsed at absolute time.
     * @param t Time in seconds.
     * @return true if `t >= period_bits * clock_period_s` or config invalid.
     */
    [[nodiscard]] constexpr bool done(T t) const {
        if (!config_.valid()) {
            return true;
        }
        const T period_duration = static_cast<T>(period_bits()) * config_.clock_period_s;
        return t >= period_duration;
    }

    /**
     * @brief Query whether one full PRBS period has been generated internally.
     * @return true once `period_bits()` chips have been emitted.
     */
    [[nodiscard]] constexpr bool done() const {
        return chips_generated_ >= period_bits();
    }

    /**
     * @brief Reset LFSR and internal chip timing to initial state.
     */
    constexpr void reset() {
        state_ = seeded_state();
        chips_generated_ = 0;
        elapsed_chip_time_ = T{0};
    }

    /**
     * @brief Return PRBS period in chips.
     * @return `2^lfsr_order - 1` for valid config, 0 otherwise.
     */
    [[nodiscard]] constexpr std::size_t period_bits() const {
        if (!config_.valid()) {
            return 0;
        }
        return (std::size_t{1} << config_.lfsr_order) - std::size_t{1};
    }

    /**
     * @brief Return current internal LFSR state.
     * @return Raw LFSR register value.
     */
    [[nodiscard]] constexpr std::uint32_t state() const {
        return state_;
    }

private:
    [[nodiscard]] constexpr std::uint32_t seeded_state() const {
        const std::uint32_t mask = detail::prbs_bit_mask(config_.lfsr_order);
        const std::uint32_t seeded = config_.seed & mask;
        if (seeded == 0u) {
            return 1u;
        }
        return seeded;
    }

    [[nodiscard]] constexpr T output_from_state(std::uint32_t state) const {
        return ((state & std::uint32_t{1u}) != 0u) ? config_.amplitude : -config_.amplitude;
    }

    design::PRBSConfig<T> config_{};
    T                     Ts_{T{0}};
    T                     elapsed_chip_time_{T{0}};
    std::size_t           chips_generated_{0};
    std::uint32_t         state_{1u};
    bool                  valid_{false};
};

} // namespace damp
