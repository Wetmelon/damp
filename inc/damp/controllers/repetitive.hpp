// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file repetitive.hpp
 * @brief Plug-in repetitive controller — internal-model rejection of a periodic
 *        reference/disturbance and *all* its harmonics with one delay loop.
 *
 * Repetitive control embeds the internal model of a periodic signal: a positive-
 * feedback delay line of one fundamental period N = round(fs / f0). That single
 * model puts (near-)infinite loop gain at the fundamental *and every harmonic*
 * simultaneously — the reason one repetitive block out-rejects a stack of
 * individual resonators. Targets: grid-tied inverters (track/clean a sinusoidal
 * reference, reject periodic load), precision motion with periodic trajectories,
 * rotating-machinery ripple.
 *
 * It is a *plug-in*: it sits alongside any base controller and adds a correction.
 *
 *     e   = r - y;                 // tracking error
 *     u   = base.control(r, y) + rc.step(e);
 *
 * Standard discrete realization (robustness filter Q(z), integer phase lead m):
 *
 *     w[k]    = Q(z)·w[k-N] + e[k]      (internal-model memory, length-N buffer)
 *     u_rc[k] = k_rc · w[k - N + m]     (the learned, phase-advanced correction)
 *
 * Robustness filter Q. A scalar Q ∈ (0,1] uniformly rolls the model gain off
 * the high harmonics (where plant phase is uncertain), trading some steady-state
 * rejection for stability. The better choice is a zero-phase low-pass FIR
 * Q(z) = Σ_{i=−M}^{M} q_i z^{−i} with symmetric taps (q_i = q_{−i}): it keeps
 * near-unity gain on the low harmonics (full rejection) yet rolls off to ~0 near
 * Nyquist (robust stability), *without* adding phase. It is realizable inside the
 * loop because Q multiplies the N-delayed signal — the "future" taps z^{+i} read
 * already-buffered samples w[k−N+i] (available for i < N). The default family is
 * the unity-DC-gain binomial: M=1 → [1,2,1]/4, M=2 → [1,4,6,4,1]/16, …
 *
 * The integer phase-lead m pre-applies the correction a few samples early to
 * offset plant transport lag.
 *
 * @see S. Hara, Y. Yamamoto, T. Omata, M. Nakano, "Repetitive Control System: A
 *      New Type Servo System for Periodic Exogenous Signals," IEEE TAC 33(7),
 *      1988, https://doi.org/10.1109/9.1274
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"              // damp::array
#include "damp/matrix/matrix_traits.hpp" // scalar_type_t

namespace damp {

namespace design {

/**
 * @brief Repetitive-controller tuning + period (with optional zero-phase FIR Q).
 *
 * The robustness filter is Q(z) = q_filter (center tap q_0) plus @p q_half
 * symmetric side taps `q_side[i-1] = q_i = q_{-i}`. With `q_half = 0` this is the
 * classic scalar Q. `MaxQHalf` bounds the FIR half-width at compile time.
 *
 * @tparam T        Scalar type
 * @tparam MaxQHalf Maximum FIR Q half-width (0 = scalar Q only)
 */
template<typename T = double, size_t MaxQHalf = 0>
struct RepetitiveConfig {

    using scalar = scalar_type_t<T>;

    size_t                        period{0};           ///< N = round(fs / f0), samples per fundamental period
    scalar                        gain{scalar{1}};     ///< k_rc, repetitive (learning) gain — 0 < k_rc ≲ 2 (plant-dependent)
    scalar                        q_filter{scalar{1}}; ///< Q center tap q_0 (== scalar Q when q_half = 0)
    size_t                        lead{0};             ///< m, phase-lead samples to offset plant lag (0 ≤ m < N)
    size_t                        q_half{0};           ///< FIR Q half-width M (0 = scalar); symmetric side taps below
    damp::array<scalar, MaxQHalf> q_side{};            ///< q_1 … q_M (each side; q_side[i-1] = q_i)

    /// DC gain of the Q-filter, q_0 + 2·Σ q_i (1.0 for a unity-DC robustness filter).
    [[nodiscard]] constexpr scalar q_dc_gain() const {
        scalar dc = q_filter;
        for (size_t i = 0; i < q_half; ++i) {
            dc += scalar{2} * q_side[i];
        }
        return dc;
    }

    [[nodiscard]] constexpr bool valid() const {
        if (period < 2) {
            return false;
        }
        if (gain <= scalar{0} || gain > scalar{2}) {
            return false;
        }
        if (q_filter <= scalar{0}) {
            return false;
        }
        if (q_half > MaxQHalf || q_half >= period) {
            return false; // window must fit inside the period delay (N > M)
        }
        const scalar dc = q_dc_gain();
        if (dc <= scalar{0} || dc > scalar{1} + default_tol<scalar>()) {
            return false; // a robustness filter must not amplify (DC gain ≤ 1)
        }
        if (lead >= period) {
            return false;
        }
        return true;
    }

    /// Convert the configuration to a different scalar precision.
    template<typename U>
    [[nodiscard]] constexpr RepetitiveConfig<std::remove_const_t<U>, MaxQHalf> as() const {
        using out_scalar = scalar_type_t<std::remove_const_t<U>>;
        RepetitiveConfig<std::remove_const_t<U>, MaxQHalf> out{};
        out.period = period;
        out.gain = static_cast<out_scalar>(gain);
        out.q_filter = static_cast<out_scalar>(q_filter);
        out.lead = lead;
        out.q_half = q_half;
        for (size_t i = 0; i < MaxQHalf; ++i) {
            out.q_side[i] = static_cast<out_scalar>(q_side[i]);
        }
        return out;
    }
};

/**
 * @brief Design result for the repetitive controller.
 */
template<typename T = double, size_t MaxQHalf = 0>
struct RepetitiveResult {

    RepetitiveConfig<T, MaxQHalf> config{};
    bool                          success{false};

    template<typename U>
    [[nodiscard]] constexpr RepetitiveResult<std::remove_const_t<U>, MaxQHalf> as() const {
        return RepetitiveResult<std::remove_const_t<U>, MaxQHalf>{config.template as<U>(), success};
    }
};

/**
 * @brief Synthesize a repetitive controller with a scalar robustness filter Q.
 *
 * Computes the period N = round(fs / f0) and validates the tuning. Returns
 * `success = false` (rather than asserting) on an invalid spec so it composes in
 * `constexpr` design code.
 *
 * @code
 * auto                             d = design::repetitive(10000.0, 50.0); // 10 kHz, 50 Hz
 * RepetitiveController<256, double> rc(d);
 * double                           u = base.control(r, y) + rc.step(r - y); // plug-in correction
 * @endcode
 *
 * MATLAB® equivalent (no single builtin — the internal model is a delay loop):
 * @code{.m}
 *   N = round(fs / f0);
 *   % w(k) = q*w(k-N) + e(k);   u_rc(k) = k_rc * w(k-N+m);
 * @endcode
 *
 * @param fs_hz Sample rate [Hz]
 * @param f0_hz Fundamental frequency to reject/track [Hz]
 * @param gain  Repetitive (learning) gain k_rc (default 1)
 * @param q     Robustness filter Q ∈ (0,1] (default 1)
 * @param lead  Phase-lead samples m (default 0)
 */
template<typename T = double>
[[nodiscard]] constexpr RepetitiveResult<T, 0> repetitive(
    T      fs_hz,
    T      f0_hz,
    T      gain = T{1},
    T      q = T{1},
    size_t lead = 0
) {
    if (fs_hz <= T{0} || f0_hz <= T{0} || f0_hz >= fs_hz) {
        return RepetitiveResult<T, 0>{}; // success stays false
    }
    const T      ratio = fs_hz / f0_hz;
    const size_t period = static_cast<size_t>(ratio + static_cast<T>(0.5)); // round to nearest

    const RepetitiveConfig<T, 0> config{
        .period = period,
        .gain = static_cast<scalar_type_t<T>>(gain),
        .q_filter = static_cast<scalar_type_t<T>>(q),
        .lead = lead,
        .q_half = 0,
    };
    return RepetitiveResult<T, 0>{
        .config = config,
        .success = config.valid(),
    };
}

namespace detail {
/// Binomial coefficient C(n, k) (small n; double is exact in range).
[[nodiscard]] constexpr double binomial(size_t n, size_t k) {
    double c = 1.0;
    for (size_t i = 0; i < k; ++i) {
        c = c * static_cast<double>(n - i) / static_cast<double>(i + 1);
    }
    return c;
}
} // namespace detail

/**
 * @brief Synthesize a repetitive controller with a binomial zero-phase FIR Q.
 *
 * Builds the unity-DC-gain binomial robustness filter of half-width M:
 * q_i = C(2M, M+i) / 4^M (M=1 → [1,2,1]/4, M=2 → [1,4,6,4,1]/16, …). This rolls
 * the model gain off the high harmonics for robust stability while preserving
 * (near-)full rejection of the low harmonics — the standard choice over a scalar Q.
 *
 * @tparam M  FIR Q half-width (the binomial order); MaxQHalf of the result
 * @param fs_hz Sample rate [Hz]
 * @param f0_hz Fundamental frequency [Hz]
 * @param gain  Repetitive gain k_rc (default 1)
 * @param lead  Phase-lead samples m (default 0)
 */
template<size_t M, typename T = double>
[[nodiscard]] constexpr RepetitiveResult<T, M> repetitive_binomial(
    T      fs_hz,
    T      f0_hz,
    T      gain = T{1},
    size_t lead = 0
) {
    if (fs_hz <= T{0} || f0_hz <= T{0} || f0_hz >= fs_hz) {
        return RepetitiveResult<T, M>{};
    }
    const size_t period = static_cast<size_t>((fs_hz / f0_hz) + static_cast<T>(0.5));

    using scalar = scalar_type_t<T>;
    // Normalize to unity DC gain: Σ row of C(2M,·) = 2^{2M} = 4^M.
    double norm = 1.0;
    for (size_t i = 0; i < 2 * M; ++i) {
        norm *= 2.0;
    }

    RepetitiveConfig<T, M> config{
        .period = period,
        .gain = static_cast<scalar>(gain),
        .q_filter = static_cast<scalar>(detail::binomial(2 * M, M) / norm),
        .lead = lead,
        .q_half = M,
    };
    for (size_t i = 1; i <= M; ++i) {
        config.q_side[i - 1] = static_cast<scalar>(detail::binomial(2 * M, M + i) / norm);
    }
    return RepetitiveResult<T, M>{
        .config = config,
        .success = config.valid(),
    };
}

} // namespace design

/**
 * @ingroup controllers
 * @brief Plug-in repetitive controller runtime (fixed-size internal model).
 *
 * `MaxPeriod` bounds the internal-model buffer at compile time (allocation-free);
 * the active period N ≤ MaxPeriod is set from the design. `step(error)` returns
 * the repetitive correction to add to the base controller's command. `MaxQHalf`
 * bounds the optional zero-phase FIR Q half-width (0 = scalar Q).
 *
 * @tparam MaxPeriod Buffer capacity (largest supported N + Q half-width)
 * @tparam T         Scalar type (float or double)
 * @tparam MaxQHalf  Maximum FIR Q half-width (0 = scalar Q only)
 */
template<size_t MaxPeriod, typename T = float, size_t MaxQHalf = 0>
class RepetitiveController {
public:
    using scalar = scalar_type_t<T>;
    static_assert(MaxPeriod >= 2, "RepetitiveController needs MaxPeriod >= 2");

    constexpr RepetitiveController() = default;

    constexpr explicit RepetitiveController(const design::RepetitiveConfig<T, MaxQHalf>& config)
        : config_(config), valid_(config.valid() && config.period + config.q_half <= MaxPeriod) {}

    constexpr explicit RepetitiveController(const design::RepetitiveResult<T, MaxQHalf>& design)
        : config_(design.config), valid_(design.success && design.config.period + design.config.q_half <= MaxPeriod) {}

    /// Convert across scalar precision, preserving the learned internal model.
    template<typename U>
    constexpr explicit RepetitiveController(const RepetitiveController<MaxPeriod, U, MaxQHalf>& other)
        : config_(other.config_.template as<T>()), write_idx_(other.write_idx_), valid_(other.valid_) {
        for (size_t i = 0; i < MaxPeriod; ++i) {
            mem_[i] = static_cast<T>(other.mem_[i]);
        }
    }

    /**
     * @brief Advance one tick with the current tracking error; returns the
     *        repetitive correction u_rc to add to the base command.
     */
    [[nodiscard]] constexpr T step(T error) {
        if (!valid_) {
            return T{0};
        }
        const size_t n = config_.period;
        const size_t mq = config_.q_half; // FIR Q half-width (distinct from the doc's phase-lead m = config_.lead)

        // Phase-advanced output tap w[k-N+lead].
        const size_t read_out = (write_idx_ + MaxPeriod - n + config_.lead) % MaxPeriod;
        const T      u_rc = static_cast<T>(config_.gain) * mem_[read_out];

        // Zero-phase FIR Q applied to the N-delayed signal:
        //   Q·w[k-N] = q_0·w[k-N] + Σ_{i=1}^{M} q_i·(w[k-N-i] + w[k-N+i]).
        T qw = static_cast<T>(config_.q_filter) * mem_[(write_idx_ + MaxPeriod - n) % MaxPeriod];
        for (size_t i = 1; i <= mq; ++i) {
            const T older = mem_[(write_idx_ + MaxPeriod - n - i) % MaxPeriod]; // w[k-N-i]
            const T newer = mem_[(write_idx_ + MaxPeriod - n + i) % MaxPeriod]; // w[k-N+i]
            qw += static_cast<T>(config_.q_side[i - 1]) * (older + newer);
        }

        mem_[write_idx_] = qw + error; // w[k] = Q·w[k-N] + e[k]
        write_idx_ = (write_idx_ + 1) % MaxPeriod;

        return u_rc;
    }

    /// Clear the learned internal model (buffer + write index); config preserved.
    constexpr void reset() {
        mem_ = {};
        write_idx_ = 0;
    }

    /// The active configuration (period, gain, Q taps, phase lead).
    [[nodiscard]] constexpr const auto& config() const { return config_; }
    /// True if the design is valid and fit the buffer (period + q_half ≤ MaxPeriod).
    [[nodiscard]] constexpr bool valid() const { return valid_; }

private:
    template<size_t, typename, size_t>
    friend class RepetitiveController; // cross-precision converting ctor

    design::RepetitiveConfig<T, MaxQHalf> config_{};
    damp::array<T, MaxPeriod>             mem_{};
    size_t                                write_idx_{0};
    bool                                  valid_{false};
};

} // namespace damp
