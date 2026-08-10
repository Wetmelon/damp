// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file harmonic_suppression.hpp
 * @brief Active harmonic suppression — a bank of proportional-resonant (PR)
 *        resonators that cancel a chosen set of harmonics in a control loop.
 *
 * The companion to the spectral *detector* (filters/spectral.hpp `Goertzel` /
 * `HarmonicAnalyzer`): once you know which harmonics to kill, drop a
 * `HarmonicSuppressor` into the current/voltage loop. Each PR resonator places
 * (near-)infinite loop gain at one harmonic, so the loop drives that harmonic's
 * error to zero. Targets: grid-tied inverter current control (reject 5th/7th/…
 * line harmonics), active filters, motor torque-ripple cancellation.
 *
 * Repetitive control (controllers/repetitive.hpp) rejects *every* harmonic of a
 * period with one delay loop; this PR bank instead targets a *specific, sparse*
 * set (e.g. {1, 5, 7, 11, 13}) with independent per-harmonic gains — cheaper and
 * more selective when only a few harmonics matter.
 *
 * @see controllers/pr.hpp for the single PR resonator and design::pr_harmonics
 * @see filters/spectral.hpp for the detection side
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp" // damp::array, damp::numbers
#include "damp/controllers/pr.hpp"

namespace damp {

namespace design {

/**
 * @brief Design result for a multi-resonant harmonic suppressor.
 * @tparam N Number of harmonics
 * @tparam T Scalar type
 */
template<size_t N, typename T = double>
struct HarmonicSuppressorResult {

    damp::array<PRResult<T>, N> gains{};  ///< one PR resonator per harmonic
    T                           w_fund{}; ///< fundamental frequency [rad/s] (for set_fundamental rescaling)
    bool                        success{false};

    template<typename U>
    [[nodiscard]] constexpr HarmonicSuppressorResult<N, std::remove_const_t<U>> as() const {
        HarmonicSuppressorResult<N, std::remove_const_t<U>> out{};
        for (size_t i = 0; i < N; ++i) {
            out.gains[i] = gains[i].template as<std::remove_const_t<U>>();
        }
        out.w_fund = static_cast<std::remove_const_t<U>>(w_fund);
        out.success = success;
        return out;
    }
};

/**
 * @brief Synthesize a multi-resonant harmonic suppressor.
 *
 * Builds one PR resonator per harmonic (via @ref pr_harmonics) and validates the
 * spec — frequencies positive, sampling time positive, and the highest harmonic
 * strictly below Nyquist (a resonator at or above fs/2 cannot be realized). The
 * proportional gain Kp is carried on the fundamental's resonator only.
 *
 * @tparam N Number of harmonics
 * @param Kp        Proportional gain (shared; placed on the fundamental)
 * @param Ki_fund   Resonant gain for the fundamental (scaled down per harmonic)
 * @param w_fund    Fundamental frequency [rad/s]
 * @param wc        Resonant-term bandwidth [rad/s] (0 = ideal PR)
 * @param Ts        Sample time [s]
 * @param harmonics Harmonic orders to suppress (e.g. {1, 5, 7, 11})
 */
template<size_t N, typename T = double>
[[nodiscard]] constexpr HarmonicSuppressorResult<N, T> harmonic_suppressor(
    T                             Kp,
    T                             Ki_fund,
    T                             w_fund,
    T                             wc,
    T                             Ts,
    const damp::array<size_t, N>& harmonics
) {
    if (Ts <= T{0} || w_fund <= T{0} || wc < T{0}) {
        return HarmonicSuppressorResult<N, T>{};
    }
    const T nyquist = damp::numbers::pi_v<T> / Ts; // [rad/s]
    for (size_t i = 0; i < N; ++i) {
        if (harmonics[i] == 0) {
            return HarmonicSuppressorResult<N, T>{}; // 0th harmonic (DC) is not a resonant target
        }
        if ((w_fund * static_cast<T>(harmonics[i])) >= nyquist) {
            return HarmonicSuppressorResult<N, T>{}; // resonator at/above Nyquist is unrealizable
        }
    }
    return HarmonicSuppressorResult<N, T>{
        .gains = pr_harmonics(Kp, Ki_fund, w_fund, wc, Ts, harmonics),
        .w_fund = w_fund,
        .success = true,
    };
}

} // namespace design

/**
 * @ingroup controllers
 * @brief Multi-resonant harmonic suppressor — a parallel bank of PR resonators.
 *
 * Sums the outputs of N PR resonators tuned to the target harmonics; drop it into
 * a current/voltage loop on the error signal. Satisfies the (r, y) controller
 * protocol so it composes like any other controller.
 *
 * @tparam N Number of harmonics
 * @tparam T Scalar type (float or double)
 */
template<size_t N, typename T = float>
class HarmonicSuppressor {
public:
    static_assert(N >= 1, "HarmonicSuppressor needs at least one harmonic");

    constexpr HarmonicSuppressor() = default;

    constexpr explicit HarmonicSuppressor(const damp::array<design::PRResult<T>, N>& gains)
        // No fundamental is carried with a raw gains array; assume resonator 0 is
        // the fundamental (order 1) so set_fundamental() can rescale proportionally.
        : w_fund_(gains[0].w0) {
        for (size_t i = 0; i < N; ++i) {
            resonators_[i] = PRController<T>(gains[i]);
        }
    }

    constexpr explicit HarmonicSuppressor(const design::HarmonicSuppressorResult<N, T>& design)
        : w_fund_(design.w_fund), valid_(design.success) {
        for (size_t i = 0; i < N; ++i) {
            resonators_[i] = PRController<T>(design.gains[i]);
        }
    }

    /// Suppression command from the loop error (sum of the resonator outputs).
    [[nodiscard]] constexpr T control(T error) {
        T u = T{0};
        for (size_t i = 0; i < N; ++i) {
            u += resonators_[i].control(error);
        }
        return u;
    }

    /// Reference-tracking overload (computes the error internally).
    [[nodiscard]] constexpr T control(T r, T y) { return control(r - y); }

    constexpr void reset() {
        for (size_t i = 0; i < N; ++i) {
            resonators_[i].reset();
        }
    }

    /// Re-tune the bank to a new fundamental (grid-frequency adaptation); the
    /// harmonic ratios are preserved (resonator i tracks its original order).
    constexpr void set_fundamental(T w_fund) {
        if (w_fund_ <= T{0} || w_fund <= T{0}) {
            return; // no valid prior fundamental to rescale from
        }
        // Every resonator was tuned to w0_i = order_i · w_fund_old; scaling all by
        // the same ratio w_fund/w_fund_old preserves each per-resonator order.
        const T ratio = w_fund / w_fund_;
        for (size_t i = 0; i < N; ++i) {
            resonators_[i].set_frequency(resonators_[i].w0 * ratio);
        }
        w_fund_ = w_fund;
    }

    [[nodiscard]] constexpr const auto& resonator(size_t i) const { return resonators_[i]; }
    [[nodiscard]] constexpr bool        valid() const { return valid_; }

private:
    damp::array<PRController<T>, N> resonators_{};
    T                               w_fund_{T{0}};
    bool                            valid_{true};
};

} // namespace damp
