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
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"

namespace damp {

namespace design {

/**
 * @brief Design result for a multi-resonant harmonic suppressor.
 * @tparam N Number of harmonics
 * @tparam T Scalar type
 */
template<size_t N, typename T = double>
struct [[nodiscard]] HarmonicSuppressorResult {

    damp::array<PRResult<T>, N> gains{};                             ///< one PR resonator per harmonic
    T                           w_fund{};                            ///< fundamental frequency [rad/s] (for set_fundamental rescaling)
    bool                        success{false};                      ///< status == DesignStatus::ok
    DesignStatus                status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr HarmonicSuppressorResult<N, std::remove_const_t<U>> as() const {
        HarmonicSuppressorResult<N, std::remove_const_t<U>> out{};
        for (size_t i = 0; i < N; ++i) {
            out.gains[i] = gains[i].template as<std::remove_const_t<U>>();
        }
        out.w_fund = static_cast<std::remove_const_t<U>>(w_fund);
        out.success = success;
        out.status = status;
        return out;
    }

    /**
     * @brief Parallel bank of PR resonators as one SISO C(s)
     *
     * State dimension 2N. Failed designs return a zero system.
     */
    [[nodiscard]] constexpr StateSpace<2 * N, 1, 1, T> to_ss() const {
        static_assert(N >= 1, "HarmonicSuppressorResult needs at least one harmonic");
        if (!success) {
            return {};
        }
        return pr_parallel_ss<N - 1>();
    }

    /**
     * @brief Transfer function of @ref to_ss (Leverrier)
     */
    [[nodiscard]] constexpr TransferFunction<2 * N + 1, 2 * N + 1, T> to_tf() const {
        return ::damp::to_transfer_function(to_ss());
    }

private:
    template<size_t I>
    [[nodiscard]] constexpr StateSpace<2 * (I + 1), 1, 1, T> pr_parallel_ss() const {
        if constexpr (I == 0) {
            return gains[0].to_ss();
        } else {
            return *::damp::parallel(pr_parallel_ss<I - 1>(), gains[I].to_ss());
        }
    }
};

/**
 * @brief Synthesize a multi-resonant harmonic suppressor.
 *
 * Builds one PR resonator per harmonic (via @ref pr_harmonics) and validates the
 * spec: a positive, finite sampling time, no 0th harmonic, and every resonator
 * accepted by @ref pr (finite gains, a positive frequency strictly below Nyquist,
 * since a resonator at or above fs/2 cannot be realized). The proportional gain Kp
 * is carried on the fundamental's resonator only.
 *
 * @code
 * const auto d = design::harmonic_suppressor({.Kp = 0.5, .Ki = 200.0, .w0 = w_grid, .wc = 5.0}, Ts,
 *                                            damp::array<size_t, 3>{1, 5, 7});
 * @endcode
 *
 * @tparam N Number of harmonics
 * @param fundamental The fundamental's resonator: Kp (shared, placed on the fundamental), Ki and
 *                    w0 [rad/s] of the fundamental (Ki scaled down per harmonic), wc [rad/s] shared
 *                    (0 = ideal PR)
 * @param Ts          Sample time [s]
 * @param harmonics   Harmonic orders to suppress (e.g. {1, 5, 7, 11})
 * @return The bank; refused with invalid_sample_time for a non-positive or non-finite Ts, with
 *         invalid_parameter for a 0th harmonic, and with the first refused resonator's status otherwise
 */
template<size_t N, typename T = double>
[[nodiscard]] constexpr HarmonicSuppressorResult<N, T>
harmonic_suppressor(const PRSpec<T>& fundamental, T Ts, const damp::array<size_t, N>& harmonics) {
    if (!damp::finite_positive(Ts)) {
        return ::damp::detail::refused<HarmonicSuppressorResult<N, T>>(DesignStatus::invalid_sample_time);
    }
    for (size_t i = 0; i < N; ++i) {
        if (harmonics[i] == 0) {
            return ::damp::detail::refused<HarmonicSuppressorResult<N, T>>(DesignStatus::invalid_parameter); // 0th harmonic (DC) is not a resonant target
        }
    }
    const auto gains = pr_harmonics(fundamental, Ts, harmonics);
    for (size_t i = 0; i < N; ++i) {
        if (!gains[i].success) {
            return ::damp::detail::refused<HarmonicSuppressorResult<N, T>>(gains[i].status);
        }
    }
    return HarmonicSuppressorResult<N, T>{
        .gains = gains,
        .w_fund = fundamental.w0,
        .success = true,
        .status = DesignStatus::ok,
    };
}

} // namespace design

/**
 * @ingroup controllers
 * @brief Multi-resonant harmonic suppressor — a parallel bank of PR resonators.
 *
 * Shared Kp plus N resonators with Kp stripped (applied once on the error). Drop
 * it into a current/voltage loop. Satisfies the (r, y) controller protocol so it
 * composes like any other controller. Built from @ref design::harmonic_suppressor,
 * which carries the fundamental that @ref set_fundamental rescales from.
 *
 * @tparam N Number of harmonics
 * @tparam T Scalar type (float or double)
 */
template<size_t N, typename T = float>
class HarmonicSuppressor {
public:
    static_assert(N >= 1, "HarmonicSuppressor needs at least one harmonic");

    constexpr HarmonicSuppressor() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!std::is_same_v<U, T>)
    constexpr HarmonicSuppressor(const design::HarmonicSuppressorResult<N, U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : HarmonicSuppressor(result.template as<T>()) {}

    constexpr HarmonicSuppressor(const design::HarmonicSuppressorResult<N, T>& design) // NOLINT(google-explicit-constructor): assign-from-design
        : w_fund_(design.w_fund), valid_(design.success) {
        if (!design.success && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(design.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        install_gains(design.gains);
        for (size_t i = 0; i < N; ++i) {
            valid_ = valid_ && resonators_[i].valid();
        }
    }

    template<typename U>
    constexpr explicit HarmonicSuppressor(const HarmonicSuppressor<N, U>& other)
        : Kp_(static_cast<T>(other.Kp())),
          w_fund_(static_cast<T>(other.fundamental())),
          valid_(other.valid()) {
        for (size_t i = 0; i < N; ++i) {
            resonators_[i] = PRController<T>(other.resonator(i));
        }
    }

    /** @brief Suppression command from the loop error (Kp once + sum of resonators). */
    [[nodiscard]] constexpr T control(T error) {
        if (!valid_) {
            return T{0};
        }
        T u = Kp_ * error;
        for (size_t i = 0; i < N; ++i) {
            u += resonators_[i].control(error);
        }
        return u;
    }

    /** @brief Reference-tracking overload (computes the error internally). */
    [[nodiscard]] constexpr T control(T r, T y) { return control(r - y); }

    constexpr void reset() {
        for (size_t i = 0; i < N; ++i) {
            resonators_[i].reset();
        }
    }

    /**
     * @brief Re-tune the bank to a new fundamental (grid-frequency adaptation); the
     * harmonic ratios are preserved (resonator i tracks its original order).
     *
     * Ignored when either fundamental is not positive and finite, or when it would move a
     * resonator to or past Nyquist (the bank is retuned whole or not at all).
     */
    constexpr void set_fundamental(T w_fund) {
        if (!damp::finite_positive(w_fund_) || !damp::finite_positive(w_fund)) {
            return; // no valid prior fundamental to rescale from
        }
        // Every resonator was tuned to w0_i = order_i · w_fund_old; scaling all by
        // the same ratio w_fund/w_fund_old preserves each per-resonator order.
        const T ratio = w_fund / w_fund_;
        for (size_t i = 0; i < N; ++i) {
            if (!(resonators_[i].w0 * ratio * resonators_[i].Ts < damp::numbers::pi_v<T>)) {
                return;
            }
        }
        for (size_t i = 0; i < N; ++i) {
            resonators_[i].set_frequency(resonators_[i].w0 * ratio);
        }
        w_fund_ = w_fund;
    }

    [[nodiscard]] constexpr const auto& resonator(size_t i) const { return resonators_[i]; }
    [[nodiscard]] constexpr T           Kp() const { return Kp_; }
    [[nodiscard]] constexpr T           fundamental() const { return w_fund_; }
    [[nodiscard]] constexpr bool        valid() const { return valid_; }

    template<typename U>
    [[nodiscard]] constexpr HarmonicSuppressor<N, U> as() const {
        return HarmonicSuppressor<N, U>(*this);
    }

private:
    constexpr void install_gains(const damp::array<design::PRResult<T>, N>& gains) {
        Kp_ = gains[0].Kp;
        for (size_t i = 0; i < N; ++i) {
            auto r = gains[i];
            r.Kp = T{0}; // Kp is applied once, not per-harmonic
            resonators_[i] = PRController<T>(r);
        }
    }

    damp::array<PRController<T>, N> resonators_{};
    T                               Kp_{T{0}};
    T                               w_fund_{T{0}};
    bool                            valid_{true};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<size_t N, typename T>
HarmonicSuppressor(const design::HarmonicSuppressorResult<N, T>&) -> HarmonicSuppressor<N, T>;

} // namespace damp
