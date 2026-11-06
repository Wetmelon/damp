// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lowpass.hpp
 * @brief Nth-order low-pass design (design::lowpass), runtime, and continuous LPF exact step
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/filters/iir_design.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"

namespace damp {
/**
 * @brief Exact step of continuous first-order LPF ẏ = −ω_c (y − u)
 *
 * With @p u held over @p h: @f$ y(h) = u + (y_0 - u)\,e^{-\omega_c h} @f$.
 * Use on event-driven SILs (PWM edges) where sample time is irregular.
 *
 * @param y       Filter state at segment start
 * @param u       Held input over the segment
 * @param omega_c Cutoff [rad/s] (2pi f_c)
 * @param h       Segment duration [s]
 * @return        State at segment end
 */
template<typename T = double>
[[nodiscard]] constexpr T continuous_lpf_exact_step(T y, T u, T omega_c, T h) {
    if (!(h > T{0}) || !(omega_c > T{0})) {
        return y;
    }
    return u + ((y - u) * damp::exp(-omega_c * h));
}

namespace design {

/**
 * @struct LowPassResult
 * @brief Direct-form coefficients of an Nth-order low-pass, with success / status
 *
 * y[k] = Σ bᵢ x[k−i] − Σ aᵢ y[k−i]: b holds b₀…b_N, a holds a₁…a_N (a₀ = 1).
 *
 * @tparam N Filter order
 * @tparam T Scalar type
 */
template<size_t N, typename T = double>
struct [[nodiscard]] LowPassResult {
    damp::array<T, N + 1> b{};                                 ///< Numerator taps b₀…b_N
    damp::array<T, N>     a{};                                 ///< Denominator taps a₁…a_N (a₀ = 1)
    bool                  success{false};                      ///< status == DesignStatus::ok
    DesignStatus          status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr LowPassResult<N, U> as() const {
        LowPassResult<N, U> out{};
        for (size_t i = 0; i <= N; ++i) {
            out.b[i] = static_cast<U>(b[i]);
        }
        for (size_t i = 0; i < N; ++i) {
            out.a[i] = static_cast<U>(a[i]);
        }
        out.success = success;
        out.status = status;
        return out;
    }
};

/**
 * @brief Discretize a continuous low-pass TF (Tustin) into direct-form coefficients
 *
 * Orders N = 1 and N = 2 use @ref to_coeffs (exact Tustin path). Higher orders realize the
 * TF, discretize it (Tustin) and convert back with @ref to_transfer_function (expanded
 * pole / zero factors, numerator checked against the Markov parameters).
 *
 * Status: invalid_sample_time if Ts ≤ 0; singular if the TF is not realizable (zero leading
 * denominator coefficient); invalid_parameter if the direct form, rounded to T, has a pole on
 * or outside the unit circle (an unstable TF, or a cutoff too close to 0 or Nyquist for direct
 * form in T: use lowpass_sections); otherwise the discretization status.
 *
 * @param tf Continuous transfer function (highest power of s first)
 * @param Ts [s] Sample time
 * @return LowPassResult with b, a
 */
template<size_t NC, typename T>
    requires(NC >= 2)
[[nodiscard]] constexpr LowPassResult<NC - 1, T> lowpass(const TransferFunction<NC, NC, T>& tf, T Ts) {
    constexpr size_t N = NC - 1;
    using Result = LowPassResult<N, T>;
    if (!(Ts > T{0}) || !damp::isfinite(Ts)) {
        return ::damp::detail::refused<Result>(DesignStatus::invalid_sample_time);
    }
    Result out{};
    if constexpr (N == 1 || N == 2) {
        const auto c = to_coeffs(tf, Ts);
        if (!c.success) {
            return ::damp::detail::refused<Result>(c.status);
        }
        if constexpr (N == 1) {
            out.b = {c.b0, c.b1};
            out.a = {c.a1};
        } else {
            out.b = {c.b0, c.b1, c.b2};
            out.a = {c.a1, c.a2};
        }
    } else {
        const auto sys_c_opt = tf.to_state_space();
        if (!sys_c_opt) {
            return ::damp::detail::refused<Result>(DesignStatus::singular);
        }
        const auto sys_d_opt = discretize(*sys_c_opt, Ts, DiscretizationMethod::Tustin);
        if (!sys_d_opt) {
            return ::damp::detail::refused<Result>(sys_d_opt.status);
        }
        // Certified SS → TF (expanded factors, numerator checked against the Markov parameters).
        const auto tf_d = ::damp::to_transfer_function(*sys_d_opt);
        const T    a0 = tf_d.den[0];
        for (size_t k = 0; k <= N; ++k) {
            out.b[k] = tf_d.num[k] / a0;
        }
        for (size_t k = 1; k <= N; ++k) {
            out.a[k - 1] = tf_d.den[k] / a0;
        }
    }
    // Certify the coefficients as rounded to T: near fc = 0 or Nyquist a float direct form
    // cannot hold the poles inside the unit circle, even from exact coefficients.
    damp::array<T, N + 1> den{T{1}};
    for (size_t k = 0; k < N; ++k) {
        den[k + 1] = out.a[k];
    }
    if (!tf_detail::schur_stable(den)) {
        return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter);
    }
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

/**
 * @brief Nth-order Butterworth low-pass at cutoff @p fc, prewarped Tustin, as direct-form coefficients
 *
 * N = 1 is @ref lowpass_1st, N = 2 @ref lowpass_2nd at Q = 1/√2; N = 3, 4 discretize the
 * Butterworth cascade (prototype prewarped so the −3 dB corner lands on @p fc). Tustin puts
 * every zero at z = −1, so the numerator is (1 + z⁻¹)^N scaled to unit DC gain on the
 * rounded denominator, as lowpass_2nd pins its taps.
 *
 * Status: invalid_sample_time if Ts ≤ 0; invalid_parameter if fc is outside (0, Nyquist), or
 * (N = 3, 4) too close to 0 or Nyquist for a stable direct form in T: use lowpass_sections.
 *
 * @tparam N Filter order (1–4)
 * @param fc [Hz] Cutoff frequency
 * @param Ts [s] Sample time
 * @return LowPassResult for LowPass\<N\>
 * @note Compare with MATLAB®'s butter(N, 2*fc*Ts).
 */
template<size_t N, typename T = double>
    requires(N >= 1 && N <= 4)
[[nodiscard]] constexpr LowPassResult<N, T> lowpass(T fc, T Ts) {
    using Result = LowPassResult<N, T>;
    if (!detail::valid_cutoff(fc, Ts)) {
        return ::damp::detail::refused<Result>(detail::cutoff_status(Ts));
    }
    if constexpr (N == 1) {
        const auto c = lowpass_1st(fc, Ts);
        Result     out{.b = {c.b0, c.b1}, .a = {c.a1}};
        return ::damp::detail::set_status(out, c.status);
    } else if constexpr (N == 2) {
        const auto c = lowpass_2nd(fc, damp::numbers::inv_sqrt2_v<T>, Ts);
        Result     out{.b = {c.b0, c.b1, c.b2}, .a = {c.a1, c.a2}};
        return ::damp::detail::set_status(out, c.status);
    } else {
        const T fa = detail::prewarped_frequency(fc, Ts);
        const T pi = damp::numbers::pi_v<T>;
        Result  out{};
        if constexpr (N == 3) {
            // B3(s) ∝ (s+1)(s² + s + 1) → Q = 1 on the quadratic
            out = lowpass(lowpass_1st(fa).tf * lowpass_2nd(fa, T{1}).tf, Ts);
        } else {
            // B4 factors: Q_k = 1/(2 sin((2k−1)π / 8)), k = 1, 2
            const T q1 = T{1} / (T{2} * damp::sin(pi / T{8}));
            const T q2 = T{1} / (T{2} * damp::sin(T{3} * pi / T{8}));
            out = lowpass(lowpass_2nd(fa, q1).tf * lowpass_2nd(fa, q2).tf, Ts);
        }
        if (!out.success) {
            return out;
        }
        // b = g·(1 + z⁻¹)^N with Σb = 1 + Σa: unit DC gain by construction, as in lowpass_2nd.
        T dc_sum{1};
        for (const T a_k : out.a) {
            dc_sum += a_k;
        }
        T binomial{1};
        for (size_t k = 0; k <= N; ++k) {
            out.b[k] = binomial * dc_sum / static_cast<T>(size_t{1} << N);
            binomial = binomial * static_cast<T>(N - k) / static_cast<T>(k + 1);
        }
        return out;
    }
}

} // namespace design

/**
 * @brief Nth-order low-pass filter (direct form I)
 *
 * Built from a design::lowpass result, `LowPass lpf{design::lowpass<1>(fc, Ts)}`, or
 * directly from the cutoff and sample time, `LowPass<1, float> lpf{fc, Ts}`.
 *
 * @tparam N filter order
 */
template<size_t N, typename T = float>
class LowPass {
    static_assert(N >= 1, "LowPass needs N >= 1");

private:
    damp::array<T, N + 1> b{};      ///< Numerator coefficients
    damp::array<T, N>     a{};      ///< Denominator coefficients
    damp::array<T, N>     x_prev{}; ///< Previous inputs
    damp::array<T, N>     y_prev{}; ///< Previous outputs

public:
    /**
     * @brief From a design::lowpass result (any scalar); converts the taps to T.
     *
     * A failed design is ill-formed under constant evaluation; at run time it leaves the
     * taps zero (zero output).
     */
    template<typename U>
    constexpr LowPass(const design::LowPassResult<N, U>& result) { // NOLINT(google-explicit-constructor): assign-from-design
        if (!result.success && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (result.success) {
            for (size_t i = 0; i <= N; ++i) {
                b[i] = static_cast<T>(result.b[i]);
            }
            for (size_t i = 0; i < N; ++i) {
                a[i] = static_cast<T>(result.a[i]);
            }
        }
    }

    /**
     * @brief Butterworth low-pass at cutoff @p fc [Hz] sampled at @p Ts [s]: design::lowpass<N>(fc, Ts).
     *
     * Refuses an invalid cutoff or sample time the same way as the design-result constructor;
     * call design::lowpass directly to read the status.
     */
    constexpr LowPass(T fc, T Ts) : LowPass(design::lowpass<N>(fc, Ts)) {}

    constexpr auto operator()(T x) {
        T y = b[0] * x;
        for (size_t i = 0; i < N; ++i) {
            y += b[i + 1] * x_prev[i];
            y -= a[i] * y_prev[i];
        }

        // Update previous states
        for (size_t i = N - 1; i > 0; --i) {
            x_prev[i] = x_prev[i - 1];
            y_prev[i] = y_prev[i - 1];
        }
        if (N > 0) {
            x_prev[0] = x;
            y_prev[0] = y;
        }

        return y;
    }

    constexpr void reset() {
        x_prev = {};
        y_prev = {};
    }

    // Default constructor
    constexpr LowPass() = default;

    // Copy constructors
    constexpr LowPass(const LowPass& other) = default;
    constexpr LowPass& operator=(const LowPass& other) = default;

    // Move constructors
    constexpr LowPass(LowPass&& other) noexcept = default;
    constexpr LowPass& operator=(LowPass&& other) noexcept = default;

    constexpr ~LowPass() = default;
};
/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<size_t N, typename T>
LowPass(const design::LowPassResult<N, T>&) -> LowPass<N, T>;

} // namespace damp
