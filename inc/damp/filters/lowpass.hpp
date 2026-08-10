// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lowpass.hpp
 * @brief Nth-order low-pass runtime and continuous LPF exact step
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/filters/iir_design.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

namespace damp {
/**
 * @brief Exact step of continuous first-order LPF @f$ \dot y = -\omega_c (y - u) @f$
 *
 * With @p u held over @p h: @f$ y(h) = u + (y_0 - u)\,e^{-\omega_c h} @f$.
 * Use on event-driven SILs (PWM edges) where sample time is irregular.
 *
 * @param y       Filter state at segment start
 * @param u       Held input over the segment
 * @param omega_c Cutoff [rad/s] (@f$ 2\pi f_c @f$)
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

/**
 * @brief Nth-order low-pass filter
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
     * @brief Simple constructor with cutoff frequency and sample time
     *
     * Fast online design of 1st order filter using ZOH method.
     *
     * @param fc Cutoff frequency [Hz]
     * @param Ts_sample Sample time [s]
     *
     */
    constexpr LowPass(T fc, T Ts_sample)
        requires(N == 1)
    {
        *this = LowPass<1, T>(design::lowpass_1st<T>(fc, Ts_sample));
    }

    /**
     * @brief Discretize a continuous TF and load IIR coefficients.
     *
     * Orders N = 1 and N = 2 use @ref design::to_coeffs (exact Tustin path).
     * Higher orders fall back to impulse-response / Hankel re-ID; if the Hankel
     * matrix is singular, coefficients are cleared to zero (safe no-op filter)
     * rather than leaving a partial / garbage denominator.
     *
     * @param tf        Continuous transfer function (ascending powers of s)
     * @param Ts_sample Sample time [s] (must be > 0 for a meaningful design)
     */
    constexpr LowPass(const TransferFunction<N + 1, N + 1, T>& tf, T Ts_sample) {
        if (!(Ts_sample > T{0})) {
            b = {};
            a = {};
            return;
        }
        if constexpr (N == 1) {
            const auto c = design::to_coeffs(tf, Ts_sample);
            b = {c.b0, c.b1};
            a = {c.a1};
        } else if constexpr (N == 2) {
            const auto c = design::to_coeffs(tf, Ts_sample);
            b = {c.b0, c.b1, c.b2};
            a = {c.a1, c.a2};
        } else {
            auto sys_c = tf.to_state_space().value();
            auto sys_d = *discretize(sys_c, Ts_sample, DiscretizationMethod::Tustin);

            // Compute impulse response
            constexpr size_t MaxK = 2 * N;

            damp::array<T, MaxK + 1> h{};
            h[0] = sys_d.D(0, 0);
            Matrix<N, N, T> A_pow{};
            for (size_t i = 0; i < N; ++i) {
                A_pow(i, i) = T{1};
            }
            for (size_t k = 1; k <= MaxK; ++k) {
                Matrix<N, 1, T> temp = A_pow * sys_d.B;
                h[k] = (sys_d.C * temp)(0, 0);
                A_pow = A_pow * sys_d.A;
            }

            // Set b coefficients
            for (size_t i = 0; i <= N; ++i) {
                b[i] = h[i];
            }

            // Compute a coefficients via Hankel / Prony
            Matrix<N, N, T> M{};
            Matrix<N, 1, T> v{};
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    M(i, j) = h[N + 1 + i - (j + 1)];
                }
                v(i, 0) = -h[N + 1 + i];
            }
            // Solve M·a = v. Singular Hankel (degenerate / non-realizable TF)
            // clears both a and b so the section is a safe zero output.
            if (const auto a_vec = mat::solve(M, v)) {
                for (size_t j = 0; j < N; ++j) {
                    a[j] = (*a_vec)(j, 0);
                }
            } else {
                b = {};
                a = {};
            }
        }
    }

    constexpr LowPass(const damp::array<T, N + 1>& b_, const damp::array<T, N>& a_)
        : b(b_), a(a_) {}

    constexpr LowPass(const design::FirstOrderCoeffs<T>& coeffs)
        requires(N == 1)
        : b{coeffs.b0, coeffs.b1}, a{coeffs.a1} {}

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
} // namespace damp
