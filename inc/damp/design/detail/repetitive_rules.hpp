// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private, include "damp/controllers/repetitive.hpp"

/**
 * @file repetitive_rules.hpp
 * @brief The period design::repetitive sizes its internal model to, and the binomial robustness filter's taps
 */

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp::design::detail {

/**
 * @brief N = round(2π/(ω₀·Ts)), the samples per period of a fundamental @p w0 [rad/s] at @p Ts [s],
 *        or 0 (an invalid period) when either is not finite and positive or the ratio overflows
 */
template<typename T>
[[nodiscard]] constexpr size_t repetitive_period(T w0, T Ts) {
    // finite_positive refuses NaN under -ffast-math too; the cap keeps the size_t cast defined.
    if (!damp::finite_positive(w0) || !damp::finite_positive(Ts)) {
        return 0;
    }
    const T ratio = (T{2} * damp::numbers::pi_v<T>) / (w0 * Ts);
    if (!damp::isfinite(ratio) || !(ratio < static_cast<T>(damp::numeric_limits<size_t>::max() / 2))) {
        return 0;
    }
    return static_cast<size_t>(ratio + static_cast<T>(0.5)); // round to nearest
}

/** @brief Binomial coefficient C(n, k) (small n; double is exact in range). */
[[nodiscard]] constexpr double binomial(size_t n, size_t k) {
    double c = 1.0;
    for (size_t i = 0; i < k; ++i) {
        c = c * static_cast<double>(n - i) / static_cast<double>(i + 1);
    }
    return c;
}

/** @brief Tap @p i of the unity-DC-gain binomial filter of half-width @p M: C(2M, M + i)/4^M (i = 0 the center). */
[[nodiscard]] constexpr double binomial_tap(size_t M, size_t i) {
    double norm = 1.0; // Σ row of C(2M,·) = 2^{2M} = 4^M
    for (size_t j = 0; j < 2 * M; ++j) {
        norm *= 2.0;
    }
    return binomial(2 * M, M + i) / norm;
}

} // namespace damp::design::detail
