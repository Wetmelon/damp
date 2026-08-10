// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file linspace.hpp
 * @brief Sample axes for host analysis: linspace, logspace, geomspace, arange
 *
 * @note Compare with NumPy's numpy.linspace / logspace / geomspace / arange.
 *
 * @note damp::analysis::logspace takes endpoint *values* (e.g. 0.01…100 rad/s),
 *       not decade exponents. That matches typical Bode grids and is closer to
 *       numpy.geomspace than to numpy.logspace(start_exp, stop_exp). Use
 *       geomspace when you want the NumPy name for value endpoints.
 */

#include <cstddef>
#include <vector>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp {
namespace analysis {

template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr std::vector<T> linspace(T start, T end, size_t num) {
    std::vector<T> result;
    result.reserve(num);
    if (num == 1) {
        result.push_back(start);
    } else {
        T step = (end - start) / static_cast<T>(num - 1);
        for (size_t i = 0; i < num; ++i) {
            result.push_back(start + (static_cast<T>(i) * step));
        }
    }
    return result;
}

template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr std::vector<T> linspace(const damp::pair<T, T>& span, size_t num) {
    return linspace(span.first, span.second, num);
}

template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr std::vector<T> logspace(T start, T end, size_t num, T base = T{10}) {
    std::vector<T> result;
    result.reserve(num);
    if (num == 1) {
        result.push_back(start);
    } else {
        // Compute logarithms in the requested base
        const T log_base = damp::log(base);
        T       log_start = damp::log(start) / log_base;
        T       log_end = damp::log(end) / log_base;
        T       step = (log_end - log_start) / static_cast<T>(num - 1);
        for (size_t i = 0; i < num; ++i) {
            result.push_back(damp::pow(base, log_start + (static_cast<T>(i) * step)));
        }
    }
    return result;
}

template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr std::vector<T> logspace(const damp::pair<T, T>& span, size_t num, T base = T{10}) {
    return logspace(span.first, span.second, num, base);
}

/**
 * @brief Geometric sequence from @p start to @p end (endpoint values, not exponents)
 *
 * @f$ x_i = \mathrm{start}\cdot (\mathrm{end}/\mathrm{start})^{i/(n-1)} @f$ for
 * @p num ≥ 2. Both endpoints must be strictly positive or both strictly negative.
 * Invalid inputs (including a zero endpoint or mixed signs) return an empty vector.
 *
 * @note Compare with NumPy's numpy.geomspace(start, stop, num).
 */
template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr std::vector<T> geomspace(T start, T end, size_t num) {
    std::vector<T> result;
    if (num == 0) {
        return result;
    }
    result.reserve(num);
    if (num == 1) {
        result.push_back(start);
        return result;
    }
    const bool both_pos = (start > T{0}) && (end > T{0});
    const bool both_neg = (start < T{0}) && (end < T{0});
    if (!(both_pos || both_neg)) {
        return result;
    }
    const T sign = both_neg ? T{-1} : T{1};
    const T log_s = damp::log(damp::abs(start));
    const T log_e = damp::log(damp::abs(end));
    const T step = (log_e - log_s) / static_cast<T>(num - 1);
    for (size_t i = 0; i < num; ++i) {
        result.push_back(sign * damp::exp(log_s + (static_cast<T>(i) * step)));
    }
    // Exact endpoints (NumPy guarantees these despite floating-point drift).
    result.front() = start;
    result.back() = end;
    return result;
}

template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr std::vector<T> geomspace(const damp::pair<T, T>& span, size_t num) {
    return geomspace(span.first, span.second, num);
}

/**
 * @brief Half-open arithmetic range @f$ [\mathrm{start},\,\mathrm{stop}) @f$ with step @p step
 *
 * Returns values @f$ \mathrm{start} + k\cdot\mathrm{step} @f$ for integer @f$ k \ge 0 @f$
 * while the result stays strictly before @p stop when @p step > 0 (strictly after
 * when @p step < 0). Zero step yields an empty vector.
 *
 * @note Compare with NumPy's numpy.arange(start, stop, step). Floating-point
 *       accumulation can differ slightly from NumPy's length formula on long ranges.
 */
template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr std::vector<T> arange(T start, T stop, T step) {
    std::vector<T> result;
    if (step == T{0}) {
        return result;
    }
    if (step > T{0}) {
        if (!(start < stop)) {
            return result;
        }
        for (size_t k = 0;; ++k) {
            const T x = start + (static_cast<T>(k) * step);
            if (!(x < stop)) {
                break;
            }
            result.push_back(x);
        }
    } else {
        if (!(start > stop)) {
            return result;
        }
        for (size_t k = 0;; ++k) {
            const T x = start + (static_cast<T>(k) * step);
            if (!(x > stop)) {
                break;
            }
            result.push_back(x);
        }
    }
    return result;
}

/// @brief @c arange(start, stop, 1)
template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr std::vector<T> arange(T start, T stop) {
    return arange(start, stop, T{1});
}

/// @brief @c arange(0, stop, 1)
template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr std::vector<T> arange(T stop) {
    return arange(T{0}, stop, T{1});
}

} // namespace analysis
} // namespace damp
