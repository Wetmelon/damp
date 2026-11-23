// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private: a public header that exports it provides its names

/**
 * @file config.hpp
 * @brief A field left out, and the size and matrix of a field whose type is its value's
 *
 * An option whose weight sizes it (IntegralAction) takes that weight's type as a template
 * parameter, so class template argument deduction reads it from its value: a scalar or a Matrix,
 * and Unset left out. The fields sized by it (its limits) have types computed from it, so they
 * take braces.
 */

#include "damp/backend.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp::design::detail {

/**
 * @brief A config field left out: that feature is off, or the weight takes its default.
 *
 * Bare braces in a field whose type is its value's reach the converting constructor, which does
 * not compile and says what to write instead.
 */
struct Unset {
    constexpr Unset() = default;

    /** @brief Bare braces reach this: it does not compile, and says to write a value or name the struct. */
    template<typename U>
    constexpr Unset(const U& value) { // NOLINT(google-explicit-constructor): implicit, so bare braces reach the message
        static_cast<void>(value);
        static_assert(
            sizeof(U) == 0,
            "This field's type comes from its value, and bare braces carry none: write a weight as a scalar or a Matrix "
            "(80.0, diag(1.0, 2.0), Matrix<2, 2>{...}), and name a config field's struct (LQRConfig{...}, KalmanNoise{...}, "
            "InputDisturbance{}, InputDelay<N>{theta}, MPCWeights{...})."
        );
    }
};

/** @brief Whether a field was left out. */
template<typename V>
inline constexpr bool is_unset = damp::is_same_v<V, Unset>;

/** @brief Whether a field's value is a Matrix. */
template<typename V>
inline constexpr bool is_matrix_value = false;

/** @brief A Matrix. */
template<size_t R, size_t C, typename T>
inline constexpr bool is_matrix_value<Matrix<R, C, T>> = true;

/** @brief The size a weight's value gives: R for an R × C Matrix, 1 for a scalar, 0 left out. */
template<typename V>
inline constexpr size_t size_of = 1;

/** @brief Left out. */
template<>
inline constexpr size_t size_of<Unset> = 0;

/** @brief A Matrix: its rows. */
template<size_t R, size_t C, typename T>
inline constexpr size_t size_of<Matrix<R, C, T>> = R;

/** @brief The scalar of a weight's value; double left out. */
template<typename V>
struct scalar_of {
    using type = V; ///< A scalar is its own
};

/** @brief Left out. */
template<>
struct scalar_of<Unset> {
    using type = double; ///< The default scalar
};

/** @brief A Matrix: its element. */
template<size_t R, size_t C, typename T>
struct scalar_of<Matrix<R, C, T>> {
    using type = T; ///< The element type
};

/** @brief The scalar of a weight's value. */
template<typename V>
using scalar_of_t = typename scalar_of<V>::type;

/** @brief Whether a weight's value is an N × N matrix: a square Matrix of that size, a floating-point scalar for N = 1, or left out. */
template<typename V, size_t N>
inline constexpr bool square_value = is_unset<V> || (is_matrix_value<V> ? false : (N == 1 && damp::is_floating_point_v<V>));

/** @brief A square Matrix of size N. */
template<size_t R, size_t C, typename T, size_t N>
inline constexpr bool square_value<Matrix<R, C, T>, N> = (R == N) && (C == N);

/** @brief A weight's value as an N × N Matrix in T (its size checked by the caller), or @p fallback when it was left out. */
template<size_t N, typename T, typename V>
[[nodiscard]] constexpr Matrix<N, N, T> square_matrix(const V& value, const Matrix<N, N, T>& fallback) {
    if constexpr (is_unset<V>) {
        return fallback;
    } else if constexpr (is_matrix_value<V>) {
        return value.template as<T>();
    } else if constexpr (N == 1) {
        return Matrix<N, N, T>{static_cast<T>(value)};
    } else {
        return fallback; // a scalar for N > 1: the caller's static_assert has refused it
    }
}

/** @brief An output list left out: every entry past any output, so only NI = NY reads it as every output. */
inline constexpr size_t kEveryOutput = static_cast<size_t>(-1);

/** @brief The default of a list of N outputs: left out. */
template<size_t N>
[[nodiscard]] constexpr damp::array<size_t, N> every_output() {
    damp::array<size_t, N> out{};
    for (size_t i = 0; i < N; ++i) {
        out[i] = kEveryOutput;
    }
    return out;
}

} // namespace damp::design::detail
