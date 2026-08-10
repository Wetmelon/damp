// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file backend.hpp
 * @brief Standard-library backend: the small set of `std`-replacement types the
 *        embeddable core needs, aliased into `damp`.
 *
 * The core does not hard-depend on the hosted C++ standard library. It uses a
 * handful of vocabulary types — `array`, `optional`, `tuple`/`pair`, a few
 * `<algorithm>`/`<utility>` helpers — through the `damp::` aliases defined here,
 * which resolve to one of two backends selected at compile time:
 *
 *   - stdlib (default): `damp::array` → `std::array`, etc. Hosted, and usable
 *     standalone (no third-party dependency).
 *   - ETL (`-DDAMP_BACKEND_ETL`): `damp::array` → `etl::array`, etc., for
 *     freestanding / embedded targets without a hosted standard library. Pairs
 *     with the [Embedded Template Library](https://www.etlcpp.com).
 *
 * We do not invent our own primitives — anything not from the stdlib comes from
 * ETL. `std::initializer_list` is intentionally *not* aliased: it is a core
 * language facility (`<initializer_list>`), available even freestanding.
 *
 * @note Selection is unified through `damp/config.hpp`: a single `damp_profile.hpp`
 *       sets `DAMP_BACKEND_ETL` (this header) alongside the math-backend macros
 *       (`damp/math/math_backend.hpp`). This header covers only the
 *       container/utility types.
 */

#include <cstddef>          // IWYU pragma: export
#include <initializer_list> // IWYU pragma: export core language facility; freestanding-safe (see @note above)

#include "damp/config.hpp" // IWYU pragma: export - pulls the profile's backend-selection macros

namespace damp::numbers {

// Mathematical constants — our own definitions (no <numbers> dependency), so the
// core stays freestanding under either backend. Mirrors damp::numbers::*_v names
// so usage reads identically. (Resolves #21's "std::numbers replacement" item.)
template<typename T>
inline constexpr T pi_v = static_cast<T>(3.141592653589793238462643383279502884L);
template<typename T>
inline constexpr T e_v = static_cast<T>(2.718281828459045235360287471352662498L);
template<typename T>
inline constexpr T sqrt2_v = static_cast<T>(1.414213562373095048801688724209698079L);
template<typename T>
inline constexpr T sqrt3_v = static_cast<T>(1.732050807568877293527446341505872367L);
template<typename T>
inline constexpr T inv_pi_v = static_cast<T>(0.318309886183790671537767526745028724L);
template<typename T>
inline constexpr T inv_sqrt2_v = static_cast<T>(0.707106781186547524400844362104849039L);
template<typename T>
inline constexpr T inv_sqrt3_v = static_cast<T>(0.577350269189625764509148780501957456L);
template<typename T>
inline constexpr T ln2_v = static_cast<T>(0.693147180559945309417232121458176568L);
template<typename T>
inline constexpr T log2e_v = static_cast<T>(1.442695040888963407359924681001892137L);

} // namespace damp::numbers

#if defined(DAMP_BACKEND_ETL)

#include <etl/algorithm.h> // IWYU pragma: export
#include <etl/array.h>     // IWYU pragma: export
#include <etl/optional.h>  // IWYU pragma: export
#include <etl/tuple.h>     // IWYU pragma: export
#include <etl/utility.h>   // IWYU pragma: export

namespace damp {
using etl::array;
using etl::clamp;
using etl::forward;
using etl::get;
using etl::index_sequence;
using etl::make_index_sequence;
using etl::make_tuple;
using etl::max;
using etl::min;
using etl::move;
using etl::nullopt;
using etl::nullopt_t;
using etl::optional;
using etl::pair;
using etl::tuple;

// ETL only exposes specialized swap overloads (array/pair/tuple/…), not a
// generic scalar swap. Provide one so LU pivots and friends can damp::swap(float&,
// float&) under the freestanding/ETL backend (matches std::swap).
template<typename T>
constexpr void swap(T& a, T& b) noexcept {
    T tmp = damp::move(a);
    a = damp::move(b);
    b = damp::move(tmp);
}

// ETL's min/max are binary-only; supply the initializer_list overloads (by
// value) so `damp::min({...})` / `damp::max({...})` work on this backend too,
// matching the std backend's std::min/max(initializer_list).
template<typename T>
[[nodiscard]] constexpr T min(std::initializer_list<T> values) {
    T m = *values.begin();
    for (const T& v : values) {
        if (v < m) {
            m = v;
        }
    }
    return m;
}
template<typename T>
[[nodiscard]] constexpr T max(std::initializer_list<T> values) {
    T m = *values.begin();
    for (const T& v : values) {
        if (m < v) {
            m = v;
        }
    }
    return m;
}
} // namespace damp

#else // stdlib backend (default)

#include <algorithm> // IWYU pragma: export
#include <array>     // IWYU pragma: export
#include <optional>  // IWYU pragma: export
#include <tuple>     // IWYU pragma: export
#include <utility>   // IWYU pragma: export

namespace damp {
using std::array;
using std::clamp;
using std::forward;
using std::get;
using std::index_sequence;
using std::make_index_sequence;
using std::make_tuple;
using std::max;
using std::min;
using std::move;
using std::nullopt;
using std::nullopt_t;
using std::optional;
using std::pair;
using std::swap;
using std::tuple;
} // namespace damp

#endif

namespace damp {

/**
 * @brief Ordered {min, max} pair returned by value.
 *
 * Unlike std::minmax (which returns a pair of *references* to its arguments and
 * dangles on temporaries, and lives in the non-freestanding @c algorithm header),
 * this returns by value and is backend-agnostic. Ties return {a, b}.
 */
template<typename T>
[[nodiscard]] constexpr pair<T, T> minmax(const T& a, const T& b) {
    return b < a ? pair<T, T>{b, a} : pair<T, T>{a, b};
}

/**
 * @brief Ordered {min, max} of an initializer list, returned by value.
 *
 * The by-value, backend-agnostic counterpart to `std::minmax(initializer_list)`
 * (which `etl::minmax` lacks entirely). Matches the standard's tie behaviour:
 * leftmost minimum, rightmost maximum. The list must be non-empty (UB otherwise,
 * as in the standard).
 */
template<typename T>
[[nodiscard]] constexpr pair<T, T> minmax(std::initializer_list<T> values) {
    T lo = *values.begin();
    T hi = *values.begin();
    for (const T& v : values) {
        if (v < lo) {
            lo = v;
        } // leftmost minimum
        if (!(v < hi)) {
            hi = v;
        } // rightmost maximum
    }
    return {lo, hi};
}

template<typename U, typename T, size_t N>
[[nodiscard]] constexpr damp::array<U, N> array_as(const damp::array<T, N>& in) {
    damp::array<U, N> out{};
    for (size_t i = 0; i < N; ++i) {
        out[i] = static_cast<U>(in[i]);
    }
    return out;
}
} // namespace damp
