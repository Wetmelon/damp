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
 * The core does not hard-depend on the C++ standard library. It uses vocabulary
 * types, type traits, limits and concepts — `array`, `optional`, `tuple`/`pair`,
 * `is_same_v`, `numeric_limits`, `floating_point`, … — through the `damp::` aliases
 * defined here, which resolve to one of two backends selected at compile time:
 *
 *   - stdlib (default): `damp::array` → `std::array`, etc. Usable standalone (no
 *     third-party dependency).
 *   - ETL (`DAMP_BACKEND_ETL` in `damp_profile.hpp`): `damp::array` → `etl::array`,
 *     etc., with the [Embedded Template Library](https://www.etlcpp.com) on the
 *     include path. Add `ETL_NO_STL` for a toolchain with no C++ standard headers.
 *
 * Anything not from the stdlib comes from ETL, with one exception: under `ETL_NO_STL`
 * `etl::optional` of anything but a fundamental type is not constexpr, so this header provides
 * the optional that constexpr designs return. `std::initializer_list` is intentionally
 * *not* aliased: it is a language name, which ETL declares itself without the STL.
 *
 * @note Selection is unified through `damp/config.hpp`: a single `damp_profile.hpp`
 *       sets `DAMP_BACKEND_ETL` (this header) alongside the math-backend macros
 *       (`damp/math/math_backend.hpp`). This header covers only the
 *       container/utility types.
 */

#include <stddef.h> // IWYU pragma: export C headers: every toolchain has them, even with no C++ headers
#include <stdint.h> // IWYU pragma: export

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

#ifdef DAMP_BACKEND_ETL

#include <etl/algorithm.h>        // IWYU pragma: export
#include <etl/array.h>            // IWYU pragma: export
#include <etl/bit.h>              // IWYU pragma: export
#include <etl/concepts.h>         // IWYU pragma: export
#include <etl/initializer_list.h> // IWYU pragma: export declares std::initializer_list without the STL
#include <etl/limits.h>           // IWYU pragma: export
#include <etl/optional.h>         // IWYU pragma: export
#include <etl/tuple.h>            // IWYU pragma: export
#include <etl/type_traits.h>      // IWYU pragma: export
#include <etl/utility.h>          // IWYU pragma: export

// etl/limits.h includes C's math.h, whose classification macros would rewrite damp::isfinite and
// friends. <cmath> removes them with the STL; do the same without it.
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit
#undef isgreater
#undef isgreaterequal
#undef isless
#undef islessequal
#undef islessgreater
#undef isunordered

namespace damp {
using etl::array;
using etl::bool_constant;
using etl::clamp;
using etl::common_type_t;
using etl::conditional_t;
using etl::constructible_from;
using etl::convertible_to;
using etl::declval;
using etl::false_type;
using etl::floating_point;
using etl::forward;
using etl::get;
using etl::index_sequence;
using etl::integral;
using etl::is_arithmetic_v;
using etl::is_base_of_v;
using etl::is_const_v;
using etl::is_constructible_v;
using etl::is_floating_point_v;
using etl::is_integral_v;
using etl::is_same_v;
using etl::is_signed_v;
using etl::is_unsigned_v;
using etl::make_index_sequence;
using etl::make_signed_t;
using etl::make_tuple;
using etl::make_unsigned_t;
using etl::max;
using etl::min;
using etl::move;
using etl::nullopt;
using etl::nullopt_t;
using etl::numeric_limits;
using etl::pair;
using etl::remove_const_t;
using etl::remove_cv_t;
using etl::remove_cvref_t;
using etl::remove_pointer_t;
using etl::true_type;
using etl::tuple;
using etl::type_identity_t;

/** @brief True during constant evaluation; the builtin needs no header. */
[[nodiscard]] constexpr bool is_constant_evaluated() noexcept {
    return __builtin_is_constant_evaluated();
}

/** @brief Reinterpret the bits of @p from as To, usable in constant evaluation (etl::bit_cast is not, without the STL). */
template<typename To, typename From>
    requires(sizeof(To) == sizeof(From))
[[nodiscard]] constexpr To bit_cast(const From& from) noexcept {
    return __builtin_bit_cast(To, from);
}

/** @brief Bits needed to represent @p value, as int like std::bit_width (etl::bit_width returns T). */
template<typename T>
[[nodiscard]] constexpr int bit_width(T value) noexcept {
    return static_cast<int>(etl::bit_width(value));
}

#if ETL_USING_STL
using etl::optional;
#else
/**
 * @brief The optional of a toolchain with no C++ standard headers (ETL_NO_STL).
 *
 * etl::optional of anything but a fundamental type is constexpr only with the STL, POD structs
 * included, because it needs std::construct_at, and constexpr designs return one. This holds a value and a flag
 * instead, so T must be default-constructible and an empty optional still holds a T{}.
 * value() on an empty optional returns that T{}: there is nothing to throw and no trap on the
 * tick path, so check has_value() as with any optional.
 */
template<typename T>
class optional {
    T    value_{};
    bool engaged_{false};

public:
    using value_type = T;

    constexpr optional() = default;
    constexpr optional(nullopt_t) {} // NOLINT(google-explicit-constructor): as std::optional

    template<typename U = T>
        requires(is_constructible_v<T, U&&> && !is_same_v<remove_cvref_t<U>, optional> && !is_same_v<remove_cvref_t<U>, nullopt_t>)
    constexpr optional(U&& value) // NOLINT(google-explicit-constructor): as std::optional
        : value_(static_cast<T>(damp::forward<U>(value))),
          engaged_(true) {}

    /** @brief From an optional of another type, as std::optional converts. */
    template<typename U>
        requires(!is_same_v<U, T> && is_constructible_v<T, const U&>)
    constexpr optional(const optional<U>& other) { // NOLINT(google-explicit-constructor): as std::optional
        if (other.has_value()) {
            value_ = static_cast<T>(*other);
            engaged_ = true;
        }
    }

    [[nodiscard]] constexpr bool     has_value() const noexcept { return engaged_; }
    [[nodiscard]] constexpr explicit operator bool() const noexcept { return engaged_; }

    [[nodiscard]] constexpr T&       value() noexcept { return value_; }
    [[nodiscard]] constexpr const T& value() const noexcept { return value_; }
    [[nodiscard]] constexpr T&       operator*() noexcept { return value_; }
    [[nodiscard]] constexpr const T& operator*() const noexcept { return value_; }
    [[nodiscard]] constexpr T*       operator->() noexcept { return &value_; }
    [[nodiscard]] constexpr const T* operator->() const noexcept { return &value_; }

    template<typename U>
    [[nodiscard]] constexpr T value_or(U&& fallback) const {
        return engaged_ ? value_ : static_cast<T>(damp::forward<U>(fallback));
    }

    constexpr void reset() noexcept {
        value_ = T{};
        engaged_ = false;
    }

    template<typename... Args>
    constexpr T& emplace(Args&&... args) {
        value_ = T(damp::forward<Args>(args)...);
        engaged_ = true;
        return value_;
    }
};
#endif

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

#include <algorithm>        // IWYU pragma: export
#include <array>            // IWYU pragma: export
#include <bit>              // IWYU pragma: export
#include <concepts>         // IWYU pragma: export
#include <initializer_list> // IWYU pragma: export core language facility
#include <limits>           // IWYU pragma: export
#include <optional>         // IWYU pragma: export
#include <tuple>            // IWYU pragma: export
#include <type_traits>      // IWYU pragma: export
#include <utility>          // IWYU pragma: export

namespace damp {
using std::array;
using std::bit_cast;
using std::bit_width;
using std::bool_constant;
using std::clamp;
using std::common_type_t;
using std::conditional_t;
using std::constructible_from;
using std::convertible_to;
using std::declval;
using std::false_type;
using std::floating_point;
using std::forward;
using std::get;
using std::index_sequence;
using std::integral;
using std::is_arithmetic_v;
using std::is_base_of_v;
using std::is_const_v;
using std::is_constant_evaluated;
using std::is_constructible_v;
using std::is_floating_point_v;
using std::is_integral_v;
using std::is_same_v;
using std::is_signed_v;
using std::is_unsigned_v;
using std::make_index_sequence;
using std::make_signed_t;
using std::make_tuple;
using std::make_unsigned_t;
using std::max;
using std::min;
using std::move;
using std::nullopt;
using std::nullopt_t;
using std::numeric_limits;
using std::optional;
using std::pair;
using std::remove_const_t;
using std::remove_cv_t;
using std::remove_cvref_t;
using std::remove_pointer_t;
using std::swap;
using std::true_type;
using std::tuple;
using std::type_identity_t;
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

/**
 * @brief Element-wise static_cast of a damp::array to scalar type U (precision rebind).
 */
template<typename U, typename T, size_t N>
[[nodiscard]] constexpr damp::array<U, N> array_as(const damp::array<T, N>& in) {
    damp::array<U, N> out{};
    for (size_t i = 0; i < N; ++i) {
        out[i] = static_cast<U>(in[i]);
    }
    return out;
}

namespace detail {

/**
 * @brief Failed design:: result passed to a runtime law. Not constexpr, so a
 * constant-evaluated constructor becomes ill-formed; at runtime it traps.
 */
[[noreturn]] inline void refuse_failed_design() noexcept {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_trap();
#elif defined(_MSC_VER)
    __fastfail(7);
#else
    for (;;) {
    }
#endif
}

/**
 * @brief A brace initializer with more rows, columns or elements than the Matrix it initializes.
 * Not constexpr, so a constant-evaluated initialization becomes ill-formed; at runtime it traps.
 */
[[noreturn]] inline void refuse_initializer_overflow() noexcept {
    refuse_failed_design();
}

} // namespace detail

} // namespace damp
