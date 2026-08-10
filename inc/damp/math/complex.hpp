// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file complex.hpp
 * @brief Complex number type for Damp math
 */

#include "math.hpp"

// std::complex interop is hosted-only: convenience conversions to/from
// std::complex. Guarded so the freestanding core (no <complex>) still provides
// the full damp::complex type. No C99 `_Complex` interop — use damp::complex or
// std::complex. Do not pull std::complex_literals into the global namespace.
#if __has_include(<complex>)
#define DAMP_HAS_STD_COMPLEX 1
#include <complex>
#else
#define DAMP_HAS_STD_COMPLEX 0
#endif

namespace damp {

/**
 * @brief Constexpr complex number class for compile-time computations
 * @ingroup linear_algebra
 *
 * Lightweight complex number implementation with full constexpr support.
 * Used for eigenvalue computations and complex matrix operations.
 *
 * @tparam T Underlying floating-point type (e.g., float, double)
 */
template<typename T>
struct complex {
    T real_, imag_;

    // Value type

    /**
     * @brief Default constructor, constructs complex number (r + i*i)
     * @param r Real part (default: 0)
     * @param i Imaginary part (default: 0)
     */
    constexpr complex(T r = T{}, T i = T{}) : real_(r), imag_(i) {}

    /**
     * @brief Type conversion constructor
     * @tparam U Source element type
     * @param other Complex number to convert from
     */
    template<typename U>
    constexpr complex(const complex<U>& other) : real_(static_cast<T>(other.real_)), imag_(static_cast<T>(other.imag_)) {}

#if DAMP_HAS_STD_COMPLEX
    /**
     * @brief Constructor from std::complex (hosted only)
     * @tparam U Source element type
     * @param other std::complex number to convert from
     */
    template<typename U>
    constexpr complex(const std::complex<U>& other) : real_(static_cast<T>(other.real())), imag_(static_cast<T>(other.imag())) {}
#endif

    /**
     * @brief Binary addition operator
     * @param other Complex number to add
     * @return Sum of this and other
     */
    constexpr complex operator+(const complex& other) const {
        return {real_ + other.real_, imag_ + other.imag_};
    }

    /**
     * @brief Binary subtraction operator
     * @param other Complex number to subtract
     * @return Difference of this and other
     */
    constexpr complex operator-(const complex& other) const {
        return {real_ - other.real_, imag_ - other.imag_};
    }

    /**
     * @brief Unary negation operator
     * @return Negated complex number
     */
    constexpr complex operator-() const {
        return {-real_, -imag_};
    }

    /**
     * @brief Binary multiplication operator
     * @param other Complex number to multiply
     * @return Product of this and other
     */
    constexpr complex operator*(const complex& other) const {
        return {(real_ * other.real_) - (imag_ * other.imag_), (real_ * other.imag_) + (imag_ * other.real_)};
    }

    /**
     * @brief Binary division operator
     * @param other Complex number to divide by
     * @return Quotient of this and other
     */
    constexpr complex operator/(const complex& other) const {
        T denom = (other.real_ * other.real_) + (other.imag_ * other.imag_);
        return {((real_ * other.real_) + (imag_ * other.imag_)) / denom, ((imag_ * other.real_) - (real_ * other.imag_)) / denom};
    }

    /**
     * @brief Compound addition operator
     * @param other Complex number to add
     * @return Reference to this
     */
    constexpr complex& operator+=(const complex& other) {
        real_ += other.real_;
        imag_ += other.imag_;
        return *this;
    }

    /**
     * @brief Compound subtraction operator
     * @param other Complex number to subtract
     * @return Reference to this
     */
    constexpr complex& operator-=(const complex& other) {
        real_ -= other.real_;
        imag_ -= other.imag_;
        return *this;
    }

    /**
     * @brief Compound multiplication operator
     * @param other Complex number to multiply
     * @return Reference to this
     */
    constexpr complex& operator*=(const complex& other) {
        T new_real = (real_ * other.real_) - (imag_ * other.imag_);
        T new_imag = (real_ * other.imag_) + (imag_ * other.real_);
        real_ = new_real;
        imag_ = new_imag;
        return *this;
    }

    /**
     * @brief Compound division operator
     * @param other Complex number to divide by
     * @return Reference to this
     */
    constexpr complex& operator/=(const complex& other) {
        T denom = (other.real_ * other.real_) + (other.imag_ * other.imag_);
        T new_real = ((real_ * other.real_) + (imag_ * other.imag_)) / denom;
        T new_imag = ((imag_ * other.real_) - (real_ * other.imag_)) / denom;
        real_ = new_real;
        imag_ = new_imag;
        return *this;
    }

    /**
     * @brief Scalar addition operator
     * @param scalar Real scalar to add
     * @return Sum of this and scalar
     */
    constexpr complex operator+(T scalar) const {
        return {real_ + scalar, imag_};
    }

    /**
     * @brief Scalar subtraction operator
     * @param scalar Real scalar to subtract
     * @return Difference of this and scalar
     */
    constexpr complex operator-(T scalar) const {
        return {real_ - scalar, imag_};
    }

    /**
     * @brief Scalar multiplication operator
     * @param scalar Real scalar to multiply
     * @return Product of this and scalar
     */
    constexpr complex operator*(T scalar) const {
        return {real_ * scalar, imag_ * scalar};
    }

    /**
     * @brief Scalar division operator
     * @param scalar Real scalar to divide by
     * @return Quotient of this and scalar
     */
    constexpr complex operator/(T scalar) const {
        return {real_ / scalar, imag_ / scalar};
    }

    /**
     * @brief Get real part
     * @return Real component
     */
    constexpr T real() const {
        return real_;
    }

    /**
     * @brief Get imaginary part
     * @return Imaginary component
     */
    constexpr T imag() const {
        return imag_;
    }

    /**
     * @brief Compute squared magnitude (norm squared)
     * @return |z|² = real² + imag²
     */
    constexpr T norm() const {
        return (real_ * real_) + (imag_ * imag_);
    }

    /**
     * @brief Compute complex conjugate
     * @return Complex conjugate (real - i*imag)
     */
    constexpr complex conj() const {
        return {real_, -imag_};
    }

    /**
     * @brief Compute magnitude (absolute value)
     *
     * @return |z| = sqrt(real² + imag²)
     */
    constexpr T abs() const {
        return sqrt(norm());
    }

    /**
     * @brief Equality comparison operator
     * @param other Complex number to compare
     * @return True if real and imaginary parts are equal
     */
    constexpr bool operator==(const complex& other) const {
        return real_ == other.real_ && imag_ == other.imag_;
    }

    /**
     * @brief Inequality comparison operator
     * @param other Complex number to compare
     * @return True if real or imaginary parts differ
     */
    constexpr bool operator!=(const complex& other) const {
        return !(*this == other);
    }
};

/**
 * @brief Scalar-complex multiplication (scalar * complex)
 * @tparam T Underlying floating-point type
 * @param scalar Real scalar multiplier
 * @param c Complex number
 * @return Product scalar * c
 */
template<typename T>
constexpr complex<T> operator*(T scalar, const complex<T>& c) {
    return {scalar * c.real_, scalar * c.imag_};
}

/**
 * @brief Scalar-complex multiplication (complex * scalar)
 * @tparam T Underlying floating-point type
 * @param c Complex number
 * @param scalar Real scalar multiplier
 * @return Product c * scalar
 */
template<typename T>
constexpr complex<T> operator*(const complex<T>& c, T scalar) {
    return {c.real_ * scalar, c.imag_ * scalar};
}

/**
 * @brief Scalar-complex addition (scalar + complex)
 * @tparam T Underlying floating-point type
 * @param scalar Real scalar to add
 * @param c Complex number
 * @return Sum scalar + c
 */
template<typename T>
constexpr complex<T> operator+(T scalar, const complex<T>& c) {
    return {scalar + c.real_, c.imag_};
}

/**
 * @brief Scalar-complex subtraction (scalar - complex)
 * @tparam T Underlying floating-point type
 * @param scalar Real scalar
 * @param c Complex number to subtract
 * @return Difference scalar - c
 */
template<typename T>
constexpr complex<T> operator-(T scalar, const complex<T>& c) {
    return {scalar - c.real_, -c.imag_};
}

/**
 * @brief Scalar-complex division (scalar / complex)
 * @tparam T Underlying floating-point type
 * @param scalar Real scalar numerator
 * @param c Complex number denominator
 * @return Quotient scalar / c
 */
template<typename T>
constexpr complex<T> operator/(T scalar, const complex<T>& c) {
    T denom = c.real_ * c.real_ + c.imag_ * c.imag_;
    return {scalar * c.real_ / denom, -scalar * c.imag_ / denom};
}

/**
 * @brief Compute complex square root (constexpr)
 *
 * Computes the principal square root of a complex number z.
 * Uses the formula: sqrt(z) = sqrt((|z| + Re(z))/2) + i*sign(Im(z))*sqrt((|z| - Re(z))/2)
 *
 * @tparam T Underlying floating-point type
 * @param z  Complex number to compute square root of
 *
 * @return Principal square root of z
 */
template<typename T>
constexpr complex<T> sqrt(const complex<T>& z) {
    T re = z.real();
    T im = z.imag();

    if (re == T{0} && im == T{0}) {
        return complex<T>{T{0}, T{0}};
    }

    T mag = damp::hypot(re, im);
    T r = sqrt((mag + re) / T{2});
    T i = sqrt((mag - re) / T{2});

    // Sign of imaginary part matches sign of input imaginary part
    if (im < T{0}) {
        i = -i;
    }

    return complex<T>{r, i};
}

/**
 * @brief Compute magnitude (absolute value) of a complex number
 * @tparam T Underlying floating-point type
 * @param z Complex number
 * @return |z| = sqrt(real² + imag²)
 */
template<typename T>
constexpr T abs(const complex<T>& z) {
    return sqrt(z.norm());
}

/**
 * @brief Compute argument (phase angle) of a complex number
 * @tparam T Underlying floating-point type
 * @param z Complex number
 * @return Angle in radians from positive real axis
 */
template<typename T>
constexpr T arg(const complex<T>& z) {
    return atan2(z.imag(), z.real());
}

template<typename T>
constexpr T conj(const T& value) {
    return value;
}

template<typename T>
constexpr complex<T> conj(const complex<T>& z) {
    return complex<T>{z.real(), -z.imag()};
}

template<typename T>
constexpr T real(const complex<T>& z) {
    return z.real();
}

template<typename T>
constexpr T imag(const T& value) {
    (void)value; // Suppress unused parameter warning
    return T{0};
}

template<typename T>
constexpr T imag(const complex<T>& z) {
    return z.imag();
}

} // namespace damp