// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file biquad.hpp
 * @brief Biquad and BiquadCascade runtime sections
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/filters/iir_design.hpp"

namespace damp {
/**
 * @brief Second-order IIR (biquad) section runtime.
 *
 * Runs a SecondOrderCoeffs section in Direct Form I:
 *   y[n] = b0·x[n] + b1·x[n-1] + b2·x[n-2] − a1·y[n-1] − a2·y[n-2]
 *
 * Pair with any design:: biquad designer (notch, bandpass, highpass_2nd,
 * peaking, lowshelf, highshelf, lowpass_2nd).
 *
 * @code
 * Biquad<float> notch{design::notch(50.0f, 5.0f, 1.0f / 1000.0f)};
 * float clean = notch(sample);
 * @endcode
 */
template<typename T = float>
class Biquad {
public:
    constexpr Biquad() = default;

    constexpr explicit Biquad(const design::SecondOrderCoeffs<T>& c)
        : b0_(c.b0), b1_(c.b1), b2_(c.b2), a1_(c.a1), a2_(c.a2) {}

    /// Convert across scalar precision, preserving coefficients and delay line.
    template<typename U>
    constexpr explicit Biquad(const Biquad<U>& o)
        : b0_(static_cast<T>(o.b0_)),
          b1_(static_cast<T>(o.b1_)),
          b2_(static_cast<T>(o.b2_)),
          a1_(static_cast<T>(o.a1_)),
          a2_(static_cast<T>(o.a2_)),
          x1_(static_cast<T>(o.x1_)),
          x2_(static_cast<T>(o.x2_)),
          y1_(static_cast<T>(o.y1_)),
          y2_(static_cast<T>(o.y2_)) {}

    /// Process one sample.
    constexpr T operator()(T x) {
        const T y = (b0_ * x) + (b1_ * x1_) + (b2_ * x2_) - (a1_ * y1_) - (a2_ * y2_);
        x2_ = x1_;
        x1_ = x;
        y2_ = y1_;
        y1_ = y;
        return y;
    }

    /// Reset the internal delay line.
    constexpr void reset() {
        x1_ = x2_ = y1_ = y2_ = T{0};
    }

    /// Replace the coefficients in place, keeping the delay line — for adaptive
    /// filters that retune (e.g. grid-frequency tracking) without losing state.
    constexpr void set_coefficients(const design::SecondOrderCoeffs<T>& c) {
        b0_ = c.b0;
        b1_ = c.b1;
        b2_ = c.b2;
        a1_ = c.a1;
        a2_ = c.a2;
    }

    /// Most recent output y[n-1]; read/adjust to seed anti-windup unwinding.
    [[nodiscard]] constexpr T last_output() const { return y1_; }
    constexpr void            set_last_output(T y1) { y1_ = y1; }

private:
    template<typename>
    friend class Biquad;

    T b0_{1}, b1_{0}, b2_{0}, a1_{0}, a2_{0};
    T x1_{0}, x2_{0}, y1_{0}, y2_{0};
};

/**
 * @brief Cascade of second-order sections (SOS) for higher-order IIR filters.
 *
 * Chains NSections biquads in series. Cascading is the numerically preferred
 * realization for higher-order IIR filters (vs. a single high-order section).
 *
 * @tparam NSections Number of biquad sections
 * @tparam T         Scalar type
 */
template<size_t NSections, typename T = float>
class BiquadCascade {
public:
    constexpr BiquadCascade() = default;

    constexpr explicit BiquadCascade(const damp::array<design::SecondOrderCoeffs<T>, NSections>& sections) {
        for (size_t i = 0; i < NSections; ++i) {
            sections_[i] = Biquad<T>(sections[i]);
        }
    }

    /// Process one sample through all sections in series.
    constexpr T operator()(T x) {
        for (size_t i = 0; i < NSections; ++i) {
            x = sections_[i](x);
        }
        return x;
    }

    /// Process only the first @p n_active sections (bypass the rest as identity).
    constexpr T operator()(T x, std::size_t n_active) {
        const std::size_t n = (n_active < NSections) ? n_active : NSections;
        for (std::size_t i = 0; i < n; ++i) {
            x = sections_[i](x);
        }
        return x;
    }

    /// Replace coefficients of section @p i (keeps that section's delay line).
    constexpr void set_section(std::size_t i, const design::SecondOrderCoeffs<T>& c) {
        if (i < NSections) {
            sections_[i].set_coefficients(c);
        }
    }

    /// Access one section (e.g. reset a single biquad after retune).
    [[nodiscard]] constexpr Biquad<T>& section(std::size_t i) { return sections_[i]; }

    [[nodiscard]] constexpr const Biquad<T>& section(std::size_t i) const { return sections_[i]; }

    /// Reset every section.
    constexpr void reset() {
        for (size_t i = 0; i < NSections; ++i) {
            sections_[i].reset();
        }
    }

    /// Reset and load identity (bypass) coefficients into every section.
    constexpr void reset_identity() {
        constexpr design::SecondOrderCoeffs<T> id{T{1}, T{0}, T{0}, T{0}, T{0}};
        for (size_t i = 0; i < NSections; ++i) {
            sections_[i] = Biquad<T>(id);
        }
    }

private:
    damp::array<Biquad<T>, NSections> sections_{};
};
} // namespace damp
