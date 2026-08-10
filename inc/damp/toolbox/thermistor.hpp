// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file thermistor.hpp
 * @brief NTC thermistor linearization (resistance → temperature).
 *
 * Splits the NTC conversion into a @ref damp::design "design" step
 * — which turns calibration data (Beta model or Steinhart-Hart coefficients)
 * into a small set of fitted coefficients — and a `constexpr`, allocation-free
 * runtime block (@ref damp::Thermistor) that evaluates them one
 * sample at a time.
 *
 * @see filters/blocks.hpp for the rest of the everyday signal-conditioning blocks.
 */

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"

namespace damp {

namespace design {

/**
 * @struct ThermistorCoeffs
 * @brief Fitted NTC coefficients in Steinhart-Hart form.
 *
 * Both the Beta model and a full Steinhart-Hart fit reduce to the same
 * evaluation, 1/T = a + b·ln(R) + c·ln(R)³, so a single coefficient triple
 * (a, b, c) describes either. The Beta model maps onto it with c = 0.
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
struct ThermistorCoeffs {
    T a{}; ///< Steinhart-Hart a.
    T b{}; ///< Steinhart-Hart b.
    T c{}; ///< Steinhart-Hart c (0 for a Beta-model fit).

    /// Re-cast the coefficients to a different scalar type.
    template<typename U>
    [[nodiscard]] constexpr ThermistorCoeffs<U> as() const {
        return ThermistorCoeffs<U>{
            static_cast<U>(a), static_cast<U>(b), static_cast<U>(c)
        };
    }
};

/**
 * @brief Fit NTC coefficients from the Beta-parameter model.
 *
 * The Beta model needs one calibration point plus the material constant β:
 *
 *     1/T = 1/T₀ + (1/β)·ln(R/R₀)
 *
 * where T₀ is a reference temperature [K] at which the thermistor reads R₀.
 * Simple and adequate over a moderate span. This maps onto the Steinhart-Hart
 * form with c = 0, a = 1/T₀ − (1/β)·ln(R₀), b = 1/β.
 *
 * @param r0   Reference resistance [Ω] at temperature t0.
 * @param t0_K Reference temperature [K] (e.g. 298.15 for 25 °C).
 * @param beta Material constant β [K].
 * @return ThermistorCoeffs in Steinhart-Hart form.
 *
 * @code
 * // 10 kΩ NTC, β = 3950, 25 °C reference:
 * const auto coeffs = design::beta(10000.0f, 298.15f, 3950.0f);
 * @endcode
 */
template<typename T = float>
[[nodiscard]] constexpr ThermistorCoeffs<T> beta(T r0, T t0_K, T beta) {
    const T b = T{1} / beta;
    return ThermistorCoeffs<T>{
        (T{1} / t0_K) - (b * damp::log(r0)), // a
        b,                                   // b
        T{0}                                 // c
    };
}

/**
 * @brief Fit the Steinhart-Hart coefficients from three calibration points.
 *
 * The standard high-accuracy NTC fit across a wide range:
 *
 *     1/T = a + b·ln(R) + c·ln(R)³
 *
 * Given three (resistance, temperature) measurements the three coefficients are
 * the unique solution of the linear system
 *
 *     ⎡1  ln(R₁)  ln(R₁)³⎤ ⎡a⎤   ⎡1/T₁⎤
 *     ⎢1  ln(R₂)  ln(R₂)³⎥ ⎢b⎥ = ⎢1/T₂⎥
 *     ⎣1  ln(R₃)  ln(R₃)³⎦ ⎣c⎦   ⎣1/T₃⎦
 *
 * solved here with the library's mat::solve. Pick the calibration points
 * spread across the operating range — e.g. low, mid, high — for the best
 * wide-range accuracy. Resistances must be distinct (singular fit → zero coeffs).
 *
 * @param p1 Calibration point {resistance [Ω], temperature [K]}.
 * @param p2 Calibration point {resistance [Ω], temperature [K]}.
 * @param p3 Calibration point {resistance [Ω], temperature [K]}.
 * @return Fitted ThermistorCoeffs.
 *
 * @code
 * // Three datasheet points: 32650 Ω @ 0 °C, 10000 Ω @ 25 °C, 3603 Ω @ 50 °C.
 * const auto coeffs = design::steinhart_hart({32650.0, 273.15},
 *                                            {10000.0, 298.15},
 *                                            { 3603.0, 323.15});
 * @endcode
 */
template<typename T = float>
[[nodiscard]] constexpr ThermistorCoeffs<T>
steinhart_hart(damp::pair<T, T> p1, damp::pair<T, T> p2, damp::pair<T, T> p3) {
    const T l1 = damp::log(p1.first);
    const T l2 = damp::log(p2.first);
    const T l3 = damp::log(p3.first);

    // Solve [a b c]ᵀ from the 3×3 system  V·x = y,  rows [1, ln(R), ln(R)³].
    const Matrix<3, 3, T> V{
        {T{1}, l1, l1 * l1 * l1},
        {T{1}, l2, l2 * l2 * l2},
        {T{1}, l3, l3 * l3 * l3},
    };

    const Matrix<3, 1, T> y{
        {T{1} / p1.second},
        {T{1} / p2.second},
        {T{1} / p3.second},
    };

    const auto x = mat::solve(V, y).value_or(Matrix<3, 1, T>{});

    return ThermistorCoeffs<T>{x(0, 0), x(1, 0), x(2, 0)};
}

} // namespace design

/**
 * @brief NTC thermistor linearization (resistance → temperature).
 *
 * Evaluates a fitted set of @ref design::ThermistorCoeffs to convert a measured
 * thermistor resistance to temperature, using the Steinhart-Hart relation
 * 1/T = a + b·ln(R) + c·ln(R)³ (the Beta model is the special case c = 0).
 *
 * Temperatures are handled in kelvin internally; @ref celsius returns Celsius.
 * Either give the Steinhart-Hart coefficients directly, or obtain them from
 * @ref design::beta (Beta model) or @ref design::steinhart_hart (three-point fit).
 *
 * @code
 * // 10 kΩ NTC, β = 3950, 25 °C reference:
 * const Thermistor<float> ntc{design::beta(10000.0f, 298.15f, 3950.0f)};
 * float temp_c = ntc.celsius(measured_ohms);
 *
 * // Or straight from datasheet Steinhart-Hart coefficients:
 * const Thermistor<float> ntc2{1.1e-3f, 2.3e-4f, 8.7e-8f};
 * @endcode
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
class Thermistor {
public:
    constexpr Thermistor() = default;

    /// Construct directly from Steinhart-Hart coefficients (c = 0 for a Beta model).
    constexpr Thermistor(T a, T b, T c) : a_(a), b_(b), c_(c) {}

    /// Construct from fitted coefficients (see @ref design::beta /
    /// @ref design::steinhart_hart).
    constexpr explicit Thermistor(const design::ThermistorCoeffs<T>& coeffs)
        : a_(coeffs.a), b_(coeffs.b), c_(coeffs.c) {}

    /// Temperature in kelvin for a measured resistance [Ω].
    [[nodiscard]] constexpr T kelvin(T resistance) const {
        const T ln_r = damp::log(resistance);
        const T inv_t = a_ + (b_ * ln_r) + (c_ * ln_r * ln_r * ln_r);
        return T{1} / inv_t;
    }

    /// Temperature in degrees Celsius for a measured resistance [Ω].
    [[nodiscard]] constexpr T celsius(T resistance) const {
        return kelvin(resistance) - static_cast<T>(273.15);
    }

private:
    T a_{0};
    T b_{0};
    T c_{0};
};

} // namespace damp
