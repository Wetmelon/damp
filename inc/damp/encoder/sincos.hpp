// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file damp/encoder/sincos.hpp
 * @brief Sine/cosine position decoding for resolvers, sin/cos encoders and linear scales.
 *
 * Resolvers, analog (1 Vpp) sin/cos encoders and sin/cos linear scales all deliver
 * two quadrature channels @f$ (\sin\phi,\ \cos\phi) @f$ whose phase @f$ \phi @f$
 * is the position within one signal period (one turn of a 1-speed resolver, one line
 * of an encoder, one pitch of a scale). Decoding is two independent steps:
 *
 * 1. Fine phase — @ref damp::decode_sincos trims the channels and takes
 *    @f$ \operatorname{atan2} @f$: the fraction of a period, plus the signal amplitude
 *    for loss-of-signal detection. Stateless.
 * 2. Period count — which period the axis is in, kept as a @ref damp::SplitReal "SplitReal":
 *    - @ref damp::combine_quadrature merges a hardware quadrature counter (the comparator
 *      square waves on a timer, 4 counts per period) with the fine phase. Works at any
 *      speed; the usual choice for multi-period encoders and scales.
 *    - @ref damp::PhaseUnwrapper (encoder/unwrap.hpp) counts periods from successive phases alone. Valid while the
 *      axis moves less than half a period between samples — fine for a resolver
 *      (one period per turn), not for a 2048-line encoder at speed.
 *
 * The result is a measurement. Filter it with an EncoderTracker (or a Kalman filter)
 * for a smooth angle and speed between samples:
 *
 * @code
 * const auto s   = decode_sincos(adc_sin, adc_cos, cal);    // fine phase
 * const auto pos = combine_quadrature(timer_count, s.frac); // absolute, in periods
 * tracker.predict(dt);
 * tracker.update(pos, dt);
 * @endcode
 *
 * A resolver's excitation carrier is assumed demodulated upstream (synchronous
 * rectification / S&H at the carrier peak), the usual RDC front end; sin/cos encoders
 * and scales are baseband already.
 *
 * `constexpr`, allocation-free, embeddable.
 *
 * @see Analog Devices AN-263 / the AD2S resolver-to-digital converter tracking loop.
 * @see Heidenhain, "Interfaces of HEIDENHAIN Encoders" — 1 Vpp sin/cos signals and
 *      interpolation.
 */

#include "damp/backend.hpp"
#include "damp/math/math.hpp"       // damp::atan2, damp::hypot, damp::floor, damp::wrap
#include "damp/math/split_real.hpp" // SplitReal

namespace damp {

/**
 * @brief Per-channel gain/offset trims for the sin/cos front end (hardware calibration).
 *
 * Corrects the channels so the trimmed pair traces a unit circle: subtract the DC
 * offset, then scale by the gain. Commission once against a slow reference move (fit
 * the offsets to the channel midpoints and the gains to equalize the amplitudes).
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
struct SinCosCalibration {
    T sin_gain{T{1}};   ///< [-] sine-channel scale (equalizes amplitude)
    T cos_gain{T{1}};   ///< [-] cosine-channel scale
    T sin_offset{T{0}}; ///< sine-channel DC offset to subtract
    T cos_offset{T{0}}; ///< cosine-channel DC offset to subtract
};

/**
 * @brief One decoded sin/cos sample: the fine phase and the signal amplitude.
 * @tparam T Scalar type.
 */
template<typename T = float>
struct SinCosSample {
    T frac{T{0}};      ///< [periods] phase within the current period, in [−0.5, 0.5)
    T amplitude{T{0}}; ///< [-] √(sin²+cos²) of the trimmed channels; ≈1 when healthy, → 0 on signal loss
};

/**
 * @brief Trim a (sin, cos) pair and decode its phase.
 *
 * @param sin_raw sine-channel reading
 * @param cos_raw cosine-channel reading
 * @param cal     channel gain/offset trims
 * @return the phase in periods (atan2 / 2π, in [−0.5, 0.5)) and the trimmed amplitude
 */
template<typename T>
[[nodiscard]] constexpr SinCosSample<T> decode_sincos(T sin_raw, T cos_raw, const SinCosCalibration<T>& cal = {}) {
    constexpr T inv_two_pi = T{1} / (T{2} * damp::numbers::pi_v<T>);
    constexpr T half = static_cast<T>(0.5);
    const T     s = (sin_raw - cal.sin_offset) * cal.sin_gain;
    const T     c = (cos_raw - cal.cos_offset) * cal.cos_gain;
    // atan2 is in (−π, π]; wrap maps the single +π point onto −0.5.
    return {.frac = damp::wrap(damp::atan2(s, c) * inv_two_pi, -half, half), .amplitude = damp::hypot(s, c)};
}

/**
 * @brief Absolute position from a quadrature edge counter and the fine phase.
 *
 * The comparator square waves of the sin and cos channels clock a quadrature counter:
 * four counts per period, at any speed. The counter alone resolves a quarter period;
 * the fine phase alone does not know which period. Rounding
 * @f$ \lfloor q/4 \rfloor + (q \bmod 4)/4 + 1/8 - \phi @f$ to an integer gives the period,
 * and the fine phase supplies the position inside it; the 1/8 is the mean lag of q/4 behind
 * the position. The rounding absorbs comparator edges that switch early or late relative to
 * the ADC sample, up to ±3/8 period.
 *
 * The counter must be zeroed so that count @c q means position ≈ q/4 periods, i.e.
 * at a phase of 0 (reset it on the index pulse, or at commissioning).
 *
 * @param quarter_count quadrature counter value (4 counts per period), sign-extended
 * @param frac          [periods] fine phase from @ref decode_sincos, in [−0.5, 0.5)
 * @return absolute position in periods
 */
template<typename T>
[[nodiscard]] constexpr SplitReal<T> combine_quadrature(int64_t quarter_count, T frac) {
    // Floor division by 4 (C++ division truncates toward zero; fix up negative counts).
    int64_t q_whole = quarter_count / 4;
    int64_t q_rem = quarter_count % 4;
    if (q_rem < 0) {
        q_rem += 4;
        --q_whole;
    }
    const T coarse_frac = static_cast<T>(q_rem) * static_cast<T>(0.25); // 0, ¼, ½, ¾
    // An ideal count is q = ⌊4p⌋, so q/4 trails the position by 1/8 period on average: the
    // extra 1/8 centres the rounding on it, leaving 3/8 period for early and for late edges.
    const T n = damp::floor(coarse_frac - frac + static_cast<T>(0.625)); // 0 or 1 for an ideal count
    return {q_whole + static_cast<int64_t>(n), frac};
}

} // namespace damp
