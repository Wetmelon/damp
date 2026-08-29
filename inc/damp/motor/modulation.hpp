// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file modulation.hpp
 * @brief Power-electronics modulation (carrier-based VSI duty maps)
 *
 * Maps a commanded voltage (αβ or per-phase) to switch duty cycles for a
 * three-phase voltage-source inverter. Domain-neutral: equally applicable to
 * motor drives, grid-tie inverters, and active front-ends. The reference-frame
 * transforms (Clarke, Park, symmetrical components) that produce the voltage
 * command live in @ref transforms.
 *
 * ## Shipped schemes
 *
 * Continuous (always-switching) carrier PWM via zero-sequence injection:
 * - SPWM — no injection
 * - SVPWM — min-max centering (existing svpwm_zero_sequence)
 * - THIPWM — classical 1/6 third-harmonic injection
 *
 * Discontinuous PWM (DPWM) family — one phase clamped to a rail for 60° or 120°
 * windows (≈33% fewer switch transitions than CPWM at the same carrier):
 * - DPWMMAX / DPWMMIN — always clamp max high / min low
 * - DPWM0 / DPWM1 / DPWM2 — 60° clamp windows (DPWM1 centered on phase peaks)
 * - DPWM3 — 30° split clamp windows
 *
 * Minimal overmodulation hooks: linear SVPWM hexagon limit, modulation index,
 * six-step endpoint.
 *
 * Dead-time (blanking) compensation and inverter SIL plants live with the full
 * motor/power packs when present — not required by this header.
 *
 * All maps are pure @c constexpr free functions (no runtime state). Optional
 * thin Modulator wrapper holds a scheme for FOC/deploy selection.
 *
 * @see A. M. Hava, R. J. Kerkman, T. A. Lipo, "Simple analytical and graphical
 *      methods for carrier-based PWM-VSI drives," IEEE Trans. Power Electron.,
 *      vol. 14, no. 1, pp. 49–61, 1999. doi:10.1109/63.737592
 * @see A. M. Hava, R. J. Kerkman, T. A. Lipo, "A high-performance generalized
 *      discontinuous PWM algorithm," IEEE Trans. Ind. Appl., vol. 34, no. 5,
 *      pp. 1059–1071, 1998. doi:10.1109/28.720446
 * @see D. G. Holmes, T. A. Lipo, "Pulse Width Modulation for Power Converters:
 *      Principles and Practice," IEEE Press / Wiley, 2003, ch. 5–6.
 */

#include <cstddef>
#include <cstdint>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/math/transforms.hpp"
#include "damp/matrix/colvec.hpp"

namespace damp::motor {

/**
 * @brief PWM and space-vector modulation helpers for inverter duty generation
 *
 * @defgroup modulation Power-Electronics Modulation
 * @brief PWM and space-vector modulation helpers for inverter duty generation
 */

/**
 * @brief Carrier-based three-phase VSI modulation scheme
 * @ingroup modulation
 *
 * Selects the zero-sequence (common-mode) injection used by
 * @ref modulation_duty_cycles. Continuous schemes (Spwm, Svpwm, Thipwm) switch
 * every half-bridge every carrier period; discontinuous schemes (Dpwm*) clamp
 * one phase to a DC rail for a fraction of the fundamental cycle.
 *
 * @see Hava et al., IEEE T-PEL 14(1), 1999 (CPWM/DPWM taxonomy)
 * @see Hava et al., IEEE T-IA 34(5), 1998 (DPWM0–3 / GDPWM)
 */
enum class PwmScheme : std::uint8_t {
    Spwm,    ///< Sinusoidal PWM — zero-sequence v₀ = 0
    Svpwm,   ///< Continuous SVPWM — min-max centering
    Thipwm,  ///< Third-harmonic injection v₀ = (V/6) cos 3θ
    DpwmMax, ///< Clamp largest phase high (120°/phase)
    DpwmMin, ///< Clamp smallest phase low (120°/phase)
    Dpwm0,   ///< 60° DPWM, windows lag DPWM1 by 30°
    Dpwm1,   ///< 60° DPWM, clamps centered on phase peaks
    Dpwm2,   ///< 60° DPWM, windows lead DPWM1 by 30°
    Dpwm3,   ///< 30° DPWM (split clamp windows)
};

/**
 * @brief Result of a duty-map: half-bridge duties plus an over-modulation flag
 * @ingroup modulation
 *
 * The duties are always a usable command (clamped to [0, 1]); is_clipped reports
 * whether that clamping engaged, i.e. the requested voltage exceeded the
 * realizable hexagon (continuous schemes) or the scheme’s rail geometry. It is
 * advisory, not a failure — hence a flag rather than a damp::optional wrapper.
 *
 * @tparam T Scalar type
 */
template<typename T = float>
struct SvmDuties {
    ColVec<3, T> duties = {};        ///< [pu] half-bridge duties {a, b, c}, each in [0, 1]
    bool         is_clipped = false; ///< true if any duty was clamped (over-modulation)
};

/**
 * @brief Peak phase voltage at the linear SVPWM hexagon limit
 * @ingroup modulation
 *
 * Continuous SVPWM (and THIPWM) synthesise up to the inscribed circle of the
 * voltage hexagon:
 * @f[
 *   V_{\mathrm{lin}} = \frac{V_{dc}}{\sqrt{3}}
 * @f]
 * (≈15.5% above SPWM’s @f$ V_{dc}/2 @f$). Beyond this, duties hit [0, 1] and
 * @c is_clipped becomes true.
 *
 * @param v_dc DC bus voltage [V]
 * @return Peak phase (αβ magnitude) at the linear limit [V]
 *
 * @see voltage_circle_radius in foc.hpp (same quantity with a modulation fraction)
 * @see Holmes & Lipo, PWM for Power Converters, ch. 5
 */
template<typename T = float>
[[nodiscard]] constexpr T linear_modulation_voltage(T v_dc) {
    return v_dc * damp::numbers::inv_sqrt3_v<T>;
}

/**
 * @brief Six-step (square-wave) fundamental peak phase voltage
 * @ingroup modulation
 *
 * Endpoint of overmodulation Mode II:
 * @f[
 *   V_{1,\mathrm{six\text{-}step}} = \frac{2}{\pi}\,V_{dc}
 * @f]
 *
 * @param v_dc DC bus voltage [V]
 * @return Fundamental peak of classical six-step [V]
 *
 * @see Holmes & Lipo, PWM for Power Converters, ch. 6
 */
template<typename T = float>
[[nodiscard]] constexpr T six_step_fundamental_voltage(T v_dc) {
    return (T{2} * v_dc) * damp::numbers::inv_pi_v<T>;
}

/**
 * @brief Modulation index relative to the linear SVPWM circle
 * @ingroup modulation
 *
 * @f[
 *   M = \frac{\|v_{\alpha\beta}\|}{V_{dc}/\sqrt{3}}
 * @f]
 * @f$ M \le 1 @f$ is the linear range for SVPWM/THIPWM; @f$ M > 1 @f$ is
 * overmodulation (duties will clamp). Returns 0 if @p v_dc is non-positive.
 *
 * @param v_ab αβ voltage command [V] (amplitude-invariant)
 * @param v_dc DC bus voltage [V]
 * @return Modulation index [-]
 */
template<typename T = float>
[[nodiscard]] constexpr T modulation_index(const AlphaBeta<T>& v_ab, T v_dc) {
    if (!(v_dc > T{0})) {
        return T{0};
    }
    const T v_lin = linear_modulation_voltage(v_dc);
    if (!(v_lin > T{0})) {
        return T{0};
    }
    return damp::hypot(v_ab.alpha, v_ab.beta) / v_lin;
}

/**
 * @brief Min-max zero-sequence injection for space-vector PWM
 * @ingroup modulation
 *
 * Returns the common-mode (zero-sequence) offset that, added equally to all
 * three phase references, centres them within the available bus and yields
 * continuous space-vector modulation (SVPWM). This is the offset that extends
 * the linear modulation range by 2/√3 (≈15.5%) over sinusoidal PWM without
 * affecting the line-to-line voltages.
 * @f[
 *   v_0 = -\frac{\max(v_a, v_b, v_c) + \min(v_a, v_b, v_c)}{2}
 * @f]
 *
 * @return Zero-sequence offset to add to every phase [V]
 *
 * @see https://en.wikipedia.org/wiki/Space_vector_modulation
 * @see A. M. Hava, R. J. Kerkman, T. A. Lipo, "Simple analytical and graphical
 *      methods for carrier-based PWM-VSI drives," IEEE Trans. Power Electron.,
 *      vol. 14, no. 1, pp. 49-61, 1999. doi:10.1109/63.737592
 *      (establishes the carrier-based min-max injection ⇔ SVPWM equivalence).
 * @see D. G. Holmes, T. A. Lipo, "Pulse Width Modulation for Power Converters:
 *      Principles and Practice," IEEE Press, 2003, ch. 3.
 */
template<typename T = float>
[[nodiscard]] constexpr T svpwm_zero_sequence(const ColVec<3, T>& v_abc) {
    const auto [v_min, v_max] = damp::minmax({v_abc[0], v_abc[1], v_abc[2]});

    return -(v_max + v_min) / T{2};
}

/**
 * @brief Third-harmonic injection (THIPWM) zero-sequence — 1/6 of fundamental
 * @ingroup modulation
 *
 * Classical THIPWM adds a third harmonic of amplitude @f$ V/6 @f$ to the phase
 * references, where @f$ V = \|v_{\alpha\beta}\| @f$ is the fundamental peak
 * (amplitude-invariant αβ) and @f$ \theta = \mathrm{atan2}(v_\beta, v_\alpha) @f$:
 * @f[
 *   v_0 = \frac{V}{6}\,\cos 3\theta
 * @f]
 * Extends the linear range toward the SVPWM hexagon while remaining continuous
 * (all legs switch every carrier period). Does not depend on @f$ V_{dc} @f$.
 *
 * @param v_ab αβ voltage command [V] (amplitude-invariant)
 * @return Zero-sequence offset [V]
 *
 * @see Holmes & Lipo, PWM for Power Converters, ch. 5 (THI)
 * @see Hava et al., IEEE T-PEL 14(1), 1999
 */
template<typename T = float>
[[nodiscard]] constexpr T thipwm_zero_sequence(const AlphaBeta<T>& v_ab) {
    const T mag = damp::hypot(v_ab.alpha, v_ab.beta);
    if (!(mag > T{0})) {
        return T{0};
    }
    const T th = damp::atan2(v_ab.beta, v_ab.alpha);
    return (mag / T{6}) * damp::cos(T{3} * th);
}

/**
 * @brief DPWM zero-sequence for a chosen discontinuous scheme
 * @ingroup modulation
 *
 * Selects a common-mode offset that clamps one phase to a DC rail:
 * @f[
 *   v_0 =
 *   \begin{cases}
 *     +V_{dc}/2 - v_{\max} & \text{high clamp (max phase duty = 1)} \\
 *     -V_{dc}/2 - v_{\min} & \text{low clamp (min phase duty = 0)}
 *   \end{cases}
 * @f]
 *
 * Scheme selection (Hava DPWM taxonomy):
 * - DpwmMax — always high clamp; DpwmMin — always low clamp
 * - Dpwm1 — high when @f$ |v_{\max}| \ge |v_{\min}| @f$ (60° windows on peaks)
 * - Dpwm0 / Dpwm2 — same as Dpwm1 on angle @f$ \theta \pm 30^\circ @f$ via
 *   @f$ \mathrm{sign}\cos 3(\theta \pm \pi/6) @f$
 * - Dpwm3 — 30° windows via @f$ \mathrm{sign}\cos 6\theta @f$
 *
 * @param v_abc Phase voltages [V] used for min/max rail clamp
 * @param v_dc  DC bus voltage [V] (must be positive for a meaningful clamp)
 * @param scheme Discontinuous scheme (DpwmMax … Dpwm3)
 * @param v_ab  αβ command used for angle-based schemes (Dpwm0/2/3); ignored for Max/Min/1
 * @return Zero-sequence offset [V]
 *
 * @see Hava et al., IEEE T-IA 34(5), 1998 (DPWM0–3)
 * @see Hava et al., IEEE T-PEL 14(1), 1999
 */
template<typename T = float>
[[nodiscard]] constexpr T dpwm_zero_sequence(const ColVec<3, T>& v_abc, T v_dc, PwmScheme scheme, const AlphaBeta<T>& v_ab = {}) {
    const auto [v_min, v_max] = damp::minmax({v_abc[0], v_abc[1], v_abc[2]});
    const T half = v_dc / T{2};

    const auto high = [&]() -> T {
        return half - v_max;
    };
    const auto low = [&]() -> T {
        return -half - v_min;
    };

    switch (scheme) {
        case PwmScheme::DpwmMax:
            return high();
        case PwmScheme::DpwmMin:
            return low();
        case PwmScheme::Dpwm1:
            // Clamp the phase of largest |voltage| to the nearer rail.
            if (damp::abs(v_max) >= damp::abs(v_min)) {
                return high();
            }
            return low();
        case PwmScheme::Dpwm0:
        case PwmScheme::Dpwm2:
        case PwmScheme::Dpwm3: {
            const T th = damp::atan2(v_ab.beta, v_ab.alpha);
            T       arg = th;
            if (scheme == PwmScheme::Dpwm0) {
                arg = th + damp::numbers::pi_v<T> / T{6};
            } else if (scheme == PwmScheme::Dpwm2) {
                arg = th - damp::numbers::pi_v<T> / T{6};
            }
            const T harm = (scheme == PwmScheme::Dpwm3) ? damp::cos(T{6} * th) : damp::cos(T{3} * arg);
            if (harm >= T{0}) {
                return high();
            }
            return low();
        }
        default:
            // Continuous schemes: fall back to SVPWM centering (callers should not hit this).
            return svpwm_zero_sequence(v_abc);
    }
}

/**
 * @brief Zero-sequence for any @ref PwmScheme
 * @ingroup modulation
 *
 * Dispatches to @ref svpwm_zero_sequence, @ref thipwm_zero_sequence, or
 * @ref dpwm_zero_sequence. SPWM is @f$ v_0 = 0 @f$ (no injection).
 *
 * @param v_abc  Phase voltages [V] (SVPWM / DPWM clamp)
 * @param v_ab   αβ command (THIPWM / angle-based DPWM)
 * @param v_dc   DC bus [V] (DPWM rail clamp; ignored by continuous schemes)
 * @param scheme Modulation scheme
 * @return Zero-sequence offset [V]
 */
template<typename T = float>
[[nodiscard]] constexpr T modulation_zero_sequence(const ColVec<3, T>& v_abc, const AlphaBeta<T>& v_ab, T v_dc, PwmScheme scheme) {
    switch (scheme) {
        case PwmScheme::Spwm:
            return T{0};
        case PwmScheme::Svpwm:
            return svpwm_zero_sequence(v_abc);
        case PwmScheme::Thipwm:
            return thipwm_zero_sequence(v_ab);
        case PwmScheme::DpwmMax:
        case PwmScheme::DpwmMin:
        case PwmScheme::Dpwm0:
        case PwmScheme::Dpwm1:
        case PwmScheme::Dpwm2:
        case PwmScheme::Dpwm3:
            return dpwm_zero_sequence(v_abc, v_dc, scheme, v_ab);
    }
    return svpwm_zero_sequence(v_abc);
}

/**
 * @brief Map phase voltages + zero-sequence to clamped half-bridge duties
 * @ingroup modulation
 *
 * @f[
 *   d_x = \mathrm{sat}_{[0,1]}\!\left(\frac{1}{2} + \frac{v_x + v_0}{V_{dc}}\right)
 * @f]
 *
 * @param v_phase Phase voltages [V]
 * @param v_0     Zero-sequence offset [V]
 * @param v_dc    DC bus [V] (must be positive)
 * @return Duties and clip flag
 */
template<typename T = float>
[[nodiscard]] constexpr SvmDuties<T>
duties_from_phase_voltages(const ColVec<3, T>& v_phase, T v_0, T v_dc) {
    const T inv_vdc = T{1} / v_dc;

    SvmDuties<T> result;
    for (size_t i = 0; i < 3; ++i) {
        const T raw = static_cast<T>(0.5) + ((v_phase[i] + v_0) * inv_vdc);
        const T sat = damp::clamp(raw, T{0}, T{1});
        result.is_clipped = result.is_clipped || (sat != raw);
        result.duties[i] = sat;
    }
    return result;
}

/**
 * @brief Carrier-based VSI duty cycles from an αβ voltage command
 * @ingroup modulation
 *
 * Inverse Clarke → scheme-selected zero-sequence → bus-midpoint duty map,
 * clamped to [0, 1]. Same geometry as @ref svm_duty_cycles for
 * @ref PwmScheme::Svpwm.
 *
 * @param v_ab   αβ voltage command [V] (amplitude-invariant)
 * @param v_dc   DC bus voltage [V]
 * @param scheme Modulation scheme (default: continuous SVPWM)
 * @return SvmDuties: {a,b,c} duties in [0, 1] and is_clipped
 *
 * @see svm_duty_cycles for the SVPWM-only entry point used by FOC
 * @see Hava et al., IEEE T-PEL 14(1), 1999
 * @see Holmes & Lipo, PWM for Power Converters, ch. 5–6
 *
 * Example (FOC deploy — swap scheme without touching the current loop):
 * @code
 * const auto Vab = inverse_park_transform(Vdq, theta);
 * const auto out = modulation_duty_cycles(Vab, Vdc, PwmScheme::Dpwm1);
 * // out.duties → timer CCR; out.is_clipped → optional overmod / FW feedback
 * @endcode
 */
template<typename T = float>
[[nodiscard]] constexpr SvmDuties<T>
modulation_duty_cycles(const AlphaBeta<T>& v_ab, T v_dc, PwmScheme scheme = PwmScheme::Svpwm) {
    // Degenerate / reverse bus: no meaningful duty map — park mid-rail (zero voltage).
    if (!(v_dc > T{0})) {
        return SvmDuties<T>{.duties = {static_cast<T>(0.5), static_cast<T>(0.5), static_cast<T>(0.5)}, .is_clipped = true};
    }

    const ColVec<3, T> v_phase = inverse_clarke_transform(v_ab);
    const T            v_0 = modulation_zero_sequence(v_phase, v_ab, v_dc, scheme);
    return duties_from_phase_voltages(v_phase, v_0, v_dc);
}

/**
 * @brief Space-vector PWM duty cycles from an αβ voltage command
 * @ingroup modulation
 *
 * Resolves the αβ voltage command to phase voltages (inverse Clarke), applies
 * min-max zero-sequence injection (svpwm_zero_sequence()) to realise SVPWM, then
 * maps each phase to a half-bridge duty cycle referenced to the Vdc/2 bus
 * midpoint:
 * @f[
 *   d_x = \frac{1}{2} + \frac{v_x + v_0}{V_{dc}}, \qquad x \in \{a, b, c\}
 * @f]
 *
 * Takes an amplitude-invariant αβ command (its magnitude is the peak phase volt,
 * which is what the duty math below assumes); a power-invariant AlphaBeta will not
 * compile here. Results are clamped to [0, 1]; clamping only engages in
 * over-modulation.
 *
 * Equivalent to @c modulation_duty_cycles(v_ab, v_dc, PwmScheme::Svpwm).
 *
 * @param v_ab αβ voltage command [V]
 * @param v_dc DC bus voltage [V]
 * @return SvmDuties: the {a, b, c} duties (each in [0, 1]) and an is_clipped flag.
 *         The duties are always valid (clamped); is_clipped just reports whether
 *         the command fell outside the realizable hexagon (over-modulation).
 *
 * @see modulation_duty_cycles for SPWM / THIPWM / DPWM variants
 * @see https://en.wikipedia.org/wiki/Space_vector_modulation
 * @see A. M. Hava, R. J. Kerkman, T. A. Lipo, "Simple analytical and graphical
 *      methods for carrier-based PWM-VSI drives," IEEE Trans. Power Electron.,
 *      vol. 14, no. 1, pp. 49-61, 1999. doi:10.1109/63.737592
 */
template<typename T = float>
[[nodiscard]] constexpr SvmDuties<T> svm_duty_cycles(const AlphaBeta<T>& v_ab, T v_dc) {
    return modulation_duty_cycles(v_ab, v_dc, PwmScheme::Svpwm);
}

/**
 * @brief Classical six-step (full-wave) duty pattern from the αβ angle
 * @ingroup modulation
 *
 * Overmodulation Mode II endpoint: each phase is held at duty 0 or 1 for 180°
 * of the fundamental, 120° apart. Magnitude of @p v_ab is ignored (direction
 * only); use when the current loop / FW path has already decided on six-step.
 * @f[
 *   d_x = \begin{cases} 1 & \cos(\theta - \phi_x) \ge 0 \\ 0 & \text{otherwise} \end{cases}
 * @f]
 * with @f$ \phi_a = 0 @f$, @f$ \phi_b = 2\pi/3 @f$, @f$ \phi_c = -2\pi/3 @f$.
 *
 * @param v_ab αβ voltage command [V] (angle only; zero vector → mid-rail)
 * @return Six-step duties; is_clipped is always true (non-linear)
 *
 * @see six_step_fundamental_voltage
 * @see Holmes & Lipo, PWM for Power Converters, ch. 6
 */
template<typename T = float>
[[nodiscard]] constexpr SvmDuties<T> six_step_duty_cycles(const AlphaBeta<T>& v_ab) {
    const T mag = damp::hypot(v_ab.alpha, v_ab.beta);
    if (!(mag > T{0})) {
        return SvmDuties<T>{.duties = {static_cast<T>(0.5), static_cast<T>(0.5), static_cast<T>(0.5)}, .is_clipped = true};
    }
    const T th = damp::atan2(v_ab.beta, v_ab.alpha);
    const T p120 = T{2} * damp::numbers::pi_v<T> / T{3};

    SvmDuties<T> result;
    result.is_clipped = true;
    result.duties[0] = (damp::cos(th) >= T{0}) ? T{1} : T{0};
    result.duties[1] = (damp::cos(th - p120) >= T{0}) ? T{1} : T{0};
    result.duties[2] = (damp::cos(th + p120) >= T{0}) ? T{1} : T{0};
    return result;
}

/**
 * @brief Thin scheme-holding wrapper for FOC / deploy paths
 * @ingroup modulation
 *
 * Pure data + @ref duties call-through. No per-tick state. Default scheme is
 * continuous SVPWM (matches @ref svm_duty_cycles / FOController).
 *
 * @tparam T Scalar type (embedded default: float)
 *
 * @code
 * Modulator<float> mod{.scheme = PwmScheme::Dpwm1};
 * const auto out = mod.duties(Vab, Vdc);
 * @endcode
 */
template<typename T = float>
struct Modulator {
    PwmScheme scheme = PwmScheme::Svpwm; ///< Active carrier scheme

    /**
     * @brief Duty map for the stored scheme
     * @param v_ab αβ voltage command [V]
     * @param v_dc DC bus [V]
     * @return Clamped duties + clip flag
     */
    [[nodiscard]] constexpr SvmDuties<T> duties(const AlphaBeta<T>& v_ab, T v_dc) const {
        return modulation_duty_cycles(v_ab, v_dc, scheme);
    }
};

} // namespace damp::motor
