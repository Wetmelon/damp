// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pll.hpp
 * @brief Phase-locked loops for grid and signal synchronization
 */

#include "damp/controllers/pid.hpp"
#include "damp/filters/sogi.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/math/transforms.hpp"

namespace damp {

/**
 * @brief Single-Phase PLL
 *
 * Uses a Second-Order Generalized Integrator (SOGI) for quadrature signal generation.
 *
 * The loop filter is a @ref ContinuousPID acting on the phase error and
 * producing a frequency *offset* about nominal; the inverter-style output clamp
 * plus back-calculation give proper anti-windup (the integrator is bled back out
 * of saturation rather than merely clipped). The phase is a wrapping integrator —
 * @ref damp::wrap_two_pi keeps it in @f$[0,2\pi)@f$ each step so the seam is a full
 * turn (lock about 0 is away from the discontinuity; 180° split-phase is not on
 * an edge). Long runs also stay inside the accurate range of the fast @ref damp
 * math backend.
 *
 * @note A phase-locked loop is a tracker, not a reference-following controller,
 *       and does not satisfy SISOController -- its entry point is
 *       `step(T input, T Ts)` returning `void`, with state read out via member
 *       accessors. Use it as a frequency / phase estimator that feeds *into* a
 *       SISOController (e.g. a current-loop PR controller riding the PLL's
 *       phase estimate), not as a block inside `Cascade<Outer, Inner>`.
 */
template<typename T>
struct SinglePhasePLL {
    /// PI loop filter (continuous gains; `step` supplies measured Ts). Output is
    /// the frequency offset [Hz] about nominal; `u_min` / `u_max` bound that
    /// offset and `Kbc` back-calculates the integrator out of saturation.
    ContinuousPID<T> loop_filter{};

    /// [1/s] Integrator bleed rate (0 = pure integrator). The leak pulls the
    /// accumulated frequency offset back toward zero with time constant 1/leak, so
    /// a transient or biased phase error doesn't park a permanent offset.
    T integrator_leak{};

    constexpr explicit SinglePhasePLL(T Fnom)
        : nominal_frequency(Fnom), frequency_estimate(Fnom) {
        // Frequency may deviate up to ±50% of nominal; that bounds the loop-filter
        // output (the offset about nominal).
        const T max_freq_deviation = Fnom * static_cast<T>(0.5);

        // Sample-rate-independent loop-filter gains: step() applies Ts, so the
        // gains themselves carry no Ts. Kbc = Kp is the standard PI tracking
        // constant (T_t = T_i) for back-calculation anti-windup.
        loop_filter.Kp = T{10};
        loop_filter.Ki = T{100};
        loop_filter.Kbc = loop_filter.Kp;
        loop_filter.u_max = max_freq_deviation;
        loop_filter.u_min = -max_freq_deviation;
    }

    constexpr void step(T input, const T Ts) {
        // Generate bandpass and quadrature signals using SOGI.
        const auto [_, quadrature] = mstogi(input, frequency_estimate, damp::numbers::sqrt2_v<T>, Ts);

        // Phase error from the SOGI quadrature mixer. qv' lags the input by 90° at
        // the SOGI centre, so for an input above the current estimate the DC of
        // input·qv' is negative — the loop must move the estimate by −(input·qv')
        // to chase it (same sign convention as the SOGI-FLL's ω̇ = −Γ·ε·qv', see
        // SogiFll). The earlier +input·qv' was positive feedback and ran the
        // estimate to a frequency rail.
        const T phase_error = -(input * quadrature);

        // Optional leaky integrator: bleed any parked offset toward zero (skip on a
        // non-positive Ts, which would bleed with the wrong sign).
        if (integrator_leak != T{0} && Ts > T{0}) {
            loop_filter.integral -= integrator_leak * loop_filter.integral * Ts;
        }

        // PI loop filter (with anti-windup) on the phase error → frequency offset.
        // r = phase_error, y = 0 ⇒ output = Kp·phase_error + Ki·∫phase_error.
        frequency_estimate = nominal_frequency + loop_filter.control(phase_error, T{0}, Ts);

        // Wrapping integrator: integrate frequency [Hz] to phase [rad] in [0, 2π).
        const T pi = damp::numbers::pi_v<T>;
        phase_estimate = damp::wrap_two_pi(phase_estimate + (T{2} * pi * frequency_estimate * Ts));
    }

    /// Estimated frequency [Hz].
    [[nodiscard]] constexpr T frequency() const {
        return frequency_estimate;
    }

    /// Estimated phase [rad], wrapped to [0, 2π).
    [[nodiscard]] constexpr T phase() const {
        return phase_estimate;
    }

    constexpr void reset() {
        mstogi.reset();
        loop_filter.reset();
        phase_estimate = T{0};
        frequency_estimate = nominal_frequency;
    }

private:
    MSTOGI<T> mstogi{}; // MSTOGI is used to reject DC offset

    T nominal_frequency{};  // Center frequency of the PLL bandpass filter
    T phase_estimate{};     // Estimated phase of the input signal
    T frequency_estimate{}; // Estimated frequency of the input signal
};

/**
 * @brief Synchronous-reference-frame (SRF) PLL for balanced three-phase input
 *
 * The textbook three-phase PLL: Clarke-Park the input onto the current phase
 * estimate, then drive the q-axis projection to zero with a PI loop filter.
 * At lock @f$ v_q = |v|\sin(\varphi - \hat\theta) \to 0 @f$ and @f$ v_d \to +|v| @f$.
 *
 * The loop filter is a @ref ContinuousPID (same proper-anti-windup and
 * wrapping-integrator treatment as @ref SinglePhasePLL).
 *
 * Assumes a balanced set (no negative- or zero-sequence). For unbalanced or
 * distorted grids use @ref DsogiPll, which adds SOGI-based sequence separation
 * ahead of the same loop filter.
 *
 * @note Like @ref SinglePhasePLL this is a tracker, not a SISOController:
 *       its entry point is `step(abc, Ts)` returning `void`, with state read out
 *       via accessors.
 */
template<typename T>
struct ThreePhasePLL {
    /// PI loop filter; output is the frequency offset [Hz] about nominal. See
    /// SinglePhasePLL::loop_filter.
    ContinuousPID<T> loop_filter{};

    /**
     * @param f_nom Nominal frequency [Hz]
     */
    constexpr explicit ThreePhasePLL(T f_nom) : nominal_frequency(f_nom), frequency_estimate(f_nom) {
        const T max_freq_deviation = f_nom * static_cast<T>(0.5);
        loop_filter.Kp = T{10};
        loop_filter.Ki = T{100};
        loop_filter.Kbc = loop_filter.Kp;
        loop_filter.u_max = max_freq_deviation;
        loop_filter.u_min = -max_freq_deviation;
    }

    /**
     * @brief One synchronization step.
     *
     * @param abc [V] Phase voltages (balanced)
     * @param Ts  [s] Sample time
     */
    constexpr void step(const damp::ColVec<3, T>& abc, const T Ts) {
        // q-axis projection on the current phase estimate is the phase error.
        const T phase_error = clarke_park_transform(abc, phase_estimate).q;

        // PI loop filter (with anti-windup) → frequency estimate [Hz].
        frequency_estimate = nominal_frequency + loop_filter.control(phase_error, T{0}, Ts);

        // Wrapping integrator: frequency [Hz] → phase [rad] in [0, 2π).
        const T pi = damp::numbers::pi_v<T>;
        phase_estimate = damp::wrap_two_pi(phase_estimate + (T{2} * pi * frequency_estimate * Ts));
    }

    /// Estimated frequency [Hz].
    [[nodiscard]] constexpr T frequency() const {
        return frequency_estimate;
    }

    /// Estimated phase [rad], wrapped to [0, 2π).
    [[nodiscard]] constexpr T phase() const {
        return phase_estimate;
    }

    constexpr void reset() {
        loop_filter.reset();
        phase_estimate = T{0};
        frequency_estimate = nominal_frequency;
    }

private:
    T nominal_frequency{};  // [Hz] Center frequency of the PLL
    T phase_estimate{};     // [rad] Estimated phase of the input signal
    T frequency_estimate{}; // [Hz] Estimated frequency of the input signal
};

/**
 * @brief Instantaneous positive-sequence αβ from a quadrature signal pair
 * @ingroup filters
 *
 * The αβ-domain analogue of the (phasor) Fortescue transform. Given the αβ vector
 * @p v and its 90°-lagged quadrature @p qv (the `q = e^{-j\pi/2}` operator,
 * supplied by a SOGI/MSTOGI quadrature-signal generator), the instantaneous
 * positive-sequence component is
 * @f[
 *   \begin{bmatrix} v_\alpha^+ \\ v_\beta^+ \end{bmatrix}
 *   = \frac{1}{2}\begin{bmatrix} 1 & -q \\ q & 1 \end{bmatrix}
 *     \begin{bmatrix} v_\alpha \\ v_\beta \end{bmatrix}
 *   = \frac{1}{2}\begin{bmatrix} v_\alpha - q v_\beta \\ q v_\alpha + v_\beta \end{bmatrix}
 * @f]
 *
 * @param v  αβ vector
 * @param qv 90°-lagged quadrature of @p v (per-axis SOGI quadrature output)
 * @return Instantaneous positive-sequence αβ
 *
 * @see https://en.wikipedia.org/wiki/Symmetrical_components
 * @see P. Rodríguez et al., "Decoupled double synchronous reference frame PLL for
 *      power converters control," IEEE Trans. Power Electron., vol. 22, no. 2,
 *      pp. 584-592, 2007. doi:10.1109/TPEL.2006.890000
 */
template<typename T = float>
[[nodiscard]] constexpr AlphaBeta<T> positive_sequence_ab(const AlphaBeta<T>& v, const AlphaBeta<T>& qv) {
    return {static_cast<T>(0.5) * (v.alpha - qv.beta), static_cast<T>(0.5) * (qv.alpha + v.beta)};
}

/**
 * @brief Instantaneous negative-sequence αβ from a quadrature signal pair
 * @ingroup filters
 *
 * Companion to positive_sequence_ab():
 * @f[
 *   \begin{bmatrix} v_\alpha^- \\ v_\beta^- \end{bmatrix}
 *   = \frac{1}{2}\begin{bmatrix} 1 & q \\ -q & 1 \end{bmatrix}
 *     \begin{bmatrix} v_\alpha \\ v_\beta \end{bmatrix}
 *   = \frac{1}{2}\begin{bmatrix} v_\alpha + q v_\beta \\ -q v_\alpha + v_\beta \end{bmatrix}
 * @f]
 *
 * @param v  αβ vector
 * @param qv 90°-lagged quadrature of @p v
 * @return Instantaneous negative-sequence αβ
 */
template<typename T = float>
[[nodiscard]] constexpr AlphaBeta<T> negative_sequence_ab(const AlphaBeta<T>& v, const AlphaBeta<T>& qv) {
    return {static_cast<T>(0.5) * (v.alpha + qv.beta), static_cast<T>(0.5) * (-qv.alpha + v.beta)};
}

/**
 * @brief Dual-SOGI three-phase positive-sequence PLL (DSOGI-PLL)
 *
 * Grid synchronization for three-phase systems under unbalance and
 * distortion. Two quadrature-signal generators (one per αβ axis) feed an
 * instantaneous positive-sequence calculator (positive_sequence_ab()); a
 * synchronous-reference-frame (SRF) PLL then locks to the extracted positive
 * sequence, driving its q-axis projection to zero.
 *
 * Pipeline per step (input is the αβ grid vector, e.g. clarke_transform(v_abc)):
 *   1. SOGI-QSG on @f$ v_\alpha @f$ and @f$ v_\beta @f$ → in-phase + quadrature.
 *   2. Positive-sequence αβ via the `q`-operator matrix.
 *   3. Park to dq on the current phase estimate; the q-axis is the phase error.
 *   4. PI loop filter → frequency, integrated to phase.
 *
 * The loop filter is a @ref ContinuousPID and the phase a wrapping
 * integrator — the same treatment as @ref SinglePhasePLL.
 *
 * The quadrature generator is a template parameter: the default SOGI is adequate
 * for clean inputs; substitute MSTOGI when the αβ signals carry DC offset (its
 * extra TOGI state rejects DC in the quadrature channel).
 *
 * @note Like SinglePhasePLL this is a tracker, not a SISOController: its entry
 *       point is `step(const AlphaBeta&, T Ts)` returning void, with state read
 *       via accessors. The negative-sequence component is also exposed for
 *       unbalance monitoring / sequence-domain control.
 *
 * @tparam T         Scalar type
 * @tparam Resonator Quadrature-signal generator (SOGI or MSTOGI)
 *
 * @see P. Rodríguez et al., "Decoupled double synchronous reference frame PLL for
 *      power converters control," IEEE Trans. Power Electron., 22(2), 2007.
 */
template<typename T = float, template<typename> class Resonator = SOGI>
class DsogiPll {
public:
    /// PI loop filter; output is the frequency offset [Hz] about nominal. See
    /// SinglePhasePLL::loop_filter.
    ContinuousPID<T> loop_filter{};

    /// SOGI damping gain (√2 ≈ unity-Q).
    T sogi_gain{damp::numbers::sqrt2_v<T>};

    /**
     * @param f_nom Nominal grid frequency [Hz]
     */
    constexpr explicit DsogiPll(T f_nom)
        : nominal_frequency(f_nom), frequency_estimate(f_nom) {
        const T max_freq_deviation = f_nom * static_cast<T>(0.5);
        loop_filter.Kp = T{10};
        loop_filter.Ki = T{100};
        loop_filter.Kbc = loop_filter.Kp;
        loop_filter.u_max = max_freq_deviation;
        loop_filter.u_min = -max_freq_deviation;
    }

    /// One synchronization step from the stationary-frame αβ grid vector.
    constexpr void step(const AlphaBeta<T>& v_ab, const T Ts) {
        // Per-axis quadrature signal generation at the current frequency estimate.
        const auto [va, qva] = sogi_alpha(v_ab.alpha, frequency_estimate, sogi_gain, Ts);
        const auto [vb, qvb] = sogi_beta(v_ab.beta, frequency_estimate, sogi_gain, Ts);

        const AlphaBeta<T> v = {va, vb};
        const AlphaBeta<T> qv = {qva, qvb};
        positive_ab_ = positive_sequence_ab(v, qv);
        negative_ab_ = negative_sequence_ab(v, qv);

        // SRF-PLL: q-axis projection of the positive sequence is the phase error
        // (v_q = |v⁺|·sin(φ − θ̂) → 0 at lock, with v_d → +|v⁺|).
        positive_dq_ = park_transform(positive_ab_, phase_estimate);
        const T phase_error = positive_dq_.q;

        // PI loop filter (with anti-windup) → frequency estimate [Hz].
        frequency_estimate = nominal_frequency + loop_filter.control(phase_error, T{0}, Ts);

        // Wrapping integrator: frequency [Hz] → phase [rad] in [0, 2π).
        const T pi = damp::numbers::pi_v<T>;
        phase_estimate = damp::wrap_two_pi(phase_estimate + (T{2} * pi * frequency_estimate * Ts));
    }

    /// Estimated frequency [Hz].
    [[nodiscard]] constexpr T frequency() const { return frequency_estimate; }
    /// Estimated phase [rad], wrapped to [0, 2π).
    [[nodiscard]] constexpr T                   phase() const { return phase_estimate; }
    [[nodiscard]] constexpr AlphaBeta<T>        positive_sequence() const { return positive_ab_; }
    [[nodiscard]] constexpr AlphaBeta<T>        negative_sequence() const { return negative_ab_; }
    [[nodiscard]] constexpr DirectQuadrature<T> positive_dq() const { return positive_dq_; }

    constexpr void reset() {
        sogi_alpha = {};
        sogi_beta = {};
        loop_filter.reset();
        phase_estimate = T{0};
        frequency_estimate = nominal_frequency;
        positive_ab_ = {};
        negative_ab_ = {};
        positive_dq_ = {};
    }

private:
    Resonator<T> sogi_alpha{};
    Resonator<T> sogi_beta{};

    T nominal_frequency{};
    T phase_estimate{};
    T frequency_estimate{};

    AlphaBeta<T>        positive_ab_{};
    AlphaBeta<T>        negative_ab_{};
    DirectQuadrature<T> positive_dq_{};
};

} // namespace damp
