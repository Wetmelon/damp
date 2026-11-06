// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file sogi.hpp
 * @brief Second-order generalized integrator (SOGI) filters
 */

#include "damp/backend.hpp"
#include "damp/design/model_result.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {

namespace design {

namespace detail {

/** @brief A resonator needs f0 > 0 and alpha > 0, both finite. */
template<typename T>
[[nodiscard]] constexpr bool valid_resonator(T f0, T alpha) {
    return (f0 > T{0}) && (alpha > T{0}) && damp::isfinite(f0) && damp::isfinite(alpha);
}

/** @brief The exact ZOH model at Ts of a continuous design result, or the reason it is refused. */
template<size_t NX, typename T>
[[nodiscard]] constexpr StateSpaceResult<NX, 1, 2, T> sampled(const StateSpaceResult<NX, 1, 2, T>& continuous, T Ts) {
    if (!(Ts > T{0}) || !damp::isfinite(Ts)) {
        return ::damp::detail::refused<StateSpaceResult<NX, 1, 2, T>>(DesignStatus::invalid_sample_time);
    }
    if (!continuous.success) {
        return continuous;
    }
    const auto sys_d = discretize(continuous.sys, Ts, DiscretizationMethod::ZOH);
    if (!sys_d) {
        return ::damp::detail::refused<StateSpaceResult<NX, 1, 2, T>>(DesignStatus::singular);
    }
    return designed(*sys_d);
}

} // namespace detail

/**
 * @brief Second-Order Generalized Integrator (SOGI) design
 *
 * SOGI provides bandpass filtering for grid synchronization.
 * Produces quadrature signals (90° phase shift).
 *
 * Transfer functions, ω₀ = 2π f0:
 * H_bp(s) = (alpha*ω₀*s) / (s² + alpha*ω₀*s + ω₀²)    [bandpass]
 * H_q(s) = (alpha*ω₀²) / (s² + alpha*ω₀*s + ω₀²)     [quadrature]
 * H_notch(s) = (s² + ω₀²) / (s² + alpha*ω₀*s + ω₀²) [notch]
 *
 * @param f0 Fundamental frequency [Hz] (> 0)
 * @param alpha Damping gain (typically 1.0-2.0, > 0)
 * @return StateSpaceResult with the SOGI system; refused (invalid_parameter) if f0 or alpha is not positive
 */
template<typename T = double>
[[nodiscard]] constexpr StateSpaceResult<2, 1, 2, T> sogi(T f0, T alpha = damp::numbers::sqrt2_v<T>) {
    if (!detail::valid_resonator(f0, alpha)) {
        return ::damp::detail::refused<StateSpaceResult<2, 1, 2, T>>(DesignStatus::invalid_parameter);
    }
    const T w0 = T{2} * damp::numbers::pi_v<T> * f0;
    return detail::designed(StateSpace<2, 1, 2, T>{
        .A = Matrix<2, 2, T>{
            {-alpha * w0, -w0},
            {w0, T{0}},
        },

        .B = Matrix<2, 1, T>{
            {alpha * w0},
            {T{0}},
        },

        .C = Matrix<2, 2, T>{
            {T{1}, T{0}},
            {T{0}, T{1}},
        },

        .D = Matrix<2, 1, T>::zeros(),
    });
}

/**
 * @brief Second-Order Generalized Integrator (SOGI) design (discrete-time)
 *
 * Exact ZOH discretization of the continuous SOGI for simulation/analysis at a
 * fixed sample rate. (The runtime @ref SOGI class rebuilds its own resonator each
 * tick from the live frequency and does not use this.)
 *
 * @param f0 Fundamental frequency [Hz]
 * @param alpha Damping gain (typically 1.0-2.0)
 * @param Ts Sample time [s]
 * @return StateSpaceResult with the discrete SOGI system; refused with invalid_sample_time if
 *         Ts ≤ 0, invalid_parameter for f0 / alpha, singular if the discretization fails
 */
template<typename T = float>
[[nodiscard]] constexpr StateSpaceResult<2, 1, 2, T> sogi(T f0, T alpha, T Ts) {
    return detail::sampled(sogi(f0, alpha), Ts);
}

/**
 * @brief Mixed Second/Third-Order Generalized Integrator (MSTOGI)
 *
 * A standard SOGI-QSG augmented with a Third-Order Generalized Integrator
 * (TOGI) that estimates and removes the DC / offset component the plain
 * SOGI-QSG quadrature output would otherwise pass. The plain `qv′` channel of a
 * SOGI-QSG is a low-pass that has unity gain at DC; a biased or offset grid
 * signal therefore corrupts the quadrature estimate. The MSTOGI subtracts a
 * co-tuned TOGI estimate so the quadrature channel rejects DC.
 *
 * States `[v′, v″, v‴]` = [in-phase, quadrature, TOGI-tracker]. With ω₀ = 2π f0 and the
 * post-gain bus `w = alpha·ω₀·(v − v′)`:
 *
 *     v̇′  = alpha·ω₀·(v − v′) − ω₀·v″   (SOGI in-phase integrator)
 *     v̇″  = ω₀·v′                   (SOGI quadrature integrator)
 *     v̇‴  = alpha·ω₀·(v − v′) − ω₀·v‴   (TOGI first-order, self-damped)
 *
 * Outputs: `v_o = v′` (band-pass, unity at ω₀) and `q·v_o = v″ − v‴`
 * (DC-rejecting quadrature). The transfer functions are
 *
 *     v_o / v   = alpha·ω₀·s / (s² + alpha·ω₀·s + ω₀²)
 *     q·v_o / v = alpha·ω₀·s·(ω₀ − s) / [ (s + ω₀)·(s² + alpha·ω₀·s + ω₀²) ]
 *
 * The q·v_o numerator has a zero at s = 0 (DC rejection) and the pair evaluates
 * to (1∠0°, 1∠−90°) at s = jω₀ — a clean unity-magnitude quadrature pair.
 *
 * @note The discrete `mstogi(omega_0, alpha, Ts)` overload is the exact ZOH
 *       discretization for simulation/analysis. The runtime @ref MSTOGI class uses
 *       a separate, lighter realization (open-loop resonator + forward-Euler
 *       washout) and does not share this state-space.
 *
 * @see Rodríguez et al., "Discrete-time implementation of second order
 *      generalized integrators for grid converters," IECON 2008,
 *      https://doi.org/10.1109/IECON.2008.4757983
 *
 * @param f0 Fundamental frequency [Hz] (> 0)
 * @param alpha Damping gain (typically √2 for ~unity-Q SOGI tuning, > 0)
 * @return StateSpaceResult with the MSTOGI system (outputs: v_o, q·v_o); refused
 *         (invalid_parameter) if f0 or alpha is not positive
 */
template<typename T = double>
[[nodiscard]] constexpr StateSpaceResult<3, 1, 2, T> mstogi(T f0, T alpha = damp::numbers::sqrt2_v<T>) {
    if (!detail::valid_resonator(f0, alpha)) {
        return ::damp::detail::refused<StateSpaceResult<3, 1, 2, T>>(DesignStatus::invalid_parameter);
    }
    const T w0 = T{2} * damp::numbers::pi_v<T> * f0;
    const T alpha_omega = alpha * w0;

    return detail::designed(StateSpace<3, 1, 2, T>{
        .A = Matrix<3, 3, T>{
            {-alpha_omega, -w0, T{0}},
            {w0, T{0}, T{0}},
            {-alpha_omega, T{0}, -w0},
        },

        .B = Matrix<3, 1, T>{
            {alpha_omega},
            {T{0}},
            {alpha_omega},
        },

        .C = Matrix<2, 3, T>{
            {T{1}, T{0}, T{0}},  // v_o   = v′      (band-pass)
            {T{0}, T{1}, T{-1}}, // q·v_o = v″ − v‴ (DC-rejecting quadrature)
        },

        .D = Matrix<2, 1, T>::zeros(),
    });
}

/**
 * @brief Mixed Second/Third-Order Generalized Integrator (MSTOGI) design (discrete-time)
 *
 * Exact ZOH discretization of the continuous MSTOGI for simulation/analysis at a
 * fixed sample rate (no Tustin frequency warping).
 *
 * @param f0 Fundamental frequency [Hz]
 * @param alpha Damping gain
 * @param Ts Sample time [s]
 * @return StateSpaceResult with the discrete MSTOGI system; refused with invalid_sample_time if
 *         Ts ≤ 0, invalid_parameter for f0 / alpha, singular if the discretization fails
 */
template<typename T = float>
[[nodiscard]] constexpr StateSpaceResult<3, 1, 2, T> mstogi(T f0, T alpha, T Ts) {
    return detail::sampled(mstogi(f0, alpha), Ts);
}

} // namespace design

/**
 * @brief Runtime SOGI quadrature-signal generator (live-frequency resonator)
 *
 * This minimal runtime object stores state and performs one-step updates each
 * call:
 *
 *   y = Cx
 *   x = Ax + Bu,   u = alpha·(in − v′)
 *
 * where A = [[cos ω₀Ts, sin ω₀Ts], [−sin ω₀Ts, cos ω₀Ts]] and B = [1 − cos ω₀Ts, sin ω₀Ts]
 * are the exact ZOH step of the undamped resonator at `w0 = 2*pi*freq`, rebuilt
 * each tick so @p freq may track a PLL/FLL estimate. The damping term enters
 * through u, held over the step (so this is not the same realization as
 * `design::sogi(freq, alpha, Ts)`, which discretizes the damped system).
 * Returns {v′ (band-pass), qv′ (quadrature, 90° lag)}.
 *
 * Feeding u back from the previous sample bounds the damping gain: the closed loop
 * z² + (α sin ω₀Ts − 2 cos ω₀Ts) z + (1 − α sin ω₀Ts) is stable only while
 * α < cot(ω₀Ts / 2), i.e. ω₀Ts < 2 atan(1/α) (f/fs < 0.196 for α = √2).
 *
 * @see "Understanding Digital Signal Processing" (Lyons, 2011), §13.36
 *
 * @tparam T Scalar type (float for embedded deployment)
 */
template<typename T = float>
class SOGI {
public:
    constexpr SOGI() = default;

    /**
     * @brief One tick: returns {v′ (band-pass), qv′ (quadrature)} and advances the resonator.
     * @param in    Input sample
     * @param freq  [Hz] Centre frequency f0 (stable while α < cot(π f0 Ts))
     * @param alpha Damping gain α (> 0)
     * @param Ts    [s] Sample time
     */
    [[nodiscard]] constexpr damp::pair<T, T> operator()(T in, T freq, T alpha, T Ts) {
        const auto wT = freq * T{2} * damp::numbers::pi_v<T> * Ts;
        const auto [s, c] = damp::sincos(wT); // sin(wT), cos(wT)

        const StateSpace sys = {
            .A = Matrix<2, 2, T>{
                {c, s},  // Quadrature
                {-s, c}, // Band-pass
            },

            .B = ColVec<2, T>{
                T{1} - c,
                s,
            },

            .C = Matrix<2, 2, T>{
                {T{0}, T{1}}, // band-pass
                {T{1}, T{0}}, // quadrature
            },
        };

        const auto u = (in - x(1)) * alpha;

        const ColVec<2, T> y = sys.C * x;
        x = (sys.A * x) + (sys.B * u);

        return {y(0), y(1)};
    }

    /** @brief Clear the resonator state. */
    constexpr void reset() { x = {}; }

private:
    ColVec<2, T> x = {};
};

/**
 * @brief Runtime MSTOGI with exact resonator and forward-Euler washout
 *
 * This minimal runtime object stores state and performs one-step updates each
 * call:
 *
 *   y = Cx
 *   x = Ax + Bu
 *
 * where `w0 = 2*pi*freq`.
 *
 * Runtime keeps the resonator on the exact sin/cos discretization and uses a
 * forward-Euler step for the TOGI washout branch:
 *
 *   x_t[k+1] = (1 - w0*Ts) * x_t[k] + alpha*w0*Ts * (u[k] - bp[k])
 *
 * This removes the per-sample exp() from the hot path while preserving DC
 * rejection behavior for practical sample rates (w0*Ts << 1). The resonator loop
 * is the one of @ref SOGI, stable while alpha < cot(w0*Ts / 2).
 *
 * @see design::mstogi() for the discrete design model and transfer functions
 * @see Rodríguez et al., IECON 2008, https://doi.org/10.1109/IECON.2008.4757983
 *
 * @tparam T Scalar type (float for embedded deployment)
 */
template<typename T = float>
class MSTOGI {
public:
    constexpr MSTOGI() = default;

    /**
     * @brief One tick: returns {v_o (band-pass), q·v_o (DC-rejecting quadrature)}.
     * @param in    Input sample
     * @param freq  [Hz] Centre frequency f0 (stable while α < cot(π f0 Ts))
     * @param alpha Damping gain α (> 0)
     * @param Ts    [s] Sample time
     */
    [[nodiscard]] constexpr damp::pair<T, T> operator()(T in, T freq, T alpha, T Ts) {
        const T wT = freq * T{2} * damp::numbers::pi_v<T> * Ts;
        const auto [s, c] = damp::sincos(wT); // sin(wT), cos(wT)

        const StateSpace sys = {
            .A = Matrix<2, 2, T>{
                {c, s},
                {-s, c},
            },

            .B = ColVec<2, T>{
                T{1} - c,
                s,
            },

            .C = Matrix<2, 2, T>{
                {T{0}, T{1}}, // band-pass
                {T{1}, T{0}}, // quadrature
            },
        };

        // Update outputs
        const T u = (in - x(1)) * alpha;

        // Washout of DC offset in quadrature output
        togi_state += (u - togi_state) * wT;

        // Update internal state
        const ColVec<2, T> y = sys.C * x;
        x = (sys.A * x) + (sys.B * u);

        return {y(0), y(1) - togi_state};
    }

    /** @brief Clear the resonator and washout state. */
    constexpr void reset() {
        x = {};
        togi_state = {};
    }

private:
    ColVec<2, T> x = {};

    T togi_state = {};
};

/**
 * @ingroup filters
 * @brief SOGI with a Frequency-Locked Loop — self-tuning single-tone tracker.
 *
 * A SOGI quadrature generator whose center frequency ω *adapts* to lock onto the
 * dominant frequency of the input, with no sweep and no prior knowledge beyond an
 * initial guess. Outputs the locked frequency, the in-phase (band-pass) and
 * quadrature signals, and the amplitude. Unlike the fixed-ω `SOGI`, you don't tell
 * it the frequency — it finds it.
 *
 * Uses:
 *  - Grid synchronization: lock to the line frequency for a single-phase PLL
 *    under frequency drift.
 *  - Online resonance tracking: lock onto a structure's vibratory mode and
 *    follow it as it drifts (temperature, payload), e.g. to re-tune an input
 *    shaper or a notch online — the streaming counterpart to an offline
 *    Goertzel/DFT sweep (filters/spectral.hpp), which needs the frequency known.
 *
 * SOGI band-pass D(s)=v′/v and quadrature Q(s)=qv′/v at center ω:
 *
 *     D(s) = k·ω·s / (s² + k·ω·s + ω²),   Q(s) = k·ω² / (s² + k·ω·s + ω²)
 *
 * The FLL drives ω from the frequency-error signal ε_f = ε·qv′ (ε = v − v′), whose
 * DC component vanishes at lock. The update is amplitude-normalized so the lock
 * dynamics are independent of input level:
 *
 *     ω̇ = −Γ·ω·(ε·qv′) / (v′² + qv′²)
 *
 * @note The SOGI step uses the exact discrete resonator at the current ω each tick,
 *       with the damping fed back from the previous sample, which is stable only for
 *       ω·Ts < 2·atan(1/k) (f/fs < 0.196 at k = √2). The upper clamp is held to 90 %
 *       of that bound; the FLL integrates ω by forward Euler. All constexpr,
 *       allocation-free.
 *
 * @see "Multiresonant Frequency-Locked Loops" (P. Rodriguez et al., IEEE TIE),
 *      https://doi.org/10.1109/TIE.2010.2042420
 *
 * @tparam T Scalar type (float for embedded deployment)
 */
template<typename T = float>
class SogiFll {
public:
    constexpr SogiFll() = default;

    /**
     * @param initial_freq_hz Starting frequency guess [Hz]
     * @param Ts              Sample time [s]
     * @param sogi_gain       SOGI gain k (√2 ≈ 0.707 damping; default)
     * @param fll_gain        Normalized FLL gain Γ (larger = faster lock, less smooth)
     * @param freq_lim_hz     Clamp on the tracked frequency [Hz]; a lower end ≤ 0 means 0.5 Hz,
     *                        an upper end ≤ 0 or past the stability cap means the cap,
     *                        0.9·2·atan(1/k) / (2π·Ts) (default: both)
     *
     * valid() is false unless Ts > 0, k > 0, the band is non-empty and the initial guess
     * is below its upper end.
     */
    constexpr SogiFll(
        T            initial_freq_hz,
        T            Ts,
        T            sogi_gain = damp::numbers::sqrt2_v<T>,
        T            fll_gain = T{2},
        Bounds<1, T> freq_lim_hz = {}
    ) : Ts_(Ts),
        k_(sogi_gain),
        gamma_(fll_gain),
        omega_(two_pi() * initial_freq_hz),
        omega_min_(two_pi() * (freq_lim_hz.lower(0) > T{0} ? freq_lim_hz.lower(0) : static_cast<T>(0.5))),
        omega_max_(two_pi() * upper_frequency(freq_lim_hz.upper(0), Ts, sogi_gain)),
        valid_(
            Ts > T{0} && damp::isfinite(Ts) && sogi_gain > T{0} && initial_freq_hz > T{0} && omega_min_ < omega_max_
            && omega_ <= omega_max_
        ) {}

    /** @brief One tick with a new input sample: advance the resonator, then the FLL frequency. */
    constexpr void step(T in) {
        if (!valid_) {
            return;
        }
        const T wt = omega_ * Ts_;
        const auto [s, c] = damp::sincos(wt);

        const T bandpass = x1_; // in-phase v′
        const T quad = x0_;     // quadrature qv′
        const T err = in - bandpass;
        const T u = err * k_;

        // Exact discrete SOGI resonator step (same A,B as the fixed-ω SOGI).
        const T x0n = (c * x0_) + (s * x1_) + ((T{1} - c) * u);
        const T x1n = (-s * x0_) + (c * x1_) + (s * u);
        x0_ = x0n;
        x1_ = x1n;

        // Amplitude-normalized FLL: ω̇ = −Γ·ω·(ε·qv′)/(v′²+qv′²). The ratio is
        // scale-free; skip only the exact-zero state (no signal yet, 0/0).
        const T amp_sq = (bandpass * bandpass) + (quad * quad);
        if (amp_sq > T{0}) {
            const T eps_f = err * quad;
            omega_ -= Ts_ * gamma_ * omega_ * eps_f / amp_sq;
            omega_ = damp::clamp(omega_, omega_min_, omega_max_);
        }
    }

    /** @brief Locked frequency [Hz]. */
    [[nodiscard]] constexpr T frequency_hz() const { return omega_ / two_pi(); }

    /** @brief Locked frequency [rad/s]. */
    [[nodiscard]] constexpr T frequency_rad() const { return omega_; }

    /** @brief In-phase (band-pass) output — the input filtered to the locked frequency. */
    [[nodiscard]] constexpr T in_phase() const { return x1_; }

    /** @brief Quadrature output (90° lag). */
    [[nodiscard]] constexpr T quadrature() const { return x0_; }

    /** @brief Estimated amplitude of the tracked tone. */
    [[nodiscard]] constexpr T amplitude() const { return damp::hypot(x0_, x1_); }

    /** @brief Estimated phase of the tracked tone [rad]. */
    [[nodiscard]] constexpr T    phase() const { return damp::atan2(x0_, x1_); }
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /** @brief Clear the resonator and re-seed the frequency, clamped to the band. */
    constexpr void reset(T initial_freq_hz) {
        x0_ = T{0};
        x1_ = T{0};
        omega_ = two_pi() * initial_freq_hz;
        if (valid_) {
            omega_ = damp::clamp(omega_, omega_min_, omega_max_);
        }
    }

private:
    static constexpr T two_pi() { return T{2} * damp::numbers::pi_v<T>; }

    /** @brief Upper clamp [Hz]: @p hi when in (0, cap), else the cap, 90 % of the loop's stability bound 2·atan(1/k) / (2π·Ts). */
    static constexpr T upper_frequency(T hi, T Ts, T k) {
        if (!(k > T{0}) || !(Ts > T{0})) {
            return T{0}; // no band: the constructor marks the tracker invalid
        }
        const T cap = static_cast<T>(0.9) * T{2} * damp::atan(T{1} / k) / (two_pi() * Ts);
        return (hi > T{0} && hi < cap) ? hi : cap;
    }

    T    Ts_{T{1}};
    T    k_{damp::numbers::sqrt2_v<T>};
    T    gamma_{T{2}};
    T    omega_{T{0}};
    T    omega_min_{T{0}};
    T    omega_max_{T{0}};
    T    x0_{T{0}}; // quadrature state
    T    x1_{T{0}}; // in-phase state
    bool valid_{false};
};

} // namespace damp
