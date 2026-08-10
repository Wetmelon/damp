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
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @brief Second-Order Generalized Integrator (SOGI) design
 *
 * SOGI provides bandpass filtering for grid synchronization.
 * Produces quadrature signals (90° phase shift).
 *
 * Transfer functions:
 * H_bp(s) = (alpha*ω₀*s) / (s² + alpha*ω₀*s + ω₀²)    [bandpass]
 * H_q(s) = (alpha*ω₀²) / (s² + alpha*ω₀*s + ω₀²)     [quadrature]
 * H_notch(s) = (s² + ω₀²) / (s² + alpha*ω₀*s + ω₀²) [notch]
 *
 * @param w0 Fundamental frequency [rad/s]
 * @param alpha Damping gain (typically 1.0-2.0)
 * @return StateSpace<2, 1, 2, T> SOGI system
 */
template<typename T = double>
[[nodiscard]] constexpr StateSpace<2, 1, 2, T> sogi(T w0, T alpha = damp::numbers::sqrt2_v<T>) {
    return {
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
    };
}

/**
 * @brief Second-Order Generalized Integrator (SOGI) design (discrete-time)
 *
 * Exact ZOH discretization of the continuous SOGI for simulation/analysis at a
 * fixed sample rate. (The runtime @ref SOGI class rebuilds its own resonator each
 * tick from the live frequency and does not use this.)
 *
 * @param w0 Fundamental frequency [rad/s]
 * @param alpha Damping gain (typically 1.0-2.0)
 * @param Ts Sample time [s]
 * @return StateSpace<2, 1, 2, T> Discrete-time SOGI system
 */
template<typename T = float>
[[nodiscard]] constexpr StateSpace<2, 1, 2, T> sogi(T w0, T alpha, T Ts) {
    return *discretize(sogi(w0, alpha), Ts, DiscretizationMethod::ZOH);
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
 * States `[v′, v″, v‴]` = [in-phase, quadrature, TOGI-tracker]. With the
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
 * @param w0 Fundamental frequency [rad/s]
 * @param alpha Damping gain (typically √2 for ~unity-Q SOGI tuning)
 * @return StateSpace<3, 1, 2, T> MSTOGI system (outputs: v_o, q·v_o)
 */
template<typename T = double>
[[nodiscard]] constexpr StateSpace<3, 1, 2, T> mstogi(T w0, T alpha = damp::numbers::sqrt2_v<T>) {
    const T alpha_omega = alpha * w0;

    return {
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
    };
}

/**
 * @brief Mixed Second/Third-Order Generalized Integrator (MSTOGI) design (discrete-time)
 *
 * Exact ZOH discretization of the continuous MSTOGI for simulation/analysis at a
 * fixed sample rate (no Tustin frequency warping).
 *
 * @param w0 Fundamental frequency [rad/s]
 * @param alpha Damping gain
 * @param Ts Sample time [s]
 * @return StateSpace<3, 1, 2, T> Discrete-time MSTOGI system
 */
template<typename T = float>
[[nodiscard]] constexpr StateSpace<3, 1, 2, T> mstogi(T w0, T alpha, T Ts) {
    return *discretize(mstogi(w0, alpha), Ts, DiscretizationMethod::ZOH);
}

} // namespace design

/**
 * @brief Runtime SOGI wrapper around design::sogi(w0, alpha, Ts)
 *
 * This minimal runtime object stores state and performs one-step updates using
 * the discrete design realization each call:
 *
 *   y = Cx
 *   x = Ax + Bu
 *
 * where A/B/C come from `design::sogi(w0, alpha, Ts)` and
 * `w0 = 2*pi*freq`.
 *
 * @see "Understanding Digital Signal Processing" (Lyons, 2011), §13.36
 *
 * @tparam T Scalar type (float for embedded deployment)
 */
template<typename T = float>
class SOGI {
public:
    constexpr SOGI() = default;

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

    /// Clear the resonator state.
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
 * rejection behavior for practical sample rates (w0*Ts << 1).
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

    /// Clear the resonator and washout state.
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
 * @note The SOGI step uses the exact discrete resonator at the current ω each tick
 *       (stable for any ω·Ts < π); the FLL integrates ω by forward Euler. All
 *       constexpr, allocation-free.
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
     * @param freq_min_hz     Lower clamp on the tracked frequency [Hz] (0 → 0.5 Hz)
     * @param freq_max_hz     Upper clamp [Hz] (0 → Nyquist)
     */
    constexpr SogiFll(
        T initial_freq_hz,
        T Ts,
        T sogi_gain = damp::numbers::sqrt2_v<T>,
        T fll_gain = T{2},
        T freq_min_hz = T{0},
        T freq_max_hz = T{0}
    ) : Ts_(Ts),
        k_(sogi_gain),
        gamma_(fll_gain),
        omega_(two_pi() * initial_freq_hz),
        omega_min_(two_pi() * (freq_min_hz > T{0} ? freq_min_hz : static_cast<T>(0.5))),
        omega_max_(two_pi() * (freq_max_hz > T{0} ? freq_max_hz : T{1} / (T{2} * Ts))),
        valid_(Ts > T{0} && initial_freq_hz > T{0}) {}

    /// Advance one step with a new input sample.
    constexpr void update(T in) {
        if (!valid_) {
            return;
        }
        const T wt = omega_ * Ts_;
        const auto [s, c] = damp::sincos(wt);

        const T bandpass = x1_; // in-phase v′
        const T quad = x0_;     // quadrature qv′
        const T eps = in - bandpass;
        const T u = eps * k_;

        // Exact discrete SOGI resonator step (same A,B as the fixed-ω SOGI).
        const T x0n = (c * x0_) + (s * x1_) + ((T{1} - c) * u);
        const T x1n = (-s * x0_) + (c * x1_) + (s * u);
        x0_ = x0n;
        x1_ = x1n;

        // Amplitude-normalized FLL: ω̇ = −Γ·ω·(ε·qv′)/(v′²+qv′²).
        const T amp_sq = (bandpass * bandpass) + (quad * quad);
        const T eps_f = eps * quad;
        omega_ -= Ts_ * gamma_ * omega_ * eps_f / (amp_sq + tiny());
        omega_ = damp::clamp(omega_, omega_min_, omega_max_);
    }

    /// Locked frequency [Hz].
    [[nodiscard]] constexpr T frequency_hz() const { return omega_ / two_pi(); }

    /// Locked frequency [rad/s].
    [[nodiscard]] constexpr T frequency_rad() const { return omega_; }

    /// In-phase (band-pass) output — the input filtered to the locked frequency.
    [[nodiscard]] constexpr T in_phase() const { return x1_; }

    /// Quadrature output (90° lag).
    [[nodiscard]] constexpr T quadrature() const { return x0_; }

    /// Estimated amplitude of the tracked tone.
    [[nodiscard]] constexpr T amplitude() const { return damp::hypot(x0_, x1_); }

    /// Estimated phase of the tracked tone [rad].
    [[nodiscard]] constexpr T    phase() const { return damp::atan2(x0_, x1_); }
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    constexpr void reset(T initial_freq_hz) {
        x0_ = T{0};
        x1_ = T{0};
        omega_ = two_pi() * initial_freq_hz;
    }

private:
    static constexpr T two_pi() { return T{2} * damp::numbers::pi_v<T>; }
    static constexpr T tiny() { return static_cast<T>(1e-12); }

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
