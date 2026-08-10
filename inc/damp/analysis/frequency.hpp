// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file frequency.hpp
 * @brief Bode, Nyquist, Nichols, sigma, margins, loop metrics, impedance / Middlebrook
 */

#include <cstddef>
#include <limits>
#include <numbers>
#include <vector>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/svd.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"


namespace damp {
namespace analysis {

// ============================================================================
// Frequency Response Data
// ============================================================================

/**
 * @brief Single-point frequency response result
 */
template<typename T = double>
struct FrequencyPoint {
    T omega{};        ///< Frequency (rad/s)
    T magnitude{};    ///< Magnitude (absolute)
    T magnitude_db{}; ///< Magnitude (dB)
    T phase_deg{};    ///< Phase (degrees)
};

/**
 * @brief Bode plot data for a SISO system
 *
 * Member gain_margin()/phase_margin() delegate to the unwrapped free functions
 * (defined below) so multi-wrap Bode data is handled correctly.
 */
template<typename T = double>
struct BodeResult; // incomplete: needed for unwrapped-margin free-function protos

template<typename T>
[[nodiscard]] constexpr damp::optional<damp::pair<T, T>> phase_margin_unwrapped(const BodeResult<T>& result);
template<typename T>
[[nodiscard]] constexpr damp::optional<damp::pair<T, T>> gain_margin_unwrapped(const BodeResult<T>& result);

template<typename T>
struct BodeResult {
    std::vector<FrequencyPoint<T>> points; ///< Frequency response data points

    /**
     * @brief Find gain margin (dB above 1.0 at -180° phase crossing)
     *
     * Gain margin is 1/|G(jω)| at the frequency where phase crosses -180°.
     * Positive values indicate stability. Uses the unwrapped phase trajectory
     * (same as @ref gain_margin_unwrapped) so multi-wrap Bode data is handled
     * correctly.
     *
     * @return {gain_margin_dB, crossover_frequency} or nullopt if no -180° crossing
     */
    [[nodiscard]] constexpr damp::optional<damp::pair<T, T>> gain_margin() const {
        return gain_margin_unwrapped(*this);
    }

    /**
     * @brief Find phase margin (degrees above -180° at 0dB gain crossing)
     *
     * Phase margin is 180° + phase(G(jω)) at the frequency where |G(jω)| = 1 (0 dB).
     * Positive values indicate stability. Uses the unwrapped phase trajectory
     * (same as @ref phase_margin_unwrapped) so multi-wrap Bode data is handled
     * correctly; the result is mapped to (-180, 180].
     *
     * @return {phase_margin_deg, crossover_frequency} or nullopt if no 0dB crossing
     */
    [[nodiscard]] constexpr damp::optional<damp::pair<T, T>> phase_margin() const {
        return phase_margin_unwrapped(*this);
    }

    /**
     * @brief Find -3dB bandwidth
     *
     * The frequency at which the magnitude drops 3 dB below the DC (or peak) value.
     *
     * @return Bandwidth in rad/s, or nullopt if not found
     */
    [[nodiscard]] constexpr damp::optional<T> bandwidth() const {
        if (points.empty()) {
            return damp::nullopt;
        }
        T dc_db = points[0].magnitude_db;
        T threshold = dc_db - T{3};
        for (size_t i = 1; i < points.size(); ++i) {
            if (points[i].magnitude_db < threshold) {
                // Interpolate
                T frac = (threshold - points[i - 1].magnitude_db) / (points[i].magnitude_db - points[i - 1].magnitude_db);
                return points[i - 1].omega + frac * (points[i].omega - points[i - 1].omega);
            }
        }
        return damp::nullopt;
    }
};

/**
 * @brief Unwrap phase data in degrees to avoid +/-180 discontinuities
 *
 * Given wrapped phase samples (typically in [-180, 180]), produces a
 * continuous phase trajectory by adding/subtracting 360 deg at jumps.
 *
 * @param phase_deg Wrapped phase samples in degrees
 * @return Unwrapped phase samples in degrees
 */
template<typename T>
[[nodiscard]] constexpr std::vector<T> unwrap_phase_deg(const std::vector<T>& phase_deg) {
    if (phase_deg.empty()) {
        return {};
    }

    std::vector<T> unwrapped = phase_deg;
    for (size_t i = 1; i < unwrapped.size(); ++i) {
        const T delta = unwrapped[i] - unwrapped[i - 1];
        // Range-reduce to nearest equivalent increment in [-180, 180).
        const T n = damp::floor(delta / T{360} + static_cast<T>(0.5));
        const T delta_reduced = delta - T{360} * n;
        unwrapped[i] = unwrapped[i - 1] + delta_reduced;
    }
    return unwrapped;
}

/**
 * @brief Normalize phase margin to (-180, 180]
 *
 * @param pm_deg Raw phase margin in degrees
 * @return Canonical phase margin in (-180, 180]
 */
template<typename T>
[[nodiscard]] constexpr T canonical_phase_margin(T pm_deg) {
    // Reuse shared wrapping helper, then remap to (-180, 180].
    const T wrapped = wrap(pm_deg, T{-180}, T{180});
    if (wrapped <= T{-180}) {
        return wrapped + T{360};
    }
    return wrapped;
}

/**
 * @brief Find phase margin using unwrapped phase trajectory
 *
 * Uses 0 dB crossing from Bode magnitude and computes PM as 180 + phase(wc),
 * then maps to canonical range (-180, 180].
 *
 * @return {phase_margin_deg, crossover_frequency} or nullopt if no 0dB crossing
 */
template<typename T>
[[nodiscard]] constexpr damp::optional<damp::pair<T, T>> phase_margin_unwrapped(const BodeResult<T>& result) {
    if (result.points.size() < 2) {
        return damp::nullopt;
    }

    std::vector<T> wrapped_phase;
    wrapped_phase.reserve(result.points.size());
    for (const auto& pt : result.points) {
        wrapped_phase.push_back(pt.phase_deg);
    }
    const auto unwrapped = unwrap_phase_deg(wrapped_phase);

    for (size_t i = 1; i < result.points.size(); ++i) {
        const T mag_prev = result.points[i - 1].magnitude_db;
        const T mag_curr = result.points[i].magnitude_db;
        if (mag_prev >= T{0} && mag_curr < T{0}) {
            const T frac = (T{0} - mag_prev) / (mag_curr - mag_prev);
            const T omega_cross = result.points[i - 1].omega + frac * (result.points[i].omega - result.points[i - 1].omega);
            const T phase_cross = unwrapped[i - 1] + frac * (unwrapped[i] - unwrapped[i - 1]);
            const T pm = canonical_phase_margin(T{180} + phase_cross);
            return damp::pair{pm, omega_cross};
        }
    }

    return damp::nullopt;
}

/**
 * @brief Find gain margin using unwrapped phase trajectory
 *
 * Finds the first -180 deg crossing on the unwrapped phase trajectory and
 * computes gain margin as -|L|_dB at that crossing.
 *
 * @return {gain_margin_dB, crossover_frequency} or nullopt if no -180 crossing
 */
template<typename T>
[[nodiscard]] constexpr damp::optional<damp::pair<T, T>> gain_margin_unwrapped(const BodeResult<T>& result) {
    if (result.points.size() < 2) {
        return damp::nullopt;
    }

    std::vector<T> wrapped_phase;
    wrapped_phase.reserve(result.points.size());
    for (const auto& pt : result.points) {
        wrapped_phase.push_back(pt.phase_deg);
    }
    const auto unwrapped = unwrap_phase_deg(wrapped_phase);

    for (size_t i = 1; i < result.points.size(); ++i) {
        const T p0 = unwrapped[i - 1] + T{180};
        const T p1 = unwrapped[i] + T{180};

        if ((p0 >= T{0} && p1 < T{0}) || (p0 <= T{0} && p1 > T{0})) {
            const T frac = (T{0} - p0) / (p1 - p0);
            const T omega_cross = result.points[i - 1].omega + frac * (result.points[i].omega - result.points[i - 1].omega);
            const T mag_cross = result.points[i - 1].magnitude_db + frac * (result.points[i].magnitude_db - result.points[i - 1].magnitude_db);
            return damp::pair{-mag_cross, omega_cross};
        }
    }

    return damp::nullopt;
}

/**
 * @brief Compute Bode plot data for a SISO state-space system
 *
 * Continuous-time systems evaluate G(jω) = C(jωI - A)^{-1}B + D.
 * Discrete-time systems (Ts > 0) evaluate G(z) on the unit circle
 * z = e^{jωTs}. Dispatch is automatic based on sys.is_discrete().
 *
 * Frequencies where @f$ sI-A @f$ (or @f$ zI-A @f$) is singular — i.e. the
 * evaluation point sits exactly on a jω-axis / unit-circle pole — are
 * omitted from the result. No infinite or NaN samples are fabricated;
 * the returned grid may be shorter than @p omega. @ref eval_frf returns
 * @c nullopt at those points and this function skips them.
 *
 * @param sys   SISO state-space system (continuous or discrete)
 * @param omega Vector of frequencies (rad/s)
 * @return BodeResult with magnitude and phase at each successfully evaluated frequency
 *
 * @note Compare with MATLAB®'s bode(sys, w) (data portion).
 * @see eval_frf(), nyquist(), sigma()
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr BodeResult<T> bode(
    const StateSpace<NX, 1, 1, T, NW, NV>& sys,
    const std::vector<T>&                  omega
) {
    using Cplx = damp::complex<T>;
    BodeResult<T> result;
    result.points.reserve(omega.size());

    for (const auto& w : omega) {
        Cplx s_or_z{T{0}, w};
        if (sys.is_discrete()) {
            const auto [s, c] = damp::sincos(w * sys.Ts);
            s_or_z = Cplx{c, s};
        }
        const auto G_frf = eval_frf(sys, s_or_z);
        if (!G_frf) {
            continue; // pole on the evaluation contour — skip honestly
        }
        const Cplx G = (*G_frf)(0, 0);
        const T    mag = damp::abs(G);
        const T    phase_rad = damp::arg(G);
        const T    mag_db = mag > T{0} ? T{20} * damp::log10(mag) : T{-300};
        const T    phase_deg = phase_rad * T{180} / damp::numbers::pi_v<T>;
        result.points.push_back({w, mag, mag_db, phase_deg});
    }
    return result;
}

/**
 * @brief Compute Bode plot data for a SISO transfer function
 *
 * Evaluates H(jω) = num(jω)/den(jω) over a vector of frequencies.
 * Samples where |den(jω)| is below default_tol (jω-axis poles) are
 * omitted — same honest skip policy as the state-space overload.
 *
 * @param num   Numerator coefficients (ascending powers of s)
 * @param den   Denominator coefficients (ascending powers of s)
 * @param omega Vector of frequencies (rad/s)
 * @return BodeResult with magnitude and phase at each successfully evaluated frequency
 */
template<size_t Nnum, size_t Nden, typename T>
[[nodiscard]] constexpr BodeResult<T> bode(
    const damp::array<T, Nnum>& num,
    const damp::array<T, Nden>& den,
    const std::vector<T>&       omega
) {
    using Cplx = damp::complex<T>;
    BodeResult<T> result;
    result.points.reserve(omega.size());

    for (const auto& w : omega) {
        Cplx jw{T{0}, w};

        // Evaluate numerator polynomial at jw (ascending powers)
        Cplx num_val{T{0}, T{0}};
        Cplx jw_power{T{1}, T{0}};
        for (size_t i = 0; i < Nnum; ++i) {
            num_val = num_val + num[i] * jw_power;
            jw_power = jw_power * jw;
        }

        // Evaluate denominator polynomial at jw (ascending powers)
        Cplx den_val{T{0}, T{0}};
        jw_power = Cplx{T{1}, T{0}};
        for (size_t i = 0; i < Nden; ++i) {
            den_val = den_val + den[i] * jw_power;
            jw_power = jw_power * jw;
        }

        // Skip jω-axis poles (den ≈ 0) — no Inf/NaN samples.
        if (!(damp::abs(den_val) > default_tol<T>())) {
            continue;
        }

        const Cplx G = num_val / den_val;
        const T    mag = damp::abs(G);
        const T    phase_rad = damp::arg(G);
        const T    mag_db = mag > T{0} ? T{20} * damp::log10(mag) : T{-300};
        const T    phase_deg = phase_rad * T{180} / damp::numbers::pi_v<T>;
        result.points.push_back({w, mag, mag_db, phase_deg});
    }
    return result;
}

/**
 * @brief Compute Bode plot data for a SISO transfer function object
 */
template<size_t Nnum, size_t Nden, typename T>
[[nodiscard]] constexpr BodeResult<T> bode(
    const TransferFunction<Nnum, Nden, T>& tf,
    const std::vector<T>&                  omega
) {
    return bode(tf.num, tf.den, omega);
}

/**
 * @brief Compute Bode plot data for a discrete-time SISO state-space system
 *
 * Evaluates G(z) on the unit circle z = e^(j*omega*Ts), where omega is in rad/s.
 *
 * @note bode() now auto-dispatches on sys.is_discrete(); this alias is retained
 *       for discoverability and call sites that want to be explicit.
 *
 * @param sys   Discrete-time SISO state-space system
 * @param omega Vector of frequencies (rad/s)
 * @return BodeResult with magnitude and phase at each frequency
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr BodeResult<T> bode_discrete(
    const StateSpace<NX, 1, 1, T, NW, NV>& sys,
    const std::vector<T>&                  omega
) {
    return bode(sys, omega);
}

/**
 * @brief Single-point Nyquist response data
 */
template<typename T = double>
struct NyquistPoint {
    T                omega{};                 ///< Frequency (rad/s)
    damp::complex<T> value{};                 ///< Complex loop value (L(jω) or L(e^{jω Ts}))
    T                real{};                  ///< Real part
    T                imag{};                  ///< Imaginary part
    T                distance_to_minus_one{}; ///< |1 + L|
};

/**
 * @brief Nyquist response data across a frequency sweep
 */
template<typename T = double>
struct NyquistResult {
    std::vector<NyquistPoint<T>> points;

    /**
     * @brief Minimum Nyquist distance to the critical point -1 + j0
     * @return {minimum_distance, frequency} or nullopt when empty
     */
    [[nodiscard]] constexpr damp::optional<damp::pair<T, T>> min_distance_to_minus_one() const {
        if (points.empty()) {
            return damp::nullopt;
        }

        T min_dist = points[0].distance_to_minus_one;
        T at_omega = points[0].omega;
        for (size_t i = 1; i < points.size(); ++i) {
            if (points[i].distance_to_minus_one < min_dist) {
                min_dist = points[i].distance_to_minus_one;
                at_omega = points[i].omega;
            }
        }
        return damp::pair{min_dist, at_omega};
    }
};

/**
 * @brief Open-loop and closed-loop frequency response package
 *
 * For a loop transfer L, includes:
 * - open_loop: L
 * - sensitivity: S = 1 / (1 + L)
 * - complementary_sensitivity: T = L / (1 + L)
 * - nyquist: complex Nyquist response of L
 */
template<typename T = double>
struct LoopResponseResult {
    BodeResult<T>    open_loop;
    BodeResult<T>    sensitivity;
    BodeResult<T>    complementary_sensitivity;
    NyquistResult<T> nyquist;

    [[nodiscard]] constexpr damp::optional<damp::pair<T, T>> phase_margin_unwrapped() const {
        return analysis::phase_margin_unwrapped(open_loop);
    }

    [[nodiscard]] constexpr damp::optional<damp::pair<T, T>> gain_margin_unwrapped() const {
        return analysis::gain_margin_unwrapped(open_loop);
    }

    [[nodiscard]] constexpr damp::optional<T> closed_loop_bandwidth() const {
        return complementary_sensitivity.bandwidth();
    }
};

template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr LoopResponseResult<T> loop_response(
    const StateSpace<NX, 1, 1, T, NW, NV>& loop,
    const std::vector<T>&                  omega
);

/**
 * @brief Compact loop summary metrics for quick stability/robustness checks
 *
 * Designed as a single result payload similar to what users expect from
 * MATLAB®/control-toolbox workflows: margins, bandwidth, Nyquist distance,
 * and peak sensitivity in one object.
 */
template<typename T = double>
struct LoopSummary {
    damp::optional<damp::pair<T, T>> phase_margin;                                        ///< {PM [deg], gain crossover omega [rad/s]}
    damp::optional<damp::pair<T, T>> gain_margin;                                         ///< {GM [dB], phase crossover omega [rad/s]}
    damp::optional<T>                bandwidth;                                           ///< Closed-loop bandwidth from T=L/(1+L), rad/s
    damp::optional<damp::pair<T, T>> min_nyquist_distance;                                ///< {min|1+L|, omega [rad/s]}
    T                                peak_sensitivity_db{-std::numeric_limits<T>::max()}; ///< max 20*log10|S|
};

/**
 * @brief Summarize loop_response() results into one compact metrics struct
 */
template<typename T>
[[nodiscard]] constexpr LoopSummary<T> summarize_loop_response(const LoopResponseResult<T>& response) {
    LoopSummary<T> summary{};
    summary.phase_margin = response.phase_margin_unwrapped();
    summary.gain_margin = response.gain_margin_unwrapped();
    summary.bandwidth = response.closed_loop_bandwidth();
    summary.min_nyquist_distance = response.nyquist.min_distance_to_minus_one();

    for (const auto& pt : response.sensitivity.points) {
        if (pt.magnitude_db > summary.peak_sensitivity_db) {
            summary.peak_sensitivity_db = pt.magnitude_db;
        }
    }

    return summary;
}

/**
 * @brief One-call loop analysis: compute L/S/T response and return compact metrics
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr LoopSummary<T> loop_metrics(
    const StateSpace<NX, 1, 1, T, NW, NV>& loop,
    const std::vector<T>&                  omega
) {
    return summarize_loop_response(loop_response(loop, omega));
}

/**
 * @brief Compute Nyquist data for a SISO state-space system
 *
 * Continuous-time systems use s = jω.
 * Discrete-time systems use z = e^{jωTs}.
 * Contour poles (@ref eval_frf fails) are omitted — same skip policy as bode().
 *
 * @note Compare with MATLAB®'s nyquist(sys, w) (data portion).
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr NyquistResult<T> nyquist(
    const StateSpace<NX, 1, 1, T, NW, NV>& sys,
    const std::vector<T>&                  omega
) {
    using Cplx = damp::complex<T>;
    NyquistResult<T> result;
    result.points.reserve(omega.size());

    for (const auto& w : omega) {
        Cplx s_or_z{T{0}, w};
        if (sys.is_discrete()) {
            const auto [s, c] = damp::sincos(w * sys.Ts);
            s_or_z = Cplx{c, s};
        }
        const auto G_frf_opt = eval_frf(sys, s_or_z);
        if (!G_frf_opt) {
            continue; // pole on contour — skip
        }
        const Cplx      G = (*G_frf_opt)(0, 0);
        const T         dist = damp::abs(Cplx{T{1}, T{0}} + G);
        NyquistPoint<T> point{};
        point.omega = w;
        point.value = G;
        point.real = G.real();
        point.imag = G.imag();
        point.distance_to_minus_one = dist;
        result.points.push_back(point);
    }

    return result;
}

/**
 * @brief Compute Nyquist data for a SISO transfer function
 *
 * Realizes @p tf then delegates to the state-space path. Returns an empty
 * result if companion realization fails (zero leading denominator coefficient).
 */
template<size_t Nnum, size_t Nden, typename T>
[[nodiscard]] constexpr NyquistResult<T> nyquist(
    const TransferFunction<Nnum, Nden, T>& tf,
    const std::vector<T>&                  omega
) {
    const auto ss = tf.to_state_space();
    if (!ss) {
        return {};
    }
    return nyquist(*ss, omega);
}

/**
 * @brief Compute open-loop L, sensitivity S, complementary sensitivity T, and Nyquist data
 *
 * For each frequency point:
 *   S = 1/(1+L),  T = L/(1+L)
 *
 * Continuous-time systems use s = jω.
 * Discrete-time systems use z = e^{jωTs}.
 * Contour poles are omitted (same skip policy as bode() / nyquist()).
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr LoopResponseResult<T> loop_response(
    const StateSpace<NX, 1, 1, T, NW, NV>& loop,
    const std::vector<T>&                  omega
) {
    using Cplx = damp::complex<T>;

    LoopResponseResult<T> result;
    result.open_loop.points.reserve(omega.size());
    result.sensitivity.points.reserve(omega.size());
    result.complementary_sensitivity.points.reserve(omega.size());
    result.nyquist.points.reserve(omega.size());

    for (const auto& w : omega) {
        Cplx s_or_z{T{0}, w};
        if (loop.is_discrete()) {
            const auto [s, c] = damp::sincos(w * loop.Ts);
            s_or_z = Cplx{c, s};
        }

        const auto L_frf_opt = eval_frf(loop, s_or_z);
        if (!L_frf_opt) {
            continue; // pole on contour — skip
        }
        const Cplx L = (*L_frf_opt)(0, 0);
        const Cplx one{T{1}, T{0}};
        const Cplx S = one / (one + L);
        const Cplx Tresp = L / (one + L);

        const T L_mag = damp::abs(L);
        const T L_mag_db = L_mag > T{0} ? T{20} * damp::log10(L_mag) : T{-300};
        const T L_phase = damp::arg(L) * T{180} / damp::numbers::pi_v<T>;
        result.open_loop.points.push_back({w, L_mag, L_mag_db, L_phase});

        const T S_mag = damp::abs(S);
        const T S_mag_db = S_mag > T{0} ? T{20} * damp::log10(S_mag) : T{-300};
        const T S_phase = damp::arg(S) * T{180} / damp::numbers::pi_v<T>;
        result.sensitivity.points.push_back({w, S_mag, S_mag_db, S_phase});

        const T T_mag = damp::abs(Tresp);
        const T T_mag_db = T_mag > T{0} ? T{20} * damp::log10(T_mag) : T{-300};
        const T T_phase = damp::arg(Tresp) * T{180} / damp::numbers::pi_v<T>;
        result.complementary_sensitivity.points.push_back({w, T_mag, T_mag_db, T_phase});

        const T         dist = damp::abs(one + L);
        NyquistPoint<T> point{};
        point.omega = w;
        point.value = L;
        point.real = L.real();
        point.imag = L.imag();
        point.distance_to_minus_one = dist;
        result.nyquist.points.push_back(point);
    }

    return result;
}

/**
 * @brief Compute open-loop L, sensitivity S, complementary sensitivity T, and Nyquist data for a transfer function
 *
 * Returns an empty package if companion realization fails.
 */
template<size_t Nnum, size_t Nden, typename T>
[[nodiscard]] constexpr LoopResponseResult<T> loop_response(
    const TransferFunction<Nnum, Nden, T>& loop,
    const std::vector<T>&                  omega
) {
    const auto ss = loop.to_state_space();
    if (!ss) {
        return {};
    }
    return loop_response(*ss, omega);
}

/**
 * @brief One-call loop analysis for transfer-function loops
 */
template<size_t Nnum, size_t Nden, typename T>
[[nodiscard]] constexpr LoopSummary<T> loop_metrics(
    const TransferFunction<Nnum, Nden, T>& loop,
    const std::vector<T>&                  omega
) {
    return summarize_loop_response(loop_response(loop, omega));
}

// ============================================================================
// Nichols Chart Data
// ============================================================================

/**
 * @brief Single-point Nichols chart sample (open-loop phase vs magnitude)
 *
 * Axes match the classical Nichols chart: phase (deg) on the abscissa and
 * open-loop magnitude (dB) on the ordinate. The critical point for unity
 * negative feedback is @f$ (-180^\circ,\,0\,\mathrm{dB}) @f$.
 */
template<typename T = double>
struct NicholsPoint {
    T omega{};        ///< Frequency (rad/s)
    T phase_deg{};    ///< Open-loop phase (degrees)
    T magnitude_db{}; ///< Open-loop magnitude (dB)
};

/**
 * @brief Nichols chart locus across a frequency sweep
 *
 * @note Compare with MATLAB®'s nichols / nicholsplot data.
 */
template<typename T = double>
struct NicholsResult {
    std::vector<NicholsPoint<T>> points;
};

/**
 * @brief Build Nichols points from existing Bode data
 *
 * Copies @f$(\phi(\omega),\,|G|_{\mathrm{dB}}(\omega))@f$ from each Bode sample.
 * Phase is left as stored on the Bode result (typically principal value in
 * @f$[-180, 180]@f$); unwrap with unwrap_phase_deg if a continuous locus is needed.
 *
 * @param bode BodeResult from analysis::bode()
 * @return NicholsResult with the same frequency samples
 *
 * @note Compare with MATLAB®'s nichols(sys) (data portion).
 */
template<typename T>
[[nodiscard]] constexpr NicholsResult<T> nichols(const BodeResult<T>& bode) {
    NicholsResult<T> result;
    result.points.reserve(bode.points.size());
    for (const auto& pt : bode.points) {
        NicholsPoint<T> np{};
        np.omega = pt.omega;
        np.phase_deg = pt.phase_deg;
        np.magnitude_db = pt.magnitude_db;
        result.points.push_back(np);
    }
    return result;
}

/**
 * @brief Compute Nichols chart data for a SISO state-space system
 *
 * Continuous-time systems use @f$ s = j\omega @f$; discrete-time systems use
 * @f$ z = e^{j\omega T_s} @f$ (same dispatch as bode()).
 *
 * @param sys   SISO state-space system
 * @param omega Frequency grid (rad/s)
 * @return NicholsResult open-loop phase vs magnitude locus
 *
 * @note Compare with MATLAB®'s nichols(sys, w).
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr NicholsResult<T> nichols(
    const StateSpace<NX, 1, 1, T, NW, NV>& sys,
    const std::vector<T>&                  omega
) {
    return nichols(bode(sys, omega));
}

/**
 * @brief Compute Nichols chart data for a SISO transfer function
 *
 * @note Compare with MATLAB®'s nichols(sys, w).
 */
template<size_t Nnum, size_t Nden, typename T>
[[nodiscard]] constexpr NicholsResult<T> nichols(
    const TransferFunction<Nnum, Nden, T>& tf,
    const std::vector<T>&                  omega
) {
    return nichols(bode(tf, omega));
}

// ============================================================================
// Singular-value frequency response (sigma)
// ============================================================================

/**
 * @brief Singular values of G(jω) at one frequency
 *
 * @tparam NS Number of singular values stored (= min(NY, NU))
 */
template<size_t NS, typename T = double>
struct SigmaPoint {
    T                  omega{};    ///< Frequency (rad/s)
    damp::array<T, NS> sigma{};    ///< Singular values σ₁ ≥ … ≥ σ_NS (≥ 0)
    damp::array<T, NS> sigma_db{}; ///< 20·log₁₀(σᵢ); −300 dB when σᵢ = 0
};

/**
 * @brief Singular-value frequency response over a grid
 *
 * @tparam NS Number of singular values per sample (= min(NY, NU))
 *
 * @note Compare with MATLAB®'s sigma / svplot data.
 */
template<size_t NS, typename T = double>
struct SigmaResult {
    std::vector<SigmaPoint<NS, T>> points; ///< One sample per frequency

    /**
     * @brief Peak of the largest singular value over the grid
     * @return {σ̄_max, ω_at_peak} or nullopt when empty
     */
    [[nodiscard]] constexpr damp::optional<damp::pair<T, T>> peak_sigma_max() const {
        if (points.empty()) {
            return damp::nullopt;
        }
        T max_s = points[0].sigma[0];
        T at_w = points[0].omega;
        for (size_t i = 1; i < points.size(); ++i) {
            if (points[i].sigma[0] > max_s) {
                max_s = points[i].sigma[0];
                at_w = points[i].omega;
            }
        }
        return damp::pair{max_s, at_w};
    }
};

/**
 * @brief Singular-value frequency response of a (possibly MIMO) state-space system
 *
 * At each @p omega sample evaluates @f$ G @f$ on the stability boundary
 * (continuous: @f$ s=j\omega @f$; discrete: @f$ z=e^{j\omega T_s} @f$) and
 * returns all singular values of @f$ G @f$ in descending order. SISO reduces
 * to @f$ \sigma_1 = |G| @f$ without forming an SVD.
 *
 * Contour poles (@ref eval_frf fails) are omitted — same honest skip policy as
 * bode() / nyquist(); the returned grid may be shorter than @p omega.
 *
 * @param sys   State-space system (any NU, NY)
 * @param omega Frequency grid (rad/s)
 * @return SigmaResult with NS = min(NY, NU) singular values per frequency
 *
 * @note Compare with MATLAB®'s sigma(sys, w).
 * @see norm_hinf() for a peak estimate over an auto-seeded grid
 * @see mat::svd() for the underlying decomposition
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr SigmaResult<(NY < NU ? NY : NU), T> sigma(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const std::vector<T>&                    omega
) {
    using Cplx = damp::complex<T>;
    constexpr size_t NS = (NY < NU) ? NY : NU;

    SigmaResult<NS, T> result;
    result.points.reserve(omega.size());

    for (const auto& w : omega) {
        // Same boundary dispatch as bode()/nyquist() (frf_s_or_z is defined later).
        Cplx s_or_z{T{0}, w};
        if (sys.is_discrete()) {
            const auto [s, c] = damp::sincos(w * sys.Ts);
            s_or_z = Cplx{c, s};
        }
        const auto G_opt = eval_frf(sys, s_or_z);
        if (!G_opt) {
            continue; // pole on contour — skip
        }
        const Matrix<NY, NU, Cplx> G = *G_opt;

        SigmaPoint<NS, T> pt{};
        pt.omega = w;

        if constexpr (NY == 1 && NU == 1) {
            const T s0 = damp::abs(G(0, 0));
            pt.sigma[0] = s0;
            pt.sigma_db[0] = s0 > T{0} ? T{20} * damp::log10(s0) : T{-300};
        } else {
            const auto svd_res = mat::svd(G);
            for (size_t k = 0; k < NS; ++k) {
                const T sk = svd_res.singular_values[k];
                pt.sigma[k] = sk;
                pt.sigma_db[k] = sk > T{0} ? T{20} * damp::log10(sk) : T{-300};
            }
        }

        result.points.push_back(pt);
    }

    return result;
}

/**
 * @brief Singular-value frequency response of a SISO transfer function
 *
 * Returns an empty result if companion realization fails.
 *
 * @note Compare with MATLAB®'s sigma(sys, w).
 */
template<size_t Nnum, size_t Nden, typename T>
[[nodiscard]] constexpr SigmaResult<1, T> sigma(
    const TransferFunction<Nnum, Nden, T>& tf,
    const std::vector<T>&                  omega
) {
    const auto ss = tf.to_state_space();
    if (!ss) {
        return {};
    }
    return sigma(*ss, omega);
}


// ============================================================================
// Impedance / Middlebrook Stability Analysis
// ============================================================================

/**
 * @brief Result of impedance frequency response evaluation
 *
 * Contains complex impedance at each frequency point, plus
 * magnitude/phase representation for Bode-style plotting.
 */
template<typename T = double>
struct ImpedanceResult {
    struct Point {
        T                omega{};        ///< Frequency (rad/s)
        damp::complex<T> Z{};            ///< Complex impedance
        T                magnitude{};    ///< |Z| (ohms)
        T                magnitude_db{}; ///< |Z| in dB
        T                phase_deg{};    ///< Phase of Z (degrees)
    };
    std::vector<Point> points;
};

/**
 * @brief Result of Middlebrook minor loop gain analysis
 *
 * The minor loop gain is T_m(s) = Z_s(s) / Z_L(s).
 * The interconnected system is stable if T_m satisfies the Nyquist criterion,
 * which for a stable T_m reduces to:
 *   - Gain margin: |T_m| should be well below 0 dB at -180° phase crossing
 *   - Phase margin: phase(T_m) should be well above -180° at 0 dB crossing
 *
 * A sufficient (conservative) condition: |T_m(jω)| < 1 for all ω,
 * i.e. |Z_s| < |Z_L| at every frequency.
 */
template<typename T = double>
struct MiddlebrookResult {
    BodeResult<T> minor_loop_gain; ///< Bode data for T_m = Z_s / Z_L

    ImpedanceResult<T> source_impedance; ///< Z_s frequency response
    ImpedanceResult<T> load_impedance;   ///< Z_L frequency response

    /**
     * @brief Check if the sufficient stability condition holds
     *
     * Returns true if |Z_s(jω)| < |Z_L(jω)| at every frequency point,
     * i.e. the minor loop gain magnitude is strictly below 0 dB everywhere.
     * This is a conservative criterion (sufficient but not necessary).
     */
    [[nodiscard]] constexpr bool is_stable_sufficient() const {
        for (const auto& pt : minor_loop_gain.points) {
            if (pt.magnitude_db >= T{0}) {
                return false;
            }
        }
        return !minor_loop_gain.points.empty();
    }

    /**
     * @brief Gain margin of the minor loop gain
     *
     * Gain margin at the -180° phase crossing of T_m.
     * Positive dB = stable.
     *
     * @return {gain_margin_dB, frequency} or nullopt
     */
    [[nodiscard]] constexpr damp::optional<damp::pair<T, T>> gain_margin() const {
        return minor_loop_gain.gain_margin();
    }

    /**
     * @brief Phase margin of the minor loop gain
     *
     * Phase margin at the 0 dB crossing of T_m.
     * Positive degrees = stable.
     *
     * @return {phase_margin_deg, frequency} or nullopt
     */
    [[nodiscard]] constexpr damp::optional<damp::pair<T, T>> phase_margin() const {
        return minor_loop_gain.phase_margin();
    }

    /**
     * @brief Find the worst-case (smallest) impedance ratio across all frequencies
     *
     * Returns the minimum |Z_L|/|Z_s| ratio, i.e. how close the system
     * is to violating the sufficient condition. Values > 1 mean the
     * sufficient condition holds at all frequencies.
     *
     * @return {min_ratio, frequency_of_worst_case}
     */
    [[nodiscard]] constexpr damp::pair<T, T> worst_case_margin() const {
        T min_ratio = std::numeric_limits<T>::max();
        T worst_freq = T{0};
        for (size_t i = 0; i < minor_loop_gain.points.size(); ++i) {
            T ratio = minor_loop_gain.points[i].magnitude;
            if (ratio > T{0}) {
                T inv_ratio = T{1} / ratio; // |Z_L|/|Z_s|
                if (inv_ratio < min_ratio) {
                    min_ratio = inv_ratio;
                    worst_freq = minor_loop_gain.points[i].omega;
                }
            }
        }
        return {min_ratio, worst_freq};
    }
};

/**
 * @brief FRF evaluation point: jω (continuous) or e^{jωTs} (discrete)
 *
 * Matches the auto-dispatch used by bode() / nyquist() / loop_response().
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr damp::complex<T> frf_s_or_z(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        omega
) {
    using Cplx = damp::complex<T>;
    if (sys.is_discrete()) {
        const auto [s, c] = damp::sincos(omega * sys.Ts);
        return Cplx{c, s};
    }
    return Cplx{T{0}, omega};
}

/**
 * @brief Compute impedance frequency response from a SISO admittance system
 *
 * Given a system G = I/V (admittance: current out per voltage in),
 * computes Z = 1/G = V/I at each frequency. Continuous systems use s = jω;
 * discrete systems (Ts > 0) evaluate on the unit circle z = e^{jωTs}, matching
 * bode().
 *
 * @param admittance_sys  SISO state-space system where output = current, input = voltage
 * @param omega           Vector of frequencies (rad/s)
 * @return ImpedanceResult with Z at each frequency
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr ImpedanceResult<T> impedance(
    const StateSpace<NX, 1, 1, T, NW, NV>& admittance_sys,
    const std::vector<T>&                  omega
) {
    using Cplx = damp::complex<T>;
    ImpedanceResult<T> result;
    result.points.reserve(omega.size());

    for (const auto& w : omega) {
        const Cplx s_or_z = frf_s_or_z(admittance_sys, w);
        auto       G_opt = eval_frf(admittance_sys, s_or_z);
        if (!G_opt) {
            continue;
        }
        auto G = *G_opt;
        Cplx Y = G(0, 0);              // Admittance Y = I/V
        Cplx Z = Cplx{T{1}, T{0}} / Y; // Impedance Z = 1/Y

        T mag = damp::abs(Z);
        T phase_rad = damp::arg(Z);
        T mag_db = mag > T{0} ? T{20} * damp::log10(mag) : T{-300};
        T phase_deg = phase_rad * T{180} / damp::numbers::pi_v<T>;

        result.points.push_back({w, Z, mag, mag_db, phase_deg});
    }
    return result;
}

/**
 * @brief Compute impedance frequency response from a SISO impedance transfer function
 *
 * Given Z directly as a state-space system (voltage out per current in),
 * evaluates the frequency response. Continuous systems use s = jω; discrete
 * systems evaluate z = e^{jωTs} (same dispatch as bode()).
 *
 * @param impedance_sys  SISO state-space system where output = voltage, input = current
 * @param omega          Vector of frequencies (rad/s)
 * @return ImpedanceResult with Z at each frequency
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr ImpedanceResult<T> impedance_direct(
    const StateSpace<NX, 1, 1, T, NW, NV>& impedance_sys,
    const std::vector<T>&                  omega
) {
    using Cplx = damp::complex<T>;
    ImpedanceResult<T> result;
    result.points.reserve(omega.size());

    for (const auto& w : omega) {
        const Cplx s_or_z = frf_s_or_z(impedance_sys, w);
        auto       G_opt = eval_frf(impedance_sys, s_or_z);
        if (!G_opt) {
            continue;
        }
        auto G = *G_opt;
        Cplx Z = G(0, 0);

        T mag = damp::abs(Z);
        T phase_rad = damp::arg(Z);
        T mag_db = mag > T{0} ? T{20} * damp::log10(mag) : T{-300};
        T phase_deg = phase_rad * T{180} / damp::numbers::pi_v<T>;

        result.points.push_back({w, Z, mag, mag_db, phase_deg});
    }
    return result;
}

/**
 * @brief Middlebrook stability analysis for cascaded source-load systems
 *
 * Evaluates the minor loop gain T_m = Z_s / Z_L over a range of frequencies
 * and returns stability margins. Continuous systems use s = jω; discrete
 * systems evaluate z = e^{jωTs} (same dispatch as bode()).
 *
 * The source and load are specified as admittance systems (current/voltage):
 *   - source_admittance: G_s = I_s/V_s, so Z_s = 1/G_s
 *   - load_admittance:   G_L = I_L/V_L, so Z_L = 1/G_L
 *
 * @param source_admittance  Source admittance SISO system
 * @param load_admittance    Load admittance SISO system
 * @param omega              Vector of frequencies (rad/s)
 * @return MiddlebrookResult with minor loop gain Bode data and margins
 */
template<size_t NX_S, size_t NW_S, size_t NV_S, size_t NX_L, size_t NW_L, size_t NV_L, typename T>
[[nodiscard]] constexpr MiddlebrookResult<T> middlebrook(
    const StateSpace<NX_S, 1, 1, T, NW_S, NV_S>& source_admittance,
    const StateSpace<NX_L, 1, 1, T, NW_L, NV_L>& load_admittance,
    const std::vector<T>&                        omega
) {
    using Cplx = damp::complex<T>;
    MiddlebrookResult<T> result;
    result.minor_loop_gain.points.reserve(omega.size());
    result.source_impedance.points.reserve(omega.size());
    result.load_impedance.points.reserve(omega.size());

    for (const auto& w : omega) {
        const Cplx s_or_z_s = frf_s_or_z(source_admittance, w);
        const Cplx s_or_z_L = frf_s_or_z(load_admittance, w);

        // Source impedance: Z_s = 1/Y_s
        auto G_s_opt = eval_frf(source_admittance, s_or_z_s);
        if (!G_s_opt) {
            continue;
        }
        auto G_s = *G_s_opt;
        Cplx Y_s = G_s(0, 0);
        Cplx Z_s = Cplx{T{1}, T{0}} / Y_s;

        // Load impedance: Z_L = 1/Y_L
        auto G_L_opt = eval_frf(load_admittance, s_or_z_L);
        if (!G_L_opt) {
            continue;
        }
        auto G_L = *G_L_opt;
        Cplx Y_L = G_L(0, 0);
        Cplx Z_L = Cplx{T{1}, T{0}} / Y_L;

        // Minor loop gain: T_m = Z_s / Z_L = Y_L / Y_s
        Cplx Tm = Z_s / Z_L;

        // Source impedance point
        T zs_mag = damp::abs(Z_s);
        T zs_phase = damp::arg(Z_s) * T{180} / damp::numbers::pi_v<T>;
        T zs_db = zs_mag > T{0} ? T{20} * damp::log10(zs_mag) : T{-300};
        result.source_impedance.points.push_back({w, Z_s, zs_mag, zs_db, zs_phase});

        // Load impedance point
        T zl_mag = damp::abs(Z_L);
        T zl_phase = damp::arg(Z_L) * T{180} / damp::numbers::pi_v<T>;
        T zl_db = zl_mag > T{0} ? T{20} * damp::log10(zl_mag) : T{-300};
        result.load_impedance.points.push_back({w, Z_L, zl_mag, zl_db, zl_phase});

        // Minor loop gain point
        T tm_mag = damp::abs(Tm);
        T tm_phase = damp::arg(Tm) * T{180} / damp::numbers::pi_v<T>;
        T tm_db = tm_mag > T{0} ? T{20} * damp::log10(tm_mag) : T{-300};
        result.minor_loop_gain.points.push_back({w, tm_mag, tm_db, tm_phase});
    }

    return result;
}

/**
 * @brief Middlebrook analysis from pre-computed impedance data
 *
 * For cases where Z_s and Z_L are already known (e.g., from measurement
 * or from separate impedance models).
 *
 * @param Z_source  Source impedance frequency response
 * @param Z_load    Load impedance frequency response
 * @return MiddlebrookResult with minor loop gain and margins
 */
template<typename T>
[[nodiscard]] constexpr MiddlebrookResult<T> middlebrook(
    const ImpedanceResult<T>& Z_source,
    const ImpedanceResult<T>& Z_load
) {
    MiddlebrookResult<T> result;
    result.source_impedance = Z_source;
    result.load_impedance = Z_load;

    size_t n = damp::min(Z_source.points.size(), Z_load.points.size());
    result.minor_loop_gain.points.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        auto& zs = Z_source.points[i];
        auto& zl = Z_load.points[i];

        // T_m = Z_s / Z_L
        auto Tm = zs.Z / zl.Z;
        T    tm_mag = damp::abs(Tm);
        T    tm_phase = damp::arg(Tm) * T{180} / damp::numbers::pi_v<T>;
        T    tm_db = tm_mag > T{0} ? T{20} * damp::log10(tm_mag) : T{-300};

        result.minor_loop_gain.points.push_back({zs.omega, tm_mag, tm_db, tm_phase});
    }

    return result;
}

} // namespace analysis
} // namespace damp
