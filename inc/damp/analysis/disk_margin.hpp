// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file disk_margin.hpp
 * @brief Disk margins: the simultaneous gain and phase variation a loop tolerates, SISO and loop-at-a-time MIMO (host-only)
 */

#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/analysis/norms.hpp"
#include "damp/backend.hpp"
#include "damp/design/stability.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace analysis {

/**
 * @brief Disk margin of a SISO loop: the disk size α and the gain and phase margins it guarantees
 *
 * The loop stays stable for every multiplicative perturbation f of L inside the disk
 * D(α, σ) = { (1 + (1−σ)δ/2) / (1 − (1+σ)δ/2) : |δ| < α }: gain variations in
 * [gain_low, gain_high] and phase variations of ± phase_margin_deg, also applied together.
 */
template<typename T = double>
struct [[nodiscard]] DiskMarginResult {
    T    alpha{0};                   ///< Disk margin α (0: the nominal closed loop is unstable or ill-posed)
    T    skew{0};                    ///< Skew σ of the disk (0: balanced, 1: from S, −1: from T)
    T    gain_low{1};                ///< Smallest gain factor tolerated (linear, ≥ 0)
    T    gain_high{1};               ///< Largest gain factor tolerated (linear); the largest finite T when gain_high_unbounded
    bool gain_high_unbounded{false}; ///< The gain can grow without bound (MATLAB®'s Inf)
    T    phase_margin_deg{0};        ///< Phase variation tolerated, ± degrees
    T    frequency{0};               ///< [rad/s] Where |S + (σ−1)/2| peaks; the largest finite T when the peak is approached as ω → ∞
    bool stable{false};              ///< Nominal closed loop (unity negative feedback) stable
    bool success{false};             ///< The margin was computed: stable, and the peak of |S + (σ−1)/2| converged
};

/**
 * @brief Loop-at-a-time disk margins of a square MIMO loop: one SISO disk margin per channel
 *
 * channel[i] is the disk margin of loop i broken with the other loops closed, so it holds for
 * a perturbation in that channel alone.
 */
template<size_t N, typename T = double>
struct [[nodiscard]] LoopDiskMarginResult {
    damp::array<DiskMarginResult<T>, N> channel{};      ///< Disk margin of each channel, the others closed
    bool                                stable{false};  ///< Nominal MIMO closed loop (unity negative feedback) stable
    bool                                success{false}; ///< Every channel's margin was computed
};

namespace detail {

/**
 * @brief Sensitivity S = (I + L)⁻¹ of the loop under unity negative feedback, as a realization
 *
 * A − B(I + D)⁻¹C, B(I + D)⁻¹, −(I + D)⁻¹C, (I + D)⁻¹, by solves with I + D; nullopt when
 * I + D is singular (the loop is not well posed).
 */
template<size_t NX, size_t N, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr damp::optional<StateSpace<NX, N, N, T>> sensitivity_realization(const StateSpace<NX, N, N, T, NW, NV>& loop) {
    const Matrix<N, N, T> I = Matrix<N, N, T>::identity();
    const Matrix<N, N, T> M = I + loop.D;
    const auto            MinvC = mat::solve(M, loop.C);
    const auto            MinvT_Bt = mat::solve(Matrix<N, N, T>(M.transpose()), Matrix<N, NX, T>(loop.B.transpose()));
    const auto            Minv = mat::solve(M, I);
    if (!MinvC || !MinvT_Bt || !Minv) {
        return damp::nullopt;
    }
    return StateSpace<NX, N, N, T>{
        .A = loop.A - (loop.B * (*MinvC)),
        .B = MinvT_Bt->transpose(),
        .C = -(*MinvC),
        .D = *Minv,
        .Ts = loop.Ts,
    };
}

/** @brief Disk margin from the peak of |S + (σ−1)/2|: α, the gain interval and the phase margin. */
template<typename T>
[[nodiscard]] constexpr DiskMarginResult<T> disk_margin_from_peak(const HinfNormResult<T>& peak, T skew) {
    DiskMarginResult<T> out{.skew = skew, .stable = true};
    if (!peak.success || !(peak.gain > T{0})) {
        return out;
    }
    out.success = true;
    out.frequency = peak.frequency;
    out.alpha = T{1} / peak.gain;

    //! Gain margins where the disk meets the real axis, phase margin where it meets the unit circle.
    const T a = out.alpha;
    out.gain_low = damp::max(T{0}, (T{2} - (a * (T{1} - skew))) / (T{2} + (a * (T{1} + skew))));
    const T den_high = T{2} - (a * (T{1} + skew));
    out.gain_high_unbounded = !(den_high > T{0});
    out.gain_high = out.gain_high_unbounded ? std::numeric_limits<T>::max() : (T{2} + (a * (T{1} - skew))) / den_high;
    const T den_phase = T{4} + (a * a * (T{1} - (skew * skew)));
    const T cos_phi = den_phase > T{0} ? damp::clamp((T{4} - (a * a * (T{1} + (skew * skew)))) / den_phase, T{-1}, T{1}) : T{-1};
    out.phase_margin_deg = damp::acos(cos_phi) * T{180} / damp::numbers::pi_v<T>;
    return out;
}

/** @brief Channel @p i of a sensitivity realization plus the disk offset (σ − 1)/2: S_ii + (σ − 1)/2. */
template<size_t NX, size_t N, typename T>
[[nodiscard]] constexpr StateSpace<NX, 1, 1, T> disk_channel(const StateSpace<NX, N, N, T>& S, size_t i, T skew) {
    StateSpace<NX, 1, 1, T> out{.A = S.A, .Ts = S.Ts};
    for (size_t k = 0; k < NX; ++k) {
        out.B(k, 0) = S.B(k, i);
        out.C(0, k) = S.C(i, k);
    }
    out.D(0, 0) = S.D(i, i) + ((skew - T{1}) / T{2});
    return out;
}

} // namespace detail

/**
 * @brief Disk margin of a SISO loop L under unity negative feedback
 *
 * With S = 1/(1 + L), the loop tolerates every perturbation in the disk D(α, σ) up to
 * @f[
 *   \alpha_{\max} = \frac{1}{\max_\omega \left| S(j\omega) + \tfrac{\sigma - 1}{2} \right|},
 * @f]
 * which for the balanced disk (σ = 0) is @f$ 2 / \lVert S - T \rVert_\infty @f$. The disk
 * meets the real axis at the gain margins
 * @f$ g_{\min} = \frac{2 - \alpha(1-\sigma)}{2 + \alpha(1+\sigma)} @f$ and
 * @f$ g_{\max} = \frac{2 + \alpha(1-\sigma)}{2 - \alpha(1+\sigma)} @f$ (unbounded once the
 * denominator reaches 0), and the unit circle at the phase margin
 * @f$ \cos\phi = \frac{4 - \alpha^2(1+\sigma^2)}{4 + \alpha^2(1-\sigma^2)} @f$; for σ = 0 that
 * is @f$ \phi = 2\arctan(\alpha/2) @f$. Unlike the classical margins, these hold for gain
 * and phase varying together.
 *
 * The peak is the H∞ norm of @f$ S + (\sigma - 1)/2 @f$, computed by norm_hinf() on the
 * realization @f$ (A - B(1+D)^{-1}C,\; B(1+D)^{-1},\; -(1+D)^{-1}C,\; (1+D)^{-1} + (\sigma-1)/2) @f$
 * to the relative tolerance default_tol(), so a narrow peak between frequency samples is not
 * missed. A nominal closed loop that is not stable, or 1 + D = 0, has no disk margin:
 * @c alpha = 0, @c stable = false.
 *
 * @note Compare with MATLAB®'s DM = diskmargin(L, sigma): DiskMargin, GainMargin, PhaseMargin
 *       and Frequency. matlab::diskmargin is the alias.
 * @see "An Introduction to Disk Margins" (Seiler, Packard & Gahinet, 2020), IEEE Control
 *      Systems Magazine 40(5), https://doi.org/10.1109/MCS.2020.3005277
 * @see loop_at_a_time_disk_margin() for a MIMO loop
 *
 * @param loop   Open-loop SISO L (continuous or discrete)
 * @param skew   Disk skew σ (default 0, the balanced disk)
 * @return Disk margin α, gain interval, phase margin and critical frequency
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr DiskMarginResult<T> disk_margin(const StateSpace<NX, 1, 1, T, NW, NV>& loop, std::type_identity_t<T> skew = T{0}) {
    const auto S = detail::sensitivity_realization(loop);
    if (!S || !stability::is_stable(*S)) {
        return DiskMarginResult<T>{.skew = skew};
    }
    return detail::disk_margin_from_peak(norm_hinf(detail::disk_channel(*S, 0, skew)), skew);
}

/**
 * @brief Loop-at-a-time disk margins of a square MIMO loop L under unity negative feedback
 *
 * Breaking loop i with the others closed leaves the SISO loop @f$ L_i @f$ with
 * @f$ 1/(1 + L_i) = S_{ii} @f$, the i-th diagonal entry of @f$ S = (I + L)^{-1} @f$. Each
 * channel's margin is therefore disk_margin()'s computation on @f$ S_{ii} @f$: the peak of
 * @f$ |S_{ii} + (\sigma - 1)/2| @f$ by norm_hinf(), and the same gain and phase margins.
 * Where the loop is broken (plant input or output) is the caller's choice of L. A diagonal
 * L gives each channel's SISO disk margin.
 *
 * Each margin holds for a perturbation in its own channel with the other channels at their
 * nominal values. The multiloop (simultaneous) disk margin, for independent perturbations in
 * every channel at once, is a structured singular value (μ) problem and is not computed: it
 * can be much smaller than every loop-at-a-time margin, as for a spinning satellite whose
 * channels each tolerate any gain increase while a small simultaneous perturbation
 * destabilizes the loop.
 *
 * @note Compare with MATLAB®'s DM = diskmargin(L, sigma) for a MIMO L (the loop-at-a-time
 *       margins; the multiloop MM is not computed). matlab::diskmargin is the alias.
 * @see "An Introduction to Disk Margins" (Seiler, Packard & Gahinet, 2020), IEEE Control
 *      Systems Magazine 40(5), https://doi.org/10.1109/MCS.2020.3005277
 * @see disk_margin()
 *
 * @param loop Open-loop square MIMO L (continuous or discrete), N inputs and N outputs
 * @param skew Disk skew σ (default 0, the balanced disk)
 * @return Per-channel disk margins, nominal stability and success
 */
template<size_t NX, size_t N, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr LoopDiskMarginResult<N, T> loop_at_a_time_disk_margin(
    const StateSpace<NX, N, N, T, NW, NV>& loop,
    std::type_identity_t<T>                skew = T{0}
) {
    LoopDiskMarginResult<N, T> out{};
    for (auto& ch : out.channel) {
        ch.skew = skew;
    }
    const auto S = detail::sensitivity_realization(loop);
    if (!S || !stability::is_stable(*S)) {
        return out;
    }
    out.stable = true;
    out.success = true;
    for (size_t i = 0; i < N; ++i) {
        out.channel[i] = detail::disk_margin_from_peak(norm_hinf(detail::disk_channel(*S, i, skew)), skew);
        out.success = out.success && out.channel[i].success;
    }
    return out;
}

} // namespace analysis
} // namespace damp
