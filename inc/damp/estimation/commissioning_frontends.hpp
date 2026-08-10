// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file commissioning_frontends.hpp
 * @brief D21 autotune front-ends over FRF + compensator + safety primitives.
 *
 * - design::resonance_autotune_from_frf — peaks→notch, valleys→band-pass
 * - design::frf_pid_autotune — margins / target crossover → PID gains
 *
 * These are pure design reducers (Tier 1): they do not own the experiment loop.
 * Pair with FrequencyResponseEstimator, SuccessiveCompensatorCommissioner, and
 * ExperimentSafety on the runtime path.
 *
 * @see frequency_response.hpp, successive_compensator.hpp, experiment_safety.hpp
 * @see design::pid_from_bandwidth
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/estimation/frequency_response.hpp"
#include "damp/estimation/successive_compensator.hpp"
#include "damp/filters/filters.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp"

namespace damp {
namespace design {

/**
 * @brief Result of resonance / anti-resonance compensator design from an FRF.
 *
 * @tparam NMax Compensator capacity
 * @tparam T    Scalar type
 */
template<std::size_t NMax, typename T = float>
struct ResonanceAutotuneResult {
    damp::array<ResonantMode<T>, NMax>      features{};
    damp::array<CompensatorKind, NMax>      kinds{};
    damp::array<SecondOrderCoeffs<T>, NMax> sections{};
    std::size_t                             count{0};
    bool                                    success{false};

    template<typename U>
    [[nodiscard]] constexpr ResonanceAutotuneResult<NMax, U> as() const {
        ResonanceAutotuneResult<NMax, U> out{};
        out.count = count;
        out.success = success;
        for (std::size_t i = 0; i < NMax; ++i) {
            out.features[i] = features[i].template as<U>();
            out.kinds[i] = kinds[i];
            out.sections[i] = sections[i].template as<U>();
        }
        return out;
    }
};

/**
 * @brief Map FRF peaks and valleys to a fixed bank of biquad compensators.
 *
 * Peaks → notch, valleys → band-pass (or peaking per @p cfg.valley_map).
 * Does not run an experiment; feed an already-measured FRF table.
 *
 * @tparam NFreq  FRF length
 * @tparam NMax   Max compensators
 * @tparam T      Scalar type
 */
template<std::size_t NFreq, std::size_t NMax, typename T = float>
[[nodiscard]] constexpr ResonanceAutotuneResult<NMax, T> resonance_autotune_from_frf(
    const damp::array<FrfPoint<T>, NFreq>& table,
    const SuccessiveCompensatorConfig<T>&  cfg = {},
    const ModeExtractorConfig<T>&          mode_cfg = {}
) {
    ResonanceAutotuneResult<NMax, T> out{};
    if (!cfg.valid()) {
        return out;
    }

    constexpr std::size_t MaxP = NMax;
    constexpr std::size_t MaxV = NMax;
    const auto            features = extract_frf_features<NFreq, MaxP, MaxV, T>(table, mode_cfg);

    SuccessiveCompensatorCommissioner<NMax, T> cm{cfg};
    (void)cm.assign_from_frf_features(features);

    out.count = cm.count();
    for (std::size_t i = 0; i < out.count; ++i) {
        out.features[i] = cm.feature(i);
        out.kinds[i] = cm.kind(i);
        out.sections[i] = compensator_from_feature(cm.feature(i), cfg);
    }
    out.success = out.count > 0;
    return out;
}

/**
 * @brief Configuration for FRF-based PID gain selection.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
struct FrfPidAutotuneConfig {
    T       target_crossover_hz{T{0}}; ///< Desired ωc / 2π [Hz]; 0 → use measured gain crossover
    T       target_phase_margin_deg{T{60}};
    T       Ts{static_cast<T>(0.001)};
    PIDType type{PIDType::PI};

    [[nodiscard]] constexpr bool valid() const {
        return (Ts > T{0}) && damp::isfinite(target_phase_margin_deg)
            && damp::isfinite(target_crossover_hz);
    }
};

/**
 * @brief FRF PID autotune hand-off (gains + success)
 */
template<typename T = float>
struct FrfPidAutotuneResult {
    PIDResult<T> pid{};
    bool         success{false};

    template<typename U>
    [[nodiscard]] constexpr FrfPidAutotuneResult<U> as() const {
        return {pid.template as<U>(), success};
    }
};

/**
 * @brief PID gains from an open-loop FRF table + margin / bandwidth targets.
 *
 * 1. Read margins via @ref margins_from_frf.
 * 2. Choose crossover frequency f_c (config target, else measured gain crossover).
 * 3. If a plant FRF bin is near f_c, scale so loop |K·G| ≈ 1 at f_c for a pure
 *    gain (stored in result as diagnostic); primary gains come from
 *    @ref pid_from_bandwidth(ωc, PM_target, Ts).
 *
 * @tparam NFreq Table length
 * @tparam T     Scalar type
 * @param table  Open-loop plant (or loop) FRF
 * @param cfg    Targets
 * @return FrfPidAutotuneResult with success when a positive crossover was found
 *         and gains are finite
 */
template<std::size_t NFreq, typename T = float>
[[nodiscard]] constexpr FrfPidAutotuneResult<T> frf_pid_autotune(
    const damp::array<FrfPoint<T>, NFreq>& table,
    const FrfPidAutotuneConfig<T>&         cfg = {}
) {
    FrfPidAutotuneResult<T> out{};
    if (!cfg.valid()) {
        return out;
    }

    const auto margins = margins_from_frf(table);

    T f_c = cfg.target_crossover_hz;
    if (!(f_c > T{0})) {
        if (margins.has_phase_margin && margins.gain_crossover_hz > T{0}) {
            f_c = margins.gain_crossover_hz;
        } else {
            // Fall back: geometric mid of first/last valid bin
            T f_lo = T{0};
            T f_hi = T{0};
            for (std::size_t i = 0; i < NFreq; ++i) {
                if (table[i].valid && table[i].freq_hz > T{0}) {
                    if (!(f_lo > T{0})) {
                        f_lo = table[i].freq_hz;
                    }
                    f_hi = table[i].freq_hz;
                }
            }
            if (f_lo > T{0} && f_hi > f_lo) {
                f_c = damp::sqrt(f_lo * f_hi);
            }
        }
    }
    if (!(f_c > T{0})) {
        return out;
    }

    const T w_c = T{2} * damp::numbers::pi_v<T> * f_c;
    out.pid = pid_from_bandwidth(w_c, cfg.target_phase_margin_deg, cfg.Ts, cfg.type);

    // Optional plant-gain scale: find nearest FRF bin to f_c and scale Kp so
    // |Kp * G| ~ 1 when the designer assumed unit plant.
    T    mag = T{1};
    T    best_df = T{0};
    bool found = false;
    for (std::size_t i = 0; i < NFreq; ++i) {
        if (!table[i].valid) {
            continue;
        }
        const T df = damp::abs(table[i].freq_hz - f_c);
        if (!found || df < best_df) {
            best_df = df;
            mag = table[i].magnitude;
            found = true;
        }
    }
    if (found && mag > default_tol<T>()) {
        const T scale = T{1} / mag;
        out.pid.Kp *= scale;
        out.pid.Ki *= scale;
        out.pid.Kd *= scale;
    }
    out.success = damp::isfinite(out.pid.Kp) && damp::isfinite(out.pid.Ki) && damp::isfinite(out.pid.Kd)
               && (out.pid.Kp > T{0} || out.pid.Ki > T{0} || out.pid.Kd > T{0});
    return out;
}

} // namespace design
} // namespace damp
