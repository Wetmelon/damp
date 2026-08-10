// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file frequency_response.hpp
 * @brief On-target FRF estimation and modal feature extraction.
 *
 * Stage 2–3 of the excite → measure → reduce → map → commit commissioning
 * pipeline ([decisions.md] D21). Allocation-free, fixed-size tables so the
 * full FRF lives on the target.
 *
 * - FrequencyResponseEstimator — lock-in demodulation via a pair of
 *   Goertzel bins on (u, y), H1 estimate G = S_yu / S_uu, coherence γ².
 *   Sequential path for SteppedSine; parallel bank for MultiSine.
 * - extract_modes / ModeExtractorResult — reduce an FRF table to
 *   resonances (ωₙ, ζ, peak gain) via half-power bandwidth. Feeds
 *   `input_shaper` and harmonic suppressors.
 *
 * Example: stepped-sine FRF of a known first-order plant, then mode pick.
 * @code
 * #include "damp/estimation/frequency_response.hpp"
 * #include "damp/estimation/excitation.hpp"
 *
 * using namespace damp;
 *
 * constexpr std::size_t Nf = 8;
 * FrequencyResponseEstimator<Nf, float> frf(1000.0f);
 * // For each frequency f: begin_bin → push(u,y) until true → table filled.
 * // Then: auto modes = design::extract_modes<Nf, 2>(frf.table());
 * @endcode
 *
 * @see excitation.hpp for SteppedSine / MultiSine stimuli
 * @see filters/spectral.hpp for the underlying Goertzel primitive
 * @see "System Identification" (Ljung, 1999), §6 (frequency-domain ID)
 * @see "Vibration Testing" (McConnell, 1995) — half-power bandwidth damping
 */

#include <cstddef>
#include <cstdint>

#include "damp/backend.hpp"
#include "damp/filters/spectral.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp"

namespace damp {

/**
 * @brief One measured frequency-response point (on-target FRF table entry).
 *
 * Stores G(j2πf) as magnitude / phase with an H1 coherence gate.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
struct FrfPoint {
    T    freq_hz{T{0}};   ///< Frequency [Hz]
    T    magnitude{T{0}}; ///< |G(jω)|
    T    phase_rad{T{0}}; ///< arg(G) [rad]
    T    coherence{T{0}}; ///< γ² ∈ [0, 1] (H1 single-input coherence)
    bool valid{false};    ///< true if |U| power was sufficient to form G

    /**
     * @brief Convert to another scalar type.
     * @tparam U Target scalar type.
     */
    template<typename U>
    [[nodiscard]] constexpr FrfPoint<U> as() const {
        return FrfPoint<U>{
            static_cast<U>(freq_hz),
            static_cast<U>(magnitude),
            static_cast<U>(phase_rad),
            static_cast<U>(coherence),
            valid,
        };
    }
};

/**
 * @brief Lock-in FRF estimator over a fixed table of NFreq bins.
 *
 * Accumulates the classical H1 estimate from simultaneous Goertzel DFTs of
 * the applied input u and measured output y:
 *
 *     S_yu = Σ Y · conj(U),   S_uu = Σ |U|²,   S_yy = Σ |Y|²
 *     G    = S_yu / S_uu
 *     γ²   = |S_yu|² / (S_uu · S_yy)
 *
 * Two usage modes:
 *
 * 1. Sequential (stepped sine) — @ref begin_bin then @ref push until it
 *    returns true; repeat for each frequency index.
 * 2. Parallel bank (multi-sine) — @ref configure_bin for each tone, then
 *    @ref push_bank; one coherent block updates every configured bin.
 *
 * O(1) memory per frequency: the table is `array<FrfPoint, NFreq>`.
 *
 * @tparam NFreq Number of frequency bins in the table.
 * @tparam T     Scalar type (default float for target).
 */
template<std::size_t NFreq, typename T = float>
class FrequencyResponseEstimator {
public:
    static_assert(NFreq >= 1, "FrequencyResponseEstimator needs at least one bin");

    constexpr FrequencyResponseEstimator() = default;

    /**
     * @brief Construct with sample rate.
     * @param fs_hz Sample rate [Hz] (> 0). Used by @ref coherent_block_size helpers.
     */
    constexpr explicit FrequencyResponseEstimator(T fs_hz) : fs_(fs_hz) {}

    /**
     * @brief Coherent Goertzel block length for an integer number of cycles.
     *
     * N = round(cycles · fs / f). Prefer frequencies that make this exact.
     *
     * @param freq_hz Target frequency [Hz]
     * @param fs_hz   Sample rate [Hz]
     * @param cycles  Periods in the block (>= 1)
     * @return Block length in samples (>= 1)
     */
    [[nodiscard]] static constexpr std::size_t coherent_block_size(T freq_hz, T fs_hz, std::size_t cycles) {
        if (!(freq_hz > T{0}) || !(fs_hz > T{0}) || cycles == 0) {
            return 1;
        }
        const T n = static_cast<T>(cycles) * fs_hz / freq_hz;
        auto    N = static_cast<std::size_t>(n + static_cast<T>(0.5));
        if (N < 1) {
            N = 1;
        }
        return N;
    }

    /**
     * @brief Coherent block length using the estimator's sample rate.
     */
    [[nodiscard]] constexpr std::size_t coherent_block_size(T freq_hz, std::size_t cycles) const {
        return coherent_block_size(freq_hz, fs_, cycles);
    }

    /**
     * @brief Start sequential measurement of table bin @p index.
     *
     * Averages @p num_blocks successive Goertzel blocks of @p block_size samples
     * (after discarding @p settle_samples). Call @ref push until it returns true.
     *
     * @param index          Table index in [0, NFreq)
     * @param freq_hz        Bin frequency [Hz]
     * @param block_size     Samples per Goertzel block
     * @param num_blocks     Blocks to average for H1 (>= 1)
     * @param settle_samples Leading samples discarded (transients)
     */
    constexpr void begin_bin(
        std::size_t index,
        T           freq_hz,
        std::size_t block_size,
        std::size_t num_blocks = 1,
        std::size_t settle_samples = 0
    ) {
        if (index >= NFreq || !(freq_hz > T{0}) || !(fs_ > T{0}) || block_size == 0
            || num_blocks == 0) {
            sequential_active_ = false;
            return;
        }
        active_index_ = index;
        active_freq_ = freq_hz;
        blocks_target_ = num_blocks;
        blocks_done_ = 0;
        settle_left_ = settle_samples;
        syu_re_ = T{0};
        syu_im_ = T{0};
        suu_ = T{0};
        syy_ = T{0};
        gu_ = Goertzel<T>(freq_hz, fs_, block_size);
        gy_ = Goertzel<T>(freq_hz, fs_, block_size);
        sequential_active_ = true;
        table_[index] = FrfPoint<T>{};
        table_[index].freq_hz = freq_hz;
    }

    /**
     * @brief Configure a parallel-bank bin (MultiSine path).
     *
     * All bank bins should use the same @p block_size so @ref push_bank
     * completes them together. Single-block H1 (no multi-block average).
     *
     * @param index      Table index in [0, NFreq)
     * @param freq_hz    Bin frequency [Hz]
     * @param block_size Samples per Goertzel block
     */
    constexpr void configure_bin(std::size_t index, T freq_hz, std::size_t block_size) {
        if (index >= NFreq || !(freq_hz > T{0}) || block_size == 0 || !(fs_ > T{0})) {
            return;
        }
        bank_u_[index] = Goertzel<T>(freq_hz, fs_, block_size);
        bank_y_[index] = Goertzel<T>(freq_hz, fs_, block_size);
        bank_active_[index] = true;
        table_[index] = FrfPoint<T>{};
        table_[index].freq_hz = freq_hz;
    }

    /**
     * @brief Sequential path: feed one (u, y) sample pair.
     * @return true when the active bin has finished all averaging blocks.
     */
    [[nodiscard]] constexpr bool push(T u, T y) {
        if (!sequential_active_) {
            return false;
        }
        if (settle_left_ > 0) {
            --settle_left_;
            return false;
        }
        const bool done_u = gu_.push(u);
        const bool done_y = gy_.push(y);
        if (!(done_u && done_y)) {
            return false;
        }
        accumulate_block(gu_.real(), gu_.imag(), gy_.real(), gy_.imag());
        ++blocks_done_;
        if (blocks_done_ < blocks_target_) {
            return false;
        }
        finalize_point(active_index_, active_freq_);
        sequential_active_ = false;
        return true;
    }

    /**
     * @brief Parallel path: feed one sample into every configured bank bin.
     * @return true when a full block completes (all active bins updated).
     */
    [[nodiscard]] constexpr bool push_bank(T u, T y) {
        bool any_active = false;
        bool completed = false;
        for (std::size_t i = 0; i < NFreq; ++i) {
            if (!bank_active_[i]) {
                continue;
            }
            any_active = true;
            const bool du = bank_u_[i].push(u);
            const bool dy = bank_y_[i].push(y);
            if (du && dy) {
                // Single-block H1 for this bin.
                syu_re_ = T{0};
                syu_im_ = T{0};
                suu_ = T{0};
                syy_ = T{0};
                accumulate_block(
                    bank_u_[i].real(), bank_u_[i].imag(), bank_y_[i].real(), bank_y_[i].imag()
                );
                finalize_point(i, table_[i].freq_hz);
                completed = true;
            }
        }
        return any_active && completed;
    }

    /// FRF table (fixed size; check @c FrfPoint::valid per entry).
    [[nodiscard]] constexpr const damp::array<FrfPoint<T>, NFreq>& table() const { return table_; }

    /// Mutable table (e.g. host post-processing).
    [[nodiscard]] constexpr damp::array<FrfPoint<T>, NFreq>& table() { return table_; }

    /// One table entry.
    [[nodiscard]] constexpr const FrfPoint<T>& point(std::size_t i) const { return table_[i]; }

    /// Number of table slots (template size).
    [[nodiscard]] static constexpr std::size_t size() { return NFreq; }

    /// Count of currently valid points.
    [[nodiscard]] constexpr std::size_t valid_count() const {
        std::size_t n = 0;
        for (std::size_t i = 0; i < NFreq; ++i) {
            if (table_[i].valid) {
                ++n;
            }
        }
        return n;
    }

    /// Sample rate [Hz] (0 if unset).
    [[nodiscard]] constexpr T sample_rate() const { return fs_; }

    /// Set sample rate [Hz] (must be > 0 before begin_bin / configure_bin).
    constexpr void set_sample_rate(T fs_hz) { fs_ = fs_hz; }

    /// Clear the table and abort any in-progress measurement.
    constexpr void reset() {
        sequential_active_ = false;
        blocks_done_ = 0;
        settle_left_ = 0;
        for (std::size_t i = 0; i < NFreq; ++i) {
            table_[i] = FrfPoint<T>{};
            bank_active_[i] = false;
            bank_u_[i].reset();
            bank_y_[i].reset();
        }
        gu_.reset();
        gy_.reset();
    }

private:
    constexpr void accumulate_block(T ur, T ui, T yr, T yi) {
        // S_yu += Y * conj(U) = (yr + j yi)(ur - j ui)
        syu_re_ += (yr * ur) + (yi * ui);
        syu_im_ += (yi * ur) - (yr * ui);
        suu_ += (ur * ur) + (ui * ui);
        syy_ += (yr * yr) + (yi * yi);
    }

    constexpr void finalize_point(std::size_t index, T freq_hz) {
        FrfPoint<T> p{};
        p.freq_hz = freq_hz;
        const T floor = default_tol<T>() * default_tol<T>();
        if (suu_ <= floor) {
            p.valid = false;
            table_[index] = p;
            return;
        }
        const T gre = syu_re_ / suu_;
        const T gim = syu_im_ / suu_;
        p.magnitude = damp::hypot(gre, gim);
        p.phase_rad = damp::atan2(gim, gre);
        const T num = (syu_re_ * syu_re_) + (syu_im_ * syu_im_);
        const T den = suu_ * syy_;
        if (den > floor) {
            T coh = num / den;
            if (coh < T{0}) {
                coh = T{0};
            }
            if (coh > T{1}) {
                coh = T{1};
            }
            p.coherence = coh;
        }
        p.valid = true;
        table_[index] = p;
    }

    T                               fs_{T{0}};
    damp::array<FrfPoint<T>, NFreq> table_{};

    // Sequential state
    Goertzel<T> gu_{};
    Goertzel<T> gy_{};
    T           active_freq_{T{0}};
    T           syu_re_{T{0}};
    T           syu_im_{T{0}};
    T           suu_{T{0}};
    T           syy_{T{0}};
    std::size_t active_index_{0};
    std::size_t blocks_target_{1};
    std::size_t blocks_done_{0};
    std::size_t settle_left_{0};
    bool        sequential_active_{false};

    // Parallel bank
    damp::array<Goertzel<T>, NFreq> bank_u_{};
    damp::array<Goertzel<T>, NFreq> bank_y_{};
    damp::array<bool, NFreq>        bank_active_{};
};

/**
 * @brief Peak (resonance) or valley (anti-resonance) on an FRF magnitude curve.
 */
enum class FrfFeatureKind : std::uint8_t {
    Peak,  ///< Local |G| maximum — mechanical / filter resonance
    Valley ///< Local |G| minimum — anti-resonance / transmission zero
};

/**
 * @brief One extracted FRF extremum (peak or valley).
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
struct ResonantMode {
    FrfFeatureKind kind{FrfFeatureKind::Peak}; ///< Peak vs valley
    T              omega_n{T{0}};              ///< Natural frequency [rad/s]
    T              frequency_hz{T{0}};         ///< Natural frequency [Hz]
    T              zeta{T{0}};                 ///< Damping ratio (half-power / double-power BW)
    T              peak_gain{T{0}};            ///< |G| at the extremum
    bool           valid{false};               ///< true if extremum + ζ estimate succeeded

    template<typename U>
    [[nodiscard]] constexpr ResonantMode<U> as() const {
        return ResonantMode<U>{
            kind,
            static_cast<U>(omega_n),
            static_cast<U>(frequency_hz),
            static_cast<U>(zeta),
            static_cast<U>(peak_gain),
            valid,
        };
    }
};

/**
 * @brief Configuration for FRF peak / valley extraction.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
struct ModeExtractorConfig {
    T    min_peak_gain{T{0}};      ///< Peaks: ignore if |G| below this
    T    max_valley_gain{T{0}};    ///< Valleys: if > 0, ignore if |G| above this (0 = no cap)
    T    min_prominence{T{0}};     ///< Extremum must differ from neighbors by this amount
    T    min_coherence{T{0}};      ///< If > 0, require point.coherence >= this
    bool require_valid_zeta{true}; ///< Drop extrema where half-power BW cannot be found
};

/**
 * @brief Result of reducing an FRF table to at most MaxModes peaks (resonances).
 *
 * Sorted by descending |G|. Feeds input shapers and notch compensators.
 *
 * @tparam MaxModes Capacity of the mode table.
 * @tparam T        Scalar type.
 */
template<std::size_t MaxModes, typename T = float>
struct ModeExtractorResult {
    damp::array<ResonantMode<T>, MaxModes> modes{};
    std::size_t                            count{0};
    bool                                   success{false};

    template<typename U>
    [[nodiscard]] constexpr ModeExtractorResult<MaxModes, U> as() const {
        ModeExtractorResult<MaxModes, U> out{};
        out.count = count;
        out.success = success;
        for (std::size_t i = 0; i < MaxModes; ++i) {
            out.modes[i] = modes[i].template as<U>();
        }
        return out;
    }
};

/**
 * @brief Combined peak + valley reduction of an FRF table.
 *
 * @tparam MaxPeaks   Capacity for resonances
 * @tparam MaxValleys Capacity for anti-resonances
 * @tparam T          Scalar type
 */
template<std::size_t MaxPeaks, std::size_t MaxValleys, typename T = float>
struct FrfFeatureResult {
    damp::array<ResonantMode<T>, MaxPeaks>   peaks{};
    damp::array<ResonantMode<T>, MaxValleys> valleys{};
    std::size_t                              peak_count{0};
    std::size_t                              valley_count{0};
    bool                                     success{false};

    template<typename U>
    [[nodiscard]] constexpr FrfFeatureResult<MaxPeaks, MaxValleys, U> as() const {
        FrfFeatureResult<MaxPeaks, MaxValleys, U> out{};
        out.peak_count = peak_count;
        out.valley_count = valley_count;
        out.success = success;
        for (std::size_t i = 0; i < MaxPeaks; ++i) {
            out.peaks[i] = peaks[i].template as<U>();
        }
        for (std::size_t i = 0; i < MaxValleys; ++i) {
            out.valleys[i] = valleys[i].template as<U>();
        }
        return out;
    }
};

/**
 * @brief Open-loop margins read from a discrete FRF table (on-target).
 *
 * Frequencies are in Hz (table units). Phase margins in degrees; gain margin
 * in dB (positive ⇒ stable for classical loop).
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
struct FrfMargins {
    T    phase_margin_deg{T{0}};   ///< PM at gain crossover
    T    gain_crossover_hz{T{0}};  ///< |G| = 1 crossing
    T    gain_margin_db{T{0}};     ///< GM at phase crossover (−180°)
    T    phase_crossover_hz{T{0}}; ///< phase = −π crossing
    bool has_phase_margin{false};
    bool has_gain_margin{false};
    bool success{false}; ///< true if at least one margin was found

    template<typename U>
    [[nodiscard]] constexpr FrfMargins<U> as() const {
        return FrfMargins<U>{
            static_cast<U>(phase_margin_deg),
            static_cast<U>(gain_crossover_hz),
            static_cast<U>(gain_margin_db),
            static_cast<U>(phase_crossover_hz),
            has_phase_margin,
            has_gain_margin,
            success,
        };
    }
};

namespace design {
namespace detail_frf {

template<std::size_t NFreq, typename T>
[[nodiscard]] constexpr bool point_ok(
    const damp::array<FrfPoint<T>, NFreq>& table, std::size_t i, const ModeExtractorConfig<T>& cfg
) {
    const auto& p = table[i];
    if (!p.valid || !(p.freq_hz > T{0})) {
        return false;
    }
    if (cfg.min_coherence > T{0} && p.coherence < cfg.min_coherence) {
        return false;
    }
    return true;
}

template<std::size_t N, typename T>
constexpr void sort_by_gain_desc(damp::array<ResonantMode<T>, N>& c, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t best = i;
        for (std::size_t j = i + 1; j < n; ++j) {
            if (c[j].peak_gain > c[best].peak_gain) {
                best = j;
            }
        }
        if (best != i) {
            const ResonantMode<T> tmp = c[i];
            c[i] = c[best];
            c[best] = tmp;
        }
    }
}

template<std::size_t N, typename T>
constexpr void sort_by_gain_asc(damp::array<ResonantMode<T>, N>& c, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t best = i;
        for (std::size_t j = i + 1; j < n; ++j) {
            if (c[j].peak_gain < c[best].peak_gain) {
                best = j;
            }
        }
        if (best != i) {
            const ResonantMode<T> tmp = c[i];
            c[i] = c[best];
            c[best] = tmp;
        }
    }
}

/// Half-power (peaks) or double-power (valleys) bandwidth → ζ.
template<std::size_t NFreq, typename T>
[[nodiscard]] constexpr ResonantMode<T> finish_extremum(
    const damp::array<FrfPoint<T>, NFreq>& table,
    std::size_t                            i,
    FrfFeatureKind                         kind,
    const ModeExtractorConfig<T>&          cfg
) {
    ResonantMode<T> mode{};
    const auto&     p = table[i];
    mode.kind = kind;
    mode.frequency_hz = p.freq_hz;
    mode.omega_n = T{2} * damp::numbers::pi_v<T> * p.freq_hz;
    mode.peak_gain = p.magnitude;

    // Peaks: −3 dB level = mag/√2. Valleys: +3 dB level = mag·√2.
    const T level = (kind == FrfFeatureKind::Peak)
                      ? (p.magnitude * (T{1} / damp::numbers::sqrt2_v<T>))
                      : (p.magnitude * damp::numbers::sqrt2_v<T>);

    bool found_lo = false;
    T    f_lo = p.freq_hz;
    for (std::size_t j = i; j > 0; --j) {
        const auto& a = table[j - 1];
        const auto& b = table[j];
        if (!a.valid || !b.valid) {
            break;
        }
        const bool cross = (kind == FrfFeatureKind::Peak)
                             ? (a.magnitude <= level && b.magnitude >= level)
                             : (a.magnitude >= level && b.magnitude <= level);
        if (cross) {
            const T denom = b.magnitude - a.magnitude;
            const T t = (damp::abs(denom) > default_tol<T>()) ? ((level - a.magnitude) / denom) : T{0};
            f_lo = a.freq_hz + t * (b.freq_hz - a.freq_hz);
            found_lo = true;
            break;
        }
    }

    bool found_hi = false;
    T    f_hi = p.freq_hz;
    for (std::size_t j = i; j + 1 < NFreq; ++j) {
        const auto& a = table[j];
        const auto& b = table[j + 1];
        if (!a.valid || !b.valid) {
            break;
        }
        const bool cross = (kind == FrfFeatureKind::Peak)
                             ? (a.magnitude >= level && b.magnitude <= level)
                             : (a.magnitude <= level && b.magnitude >= level);
        if (cross) {
            const T denom = b.magnitude - a.magnitude;
            const T t = (damp::abs(denom) > default_tol<T>()) ? ((level - a.magnitude) / denom) : T{0};
            f_hi = a.freq_hz + t * (b.freq_hz - a.freq_hz);
            found_hi = true;
            break;
        }
    }

    if (found_lo && found_hi && f_hi > f_lo) {
        mode.zeta = (f_hi - f_lo) / (T{2} * p.freq_hz);
        if (mode.zeta < T{0}) {
            mode.zeta = T{0};
        }
        if (mode.zeta > T{1}) {
            mode.zeta = T{1};
        }
        mode.valid = true;
    } else if (!cfg.require_valid_zeta) {
        mode.zeta = T{0};
        mode.valid = true;
    }
    return mode;
}

template<std::size_t NFreq, std::size_t MaxModes, typename T>
[[nodiscard]] constexpr ModeExtractorResult<MaxModes, T> extract_extrema(
    const damp::array<FrfPoint<T>, NFreq>& table,
    const ModeExtractorConfig<T>&          cfg,
    FrfFeatureKind                         kind
) {
    ModeExtractorResult<MaxModes, T> result{};
    if constexpr (NFreq < 3 || MaxModes == 0) {
        return result;
    }

    damp::array<ResonantMode<T>, NFreq> candidates{};
    std::size_t                         n_cand = 0;

    for (std::size_t i = 1; i + 1 < NFreq; ++i) {
        if (!point_ok(table, i, cfg) || !point_ok(table, i - 1, cfg) || !point_ok(table, i + 1, cfg)) {
            continue;
        }
        const auto& p = table[i];
        const auto& pl = table[i - 1];
        const auto& pr = table[i + 1];

        if (kind == FrfFeatureKind::Peak) {
            if (p.magnitude < cfg.min_peak_gain) {
                continue;
            }
            if (!(p.magnitude > pl.magnitude && p.magnitude >= pr.magnitude)) {
                continue;
            }
            if ((p.magnitude - pl.magnitude) < cfg.min_prominence
                || (p.magnitude - pr.magnitude) < cfg.min_prominence) {
                continue;
            }
        } else {
            if (cfg.max_valley_gain > T{0} && p.magnitude > cfg.max_valley_gain) {
                continue;
            }
            if (!(p.magnitude < pl.magnitude && p.magnitude <= pr.magnitude)) {
                continue;
            }
            if ((pl.magnitude - p.magnitude) < cfg.min_prominence
                || (pr.magnitude - p.magnitude) < cfg.min_prominence) {
                continue;
            }
        }

        const auto mode = finish_extremum(table, i, kind, cfg);
        if (!mode.valid) {
            continue;
        }
        if (n_cand < NFreq) {
            candidates[n_cand++] = mode;
        }
    }

    if (kind == FrfFeatureKind::Peak) {
        sort_by_gain_desc(candidates, n_cand);
    } else {
        sort_by_gain_asc(candidates, n_cand); // deepest valleys first
    }

    const std::size_t take = (n_cand < MaxModes) ? n_cand : MaxModes;
    for (std::size_t i = 0; i < take; ++i) {
        result.modes[i] = candidates[i];
    }
    result.count = take;
    result.success = take > 0;
    return result;
}

} // namespace detail_frf

/**
 * @brief Extract resonant peaks (local |G| maxima) from an FRF table.
 *
 * Half-power bandwidth → ζ ≈ Δf / (2 fₙ). Strongest peaks first.
 *
 * @see "Vibration Testing: Theory and Practice" (McConnell, 1995)
 */
template<std::size_t NFreq, std::size_t MaxModes, typename T = float>
[[nodiscard]] constexpr ModeExtractorResult<MaxModes, T> extract_modes(
    const damp::array<FrfPoint<T>, NFreq>& table,
    const ModeExtractorConfig<T>&          cfg = {}
) {
    return detail_frf::extract_extrema<NFreq, MaxModes, T>(table, cfg, FrfFeatureKind::Peak);
}

/**
 * @brief Extract anti-resonance valleys (local |G| minima) from an FRF table.
 *
 * Double-power (+3 dB) bandwidth → ζ. Deepest valleys first.
 * Feeds band-pass / peaking compensators and two-mass zero placement.
 */
template<std::size_t NFreq, std::size_t MaxModes, typename T = float>
[[nodiscard]] constexpr ModeExtractorResult<MaxModes, T> extract_valleys(
    const damp::array<FrfPoint<T>, NFreq>& table,
    const ModeExtractorConfig<T>&          cfg = {}
) {
    return detail_frf::extract_extrema<NFreq, MaxModes, T>(table, cfg, FrfFeatureKind::Valley);
}

/**
 * @brief Extract peaks and valleys in one pass over the FRF table.
 */
template<std::size_t NFreq, std::size_t MaxPeaks, std::size_t MaxValleys, typename T = float>
[[nodiscard]] constexpr FrfFeatureResult<MaxPeaks, MaxValleys, T> extract_frf_features(
    const damp::array<FrfPoint<T>, NFreq>& table,
    const ModeExtractorConfig<T>&          cfg = {}
) {
    FrfFeatureResult<MaxPeaks, MaxValleys, T> out{};
    const auto                                peaks = extract_modes<NFreq, MaxPeaks, T>(table, cfg);
    const auto                                valleys = extract_valleys<NFreq, MaxValleys, T>(table, cfg);
    out.peaks = peaks.modes;
    out.valleys = valleys.modes;
    out.peak_count = peaks.count;
    out.valley_count = valleys.count;
    out.success = peaks.success || valleys.success;
    return out;
}

/**
 * @brief Gain / phase margins from a fixed on-target FRF table.
 *
 * Treats the table as open-loop L(jω) (or plant G if that is what was measured).
 * - Gain crossover: |G| crosses 1 (descending). PM = 180° + ∠G at that f.
 * - Phase crossover: ∠G crosses −π. GM_dB = −20 log₁₀(|G|) at that f.
 *
 * Phase is used as stored (rad); unwrap is not applied — prefer stepped-sine
 * tables that do not jump ±π between bins.
 *
 * @tparam NFreq Table length
 * @tparam T     Scalar type
 * @param table  FRF points
 * @return FrfMargins
 */
template<std::size_t NFreq, typename T = float>
[[nodiscard]] constexpr FrfMargins<T> margins_from_frf(const damp::array<FrfPoint<T>, NFreq>& table) {
    FrfMargins<T> m{};
    if constexpr (NFreq < 2) {
        return m;
    }

    const T pi = damp::numbers::pi_v<T>;

    for (std::size_t i = 1; i < NFreq; ++i) {
        const auto& a = table[i - 1];
        const auto& b = table[i];
        if (!a.valid || !b.valid) {
            continue;
        }

        // Gain crossover: |G| from ≥1 to <1
        if (!m.has_phase_margin && a.magnitude >= T{1} && b.magnitude < T{1}) {
            const T denom = a.magnitude - b.magnitude;
            const T t = (damp::abs(denom) > default_tol<T>()) ? ((a.magnitude - T{1}) / denom) : T{0};
            const T f = a.freq_hz + t * (b.freq_hz - a.freq_hz);
            const T ph = a.phase_rad + t * (b.phase_rad - a.phase_rad);
            m.gain_crossover_hz = f;
            m.phase_margin_deg = (T{180} / pi) * (pi + ph);
            m.has_phase_margin = true;
        }

        // Phase crossover: phase from > −π to ≤ −π (or reverse)
        if (!m.has_gain_margin) {
            const bool cross = (a.phase_rad > -pi && b.phase_rad <= -pi) || (a.phase_rad < -pi && b.phase_rad >= -pi);
            if (cross) {
                const T denom = b.phase_rad - a.phase_rad;
                const T t = (damp::abs(denom) > default_tol<T>()) ? ((-pi - a.phase_rad) / denom) : T{0};
                const T f = a.freq_hz + t * (b.freq_hz - a.freq_hz);
                const T mag = a.magnitude + t * (b.magnitude - a.magnitude);
                m.phase_crossover_hz = f;
                if (mag > default_tol<T>()) {
                    m.gain_margin_db = -T{20} * damp::log10(mag);
                }
                m.has_gain_margin = true;
            }
        }
    }

    m.success = m.has_phase_margin || m.has_gain_margin;
    return m;
}

} // namespace design

} // namespace damp
