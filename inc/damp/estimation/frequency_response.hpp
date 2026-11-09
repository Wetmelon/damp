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
 * - FrequencyResponseEstimator — lock-in demodulation of (u, y) against a
 *   unit phasor at each bin frequency, H1 estimate G = S_yu / S_uu,
 *   coherence γ². Sequential path for SteppedSine; parallel bank for MultiSine.
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
 * // For each frequency f: begin_bin (true) → push(u,y) until true → table filled.
 * // Then: auto modes = design::extract_modes<Nf, 2>(frf.table());
 * @endcode
 *
 * @see excitation.hpp for SteppedSine / MultiSine stimuli
 * @see "System Identification" (Ljung, 1999), §6 (frequency-domain ID)
 * @see "Vibration Testing" (McConnell, 1995) — half-power bandwidth damping
 */


#include "damp/backend.hpp"
#include "damp/math/math.hpp"

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
    T    coherence{T{0}}; ///< γ² ∈ [0, 1] (H1 single-input coherence); 0 = not measured (fewer than two blocks)
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

namespace detail {

/**
 * @brief Lock-in sums for one frequency: U = Σ u·e^{−jωn}, Y = Σ y·e^{−jωn} over a block.
 *
 * The phasor e^{−jωn} is rotated by e^{−jω} each sample and pulled back to unit length, so
 * its angle error comes from the rounding of sin ω alone. A Goertzel recursion instead
 * differences terms scaled by 2cos ω ≈ 2, which loses tens of percent in float at small ω.
 */
template<typename T>
struct FrfLockIn {
    T      rot_re{T{1}}; ///< cos ω
    T      rot_im{T{0}}; ///< −sin ω
    T      p_re{T{1}};   ///< Phasor e^{−jωn}, real part
    T      p_im{T{0}};   ///< Phasor e^{−jωn}, imaginary part
    T      u_re{T{0}};   ///< Re U
    T      u_im{T{0}};   ///< Im U
    T      y_re{T{0}};   ///< Re Y
    T      y_im{T{0}};   ///< Im Y
    size_t block{0};     ///< Samples per block (0 = inert)
    size_t count{0};     ///< Samples in the current block

    constexpr FrfLockIn() = default;

    /** @brief Bin at ω = 2πf/fs [rad/sample] over blocks of @p block_size samples. */
    constexpr FrfLockIn(T omega, size_t block_size) : block(block_size) {
        const auto [s, c] = damp::sincos(omega);
        rot_re = c;
        rot_im = -s;
    }

    /** @brief Feed one (u, y) pair; true when the block just completed (sums valid until the next push). */
    constexpr bool push(T u, T y) {
        if (block == 0) {
            return false;
        }
        if (count == block) {
            p_re = T{1};
            p_im = T{0};
            u_re = T{0};
            u_im = T{0};
            y_re = T{0};
            y_im = T{0};
            count = 0;
        }
        u_re += u * p_re;
        u_im += u * p_im;
        y_re += y * p_re;
        y_im += y * p_im;
        const T re = (p_re * rot_re) - (p_im * rot_im);
        const T im = (p_re * rot_im) + (p_im * rot_re);
        const T g = T{1.5} - (T{0.5} * ((re * re) + (im * im))); // one Newton step toward |p| = 1
        p_re = re * g;
        p_im = im * g;
        ++count;
        return count == block;
    }
};

} // namespace detail

/**
 * @brief Lock-in FRF estimator over a fixed table of NFreq bins.
 *
 * Accumulates the classical H1 estimate from single-bin DFTs of the applied
 * input u and measured output y over blocks of samples (both multiplied by the
 * same unit phasor e^{−jωn}, rotated each sample, and summed):
 *
 *     S_yu = Σ Y · conj(U),   S_uu = Σ |U|²,   S_yy = Σ |Y|²
 *     G    = S_yu / S_uu
 *     γ²   = |S_yu|² / (S_uu · S_yy)
 *
 * γ² needs two or more blocks: over one block |Y·conj(U)|² = |U|²·|Y|² for any
 * data, so a single-block point reports coherence 0 (not measured).
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
template<size_t NFreq, typename T = float>
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
     * @brief Coherent block length for an integer number of cycles.
     *
     * N = round(cycles · fs / f). Prefer frequencies that make this exact.
     *
     * @param freq_hz Target frequency [Hz]
     * @param fs_hz   Sample rate [Hz]
     * @param cycles  Periods in the block (>= 1)
     * @return Block length in samples (>= 1)
     */
    [[nodiscard]] static constexpr size_t coherent_block_size(T freq_hz, T fs_hz, size_t cycles) {
        if (!(freq_hz > T{0}) || !(fs_hz > T{0}) || cycles == 0) {
            return 1;
        }
        const T n = static_cast<T>(cycles) * fs_hz / freq_hz;
        auto    N = static_cast<size_t>(n + static_cast<T>(0.5));
        if (N < 1) {
            N = 1;
        }
        return N;
    }

    /**
     * @brief Coherent block length using the estimator's sample rate.
     */
    [[nodiscard]] constexpr size_t coherent_block_size(T freq_hz, size_t cycles) const {
        return coherent_block_size(freq_hz, fs_, cycles);
    }

    /**
     * @brief Start sequential measurement of table bin @p index.
     *
     * Averages @p num_blocks successive blocks of @p block_size samples (after
     * discarding @p settle_samples). Call @ref push until it returns true. One
     * block gives no coherence (0); two or more measure it.
     *
     * @param index          Table index in [0, NFreq)
     * @param freq_hz        Bin frequency [Hz], in (0, fs/2)
     * @param block_size     Samples per block
     * @param num_blocks     Blocks to average for H1 (>= 1; >= 2 for coherence)
     * @param settle_samples Leading samples discarded (transients)
     * @return false if refused (bad index, frequency at or above Nyquist, sample rate unset,
     *         zero block size or count); no measurement is then active.
     */
    [[nodiscard]] constexpr bool begin_bin(
        size_t index,
        T      freq_hz,
        size_t block_size,
        size_t num_blocks = 1,
        size_t settle_samples = 0
    ) {
        if (!bin_ok(index, freq_hz, block_size) || num_blocks == 0) {
            sequential_active_ = false;
            return false;
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
        lock_ = detail::FrfLockIn<T>(omega_of(freq_hz), block_size);
        sequential_active_ = true;
        table_[index] = FrfPoint<T>{};
        table_[index].freq_hz = freq_hz;
        return true;
    }

    /**
     * @brief Configure a parallel-bank bin (MultiSine path).
     *
     * All bank bins should use the same @p block_size so @ref push_bank
     * completes them together. Single-block H1 (no multi-block average), so
     * bank points report no coherence (0).
     *
     * @param index      Table index in [0, NFreq)
     * @param freq_hz    Bin frequency [Hz], in (0, fs/2)
     * @param block_size Samples per block
     * @return false if refused (bad index, frequency at or above Nyquist, sample rate unset,
     *         zero block size); the bin is then left as it was.
     */
    [[nodiscard]] constexpr bool configure_bin(size_t index, T freq_hz, size_t block_size) {
        if (!bin_ok(index, freq_hz, block_size)) {
            return false;
        }
        bank_[index] = detail::FrfLockIn<T>(omega_of(freq_hz), block_size);
        bank_active_[index] = true;
        table_[index] = FrfPoint<T>{};
        table_[index].freq_hz = freq_hz;
        return true;
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
        if (!lock_.push(u, y)) {
            return false;
        }
        accumulate_block(lock_);
        ++blocks_done_;
        if (blocks_done_ < blocks_target_) {
            return false;
        }
        finalize_point(active_index_, active_freq_, blocks_done_);
        sequential_active_ = false;
        return true;
    }

    /**
     * @brief Parallel path: feed one sample into every configured bank bin.
     * @return true when at least one bin completed a block on this sample (its table entry
     *         was updated); bins with equal block sizes complete together.
     */
    [[nodiscard]] constexpr bool push_bank(T u, T y) {
        bool any_active = false;
        bool completed = false;
        for (size_t i = 0; i < NFreq; ++i) {
            if (!bank_active_[i]) {
                continue;
            }
            any_active = true;
            if (bank_[i].push(u, y)) {
                // Single-block H1 for this bin.
                syu_re_ = T{0};
                syu_im_ = T{0};
                suu_ = T{0};
                syy_ = T{0};
                accumulate_block(bank_[i]);
                finalize_point(i, table_[i].freq_hz, 1);
                completed = true;
            }
        }
        return any_active && completed;
    }

    /** @brief FRF table (fixed size; check @c FrfPoint::valid per entry). */
    [[nodiscard]] constexpr const damp::array<FrfPoint<T>, NFreq>& table() const { return table_; }

    /** @brief Mutable table (e.g. host post-processing). */
    [[nodiscard]] constexpr damp::array<FrfPoint<T>, NFreq>& table() { return table_; }

    /** @brief One table entry. */
    [[nodiscard]] constexpr const FrfPoint<T>& point(size_t i) const { return table_[i]; }

    /** @brief Number of table slots (template size). */
    [[nodiscard]] static constexpr size_t size() { return NFreq; }

    /** @brief Count of currently valid points. */
    [[nodiscard]] constexpr size_t valid_count() const {
        size_t n = 0;
        for (size_t i = 0; i < NFreq; ++i) {
            if (table_[i].valid) {
                ++n;
            }
        }
        return n;
    }

    /** @brief Sample rate [Hz] (0 if unset). */
    [[nodiscard]] constexpr T sample_rate() const { return fs_; }

    /** @brief Set sample rate [Hz] (must be > 0 before begin_bin / configure_bin). */
    constexpr void set_sample_rate(T fs_hz) { fs_ = fs_hz; }

    /** @brief Clear the table and abort any in-progress measurement. */
    constexpr void reset() {
        sequential_active_ = false;
        blocks_done_ = 0;
        settle_left_ = 0;
        for (size_t i = 0; i < NFreq; ++i) {
            table_[i] = FrfPoint<T>{};
            bank_active_[i] = false;
            bank_[i] = detail::FrfLockIn<T>{};
        }
        lock_ = detail::FrfLockIn<T>{};
    }

private:
    [[nodiscard]] constexpr bool bin_ok(size_t index, T freq_hz, size_t block_size) const {
        return index < NFreq && freq_hz > T{0} && freq_hz < fs_ / T{2} && block_size > 0;
    }

    [[nodiscard]] constexpr T omega_of(T freq_hz) const { return (T{2} * damp::numbers::pi_v<T> * freq_hz) / fs_; }

    constexpr void accumulate_block(const detail::FrfLockIn<T>& b) {
        // S_yu += Y * conj(U) = (yr + j yi)(ur - j ui)
        syu_re_ += (b.y_re * b.u_re) + (b.y_im * b.u_im);
        syu_im_ += (b.y_im * b.u_re) - (b.y_re * b.u_im);
        suu_ += (b.u_re * b.u_re) + (b.u_im * b.u_im);
        syy_ += (b.y_re * b.y_re) + (b.y_im * b.y_im);
    }

    constexpr void finalize_point(size_t index, T freq_hz, size_t blocks) {
        FrfPoint<T> p{};
        p.freq_hz = freq_hz;
        // Exact-zero guards only: an absolute floor would depend on the units of u, y.
        if (!(suu_ > T{0})) {
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
        if (blocks >= 2 && den > T{0}) { // one block: γ² ≡ 1 by Cauchy–Schwarz, not a measurement
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
    detail::FrfLockIn<T> lock_{};
    T                    active_freq_{T{0}};
    T                    syu_re_{T{0}};
    T                    syu_im_{T{0}};
    T                    suu_{T{0}};
    T                    syy_{T{0}};
    size_t               active_index_{0};
    size_t               blocks_target_{1};
    size_t               blocks_done_{0};
    size_t               settle_left_{0};
    bool                 sequential_active_{false};

    // Parallel bank
    damp::array<detail::FrfLockIn<T>, NFreq> bank_{};
    damp::array<bool, NFreq>                 bank_active_{};
};

/**
 * @brief Peak (resonance) or valley (anti-resonance) on an FRF magnitude curve.
 */
enum class FrfFeatureKind : uint8_t {
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
    T              zeta{T{0}};                 ///< Damping ratio (half-power / double-power BW); 0 = not estimated
    T              peak_gain{T{0}};            ///< |G| at the extremum
    bool           valid{false};               ///< Extremum found, with ζ estimated (or with ζ = 0 when require_valid_zeta is false)

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
    bool require_valid_zeta{true}; ///< Drop extrema without a half-power BW; if false they are kept, valid, with ζ = 0 (not estimated)
};

/**
 * @brief Result of reducing an FRF table to at most MaxModes peaks (resonances).
 *
 * Sorted by descending |G|. Feeds input shapers and notch compensators.
 *
 * @tparam MaxModes Capacity of the mode table.
 * @tparam T        Scalar type.
 */
template<size_t MaxModes, typename T = float>
struct [[nodiscard]] ModeExtractorResult {
    damp::array<ResonantMode<T>, MaxModes> modes{};
    size_t                                 count{0};
    bool                                   success{false};

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr ModeExtractorResult<MaxModes, U> as() const {
        ModeExtractorResult<MaxModes, U> out{};
        out.count = count;
        out.success = success;
        for (size_t i = 0; i < MaxModes; ++i) {
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
template<size_t MaxPeaks, size_t MaxValleys, typename T = float>
struct [[nodiscard]] FrfFeatureResult {
    damp::array<ResonantMode<T>, MaxPeaks>   peaks{};
    damp::array<ResonantMode<T>, MaxValleys> valleys{};
    size_t                                   peak_count{0};
    size_t                                   valley_count{0};
    bool                                     success{false};

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr FrfFeatureResult<MaxPeaks, MaxValleys, U> as() const {
        FrfFeatureResult<MaxPeaks, MaxValleys, U> out{};
        out.peak_count = peak_count;
        out.valley_count = valley_count;
        out.success = success;
        for (size_t i = 0; i < MaxPeaks; ++i) {
            out.peaks[i] = peaks[i].template as<U>();
        }
        for (size_t i = 0; i < MaxValleys; ++i) {
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
    T    phase_margin_deg{T{0}};   ///< PM at gain crossover, in (−180°, 180°] (negative: closed loop unstable)
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

template<size_t NFreq, typename T>
[[nodiscard]] constexpr bool point_ok(
    const damp::array<FrfPoint<T>, NFreq>& table, size_t i, const ModeExtractorConfig<T>& cfg
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

template<size_t N, typename T>
constexpr void sort_by_gain_desc(damp::array<ResonantMode<T>, N>& c, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        size_t best = i;
        for (size_t j = i + 1; j < n; ++j) {
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

template<size_t N, typename T>
constexpr void sort_by_gain_asc(damp::array<ResonantMode<T>, N>& c, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        size_t best = i;
        for (size_t j = i + 1; j < n; ++j) {
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

/** @brief Half-power (peaks) or double-power (valleys) bandwidth → ζ. */
template<size_t NFreq, typename T>
[[nodiscard]] constexpr ResonantMode<T> finish_extremum(
    const damp::array<FrfPoint<T>, NFreq>& table,
    size_t                                 i,
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
    for (size_t j = i; j > 0; --j) {
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
            const T t = (denom != T{0}) ? ((level - a.magnitude) / denom) : T{0};
            f_lo = a.freq_hz + t * (b.freq_hz - a.freq_hz);
            found_lo = true;
            break;
        }
    }

    bool found_hi = false;
    T    f_hi = p.freq_hz;
    for (size_t j = i; j + 1 < NFreq; ++j) {
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
            const T t = (denom != T{0}) ? ((level - a.magnitude) / denom) : T{0};
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

template<size_t NFreq, size_t MaxModes, typename T>
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
    size_t                              n_cand = 0;

    for (size_t i = 1; i + 1 < NFreq; ++i) {
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

    const size_t take = (n_cand < MaxModes) ? n_cand : MaxModes;
    for (size_t i = 0; i < take; ++i) {
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
template<size_t NFreq, size_t MaxModes, typename T = float>
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
template<size_t NFreq, size_t MaxModes, typename T = float>
[[nodiscard]] constexpr ModeExtractorResult<MaxModes, T> extract_valleys(
    const damp::array<FrfPoint<T>, NFreq>& table,
    const ModeExtractorConfig<T>&          cfg = {}
) {
    return detail_frf::extract_extrema<NFreq, MaxModes, T>(table, cfg, FrfFeatureKind::Valley);
}

/**
 * @brief Extract peaks and valleys in one pass over the FRF table.
 */
template<size_t NFreq, size_t MaxPeaks, size_t MaxValleys, typename T = float>
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
 * - Gain crossover: |G| crosses 1 (descending). PM = 180° + ∠G at that f, wrapped to (−180°, 180°].
 * - Phase crossover: ∠G crosses −π. GM_dB = −20 log₁₀(|G|) at that f.
 *
 * Stored phases may be wrapped to (−π, π]: consecutive valid bins are unwrapped
 * on the fly (each step taken as the shortest arc), so bins must be dense
 * enough that the true phase moves less than π between neighbours.
 *
 * @tparam NFreq Table length
 * @tparam T     Scalar type
 * @param table  FRF points
 * @return FrfMargins
 */
template<size_t NFreq, typename T = float>
[[nodiscard]] constexpr FrfMargins<T> margins_from_frf(const damp::array<FrfPoint<T>, NFreq>& table) {
    FrfMargins<T> m{};
    if constexpr (NFreq < 2) {
        return m;
    }

    const T pi = damp::numbers::pi_v<T>;
    // Linear interpolation parameter of `level` between a and b (a ≠ b by the crossing test).
    const auto lerp_t = [](T a, T b, T level) { return (b != a) ? ((level - a) / (b - a)) : T{0}; };

    // Unwrapped phase of the previous valid bin; atan2-wrapped measurements jump
    // ±2π near −180°, exactly where the phase crossover sits.
    T ph_a = T{0};
    for (size_t i = 1; i < NFreq; ++i) {
        const auto& a = table[i - 1];
        const auto& b = table[i];
        if (!a.valid || !b.valid) {
            continue;
        }
        if (i == 1 || !table[i - 2].valid) {
            ph_a = a.phase_rad;
        }
        const T ph_b = ph_a + damp::wrap_pi(b.phase_rad - a.phase_rad);

        // Gain crossover: |G| from ≥1 to <1
        if (!m.has_phase_margin && a.magnitude >= T{1} && b.magnitude < T{1}) {
            const T t = lerp_t(a.magnitude, b.magnitude, T{1});
            const T ph = damp::wrap_pi(ph_a + t * (ph_b - ph_a));
            m.gain_crossover_hz = a.freq_hz + t * (b.freq_hz - a.freq_hz);
            // 180° + ∠L folded into (−180°, 180°] (wrap_pi folds into [−π, π), hence the negations).
            m.phase_margin_deg = -(T{180} / pi) * damp::wrap_pi(-(pi + ph));
            m.has_phase_margin = true;
        }

        // Phase crossover: unwrapped phase crosses −π (modulo 2π) between a and b.
        if (!m.has_gain_margin) {
            // Shift both so a lies in (−π, π]: the step (< π) can then only
            // cross −π going down or +π (≡ −π) going up.
            const T    shift = ph_a - damp::wrap_pi(ph_a);
            const T    pa = ph_a - shift;
            const T    pb = ph_b - shift;
            const bool down = (pb <= -pi);
            const bool up = (pa < pi && pb >= pi);
            if (down || up) {
                const T t = lerp_t(pa, pb, down ? -pi : pi);
                const T mag = a.magnitude + t * (b.magnitude - a.magnitude);
                m.phase_crossover_hz = a.freq_hz + t * (b.freq_hz - a.freq_hz);
                if (mag > T{0}) {
                    m.gain_margin_db = -T{20} * damp::log10(mag);
                }
                m.has_gain_margin = true;
            }
        }
        ph_a = ph_b;
    }

    m.success = m.has_phase_margin || m.has_gain_margin;
    return m;
}

} // namespace design

} // namespace damp
