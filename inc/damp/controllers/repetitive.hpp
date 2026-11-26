// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file repetitive.hpp
 * @brief Plug-in repetitive controller — internal-model rejection of a periodic
 *        reference/disturbance and *all* its harmonics with one delay loop.
 *
 * Repetitive control embeds the internal model of a periodic signal: a positive-
 * feedback delay line of one fundamental period N = round(2π/(ω₀·Ts)). That single
 * model puts (near-)infinite loop gain at the fundamental *and every harmonic*
 * simultaneously — the reason one repetitive block out-rejects a stack of
 * individual resonators. Targets: grid-tied inverters (track/clean a sinusoidal
 * reference, reject periodic load), precision motion with periodic trajectories,
 * rotating-machinery ripple.
 *
 * It is a *plug-in*: it sits alongside any base controller and adds a correction.
 *
 *     e   = r - y;                 // tracking error
 *     u   = base.step(r, y) + rc.step(e);
 *
 * Standard discrete realization (robustness filter Q(z), integer phase lead m):
 *
 *     w[k]    = Q(z)·w[k-N] + e[k]          (internal-model memory, length-N buffer)
 *     u_rc[k] = k_rc · Q(z)·w[k - N + m]    (the learned, phase-advanced correction)
 *
 * Robustness filter Q. A scalar Q ∈ (0,1] uniformly rolls the model gain off
 * the high harmonics (where plant phase is uncertain), trading some steady-state
 * rejection for stability. The better choice is a zero-phase low-pass FIR
 * Q(z) = Σ_{i=−M}^{M} q_i z^{−i} with symmetric taps (q_i = q_{−i}): it keeps
 * near-unity gain on the low harmonics (full rejection) yet rolls off to ~0 near
 * Nyquist (robust stability), *without* adding phase. It is realizable inside the
 * loop because Q multiplies the N-delayed signal — the "future" taps z^{+i} read
 * already-buffered samples w[k−N+m+i] (available for m + M < N). The default family
 * is the unity-DC-gain binomial: M=1 → [1,2,1]/4, M=2 → [1,4,6,4,1]/16, …
 *
 * Q filters the output tap as well as the memory, so it multiplies the whole learning
 * loop and the stability condition is ‖Q(1 − k_rc zᵐ G_o)‖∞ < 1, with G_o the base
 * loop's map from u_rc to y: where Q rolls off, the plant-uncertainty term drops out.
 *
 * The integer phase-lead m pre-applies the correction a few samples early to
 * offset plant transport lag.
 *
 * @see S. Hara, Y. Yamamoto, T. Omata, M. Nakano, "Repetitive Control System: A
 *      New Type Servo System for Periodic Exogenous Signals," IEEE TAC 33(7),
 *      1988, https://doi.org/10.1109/9.1274
 */


#include "damp/backend.hpp" // damp::array
#include "damp/controllers/detail/repetitive_gains.hpp"
#include "damp/controllers/law.hpp" // IWYU pragma: export
#include "damp/design/detail/repetitive_rules.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"            // damp::finite_positive
#include "damp/matrix/matrix_traits.hpp" // scalar_type_t

namespace damp {

namespace design {

/**
 * @brief What design::repetitive takes: the fundamental of the periodic signal, the learning gain,
 *        the robustness filter's DC gain and the phase lead; Ts is passed beside it
 *
 * `design::repetitive({.w0 = 2.0 * pi * 50.0, .gain = 1.0, .q = 0.98, .lead = 1}, 1e-3)`: a 50 Hz
 * period at 1 kHz, N = 20 samples.
 */
template<typename T = double>
struct RepetitiveSpec {
    T      w0{};       ///< [rad/s] Fundamental of the periodic reference or disturbance (> 0): N = round(2π/(ω₀·Ts)) samples per period
    T      gain{T{1}}; ///< k_rc, the repetitive (learning) gain, in (0, 2] (plant-dependent)
    T      q{T{1}};    ///< The robustness filter's DC gain, in (0, 1]: the scalar Q, or the scale of a binomial filter's taps
    size_t lead{0};    ///< m, phase-lead samples to offset plant lag (m + the filter's half-width < N)
};

/**
 * @brief Repetitive-controller tuning + period (with optional zero-phase FIR Q).
 *
 * The robustness filter is Q(z) = q_filter (center tap q_0) plus @p q_half
 * symmetric side taps `q_side[i-1] = q_i = q_{-i}`. With `q_half = 0` this is the
 * classic scalar Q. `MaxQHalf` bounds the FIR half-width at compile time.
 *
 * @tparam T        Scalar type
 * @tparam MaxQHalf Maximum FIR Q half-width (0 = scalar Q only)
 */
template<typename T = double, size_t MaxQHalf = 0>
struct RepetitiveConfig {

    using scalar = scalar_type_t<T>;

    size_t                        period{0};           ///< N = round(2π/(ω₀·Ts)), samples per fundamental period
    scalar                        gain{scalar{1}};     ///< k_rc, repetitive (learning) gain — 0 < k_rc ≲ 2 (plant-dependent)
    scalar                        q_filter{scalar{1}}; ///< Q center tap q_0 (== scalar Q when q_half = 0)
    size_t                        lead{0};             ///< m, phase-lead samples to offset plant lag (m + q_half < N)
    size_t                        q_half{0};           ///< FIR Q half-width M (0 = scalar); symmetric side taps below
    damp::array<scalar, MaxQHalf> q_side{};            ///< q_1 … q_M (each side; q_side[i-1] = q_i)

    /** @brief DC gain of the Q-filter, q_0 + 2·Σ q_i (1.0 for a unity-DC robustness filter). */
    [[nodiscard]] constexpr scalar q_dc_gain() const {
        scalar dc = q_filter;
        for (size_t i = 0; i < q_half; ++i) {
            dc += scalar{2} * q_side[i];
        }
        return dc;
    }

    /** @brief Whether the tuning runs: N ≥ 2, 0 < k_rc ≤ 2, q₀ > 0, the output tap's window in the past, and a filter that does not amplify. */
    [[nodiscard]] constexpr bool valid() const {
        if (period < 2) {
            return false;
        }
        if (!damp::finite_positive(gain) || gain > scalar{2}) {
            return false;
        }
        if (!damp::finite_positive(q_filter)) {
            return false;
        }
        if (q_half > MaxQHalf || lead + q_half >= period) {
            return false; // the output tap's window w[k−N+m−M … k−N+m+M] must be in the past (m + M < N)
        }
        const scalar dc = q_dc_gain();
        // Allow the rounding of summing 2M+1 taps (each ≤ 1) to unity, no more.
        const scalar dc_slack = scalar{4} * static_cast<scalar>((2 * q_half) + 1) * eps<scalar>();
        if (!damp::finite_positive(dc) || dc > scalar{1} + dc_slack) {
            return false; // a robustness filter must not amplify (DC gain ≤ 1)
        }
        return true;
    }

    /** @brief Convert the configuration to a different scalar precision. */
    template<typename U>
    [[nodiscard]] constexpr RepetitiveConfig<damp::remove_const_t<U>, MaxQHalf> as() const {
        using out_scalar = scalar_type_t<damp::remove_const_t<U>>;
        RepetitiveConfig<damp::remove_const_t<U>, MaxQHalf> out{};
        out.period = period;
        out.gain = static_cast<out_scalar>(gain);
        out.q_filter = static_cast<out_scalar>(q_filter);
        out.lead = lead;
        out.q_half = q_half;
        for (size_t i = 0; i < MaxQHalf; ++i) {
            out.q_side[i] = static_cast<out_scalar>(q_side[i]);
        }
        return out;
    }
};

/**
 * @brief Design result for the repetitive controller.
 */
template<typename T = double, size_t MaxQHalf = 0>
struct [[nodiscard]] RepetitiveResult {
    static constexpr size_t q_half_max = MaxQHalf; ///< The robustness filter's largest half-width

    RepetitiveConfig<T, MaxQHalf> config{};                            ///< The period, the gain, the robustness filter and the lead
    bool                          success{false};                      ///< status == DesignStatus::ok
    DesignStatus                  status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same design in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr RepetitiveResult<damp::remove_const_t<U>, MaxQHalf> as() const {
        return RepetitiveResult<damp::remove_const_t<U>, MaxQHalf>{config.template as<U>(), success, status};
    }
};

namespace detail {

/** @brief The result for @p config at @p Ts: invalid_sample_time for a Ts that is not finite and positive, invalid_parameter for a tuning that does not run, else ok. */
template<typename T, size_t MaxQHalf>
[[nodiscard]] constexpr RepetitiveResult<T, MaxQHalf> repetitive_result(const RepetitiveConfig<T, MaxQHalf>& config, T Ts) {
    if (!damp::finite_positive(Ts)) {
        return refused<RepetitiveResult<T, MaxQHalf>>(DesignStatus::invalid_sample_time);
    }
    if (!config.valid()) {
        return refused<RepetitiveResult<T, MaxQHalf>>(DesignStatus::invalid_parameter);
    }
    RepetitiveResult<T, MaxQHalf> out{.config = config};
    set_status(out, DesignStatus::ok);
    return out;
}

} // namespace detail

/**
 * @brief Synthesize a repetitive controller with a scalar robustness filter Q = q.
 *
 * Computes the period N = round(2π/(ω₀·Ts)) and validates the tuning. Refuses an invalid spec
 * with a status (rather than asserting) so it composes in `constexpr` design code.
 *
 * @code
 * constexpr auto kRep = design::repetitive({.w0 = 2.0 * pi * 50.0}, 1e-4); // 50 Hz at 10 kHz: N = 200
 * static constinit Law<kRep, float> rc;                                      // a buffer of N samples
 * float u = base.step(r, y) + rc.step(r - y);                                // plug-in correction
 * @endcode
 *
 * MATLAB® equivalent (no single builtin — the internal model is a delay loop):
 * @code{.m}
 *   N = round(2*pi / (w0 * Ts));
 *   % w(k) = q*w(k-N) + e(k);   u_rc(k) = k_rc * q * w(k-N+m);
 * @endcode
 *
 * @param spec ω₀ > 0, k_rc in (0, 2], q in (0, 1], the lead m < N
 * @param Ts   [s] The sample period (finite and positive)
 * @return The tuning; refused with invalid_sample_time for a Ts that is not finite and positive,
 *         and with invalid_parameter for a spec out of its domain or a period below 2 samples
 */
template<typename T = double>
[[nodiscard]] constexpr RepetitiveResult<T, 0> repetitive(const RepetitiveSpec<T>& spec, T Ts) {
    const RepetitiveConfig<T, 0> config{
        .period = detail::repetitive_period(spec.w0, Ts),
        .gain = static_cast<scalar_type_t<T>>(spec.gain),
        .q_filter = static_cast<scalar_type_t<T>>(spec.q),
        .lead = spec.lead,
        .q_half = 0,
    };
    return detail::repetitive_result(config, Ts);
}

/**
 * @brief Synthesize a repetitive controller with a binomial zero-phase FIR Q.
 *
 * Builds the unity-DC-gain binomial robustness filter of half-width M, scaled by the spec's q:
 * q_i = q·C(2M, M+i)/4^M (M=1 → q·[1,2,1]/4, M=2 → q·[1,4,6,4,1]/16, …). This rolls the model
 * gain off the high harmonics for robust stability while preserving (near-)full rejection of the
 * low harmonics — the standard choice over a scalar Q.
 *
 * @tparam M   FIR Q half-width (the binomial order); MaxQHalf of the result
 * @param spec ω₀ > 0, k_rc in (0, 2], q in (0, 1] (1: the binomial filter itself), the lead m with m + M < N
 * @param Ts   [s] The sample period (finite and positive)
 * @return The tuning; refused as design::repetitive's
 */
template<size_t M, typename T = double>
[[nodiscard]] constexpr RepetitiveResult<T, M> repetitive_binomial(const RepetitiveSpec<T>& spec, T Ts) {
    using scalar = scalar_type_t<T>;
    const auto             q = static_cast<scalar>(spec.q);
    RepetitiveConfig<T, M> config{
        .period = detail::repetitive_period(spec.w0, Ts),
        .gain = static_cast<scalar>(spec.gain),
        .q_filter = q * static_cast<scalar>(detail::binomial_tap(M, 0)),
        .lead = spec.lead,
        .q_half = M,
    };
    for (size_t i = 1; i <= M; ++i) {
        config.q_side[i - 1] = q * static_cast<scalar>(detail::binomial_tap(M, i));
    }
    return detail::repetitive_result(config, Ts);
}

} // namespace design

/**
 * @ingroup controllers
 * @brief Plug-in repetitive controller runtime (fixed-size internal model).
 *
 * `MaxPeriod` bounds the internal-model buffer at compile time (allocation-free);
 * the active period N ≤ MaxPeriod is set from the design. `step(e)` returns
 * the repetitive correction to add to the base controller's command. `MaxQHalf`
 * bounds the optional zero-phase FIR Q half-width (0 = scalar Q). A failed design, or
 * one whose period + q_half exceeds MaxPeriod, does not compile under constant
 * evaluation; at run time it leaves the controller invalid (zero correction).
 * `Law<kRes, T>` of a `constexpr` design sizes the buffer to N + q_half itself.
 *
 * It is an add-on, clamped with the law it adds to: it has no limits of its own, and
 * back_calculate takes that law's command as the feedforward u_ff. While another source drives
 * the plant, track(u, e) takes the place of step(e).
 *
 * @tparam MaxPeriod Buffer capacity (largest supported N + Q half-width)
 * @tparam T         Scalar type (float or double)
 * @tparam MaxQHalf  Maximum FIR Q half-width (0 = scalar Q only)
 * @tparam Gains     Where the tuning lives: void, the default, for a member (detail::RepetitiveGainSet), or the
 *                   constants of a `constexpr` design's Law (detail::ConstantRepetitiveGains)
 */
template<size_t MaxPeriod, typename T = float, size_t MaxQHalf = 0, typename Gains = void>
class RepetitiveController {
    using GainSet = damp::conditional_t<damp::is_same_v<Gains, void>, detail::RepetitiveGainSet<design::RepetitiveConfig<T, MaxQHalf>>, Gains>; ///< The tuning's storage
    static constexpr bool kHeld = !GainSet::constant;                                                                                           ///< The tuning is a member

    /** @brief The tuning's storage for @p config: the member, or nothing on a constexpr design's Law. */
    [[nodiscard]] static constexpr GainSet stored(const design::RepetitiveConfig<T, MaxQHalf>& config) {
        if constexpr (kHeld) {
            return GainSet{config};
        } else {
            static_cast<void>(config);
            return GainSet{};
        }
    }

public:
    using scalar = scalar_type_t<T>;
    static_assert(MaxPeriod >= 2, "RepetitiveController needs MaxPeriod >= 2");

    constexpr RepetitiveController() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr RepetitiveController(const design::RepetitiveResult<U, MaxQHalf>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : RepetitiveController(result.template as<T>()) {}

    constexpr RepetitiveController(const design::RepetitiveResult<T, MaxQHalf>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : k_(stored(result.config)), valid_(result.success && result.config.period + result.config.q_half <= MaxPeriod) {
        if (!valid_ && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_parameter : result.status); // constant-evaluated (constinit / constexpr): a failed design, or one MaxPeriod cannot hold, does not compile
        }
    }

    /** @brief Convert across scalar precision, preserving the learned internal model (or off a Law's constant tuning, onto a member). */
    template<typename U, typename G>
        requires kHeld
    constexpr explicit RepetitiveController(const RepetitiveController<MaxPeriod, U, MaxQHalf, G>& other)
        : k_{other.config().template as<T>()},
          clamp_(other.clamp_),
          write_idx_(other.write_idx_),
          e_last_(static_cast<T>(other.e_last_)),
          u_cmd_(static_cast<T>(other.u_cmd_)),
          valid_(other.valid_) {
        for (size_t i = 0; i < MaxPeriod; ++i) {
            mem_[i] = static_cast<T>(other.mem_[i]);
        }
    }

    /**
     * @brief Advance one tick with the tracking error @p e; returns the
     *        repetitive correction u_rc to add to the base command.
     */
    [[nodiscard]] constexpr T step(T e) {
        if (!valid_) {
            return T{0};
        }
        const size_t n = k_.config.period;

        // Phase-advanced output through the same zero-phase FIR: u_rc = k_rc·Q·w[k-N+lead].
        const T u_rc = static_cast<T>(k_.config.gain) * filtered(n - k_.config.lead);

        //! w[k] = Q·w[k-N] + e[k], unless the command that read w[k-N] (m ticks ago, the lead) was
        //! cut off by a clamp and this error would push the next period's further past it. With no
        //! lead that command is this tick's, and back_calculate takes the error back: the slot here
        //! still holds the flag of the tick a ring ago.
        e_last_ = (k_.config.lead > 0 && winds_up(e, clamp_[(write_idx_ + MaxPeriod - k_.config.lead) % MaxPeriod])) ? T{0} : e;
        mem_[write_idx_] = filtered(n) + e_last_;
        clamp_[write_idx_] = 0;
        write_idx_ = (write_idx_ + 1) % MaxPeriod;
        u_cmd_ = u_rc;

        return u_rc;
    }

    /**
     * @brief Record the input the plant received, @p u_sat, when the summed command was
     *        clamped, and what was added to this correction before that clamp, @p u_ff (the
     *        base controller's command)
     *
     * The memory does not learn from a sample the actuator could not follow. The error at tick k
     * answers the command sent m ticks before it (the phase lead m is the plant's delay), which
     * read the memory sample one period back. When a clamp cut that command off, u_sat ≠ u_rc +
     * u_ff, and the error would push the next period's command further past the limit, the
     * sample is stored without it: w[k] = Q·w[k−N]. An error that pulls the command back from
     * the limit is learned as ever. With Q = 1 the memory integrates the error period over
     * period, so without this it grows for as long as the clamp lasts.
     *
     * Costs one byte per memory sample: which way each of the last commands was cut.
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) {
        if (!valid_) {
            return;
        }
        const T      asked = u_cmd_ + u_ff;
        const T      diff = u_sat - asked;
        const T      tol = T{4} * eps<T>() * (damp::abs(u_sat) + damp::abs(asked));
        const size_t last = (write_idx_ + MaxPeriod - 1) % MaxPeriod; // this tick's slot
        clamp_[last] = (damp::abs(diff) > tol) ? ((diff < T{0}) ? static_cast<signed char>(-1) : static_cast<signed char>(1)) : static_cast<signed char>(0);
        if (k_.config.lead == 0 && winds_up(e_last_, clamp_[last])) {
            mem_[last] -= e_last_; // no lead: this tick's own command is the one the error answers
            e_last_ = T{0};
        }
    }

    /**
     * @brief In place of step(@p e) while another source drives the plant: the internal model held
     *        and kept in phase, @p u recorded as the correction applied
     *
     * The ring advances one tick, each sample carried one period on through Q, w[k] = Q·w[k−N],
     * without learning @p e: the correction is not acting on the plant, so the loop is not closed
     * through it and learning the error would only wind the model up. The model keeps what it
     * learned and its phase in the period, so the first step after a track emits the correction
     * for that point of the period. @p u is recorded as the last command, so a back_calculate of it
     * finds no clamp.
     */
    constexpr void track(T u, T e) {
        static_cast<void>(e); // not learned: the add-on is not driving the plant
        if (!valid_) {
            return;
        }
        mem_[write_idx_] = filtered(k_.config.period);
        clamp_[write_idx_] = 0;
        write_idx_ = (write_idx_ + 1) % MaxPeriod;
        e_last_ = T{0};
        u_cmd_ = u;
    }

    /** @brief The correction the last step() returned, or the one the last track() recorded. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd_; }

    /** @brief Clear the learned internal model (buffer + write index); config preserved. */
    constexpr void reset() {
        mem_ = {};
        clamp_ = {};
        write_idx_ = 0;
        e_last_ = T{0};
        u_cmd_ = T{0};
    }

    /** @brief The active configuration (period, gain, Q taps, phase lead). */
    [[nodiscard]] constexpr const auto& config() const { return k_.config; }
    /** @brief True if the design is valid and fit the buffer (period + q_half ≤ MaxPeriod). */
    [[nodiscard]] constexpr bool valid() const { return valid_; }

private:
    template<size_t, typename, size_t, typename>
    friend class RepetitiveController; // cross-precision converting ctor

    /**
     * @brief Q·w[k−d] = q₀·w[k−d] + Σᵢ qᵢ·(w[k−d−i] + w[k−d+i]), the zero-phase FIR on the memory d samples back
     *
     * Needs d > q_half (the newest tap is in the past) and d + q_half ≤ MaxPeriod (the oldest is
     * still buffered), which valid() and the constructor's capacity check guarantee for d = N − m and d = N.
     */
    [[nodiscard]] constexpr T filtered(size_t d) const {
        T qw = static_cast<T>(k_.config.q_filter) * mem_[(write_idx_ + MaxPeriod - d) % MaxPeriod];
        for (size_t i = 1; i <= k_.config.q_half; ++i) {
            const T older = mem_[(write_idx_ + MaxPeriod - d - i) % MaxPeriod]; // w[k-d-i]
            const T newer = mem_[(write_idx_ + MaxPeriod - d + i) % MaxPeriod]; // w[k-d+i]
            qw += static_cast<T>(k_.config.q_side[i - 1]) * (older + newer);
        }
        return qw;
    }

    /** @brief Whether learning the error @p e would push a command cut off towards @p cut (−1 from above, +1 from below, 0 not cut) further past the limit. The learning gain is positive. */
    [[nodiscard]] static constexpr bool winds_up(T e, signed char cut) { return (cut < 0 && e > T{0}) || (cut > 0 && e < T{0}); }

    [[no_unique_address]] GainSet       k_{};     ///< The tuning: a member, or no storage on a constexpr design's Law
    damp::array<T, MaxPeriod>           mem_{};   ///< The internal model w, a ring of MaxPeriod samples
    damp::array<signed char, MaxPeriod> clamp_{}; ///< Per tick, which way a clamp cut the command off (back_calculate): −1, 0, +1
    size_t                              write_idx_{0};
    T                                   e_last_{}; ///< The error the last step() stored
    T                                   u_cmd_{};  ///< The correction the last step() returned
    bool                                valid_{false};
};

namespace detail {

/** @brief Whether @p R is a repetitive design result. */
template<typename R>
inline constexpr bool is_repetitive_result = false;

/** @brief design::repetitive's and repetitive_binomial's result. */
template<typename U, size_t M>
inline constexpr bool is_repetitive_result<design::RepetitiveResult<U, M>> = true;

/** @brief A repetitive design. */
template<const auto& Res>
concept repetitive_design = is_repetitive_result<design_t<Res>>;

/** @brief The buffer a `constexpr` repetitive design needs: its period plus its filter's half-width (2 at least, for a failed one, which its Law refuses). */
template<const auto& Res>
[[nodiscard]] consteval size_t repetitive_capacity() {
    if constexpr (constant_design<Res>) {
        return damp::max(size_t{2}, Res.config.period + Res.config.q_half);
    } else {
        static_assert(
            constant_design<Res>,
            "Law<kRes, T> of a repetitive design sizes its buffer from the design's period, so kRes is constexpr. A design made at start-up runs as RepetitiveController<MaxPeriod, T, M> rc{res}, MaxPeriod the largest period it may have."
        );
        return 2;
    }
}

} // namespace detail

/**
 * @brief Law<kRes, T> of a `constexpr` repetitive design: the RepetitiveController whose buffer is
 *        the design's period plus its filter's half-width, no larger, its tuning constants of the
 *        type (detail::ConstantRepetitiveGains)
 *
 * The period N, the gain k_rc, the taps of Q and the lead m are in flash and fold into the tick;
 * the internal model stays a member. A failed design does not compile wherever the Law is
 * declared. A law built at start-up is `Law<kRes, T>::Runtime rc{kRes.as<T>()};`.
 */
template<const auto& Res, typename T>
    requires detail::repetitive_design<Res>
struct Law<Res, T> : RepetitiveController<detail::repetitive_capacity<Res>(), T, detail::design_t<Res>::q_half_max, detail::ConstantRepetitiveGains<Res, T>> {
    using Runtime = RepetitiveController<detail::repetitive_capacity<Res>(), T, detail::design_t<Res>::q_half_max>; ///< The runtime with a member: the type of a member a redesign assigns
    using Base = RepetitiveController<detail::repetitive_capacity<Res>(), T, detail::design_t<Res>::q_half_max, detail::ConstantRepetitiveGains<Res, T>>;

    static_assert(detail::constant_design_checked<Res>());

    /** @brief The design's law: its tuning the type's, its internal model a member. */
    constexpr Law() : Base{Res} {}

    /** @brief There is nothing to replace: the tuning is in the type. */
    template<typename Other>
        requires(!damp::is_same_v<damp::remove_cvref_t<Other>, Law>)
    constexpr Law& operator=(const Other& other) {
        static_assert(sizeof(Other) == 0, "This Law is a constexpr repetitive design's: its tuning is constants of the type, so there is nothing to replace. A member a redesign assigns is Law<kRes, T>::Runtime rc{kRes.as<T>()};");
        static_cast<void>(other);
        return *this;
    }
};

} // namespace damp
