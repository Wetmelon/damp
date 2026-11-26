// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file smith_predictor.hpp
 * @brief Smith predictor dead-time compensator: design::smith_predictor(model, Ts, InputDelay<N>{θ},
 *        primary, opts...) → SmithPredictorResult → SmithPredictor, any SISO law as the primary
 *
 * Classic structure (Smith 1957): a primary controller @f$C@f$ designed for the
 * delay-free plant model @f$\hat{G}@f$ is wrapped with a parallel model path so the
 * feedback signal seen by @f$C@f$ is free of transport delay when the model matches.
 *
 * Plant: @f$ G_p(s) = G(s)\,e^{-s\theta} @f$. Internal model: the first-order
 * @f$ \hat{G}(s) = K/(\tau s + 1) @f$ under ZOH at @f$T_s@f$, twice: delay-free,
 * @f$\hat{y}^{+} = a\hat{y} + b\,u_k@f$, and with the input delay
 * @f$\theta = d\,T_s + \tau_f@f$ (@f$0 \le \tau_f < T_s@f$) the plant has, exact for a fraction
 * of a sample as discretize's InputDelay is,
 * @f$\hat{y}_d^{+} = a\hat{y}_d + b_{\mathrm{fresh}}u_{k-d} + b_{\mathrm{stale}}u_{k-d-1}@f$.
 * Per tick:
 *
 * @f[
 *   y_{\mathrm{pred}} = y + \hat{y} - \hat{y}_d,\quad
 *   u_k = C(r,\, y_{\mathrm{pred}})
 * @f]
 *
 * then both models step on the command. With @f$a = e^{-T_s/\tau}@f$,
 * @f$b = b_{\mathrm{fresh}} + b_{\mathrm{stale}} = K(1-a)@f$, so the two settle on the same
 * output, and @f$b_{\mathrm{stale}} = 0@f$ for a delay of whole samples.
 *
 * It is the dead-time compensator for a single loop whose delay spans many samples. The
 * state-space laws (LQR, LQI, LQG, LQGI, OFLQG, MPC, OFMPC) take an input delay of N samples in the
 * model `discretize(sys, Ts, InputDelay<N>{θ})` returns, one state per delay state. The classic Smith structure is
 * not internally stable for an open-loop unstable plant: the model path ŷ − ŷ_d runs open
 * loop, so the plant's unstable pole stays in the map from a load disturbance to the output.
 *
 * The primary @f$C@f$ is any SisoLaw (default PIController). The design packs a PID designed on
 * the delay-free model for the common process-control path; pass another primary with the
 * design's model when pairing Smith with ADRC, SMC, etc.
 *
 * Example: FOPDT process with a SIMC PI primary
 * @code
 * #include "damp/controllers/smith_predictor.hpp"
 *
 * using namespace damp;
 *
 * // Plant G(s) = 1/(5s+1) · e^{−2s}, sampled at 10 Hz: 20 delay states
 * constexpr auto sp = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 5.0}, 0.1, InputDelay<20>{2.0},
 *                                             design::SIMC{.tau_c = 2.0});
 * static constinit SmithPredictor ctrl{sp.as<float>()}; // SmithPredictor<20, float>; a failed design does not compile
 * // In loop: float u = ctrl.step(r, y);
 * @endcode
 *
 * Example: another primary (any SisoLaw) with the design's model
 * @code
 * SmithPredictor smith{ADRCController<1, float>{adrc_design}, sp.as<float>()}; // sp's PID is not used
 * @endcode
 *
 * @see O.J.M. Smith, "Closer control of loops with dead time," Chem. Eng. Prog., 1957
 * @see Åström & Hägglund, "Advanced PID Control" (2006), §7.5 (Smith predictor)
 * @see "Control of Dead-time Processes" (Normey-Rico & Camacho, 2007): internal stability of the
 *      Smith structure on unstable plants
 * @see design::SIMC for the delay-free PI tuning used by the FOPDT designer
 * @see SisoLaw, Cascade — same primary/composition shape
 */

#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/controllers/detail/law_detail.hpp"
#include "damp/controllers/detail/smith_model.hpp"
#include "damp/controllers/law.hpp" // IWYU pragma: export
#include "damp/controllers/pid.hpp"
#include "damp/design/pid_design.hpp" // IWYU pragma: export
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

// Arduino.h defines PI as a macro, which would rewrite PIDMode::PI: set aside for this header.
#pragma push_macro("PI")
#undef PI

namespace damp {
namespace design {

/**
 * @brief Smith predictor design result: the primary PID at Ts, the first-order model and the
 *        input delay it is sampled with
 *
 * Produced by @ref smith_predictor. The packed @ref PIDResult is the primary designed on the
 * delay-free model; construct @ref SmithPredictor from the result (a PIController primary), or
 * from another primary and the result, which then gives the model alone.
 *
 * @tparam T Scalar type (design default: double)
 * @tparam N Delay states: the longest delay the runtime holds, in samples (the InputDelay's)
 */
template<typename T = double, size_t N = 1>
struct [[nodiscard]] SmithPredictorResult {
    PIDResult<T> pid{}; ///< The primary, designed on the delay-free model at Ts

    T            a{T{0}};                             ///< Discrete FO model pole a = e^{−T_s/τ}
    T            b{T{0}};                             ///< The delay-free model's gain on u_k: K(1 − a) (ZOH)
    T            b_fresh{T{0}};                       ///< The delayed model's gain on u_{k−d}, held after the switch at τ_f
    T            b_stale{T{0}};                       ///< Its gain on u_{k−d−1}, held before the switch (0 for a delay of whole samples)
    size_t       delay_samples{0};                    ///< d, the delay's whole periods (≤ N)
    T            Ts{T{0}};                            ///< Sample period [s]
    T            K{T{0}};                             ///< Continuous plant static gain
    T            tau{T{0}};                           ///< Continuous plant time constant [s]
    T            delay{T{0}};                         ///< [s] The input delay θ = d·T_s + τ_f, the InputDelay's
    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /**
     * @brief Convert all scalar fields to another precision
     * @tparam U Target scalar type
     */
    template<typename U>
    [[nodiscard]] constexpr SmithPredictorResult<U, N> as() const {
        return SmithPredictorResult<U, N>{
            .pid = pid.template as<U>(),
            .a = static_cast<U>(a),
            .b = static_cast<U>(b),
            .b_fresh = static_cast<U>(b_fresh),
            .b_stale = static_cast<U>(b_stale),
            .delay_samples = delay_samples,
            .Ts = static_cast<U>(Ts),
            .K = static_cast<U>(K),
            .tau = static_cast<U>(tau),
            .delay = static_cast<U>(delay),
            .success = success,
            .status = status,
        };
    }

    /**
     * @brief The continuous primary C(s): the packed PID's gains
     *
     * The delay is not in this rational map (Smith hides θ from C). Failed
     * designs, and a result whose Ts is not finite and positive, return a zero PID.
     */
    [[nodiscard]] constexpr PIDGains<T> primary_pid() const { return (success && damp::finite_positive(Ts)) ? pid.gains : PIDGains<T>{}; }

    /** @brief Delay-free primary transfer function (the packed PID's). */
    [[nodiscard]] constexpr TransferFunction<3, 3, T> to_tf() const { return primary_pid().to_tf(); }

    /**
     * @brief Delay-free primary state-space (same map as PIDGains::to_ss)
     *
     * @return The realization, or nullopt when the primary has a derivative without a filter
     */
    [[nodiscard]] constexpr damp::optional<StateSpace<2, 1, 1, T>> to_ss() const { return primary_pid().to_ss(); }

    /** @brief Delay-free FO model Ĝ(s) = K / (τ s + 1). */
    [[nodiscard]] constexpr TransferFunction<1, 2, T> plant_tf() const {
        return TransferFunction<1, 2, T>{
            .num = {K},
            .den = {tau, T{1}},
        };
    }
};

namespace detail {

/** @brief Whether @p P is a PID in gains you have (Parallel, Standard), which design::pid takes in place of a model and a rule. */
template<typename P>
inline constexpr bool is_pid_gain_form = false;

/** @brief Parallel is. */
template<typename U>
inline constexpr bool is_pid_gain_form<Parallel<U>> = true;

/** @brief Standard is. */
template<typename U>
inline constexpr bool is_pid_gain_form<Standard<U>> = true;

/** @brief The primary @p pid packed with the model of @p model and @p seconds at @p Ts in N delay states; the first refusal's status. */
template<size_t N, typename T>
[[nodiscard]] constexpr SmithPredictorResult<T, N> smith_result(const PIDResult<T>& pid, const FOPDT<T>& model, T seconds, T Ts) {
    const auto m = smith_model(model.K, model.tau, seconds, Ts, N);
    if (m.status != DesignStatus::ok) {
        return refused<SmithPredictorResult<T, N>>(m.status);
    }
    if (model.L != T{0} || !damp::isfinite(model.L)) {
        return refused<SmithPredictorResult<T, N>>(DesignStatus::invalid_parameter); // the dead time is the InputDelay's
    }
    if (!pid.success) {
        return refused<SmithPredictorResult<T, N>>(pid.status);
    }
    SmithPredictorResult<T, N> out{
        .pid = pid,
        .a = m.a,
        .b = m.b,
        .b_fresh = m.b_fresh,
        .b_stale = m.b_stale,
        .delay_samples = m.whole,
        .Ts = Ts,
        .K = model.K,
        .tau = model.tau,
        .delay = seconds,
    };
    set_status(out, DesignStatus::ok);
    return out;
}

} // namespace detail

/**
 * @brief A Smith predictor for a first-order plant whose input arrives @p delay late, its primary
 *        a PID on the delay-free model: a rule, or the gains you have
 *
 * The model is the delay-free plant K/(τs + 1), `FOPDT{.K, .tau}` with L = 0: the dead time is the
 * InputDelay's, θ ≤ N·T_s, which sizes the runtime (N delay states, deduced by SmithPredictor) and
 * is exact for a fraction of a sample. The primary is `design::pid` on the delay-free model at
 * @p Ts: a rule (SIMC, AMIGO needs L > 0 and so refuses here) with the PID's options
 * (InputLimits, DerivativeFilter, SetpointWeights, AntiWindup), or the gains you have, Parallel
 * or Standard, with the same options.
 *
 * @code
 * constexpr auto sp = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 5.0}, 0.1, InputDelay<20>{2.0},
 *                                             design::SIMC{.tau_c = 2.0}, design::InputLimits{.lower = {0.0}, .upper = {10.0}});
 * constexpr auto sp_pi = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 4.0}, 0.1, InputDelay<10>{1.0},
 *                                                design::Parallel{.Kp = 2.0, .Ki = 0.4});
 * @endcode
 *
 * @note Compare with the textbook Smith structure: primary controller designed
 *       for @f$G(s)@f$ only; the delay is restored by the parallel model path.
 *
 * @param model   The delay-free plant K/(τs + 1): K finite and non-zero, τ finite and positive, L = 0
 * @param Ts      Sample period [s] (finite, > 0)
 * @param delay   The input delay θ in N delay states, 0 ≤ θ ≤ N·T_s
 * @param primary The PID's rule on the delay-free model, or its gains (Parallel, Standard)
 * @param opts    The PID's options
 * @return SmithPredictorResult with the PID design, the model's coefficients and the delay's
 *         whole samples; refused with invalid_sample_time for a bad Ts, invalid_parameter for a
 *         bad plant, an L other than 0 or a delay outside [0, N·T_s], and the PID design's status
 *         when it fails
 *
 * @see design::SIMC — delay-free SIMC tuning (L = 0)
 * @see O.J.M. Smith (1957); Åström & Hägglund (2006), §7.5
 */
template<typename T, size_t N, typename S, typename Primary, typename... Opts>
[[nodiscard]] constexpr SmithPredictorResult<T, N> smith_predictor(
    const FOPDT<T>&          model,
    damp::type_identity_t<T> Ts,
    const InputDelay<N, S>&  delay,
    const Primary&           primary,
    const Opts&... opts
) {
    static_assert(detail::channel_delays<S, 1>, "A Smith predictor has one input: InputDelay<N>{theta}, a scalar delay.");
    const T seconds = detail::channel_seconds<T>(delay.seconds, 0);
    if constexpr (detail::is_pid_gain_form<Primary>) {
        return detail::smith_result<N>(pid(primary, Ts, opts...), model, seconds, Ts);
    } else {
        return detail::smith_result<N>(pid(FOPDT<T>{.K = model.K, .L = T{0}, .tau = model.tau}, Ts, primary, opts...), model, seconds, Ts);
    }
}

/**
 * @brief A Smith predictor for a first-order plant whose input arrives @p delay late, its primary
 *        SIMC's PI on the delay-free model with τ_c = τ
 *
 * smith_predictor(model, Ts, delay, design::SIMC{.tau_c = model.tau}).
 */
template<typename T, size_t N, typename S>
[[nodiscard]] constexpr SmithPredictorResult<T, N> smith_predictor(const FOPDT<T>& model, damp::type_identity_t<T> Ts, const InputDelay<N, S>& delay) {
    return smith_predictor(model, Ts, delay, SIMC<T>{.tau_c = model.tau, .mode = PIDMode::PI});
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Smith predictor runtime: a SISO primary law + a first-order model run delay-free and delayed
 *
 * Per-tick dead-time compensator. Stores any SisoLaw primary (default PIController) and the
 * design's first-order ZOH model, run twice on the commands: delay-free, ŷ, and through the
 * input delay, ŷ_d, which reads the commands in flight from a ring of the last N + 1.
 *
 * | Member | What it is |
 * | --- | --- |
 * | step(r, y) | one tick: the primary on y + ŷ − ŷ_d, then both models on its command |
 * | track(u, r, y) | in place of step while another source drives the plant |
 * | back_calculate(u_sat, u_ff) | a clamp after the law: the models and the primary take the input applied |
 * | limits(), set_limits(lim) | the primary's: its command is the predictor's |
 * | last_command(), valid(), reset() | the last command; whether the design ran and the primary runs; clear the models and the primary |
 * | model_output(), delayed_model_output(), primary() | ŷ; ŷ_d; the primary law |
 *
 * Matches the Cascade composition style: the Smith residual feeds
 * @c primary.step(r, y_pred); the primary type is a template parameter.
 *
 * The limits are the primary's (limits() and set_limits() forward to it): the predictor's
 * command is the primary's, clamped where the primary clamps, and the models run on that
 * clamped command. While another source drives the plant, track(u, r, y) takes the place of
 * step(r, y): both models step on u, then the primary tracks u at the prediction the next
 * step forms from y, so step(r, y) at this r and y returns what the primary's own track lands on
 * (u clamped, for a PID primary). u is recorded as the command applied.
 *
 * @tparam N        Delay states: the longest delay it holds, in samples (the design's InputDelay<N>)
 * @tparam T        Scalar type (embedded default: float)
 * @tparam Primary  Primary controller (SisoLaw). Default:
 *                  PIController for the FOPDT/SIMC process path. Use
 *                  PIDController (full PID), ADRC, SMC, etc. when needed.
 *
 * @see design::smith_predictor
 * @see SisoLaw
 * @see O.J.M. Smith (1957); Åström & Hägglund (2006), §7.5
 */
template<size_t N, typename T = float, typename Primary = PIController<T>>
    requires SisoLaw<Primary, T>
class SmithPredictor {
    static_assert(N >= 1, "SmithPredictor holds at least one delay state, as InputDelay<N> does.");

    Primary               primary_{};
    damp::array<T, N + 1> ring_{};           ///< The last N + 1 commands, the newest at head_
    size_t                head_{0};          ///< Slot of the newest command
    T                     a_{T{0}};          ///< Discrete FO pole
    T                     b_{T{0}};          ///< The delay-free model's gain
    T                     b_fresh_{T{0}};    ///< The delayed model's gain on u_{k−d}
    T                     b_stale_{T{0}};    ///< Its gain on u_{k−d−1}
    T                     y_hat_{T{0}};      ///< ŷ: the delay-free model's output
    T                     y_delayed_{T{0}};  ///< ŷ_d: the delayed model's output, what the measurement matches when the model is exact
    T                     u_cmd_{T{0}};      ///< The command the last step() returned, or the input track() recorded
    size_t                delay_samples_{0}; ///< d, the delay's whole periods
    bool                  valid_{false};     ///< Design succeeded, d ≤ N, the model finite and b ≠ 0

    /** @brief Push @p u as the newest command and step both models on it. */
    constexpr void advance(T u) {
        head_ = ::damp::detail::ring_next(head_, N);
        ring_[head_] = u;
        y_hat_ = (a_ * y_hat_) + (b_ * u);
        T delayed = (a_ * y_delayed_) + (b_fresh_ * ring_[::damp::detail::ring_back(head_, delay_samples_, N)]);
        if (delay_samples_ < N) {
            delayed += b_stale_ * ring_[::damp::detail::ring_back(head_, delay_samples_ + 1, N)];
        }
        y_delayed_ = delayed;
    }

    /** @brief Whether a model of pole @p a, gains @p b, @p b_fresh, @p b_stale and @p d whole periods runs in N delay states. */
    [[nodiscard]] static constexpr bool runs(T a, T b, T b_fresh, T b_stale, size_t d) {
        const bool finite = damp::isfinite(a) && damp::isfinite(b) && damp::isfinite(b_fresh) && damp::isfinite(b_stale);
        return finite && (b != T{0}) && (d <= N) && (d < N || b_stale == T{0});
    }

public:
    constexpr SmithPredictor() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U, size_t M>
        requires(!damp::is_same_v<U, T>) && (M <= N) && damp::constructible_from<Primary, const design::PIDResult<T>&>
    constexpr SmithPredictor(const design::SmithPredictorResult<U, M>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : SmithPredictor(result.template as<T>()) {}

    /**
     * @brief From a design result: its PID as the primary and its model
     *
     * Requires @p Primary to be constructible from @ref design::PIDResult (true for PIController,
     * PIDController and @ref PController). A design of M ≤ N delay states fits.
     *
     * A failed design does not compile under constant evaluation; at run time the controller is
     * left invalid and @ref step returns zero.
     */
    template<size_t M>
        requires(M <= N) && damp::constructible_from<Primary, const design::PIDResult<T>&>
    constexpr SmithPredictor(const design::SmithPredictorResult<T, M>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : SmithPredictor(Primary(result.pid), result) {}

    /**
     * @brief From any primary law and a design's model: the design's PID is not used
     *
     * @param primary Primary controller instance (any SisoLaw: ADRC, SMC, a custom law, …)
     * @param model   The design whose model, delay and T_s the predictor runs on
     */
    template<size_t M>
        requires(M <= N)
    constexpr SmithPredictor(const Primary& primary, const design::SmithPredictorResult<T, M>& model)
        : primary_(primary),
          a_(model.a),
          b_(model.b),
          b_fresh_(model.b_fresh),
          b_stale_(model.b_stale),
          delay_samples_(model.delay_samples),
          valid_(model.success && runs(model.a, model.b, model.b_fresh, model.b_stale, model.delay_samples)) {
        if (!model.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(model.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (!valid_ && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // a model built by hand that is not finite, or b = 0
        }
    }

    /**
     * @brief Cross-precision converting constructor
     *
     * Copies the primary, the models with the commands in flight, via public accessors (no
     * friend injection: a constrained friend redeclares this class without the @c requires
     * clause).
     *
     * @tparam U         Source scalar type
     * @tparam PrimaryU  Source primary type (must convert into @p Primary)
     */
    template<typename U, typename PrimaryU>
    constexpr explicit SmithPredictor(const SmithPredictor<N, U, PrimaryU>& other)
        requires damp::constructible_from<Primary, const PrimaryU&>
        : primary_(other.primary()),
          head_(other.ring_head()),
          a_(static_cast<T>(other.model_pole())),
          b_(static_cast<T>(other.model_gain())),
          b_fresh_(static_cast<T>(other.model_gains().fresh)),
          b_stale_(static_cast<T>(other.model_gains().stale)),
          y_hat_(static_cast<T>(other.model_output())),
          y_delayed_(static_cast<T>(other.delayed_model_output())),
          u_cmd_(static_cast<T>(other.last_command())),
          delay_samples_(other.delay_samples()),
          valid_(other.valid()) {
        for (size_t i = 0; i < N + 1; ++i) {
            ring_[i] = static_cast<T>(other.commands_in_flight()[i]);
        }
    }

    /**
     * @brief One fixed-rate Smith predictor tick
     *
     * Forms the delay-compensated measurement y_pred = y + ŷ − ŷ_d, runs the primary
     * controller, then queues the command and steps both models on it.
     *
     * @param r Reference
     * @param y Plant measurement (includes true dead time)
     * @return Control command u (zero if the design was invalid)
     */
    [[nodiscard]] constexpr T step(T r, T y) {
        if (!valid_) {
            return T{0};
        }
        const T u = primary_.step(r, y + y_hat_ - y_delayed_);
        advance(u);
        u_cmd_ = u;
        return u;
    }

    /**
     * @brief In place of step(@p r, @p y) while another source drives the plant, @p u the input applied
     *
     * Both models step on @p u, then the primary tracks @p u at the prediction the next step forms
     * from @p y, y + ŷ − ŷ_d after the models moved, so step(@p r, @p y) at this r and y returns
     * what the primary's track lands on (u clamped into its limits, for a PID primary). @p u is
     * recorded as the command applied.
     */
    constexpr void track(T u, T r, T y) {
        if (!valid_) {
            return;
        }
        advance(u);
        primary_.track(u, r, y + y_hat_ - y_delayed_);
        u_cmd_ = u;
    }

    /**
     * @brief Record the input the plant received, @p u_sat, when a stage after this law clamped
     *        its command, and the feedforward @p u_ff added to the command before that clamp
     *
     * The internal models run on the input that reached the plant: the command just queued is
     * replaced by @p u_sat, and a model that has already stepped on it (the delay-free one, and
     * the delayed one with no whole sample of delay) is corrected, so the prediction ŷ − ŷ_d
     * stays the plant's, where a clamp the model did not see left it predicting an output the
     * plant never produced. The primary controller is told the same. last_command() keeps the
     * command step() returned.
     *
     * @see "Advanced PID Control" (Åström & Hägglund, 2006), §7.5
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) {
        if (!valid_) {
            return;
        }
        const T excess = u_sat - ring_[head_];
        y_hat_ += b_ * excess;
        if (delay_samples_ == 0) {
            y_delayed_ += b_fresh_ * excess; // no whole sample of delay: the delayed model has already stepped on this command
        }
        ring_[head_] = u_sat;
        primary_.back_calculate(u_sat, u_ff);
    }

    /** @brief The command the last step() returned, or the input the last track() recorded. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd_; }

    /** @brief The primary's limits: the predictor's command is the primary's. */
    [[nodiscard]] constexpr auto limits() const { return primary_.limits(); }

    /** @brief Set the primary's limits at run time; refused (false, nothing changed) as the primary refuses them. */
    template<typename Limits>
        requires requires(Primary& p, const Limits& l) { p.set_limits(l); }
    [[nodiscard]] constexpr bool set_limits(const Limits& lim) {
        return primary_.set_limits(lim);
    }

    /** @brief Clear the primary, the commands in flight and both models. */
    constexpr void reset() {
        primary_.reset();
        for (auto& u : ring_) {
            u = T{0};
        }
        head_ = 0;
        y_hat_ = T{0};
        y_delayed_ = T{0};
        u_cmd_ = T{0};
    }

    /** @brief Whether the design succeeded, its delay fits N and the primary runs. */
    [[nodiscard]] constexpr bool valid() const { return valid_ && primary_.valid(); }
    /** @brief ŷ: the delay-free model's output, which the next tick adds to y. */
    [[nodiscard]] constexpr T model_output() const { return y_hat_; }
    /** @brief ŷ_d: the delayed model's output, what the measurement matches when the model is exact. */
    [[nodiscard]] constexpr T delayed_model_output() const { return y_delayed_; }
    /** @brief Discrete FO pole a. */
    [[nodiscard]] constexpr T model_pole() const { return a_; }
    /** @brief The delay-free model's gain b = K(1 − a). */
    [[nodiscard]] constexpr T model_gain() const { return b_; }

    /** @brief The delayed model's gains: on u_{k−d} and on u_{k−d−1}. */
    struct DelayedGains {
        T fresh{T{0}}; ///< On u_{k−d}
        T stale{T{0}}; ///< On u_{k−d−1}
    };

    /** @brief The delayed model's gains on u_{k−d} and u_{k−d−1}. */
    [[nodiscard]] constexpr DelayedGains model_gains() const { return {b_fresh_, b_stale_}; }
    /** @brief Delay d in whole samples. */
    [[nodiscard]] constexpr size_t delay_samples() const { return delay_samples_; }
    /** @brief The commands in flight: the last N + 1, the newest at ring_head(). */
    [[nodiscard]] constexpr const damp::array<T, N + 1>& commands_in_flight() const { return ring_; }
    /** @brief The slot of the newest command in flight. */
    [[nodiscard]] constexpr size_t ring_head() const { return head_; }
    /** @brief The primary controller. */
    [[nodiscard]] constexpr Primary& primary() { return primary_; }
    /** @brief The primary controller. */
    [[nodiscard]] constexpr const Primary& primary() const { return primary_; }
};

/** @brief Deduce the runtime from its design: its delay states and scalar, a PIController primary. */
template<typename T, size_t N>
SmithPredictor(const design::SmithPredictorResult<T, N>&) -> SmithPredictor<N, T>;

/** @brief Deduce the runtime from a primary and a design's model: its delay states and scalar, that primary. */
template<typename Primary, typename T, size_t N>
SmithPredictor(const Primary&, const design::SmithPredictorResult<T, N>&) -> SmithPredictor<N, T, Primary>;

namespace detail {

/** @brief Whether @p R is a Smith predictor design result. */
template<typename R>
inline constexpr bool is_smith_result = false;

/** @brief design::smith_predictor's result. */
template<typename U, size_t N>
inline constexpr bool is_smith_result<design::SmithPredictorResult<U, N>> = true;

/** @brief The delay states of a Smith predictor design. */
template<typename R>
inline constexpr size_t smith_delay_states = 0;

/** @brief design::smith_predictor's: its InputDelay's N. */
template<typename U, size_t N>
inline constexpr size_t smith_delay_states<design::SmithPredictorResult<U, N>> = N;

/** @brief A Smith predictor design (design::smith_predictor). */
template<const auto& Res>
concept smith_design = is_smith_result<design_t<Res>>;

} // namespace detail

/** @brief Law<kRes, T> of a Smith predictor design initialized at start-up: the SmithPredictor it runs as, a PIController primary. A runtime assigned replaces it. */
template<const auto& Res, typename T>
    requires detail::smith_design<Res>
struct Law<Res, T> : SmithPredictor<detail::smith_delay_states<detail::design_t<Res>>, T> {
    using Runtime = SmithPredictor<detail::smith_delay_states<detail::design_t<Res>>, T>; ///< The runtime the design runs as

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res} {}

    using Runtime::operator=;
};

/**
 * @brief Law<kRes, T> of a `constexpr` Smith predictor design: the same runtime, a failed design
 *        refused wherever the Law is declared
 */
// ponytail: a constexpr design's model and its PID's gains could be constants of the type, as a PID Law's are; they stay members here.
template<const auto& Res, typename T>
    requires(detail::smith_design<Res> && detail::constant_design<Res>)
struct Law<Res, T> : SmithPredictor<detail::smith_delay_states<detail::design_t<Res>>, T> {
    using Runtime = SmithPredictor<detail::smith_delay_states<detail::design_t<Res>>, T>; ///< The runtime the design runs as

    static_assert(detail::constant_design_runs(Res.success, Res.status, Res.Ts));

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res} {}

    using Runtime::operator=;
};

} // namespace damp

#pragma pop_macro("PI")
