// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file smith_predictor.hpp
 * @brief Smith predictor dead-time compensator (SISO primary + FOPDT model)
 *
 * Classic structure (Smith 1957): a primary controller @f$C@f$ designed for the
 * delay-free plant model @f$\hat{G}@f$ is wrapped with a parallel model path so the
 * feedback signal seen by @f$C@f$ is free of transport delay when the model matches.
 *
 * Plant: @f$ G_p(s) = G(s)\,e^{-sL} @f$. Internal model: delay-free first-order
 * @f$ \hat{G}(s) = K/(\tau s + 1) @f$ plus a pure delay of @f$d = \mathrm{round}(L/T_s)@f$
 * samples. The model state @f$x_k@f$ runs on the command that reaches the plant, so its
 * output is the delayed model output @f$\hat{y}_d@f$, and the delay-free output is the same
 * model carried over the @f$d@f$ commands in flight. Per tick:
 *
 * @f[
 *   \hat{y}_d = x_k,\quad
 *   \hat{y} = x_{k+d} = a^d x_k + \sum_{j=0}^{d-1} a^{d-1-j} b\, u_{k-d+j},\quad
 *   y_{\mathrm{pred}} = y + \hat{y} - \hat{y}_d,\quad
 *   u_k = C(r,\, y_{\mathrm{pred}}),\quad
 *   x_{k+1} = a\,x_k + b\,u_{k-d}
 * @f]
 *
 * with ZOH coeffs @f$ a = e^{-T_s/\tau} @f$, @f$ b = K(1-a) @f$.
 *
 * It is the dead-time compensator for a single loop whose delay spans many samples. The
 * state-space laws (LQR, LQI, LQG, LQGI, OFLQG, MPC, OFMPC) take an input delay of N samples in the
 * model `discretize(sys, Ts, InputDelay<N>{θ})` returns, one state per register. The classic Smith structure is
 * not internally stable for an open-loop unstable plant: the model path ŷ − ŷ_d runs open
 * loop, so the plant's unstable pole stays in the map from a load disturbance to the output.
 *
 * The primary @f$C@f$ is any SisoController (default PIController). The
 * FOPDT designer packs a discrete PI via SIMC for the common process-control
 * path; pass a custom primary when pairing Smith with ADRC, lead-lag, etc.
 *
 * Example: FOPDT process with default Smith PI
 * @code
 * #include "damp/controllers/smith_predictor.hpp"
 *
 * using namespace damp;
 *
 * // Plant G(s) = 1/(5s+1) * e^{-2s}, sample at 10 Hz
 * constexpr auto sp = design::smith_predictor_from_fopdt(
 *     1.0, 5.0, 2.0, 0.1, 2.0); // K, tau, L, Ts, tau_c
 * static constinit SmithPredictor<64, float> ctrl{sp.as<float>()}; // a failed design does not compile
 * // In loop: float u = ctrl.control(r, y);
 * @endcode
 *
 * Example: custom primary (any SisoController)
 * @code
 * LeadLagController<float> c{...}; // design + discretize
 * SmithPredictor<64, float, LeadLagController<float>> smith(c, a, b, d);
 * @endcode
 *
 * @see O.J.M. Smith, "Closer control of loops with dead time," Chem. Eng. Prog., 1957
 * @see Åström & Hägglund, "Advanced PID Control" (2006), §7.5 (Smith predictor)
 * @see "Control of Dead-time Processes" (Normey-Rico & Camacho, 2007): internal stability of the
 *      Smith structure on unstable plants
 * @see design::simc for the delay-free PI/PID tuning used by the FOPDT designer
 * @see SisoController, Cascade — same primary/composition shape
 */

#include <concepts>
#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

namespace damp {
namespace design {

/**
 * @brief Smith predictor design result (discrete PID + FO model + delay samples)
 *
 * Produced by @ref smith_predictor_from_fopdt / @ref smith_predictor_from_pid.
 * The packed @ref DiscretePIDGains is the convenience primary for the process
 * path; construct @ref SmithPredictor from it (default PIController) or
 * build a custom SisoController and pass model coeffs separately.
 *
 * @tparam T Scalar type (design default: double)
 */
template<typename T = double>
struct [[nodiscard]] SmithPredictorResult {
    DiscretePIDGains<T> pid{}; ///< Convenience primary (delay-free design, discretized)

    T            a{T{0}};                             ///< Discrete FO model pole a = e^{−T_s/τ}
    T            b{T{0}};                             ///< Discrete FO model gain b = K(1−a) (ZOH)
    size_t       delay_samples{0};                    ///< Integer delay d = round(L/T_s)
    T            Ts{T{0}};                            ///< Sample period [s]
    T            K{T{0}};                             ///< Continuous plant static gain
    T            tau{T{0}};                           ///< Continuous plant time constant [s]
    T            L{T{0}};                             ///< Continuous dead time [s]
    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /**
     * @brief Convert all scalar fields to another precision
     * @tparam U Target scalar type
     */
    template<typename U>
    [[nodiscard]] constexpr SmithPredictorResult<U> as() const {
        return SmithPredictorResult<U>{
            pid.template as<U>(),
            static_cast<U>(a),
            static_cast<U>(b),
            delay_samples,
            static_cast<U>(Ts),
            static_cast<U>(K),
            static_cast<U>(tau),
            static_cast<U>(L),
            success,
            status,
        };
    }

    /**
     * @brief Continuous primary C(s) recovered from the packed discrete PID
     *
     * The delay is not in this rational map (Smith hides L from C). Failed
     * designs return a zero PID.
     */
    [[nodiscard]] constexpr PIDGains<T> primary_pid() const {
        PIDGains<T> c{};
        if (!success || !(Ts > T{0})) {
            return c;
        }
        c.Kp = pid.Kp;
        c.Ki = pid.Ki / Ts;
        c.Kd = pid.Kd;
        c.Tf = pid.Tf;
        c.b = pid.b;
        c.c = pid.c;
        return c;
    }

    /**
     * @brief Delay-free primary transfer function (SIMC / packed PID)
     */
    [[nodiscard]] constexpr TransferFunction<3, 3, T> to_tf() const { return primary_pid().to_tf(); }

    /**
     * @brief Delay-free primary state-space (same map as PIDGains::to_ss)
     *
     * @return The realization, or nullopt when the primary has a derivative without a filter
     */
    [[nodiscard]] constexpr damp::optional<StateSpace<2, 1, 1, T>> to_ss() const { return primary_pid().to_ss(); }

    /**
     * @brief Delay-free FOPDT model Ĝ(s) = K / (τ s + 1)
     */
    [[nodiscard]] constexpr TransferFunction<1, 2, T> plant_tf() const {
        return TransferFunction<1, 2, T>{
            .num = {K},
            .den = {tau, T{1}},
        };
    }
};

namespace detail {

/**
 * @brief Why the FOPDT plant (K, τ, L) sampled at Ts has no Smith model, or ok
 *
 * Ts must be positive and finite (invalid_sample_time); K finite and non-zero, τ positive and
 * finite, L non-negative and finite, and L/Ts small enough that round(L/Ts) fits size_t
 * (invalid_parameter).
 */
template<typename T>
[[nodiscard]] constexpr DesignStatus smith_plant_status(T K, T tau, T L, T Ts) {
    if (!damp::finite_positive(Ts)) {
        return DesignStatus::invalid_sample_time;
    }
    const bool ok = damp::isfinite(K) && (K != T{0}) && damp::finite_positive(tau) && damp::finite_non_negative(L)
                 && (L / Ts < static_cast<T>(std::numeric_limits<size_t>::max() / 2));
    return ok ? DesignStatus::ok : DesignStatus::invalid_parameter;
}

/**
 * @brief Pack the discrete primary @p pid with the ZOH model of a plant smith_plant_status accepted
 *
 * Carries a failed @p pid's status, and refuses a model gain that is zero or not finite.
 */
template<typename T>
[[nodiscard]] constexpr SmithPredictorResult<T> smith_result(const DiscretePIDGains<T>& pid, T K, T tau, T L, T Ts) {
    if (!pid.success) {
        return ::damp::detail::refused<SmithPredictorResult<T>>(pid.status);
    }
    // ZOH discrete FO: G(z) = b / (z - a)
    const T a = damp::exp(-Ts / tau);
    const T b = K * (T{1} - a);
    if (!damp::isfinite(b) || b == T{0}) {
        return ::damp::detail::refused<SmithPredictorResult<T>>(DesignStatus::invalid_parameter);
    }
    SmithPredictorResult<T> out{
        .pid = pid,
        .a = a,
        .b = b,
        // Integer delay samples (round L/Ts). L = 0 → d = 0 (delay buffer is pass-through).
        .delay_samples = static_cast<size_t>((L / Ts) + static_cast<T>(0.5)),
        .Ts = Ts,
        .K = K,
        .tau = tau,
        .L = L,
    };
    ::damp::detail::set_status(out, DesignStatus::ok);
    return out;
}

} // namespace detail

/**
 * @brief Design a Smith predictor from FOPDT parameters
 *
 * Builds a ZOH-discretized first-order model of the delay-free dynamics and
 * tunes a continuous PI on that delay-free plant with SIMC (@p L treated as
 * zero for the primary controller), then discretizes it at @p Ts. SIMC on the
 * delay-free first-order plant is a PI (its derivative heuristic is zero at L = 0);
 * pack another primary with @ref smith_predictor_from_pid.
 *
 * @note Compare with the textbook Smith structure: primary controller designed
 *       for @f$G(s)@f$ only; the delay is restored by the parallel model path.
 *
 * @see simc — delay-free SIMC tuning (L = 0)
 * @see smith_predictor_from_pid
 * @see O.J.M. Smith (1957); Åström & Hägglund (2006), §7.5
 *
 * @param K     Static gain of the FOPDT plant (finite, non-zero)
 * @param tau   Time constant [s] (finite, > 0)
 * @param L     Dead time [s] (finite, ≥ 0)
 * @param Ts    Sample period [s] (finite, > 0)
 * @param tau_c Desired closed-loop time constant [s] for SIMC on the delay-free
 *              plant. If ≤ 0, defaults to @p tau.
 * @return SmithPredictorResult with discrete PI, model coeffs, and delay samples; refused with
 *         invalid_sample_time for a bad Ts and invalid_parameter for a bad plant or tau_c
 */
template<typename T = double>
[[nodiscard]] constexpr SmithPredictorResult<T> smith_predictor_from_fopdt(T K, T tau, T L, T Ts, T tau_c = T{0}) {
    const DesignStatus plant = detail::smith_plant_status(K, tau, L, Ts);
    if (plant != DesignStatus::ok) {
        return ::damp::detail::refused<SmithPredictorResult<T>>(plant);
    }

    const T tau_c_eff = (tau_c > T{0}) ? tau_c : tau;

    // Primary controller on delay-free FO plant (Smith: hide L from C)
    const PIDGains<T> cont = simc(K, T{0}, tau, tau_c_eff, PIDMode::PI);
    if (!cont.success) {
        return ::damp::detail::refused<SmithPredictorResult<T>>(cont.status);
    }
    return detail::smith_result(cont.discretize(Ts), K, tau, L, Ts);
}

/**
 * @brief Pack a Smith predictor from an existing continuous PID and FOPDT plant
 *
 * Use when gains come from another designer (pole placement, manual tuning).
 * The continuous PID is discretized at @p Ts; the FO model and delay samples
 * are built from @p K, @p tau, @p L exactly as in @ref smith_predictor_from_fopdt.
 *
 * @param pid Continuous primary-controller design (for the delay-free plant)
 * @param K   Static gain
 * @param tau Time constant [s]
 * @param L   Dead time [s]
 * @param Ts  Sample period [s]
 * @return SmithPredictorResult ready for @ref SmithPredictor with a PID primary; refused as
 *         @ref smith_predictor_from_fopdt refuses the plant, and with the status of a failed
 *         @p pid or of its discretization
 */
template<typename T = double>
[[nodiscard]] constexpr SmithPredictorResult<T> smith_predictor_from_pid(
    const PIDGains<T>& pid,
    T                  K,
    T                  tau,
    T                  L,
    T                  Ts
) {
    const DesignStatus plant = detail::smith_plant_status(K, tau, L, Ts);
    if (plant != DesignStatus::ok) {
        return ::damp::detail::refused<SmithPredictorResult<T>>(plant);
    }
    return detail::smith_result(pid.discretize(Ts), K, tau, L, Ts);
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Smith predictor runtime: SISO primary + FO model + pure delay
 *
 * Per-tick dead-time compensator. Stores any SisoController primary (default PIController)
 * and a scalar first-order ZOH model whose dead time is a ring of the last d + 1 commands (at
 * most @p MaxDelay + 1): the model state x runs on the command that reaches the plant, so its
 * output is the delayed model output ŷ_d, and the model carried over the d commands in flight
 * gives the delay-free output ŷ = x_{k+d}.
 *
 * Matches the Cascade composition style: the Smith residual feeds
 * @c primary.control(r, y_pred); the primary type is a template parameter.
 *
 * @tparam MaxDelay Capacity in whole samples (must be ≥ 1): the realizable delay is
 *                  0 … MaxDelay samples.
 * @tparam T        Scalar type (embedded default: float)
 * @tparam Primary  Primary controller (SisoController). Default:
 *                  PIController for the FOPDT/SIMC process path. Use
 *                  PIDController (full PID), ADRC, lead-lag, etc. when needed.
 *
 * @see design::smith_predictor_from_fopdt
 * @see SisoController
 * @see O.J.M. Smith (1957); Åström & Hägglund (2006), §7.5
 */
template<size_t MaxDelay, typename T = float, typename Primary = PIController<T>>
    requires SisoController<Primary, T>
class SmithPredictor {
    static_assert(MaxDelay >= 1, "SmithPredictor requires MaxDelay >= 1");

    Primary                      primary_{};
    damp::array<T, MaxDelay + 1> ring_{};           ///< The last d + 1 commands in the first d + 1 slots, oldest at head_
    size_t                       head_{0};          ///< Index of the oldest command
    T                            a_{T{0}};          ///< Discrete FO pole
    T                            b_{T{0}};          ///< Discrete FO gain
    T                            x_{T{0}};          ///< Model state on the delayed command: x_{k+1} = a x_k + b u_{k−d}, output ŷ_d
    bool                         valid_{false};     ///< Design succeeded, d ≤ MaxDelay and b ≠ 0
    size_t                       delay_samples_{0}; ///< d, as requested

    /** @brief d in the ring: the delay, held at MaxDelay. */
    [[nodiscard]] constexpr size_t periods() const { return (delay_samples_ < MaxDelay) ? delay_samples_ : MaxDelay; }

    /** @brief The slot after @p i in the ring of d + 1 (i ≤ d). */
    [[nodiscard]] constexpr size_t next(size_t i) const { return ::damp::detail::ring_next(i, periods()); }

public:
    constexpr SmithPredictor() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!std::is_same_v<U, T>) && std::constructible_from<Primary, const design::DiscretePIDGains<T>&>
    constexpr SmithPredictor(const design::SmithPredictorResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : SmithPredictor(result.template as<T>()) {}

    /**
     * @brief Construct from a design result (PID-family primary)
     *
     * @param result Design from @ref design::smith_predictor_from_fopdt (or from_pid)
     *
     * Requires @p Primary to be constructible from @ref design::DiscretePIDGains
     * (true for PIController / PIDController / @ref PController).
     *
     * If @p result.success is false or @p delay_samples exceeds @p MaxDelay, the design does not
     * compile under constant evaluation; at run time the controller is left invalid and
     * @ref control returns zero.
     */
    constexpr SmithPredictor(const design::SmithPredictorResult<T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        requires std::constructible_from<Primary, const design::DiscretePIDGains<T>&>
        : SmithPredictor(Primary(result.pid), result.a, result.b, result.delay_samples) {
        valid_ = valid_ && result.success;
        if (!valid_ && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_parameter : result.status); // constant-evaluated (constinit / constexpr): a failed design, or a delay past MaxDelay, does not compile
        }
    }

    /**
     * @brief Construct from discrete PID gains, FO model coeffs, and delay samples
     *
     * Convenience for PID-family @p Primary (same constructibility requirement as
     * the design-result constructor).
     *
     * @param pid            Discretized primary controller gains
     * @param a              Discrete FO pole
     * @param b              Discrete FO gain
     * @param delay_samples  Integer pure delay (≤ MaxDelay)
     */
    constexpr SmithPredictor(
        const design::DiscretePIDGains<T>& pid,
        T                                  a,
        T                                  b,
        size_t                             delay_samples
    )
        requires std::constructible_from<Primary, const design::DiscretePIDGains<T>&>
        : SmithPredictor(Primary(pid), a, b, delay_samples) {}

    /**
     * @brief Construct from an arbitrary SISO primary and FOPDT model coeffs
     *
     * Use when @p Primary is not built from @ref design::DiscretePIDGains
     * (ADRC, lead-lag, custom law, etc.).
     *
     * @param primary        Primary controller instance
     * @param a              Discrete FO pole
     * @param b              Discrete FO gain
     * @param delay_samples  Integer pure delay (≤ MaxDelay)
     */
    constexpr SmithPredictor(const Primary& primary, T a, T b, size_t delay_samples)
        : primary_(primary),
          a_(a),
          b_(b),
          valid_(delay_samples <= MaxDelay && b != T{0}),
          delay_samples_(delay_samples) {}

    /**
     * @brief Cross-precision converting constructor
     *
     * Copies the primary, the model with its commands in flight, and the model state via
     * public accessors (no friend injection: a constrained friend redeclares this class
     * without the @c requires clause).
     *
     * @tparam U         Source scalar type
     * @tparam PrimaryU  Source primary type (must convert into @p Primary)
     */
    template<typename U, typename PrimaryU>
    constexpr explicit SmithPredictor(const SmithPredictor<MaxDelay, U, PrimaryU>& other)
        requires std::constructible_from<Primary, const PrimaryU&>
        : primary_(other.primary()),
          head_(other.ring_head()),
          a_(static_cast<T>(other.model_pole())),
          b_(static_cast<T>(other.model_gain())),
          x_(static_cast<T>(other.delayed_model_output())),
          valid_(other.valid()),
          delay_samples_(other.delay_samples()) {
        for (size_t i = 0; i < MaxDelay + 1; ++i) {
            ring_[i] = static_cast<T>(other.commands_in_flight()[i]);
        }
    }

    /**
     * @brief One fixed-rate Smith predictor tick
     *
     * Forms the delay-compensated measurement y_pred = y + ŷ − ŷ_d, runs the primary
     * controller, then queues the command and steps the model on the command that reaches
     * the plant this period.
     *
     * @param r Reference
     * @param y Plant measurement (includes true dead time)
     * @return Control command u (zero if the design was invalid)
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        if (!valid_) {
            return T{0};
        }

        const T y_delayed = x_;         // ŷ_d: the model output on the delayed command
        const T y_hat = model_output(); // ŷ = x_{k+d}: the same model past the commands in flight
        const T y_pred = y + y_hat - y_delayed;

        const T u = primary_.control(r, y_pred);

        ring_[head_] = u; // the oldest slot takes u_k; the oldest is now u_{k−d}
        head_ = next(head_);
        x_ = (a_ * x_) + (b_ * ring_[head_]);

        return u;
    }

    /** @brief Clear primary, commands in flight, and model state */
    constexpr void reset() {
        primary_.reset();
        for (auto& u : ring_) {
            u = T{0};
        }
        head_ = 0;
        x_ = T{0};
    }

    /** @brief Whether the design succeeded and its delay fits MaxDelay. */
    [[nodiscard]] constexpr bool valid() const { return valid_; }
    /** @brief ŷ: the delay-free model output the next tick uses (the model state predicted over the d commands in flight). */
    [[nodiscard]] constexpr T model_output() const {
        T      x = x_;
        size_t i = head_;
        for (size_t j = 0; j < periods(); ++j) {
            i = next(i);
            x = (a_ * x) + (b_ * ring_[i]);
        }
        return x;
    }
    /** @brief ŷ_d: the model output on the delayed command, what the measurement matches when the model is exact. */
    [[nodiscard]] constexpr T delayed_model_output() const { return x_; }
    /** @brief Discrete FO pole a. */
    [[nodiscard]] constexpr T model_pole() const { return a_; }
    /** @brief Discrete FO gain b. */
    [[nodiscard]] constexpr T model_gain() const { return b_; }
    /** @brief Delay d in samples. */
    [[nodiscard]] constexpr size_t delay_samples() const { return delay_samples_; }
    /** @brief The commands in flight: the last d + 1 in the first d + 1 slots, oldest at ring_head(). */
    [[nodiscard]] constexpr const damp::array<T, MaxDelay + 1>& commands_in_flight() const { return ring_; }
    /** @brief The slot of the oldest command in flight. */
    [[nodiscard]] constexpr size_t ring_head() const { return head_; }
    /** @brief The primary controller. */
    [[nodiscard]] constexpr Primary& primary() { return primary_; }
    /** @brief The primary controller. */
    [[nodiscard]] constexpr const Primary& primary() const { return primary_; }
};

} // namespace damp
