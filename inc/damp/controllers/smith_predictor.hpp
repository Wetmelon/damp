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
 * samples. Per tick:
 *
 * @f[
 *   \hat{y} = x,\quad
 *   y_{\mathrm{pred}} = y + \hat{y} - \hat{y}_d,\quad
 *   u = C(r,\, y_{\mathrm{pred}}),\quad
 *   x \leftarrow a\,x + b\,u
 * @f]
 *
 * with ZOH coeffs @f$ a = e^{-T_s/\tau} @f$, @f$ b = K(1-a) @f$, and
 * @f$ \hat{y}_d @f$ the @f$d@f$-sample delayed model output.
 *
 * The primary @f$C@f$ is any SisoController (default PIController). The
 * FOPDT designer packs a discrete PI/PID via SIMC for the common process-control
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
 * static_assert(sp.success);
 *
 * SmithPredictor<64, float> ctrl(sp.as<float>());
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
 * @see design::simc for the delay-free PI/PID tuning used by the FOPDT designer
 * @see SisoController, Cascade — same primary/composition shape
 */

#include <concepts>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/filters/delay.hpp"
#include "damp/math/math.hpp"

namespace damp {
namespace design {

/**
 * @brief Smith predictor design result (discrete PID + FO model + delay samples)
 *
 * Produced by @ref smith_predictor_from_fopdt / @ref smith_predictor_from_pid.
 * The packed @ref DiscretePIDResult is the convenience primary for the process
 * path; construct @ref SmithPredictor from it (default PIController) or
 * build a custom SisoController and pass model coeffs separately.
 *
 * @tparam T Scalar type (design default: double)
 */
template<typename T = double>
struct SmithPredictorResult {
    DiscretePIDResult<T> pid{}; ///< Convenience primary (delay-free design, discretized)

    T      a{T{0}};          ///< Discrete FO model pole a = e^{−T_s/τ}
    T      b{T{0}};          ///< Discrete FO model gain b = K(1−a) (ZOH)
    size_t delay_samples{0}; ///< Integer delay d = round(L/T_s)
    T      Ts{T{0}};         ///< Sample period [s]
    T      K{T{0}};          ///< Continuous plant static gain
    T      tau{T{0}};        ///< Continuous plant time constant [s]
    T      L{T{0}};          ///< Continuous dead time [s]
    bool   success{false};   ///< true if parameters produced a valid design

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
        };
    }

    /**
     * @brief Continuous primary C(s) recovered from the packed discrete PID
     *
     * The delay is not in this rational map (Smith hides L from C). Failed
     * designs return a zero PID.
     */
    [[nodiscard]] constexpr PIDResult<T> primary_pid() const {
        PIDResult<T> c{};
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
     * @brief PI state-space of the primary (same map as PIDResult::to_ss)
     */
    [[nodiscard]] constexpr StateSpace<1, 1, 1, T> to_ss() const { return primary_pid().to_ss(); }

    /**
     * @brief Delay-free FOPDT model Ĝ(s) = K / (τ s + 1)
     */
    [[nodiscard]] constexpr TransferFunction<1, 2, T> plant_tf() const {
        return TransferFunction<1, 2, T>{
            .num = {K},
            .den = {T{1}, tau},
        };
    }
};

/**
 * @brief Design a Smith predictor from FOPDT parameters
 *
 * Builds a ZOH-discretized first-order model of the delay-free dynamics and
 * tunes a continuous PI/PID on that delay-free plant with SIMC (@p L treated as
 * zero for the primary controller), then discretizes the PID at @p Ts.
 *
 * @note Compare with the textbook Smith structure: primary controller designed
 *       for @f$G(s)@f$ only; the delay is restored by the parallel model path.
 *
 * @see simc — delay-free SIMC tuning (L = 0)
 * @see smith_predictor_from_pid
 * @see O.J.M. Smith (1957); Åström & Hägglund (2006), §7.5
 *
 * @param K     Static gain of the FOPDT plant (must be non-zero)
 * @param tau   Time constant [s] (must be > 0)
 * @param L     Dead time [s] (must be ≥ 0)
 * @param Ts    Sample period [s] (must be > 0)
 * @param tau_c Desired closed-loop time constant [s] for SIMC on the delay-free
 *              plant. If ≤ 0, defaults to @p tau.
 * @param type  Convenience primary type (default PI; P/PD/PID also supported)
 * @return SmithPredictorResult with discrete PID, model coeffs, and delay samples
 */
template<typename T = double>
[[nodiscard]] constexpr SmithPredictorResult<T> smith_predictor_from_fopdt(
    T       K,
    T       tau,
    T       L,
    T       Ts,
    T       tau_c = T{0},
    PIDType type = PIDType::PI
) {
    if (!(tau > T{0}) || !(Ts > T{0}) || !(L >= T{0}) || K == T{0}) {
        return SmithPredictorResult<T>{
            .Ts = Ts,
            .K = K,
            .tau = tau,
            .L = L,
        };
    }

    const T tau_c_eff = (tau_c > T{0}) ? tau_c : tau;

    // Primary controller on delay-free FO plant (Smith: hide L from C)
    const PIDResult<T> cont = simc(K, T{0}, tau, tau_c_eff, type);
    if (cont.Kp == T{0} && cont.Ki == T{0} && cont.Kd == T{0}) {
        // simc rejected the spec
        return SmithPredictorResult<T>{
            .Ts = Ts,
            .K = K,
            .tau = tau,
            .L = L,
        };
    }

    // ZOH discrete FO: G(z) = b / (z - a)
    const T a = damp::exp(-Ts / tau);
    const T b = K * (T{1} - a);
    // Integer delay samples (round L/Ts). L = 0 → d = 0 (delay buffer is pass-through).
    const size_t delay_samples = static_cast<size_t>((L / Ts) + static_cast<T>(0.5));
    const auto   pid = cont.discretize(Ts);

    return SmithPredictorResult<T>{
        .pid = pid,
        .a = a,
        .b = b,
        .delay_samples = delay_samples,
        .Ts = Ts,
        .K = K,
        .tau = tau,
        .L = L,
        .success = damp::isfinite(a) && damp::isfinite(b) && damp::isfinite(pid.Kp) && (b != T{0}),
    };
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
 * @return SmithPredictorResult ready for @ref SmithPredictor with a PID primary
 */
template<typename T = double>
[[nodiscard]] constexpr SmithPredictorResult<T> smith_predictor_from_pid(
    const PIDResult<T>& pid,
    T                   K,
    T                   tau,
    T                   L,
    T                   Ts
) {
    if (!(tau > T{0}) || !(Ts > T{0}) || !(L >= T{0}) || K == T{0}) {
        return SmithPredictorResult<T>{
            .Ts = Ts,
            .K = K,
            .tau = tau,
            .L = L,
        };
    }

    const T      a = damp::exp(-Ts / tau);
    const T      b = K * (T{1} - a);
    const size_t delay_samples = static_cast<size_t>((L / Ts) + static_cast<T>(0.5));
    const auto   disc_pid = pid.discretize(Ts);

    return SmithPredictorResult<T>{
        .pid = disc_pid,
        .a = a,
        .b = b,
        .delay_samples = delay_samples,
        .Ts = Ts,
        .K = K,
        .tau = tau,
        .L = L,
        .success = damp::isfinite(a) && damp::isfinite(b) && (b != T{0}),
    };
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Smith predictor runtime: SISO primary + FO model + pure delay
 *
 * Lightweight per-tick dead-time compensator. Stores any SisoController
 * primary (default PIController), a scalar first-order ZOH plant model, and
 * a @ref Delay line of length @p MaxDelay for @f$\hat{G}\,e^{-sL}@f$.
 *
 * Matches the Cascade composition style: the Smith residual feeds
 * @c primary.control(r, y_pred); the primary type is a template parameter.
 *
 * @tparam MaxDelay Delay buffer capacity (must be ≥ 1). Realizable delay is
 *                  @f$0\ldots\texttt{MaxDelay}-1@f$ samples (see @ref Delay).
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

    Primary            primary_{};
    Delay<MaxDelay, T> delay_{};

    T      a_{T{0}}; ///< Discrete FO pole
    T      b_{T{0}}; ///< Discrete FO input gain
    T      x_{T{0}}; ///< Delay-free model state (= ŷ)
    size_t delay_samples_{0};
    bool   valid_{false};

public:
    constexpr SmithPredictor() = default;

    /**
     * @brief Construct from a design result (PID-family primary)
     *
     * @param result Design from @ref design::smith_predictor_from_fopdt (or from_pid)
     *
     * Requires @p Primary to be constructible from @ref design::DiscretePIDResult
     * (true for PIController / PIDController / @ref PController).
     *
     * If @p result.success is false or @p delay_samples is not in
     * @f$[0,\texttt{MaxDelay}-1]@f$, the controller is left invalid and
     * @ref control returns zero.
     */
    constexpr explicit SmithPredictor(const design::SmithPredictorResult<T>& result)
        requires std::constructible_from<Primary, const design::DiscretePIDResult<T>&>
        : primary_(result.pid),
          a_(result.a),
          b_(result.b),
          delay_samples_(result.delay_samples),
          valid_(result.success && result.delay_samples < MaxDelay) {
        if (valid_) {
            delay_.init(delay_samples_);
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
     * @param delay_samples  Integer pure delay (< MaxDelay)
     */
    constexpr SmithPredictor(
        const design::DiscretePIDResult<T>& pid,
        T                                   a,
        T                                   b,
        size_t                              delay_samples
    )
        requires std::constructible_from<Primary, const design::DiscretePIDResult<T>&>
        : primary_(pid),
          a_(a),
          b_(b),
          delay_samples_(delay_samples),
          valid_(delay_samples < MaxDelay && b != T{0}) {
        if (valid_) {
            delay_.init(delay_samples_);
        }
    }

    /**
     * @brief Construct from an arbitrary SISO primary and FOPDT model coeffs
     *
     * Use when @p Primary is not built from @ref design::DiscretePIDResult
     * (ADRC, lead-lag, custom law, etc.).
     *
     * @param primary        Primary controller instance
     * @param a              Discrete FO pole
     * @param b              Discrete FO gain
     * @param delay_samples  Integer pure delay (< MaxDelay)
     */
    constexpr SmithPredictor(const Primary& primary, T a, T b, size_t delay_samples)
        : primary_(primary),
          a_(a),
          b_(b),
          delay_samples_(delay_samples),
          valid_(delay_samples < MaxDelay && b != T{0}) {
        if (valid_) {
            delay_.init(delay_samples_);
        }
    }

    /**
     * @brief Cross-precision converting constructor
     *
     * Copies primary and model state via public accessors (no friend injection —
     * a constrained friend redeclares this class without the @c requires clause).
     * The delay buffer is re-initialized empty at the same length.
     *
     * @tparam U         Source scalar type
     * @tparam PrimaryU  Source primary type (must convert into @p Primary)
     */
    template<typename U, typename PrimaryU>
    constexpr explicit SmithPredictor(const SmithPredictor<MaxDelay, U, PrimaryU>& other)
        requires std::constructible_from<Primary, const PrimaryU&>
        : primary_(other.primary()),
          a_(static_cast<T>(other.model_pole())),
          b_(static_cast<T>(other.model_gain())),
          x_(static_cast<T>(other.model_output())),
          delay_samples_(other.delay_samples()),
          valid_(other.valid()) {
        if (valid_) {
            delay_.init(delay_samples_);
        }
    }

    /**
     * @brief One fixed-rate Smith predictor tick
     *
     * Forms the delay-compensated measurement @f$y_{\mathrm{pred}} = y + \hat{y}
     * - \hat{y}_d@f$, runs the primary controller, then steps the FO model.
     *
     * @param r Reference
     * @param y Plant measurement (includes true dead time)
     * @return Control command u (zero if the design was invalid)
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        if (!valid_) {
            return T{0};
        }

        const T y_hat = x_;
        const T y_delayed = delay_(y_hat);
        const T y_pred = y + y_hat - y_delayed;

        const T u = primary_.control(r, y_pred);

        // Step delay-free model for the next tick
        x_ = (a_ * x_) + (b_ * u);

        return u;
    }

    /// Clear primary, delay buffer, and model state
    constexpr void reset() {
        primary_.reset();
        delay_.reset();
        x_ = T{0};
    }

    [[nodiscard]] constexpr bool           valid() const { return valid_; }
    [[nodiscard]] constexpr T              model_output() const { return x_; }
    [[nodiscard]] constexpr T              model_pole() const { return a_; }
    [[nodiscard]] constexpr T              model_gain() const { return b_; }
    [[nodiscard]] constexpr size_t         delay_samples() const { return delay_samples_; }
    [[nodiscard]] constexpr Primary&       primary() { return primary_; }
    [[nodiscard]] constexpr const Primary& primary() const { return primary_; }
};

} // namespace damp
