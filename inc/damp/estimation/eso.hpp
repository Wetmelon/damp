// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file eso.hpp
 * @brief Extended-state observer: design::eso → ESOResult → ESO
 *
 * The Han–Gao ESO for the SISO chain y⁽ⁿ⁾ = f + b₀ u. It estimates the output, its
 * derivatives, and the lumped disturbance f (model error, load, unknown input). ADRC is
 * this observer plus a PD law that cancels f̂; the observer is the part people run on its
 * own (a robust differentiator / disturbance estimate without a plant A, B, C).
 *
 * The runtime is the discrete current estimator of Miklosovic, Radke and Gao: the ESO's
 * linear dynamics discretized by ZOH, poles at e^{−ω_o T_s}. predict advances the state,
 * update corrects it with the measurement; step does both.
 *
 * @code
 * constexpr auto res = design::eso<2>(40.0, 20.0, 1e-3); // ω_o, b₀, T_s
 * static constinit ESO obs{res.as<float>()};
 * // ISR: z = [ŷ, ẏ, f̂]
 * const auto z = obs.step(y, u_applied);
 * const float fhat = obs.disturbance();
 * @endcode
 */

#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/const_matrix.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {

/**
 * @brief Extended-state observer design: continuous β and, for Ts > 0, the discrete ESO
 *
 * z̄_k = Φ ẑ_{k−1} + Γ u_{k−1}, ẑ_k = z̄_k + L (y_k − z̄_k,₀). Use .as\<U\>() to convert (e.g.
 * double to float).
 *
 * @see "Discrete implementation and generalization of the extended state observer"
 *      (Miklosovic, Radke & Gao, 2006), ACC
 */
template<size_t NX, typename T = double>
struct [[nodiscard]] ESOResult {
    static_assert(NX == 1 || NX == 2, "ESO supports 1st- and 2nd-order plants only");

    T wo{}; ///< Observer bandwidth (all poles at −ω_o)
    T b0{}; ///< Plant high-frequency gain: y⁽ⁿ⁾ = f + b₀ u

    ColVec<NX + 1, T> beta{}; ///< Continuous ESO gains (poles at −ω_o)

    Matrix<NX + 1, NX + 1, T> Phi{};   ///< ZOH state map Φ of the chain of integrators (zero when Ts = 0)
    ColVec<NX + 1, T>         Gamma{}; ///< ZOH input map Γ on u_{k−1}
    ColVec<NX + 1, T>         L{};     ///< Current-estimator gain on the innovation y_k − z̄_k,₀

    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the design failed
    T            Ts{0};                               ///< [s] Sample time of the runtime: 0 continuous (analysis), > 0 the tick

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief Estimation-error dynamics (I − L C)Φ; eigenvalues e^{−ω_o T_s} (zero when Ts = 0). */
    [[nodiscard]] constexpr Matrix<NX + 1, NX + 1, T> error_dynamics() const {
        Matrix<NX + 1, NX + 1, T> F{};
        for (size_t i = 0; i <= NX; ++i) {
            for (size_t j = 0; j <= NX; ++j) {
                F(i, j) = Phi(i, j) - (L(i) * Phi(0, j)); // C picks row 0
            }
        }
        return F;
    }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return ESOResult<NX, U>{
            static_cast<U>(wo),
            static_cast<U>(b0),
            beta.template as<U>(),
            Phi.template as<U>(),
            Gamma.template as<U>(),
            L.template as<U>(),
            success,
            status,
            static_cast<U>(Ts),
        };
    }
};

namespace detail {

/**
 * @brief Fill the discrete ESO of @p r (Φ, Γ, L) from its b0, wo and Ts > 0: the ZOH model of
 *        the chain of integrators and the current-estimator gain with every pole at e^{−ω_o T_s}.
 *
 * Refuses (invalid_parameter) an ω_o T_s below the resolution of T, where λ = 1 and L = 0 would
 * never correct, and gains that do not fit in T.
 *
 * @see "Discrete implementation and generalization of the extended state observer"
 *      (Miklosovic, Radke & Gao, 2006), ACC
 */
template<size_t NX, typename T>
constexpr void discrete_eso(ESOResult<NX, T>& r) {
    constexpr size_t       n = NX + 1;
    StateSpace<n, 1, 1, T> eso{};
    for (size_t i = 0; i < NX; ++i) {
        eso.A(i, i + 1) = T{1}; // ż_i = ẑ_{i+1}; f̂ enters the highest derivative
    }
    eso.B(NX - 1, 0) = r.b0;
    eso.C(0, 0) = T{1};
    const auto d_opt = ::damp::detail::discretize_zoh_impl(eso, r.Ts); // exact: A is nilpotent
    if (!d_opt) {
        ::damp::detail::set_status(r, DesignStatus::invalid_parameter); // b0 and Ts are checked finite first
        return;
    }
    const auto& d = *d_opt;

    const T h = r.Ts;
    const T lambda = damp::exp(-r.wo * h);
    const T m = T{1} - lambda;
    if (!(m > T{0})) {
        ::damp::detail::set_status(r, DesignStatus::invalid_parameter);
        return;
    }
    if constexpr (NX == 1) {
        r.L(0) = m * (T{1} + lambda);
        r.L(1) = (m * m) / h;
    } else {
        r.L(0) = m * (T{1} + lambda + (lambda * lambda));
        r.L(1) = (T{3} * m * m * (T{1} + lambda)) / (T{2} * h);
        r.L(2) = (m * m * m) / (h * h);
    }

    r.Phi = d.A;
    for (size_t i = 0; i < n; ++i) {
        r.Gamma(i) = d.B(i, 0);
        bool finite = damp::isfinite(r.L(i)) && damp::isfinite(r.Gamma(i));
        for (size_t j = 0; j < n; ++j) {
            finite = finite && damp::isfinite(r.Phi(i, j));
        }
        if (!finite) {
            ::damp::detail::set_status(r, DesignStatus::invalid_parameter);
            return;
        }
    }
}

} // namespace detail

/**
 * @brief Extended-state observer: poles at −ω_o, plant gain b₀, optional sample time
 *
 * Estimates [ŷ, (ẏ,) f̂] for y⁽ⁿ⁾ = f + b₀ u. All observer poles at −ω_o (binomial β).
 * Only plant orders NX = 1 and NX = 2 (ESO dimension NX+1), same as linear ADRC.
 *
 * For Ts > 0 the design also builds the discrete ESO the runtime ticks (Miklosovic, Radke &
 * Gao, 2006). The ESO's linear dynamics @f$\dot z = A z + B u + \beta (y - C z)@f$, with
 * @f$A@f$ the chain of integrators, are discretized by ZOH (@ref discretize), exact for this
 * model, to @f$(\Phi, \Gamma)@f$. The current estimator
 * @f[
 *   \bar z_k = \Phi \hat z_{k-1} + \Gamma u_{k-1},\qquad
 *   \hat z_k = \bar z_k + L\,(y_k - C \bar z_k)
 * @f]
 * has error dynamics @f$(I - L C)\Phi@f$, and @f$L@f$ places all its poles at
 * @f$\lambda = e^{-\omega_o T_s}@f$, the image of the continuous poles at @f$-\omega_o@f$:
 * @f[
 *   n = 2:\ L = \begin{bmatrix} 1-\lambda^2 \\ (1-\lambda)^2/T_s \end{bmatrix},\qquad
 *   n = 3:\ L = \begin{bmatrix} 1-\lambda^3 \\ \tfrac{3}{2}(1-\lambda)^2(1+\lambda)/T_s \\
 *              (1-\lambda)^3/T_s^2 \end{bmatrix}.
 * @f]
 * The result carries @f$\Phi@f$, @f$\Gamma@f$ and @f$L@f$, and the runtime ticks those two
 * equations directly. The innovation @f$y_k - C \bar z_k@f$ is formed first, so a large
 * @f$|y|@f$ does not cancel inside the disturbance estimate. Unlike a forward-Euler step of
 * the continuous ESO, whose poles 1 − ω_o T_s leave the unit circle for ω_o T_s > 2, these
 * poles stay at e^{−ω_o T_s} for any ω_o T_s.
 * Ts = 0 is the continuous design, for analysis; Φ, Γ and L stay zero.
 *
 * Requires wo > 0 and a finite b0, else refused with invalid_parameter. The observer never
 * divides by b0: a negative b0 is a reverse-acting plant, and b0 = 0 with u = 0 is a plain
 * differentiator of y. Also refused with invalid_parameter: an ω_o T_s below the resolution of
 * T (the observer would never correct) and gains that do not fit in T. A negative or non-finite
 * Ts is refused with invalid_sample_time.
 *
 * @param wo Observer bandwidth [rad/s]
 * @param b0 Plant high-frequency gain
 * @param Ts [s] Sample time of the runtime (0 = continuous)
 *
 * @return ESOResult with β and, for Ts > 0, Φ, Γ, L
 *
 * @see "From PID to Active Disturbance Rejection Control" (Han, 1998)
 * @see "Discrete implementation and generalization of the extended state observer"
 *      (Miklosovic, Radke & Gao, 2006), ACC: the current discrete ESO and its gains
 * @see design::adrc — this observer plus a PD law that cancels f̂
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr ESOResult<NX, T> eso(T wo, T b0, T Ts = T{0}) {
    static_assert(NX == 1 || NX == 2, "ESO supports 1st- and 2nd-order plants only");

    const bool ok = damp::isfinite(wo) && damp::isfinite(b0) && (wo > T{0});
    if (!ok) {
        return ::damp::detail::refused<ESOResult<NX, T>>(DesignStatus::invalid_parameter);
    }
    if (!(Ts >= T{0}) || !damp::isfinite(Ts)) {
        return ::damp::detail::refused<ESOResult<NX, T>>(DesignStatus::invalid_sample_time);
    }

    ColVec<NX + 1, T> beta{};
    constexpr size_t  n = NX + 1;
    for (size_t i = 1; i <= n; ++i) {
        T binom = T{1};
        for (size_t k = 1; k <= i; ++k) {
            binom *= static_cast<T>(n - k + 1);
            binom /= static_cast<T>(k);
        }
        const T power = damp::pow(wo, static_cast<int>(i)); // saturates at max()
        if (!(power < (std::numeric_limits<T>::max() / binom))) {
            return ::damp::detail::refused<ESOResult<NX, T>>(DesignStatus::invalid_parameter);
        }
        beta[i - 1] = binom * power;
    }

    ESOResult<NX, T> out{
        .wo = wo,
        .b0 = b0,
        .beta = beta,
        .success = true,
        .status = DesignStatus::ok,
        .Ts = Ts,
    };
    if (Ts > T{0}) {
        detail::discrete_eso(out);
        if (!out.success) {
            return ::damp::detail::refused<ESOResult<NX, T>>(out.status);
        }
    }
    return out;
}

} // namespace design

namespace detail {

/** @brief Where an ESO keeps its maps: in the object, any value at run time. */
template<size_t NX, typename T>
struct HeldEso {
    static constexpr bool constant = false; ///< The matrices are members
    using Phi = Matrix<NX + 1, NX + 1, T>;  ///< Φ
    using Gamma = ColVec<NX + 1, T>;        ///< Γ
    using L = ColVec<NX + 1, T>;            ///< L
};

/**
 * @brief An ESO's maps as constants of the type: the chain of integrators' zeros and ones cost nothing.
 *
 * @tparam Res design::eso's result, or a design that has one in `.as_eso()` (design::adrc), `constexpr`
 */
template<const auto& Res, typename T>
struct ConstantEso {
    static constexpr bool constant = true;               ///< The matrices are constants of the type
    static constexpr auto design = Res.template as<T>(); ///< The design, in the runtime's scalar
    static constexpr auto filter = [] {
        if constexpr (requires { design.as_eso(); }) {
            return design.as_eso();
        } else {
            return design;
        }
    }(); ///< The observer's half of it
    static constexpr auto kPhi = filter.Phi;     ///< Φ
    static constexpr auto kGamma = filter.Gamma; ///< Γ
    static constexpr auto kL = filter.L;         ///< L
    using Phi = ConstMatrix<kPhi>;               ///< Φ
    using Gamma = ConstMatrix<kGamma>;           ///< Γ
    using L = ConstMatrix<kL>;                   ///< L
};

} // namespace detail

/**
 * @ingroup discrete_estimators
 * @brief Extended-state observer (ESO)
 *
 * Discrete Han–Gao ESO for y⁽ⁿ⁾ = f + b₀ u. State is [ŷ, (ẏ,) f̂]. The tick is the
 * ZOH current estimator design::eso builds, in innovation form.
 *
 * Same verbs as the other estimators. predict(u) advances the state one period with the
 * input the plant holds over it; call it alone through a missed sample. update(y) corrects
 * the state without advancing time. step(y, u) is the fused tick: u is the input applied over
 * the last period, then y corrects.
 *
 * Low-latency ISR: update(y), compute and write the command, back_calculate(u_sat), predict().
 * step(y) is predict() then update(y).
 *
 * @see "Discrete implementation and generalization of the extended state observer"
 *      (Miklosovic, Radke & Gao, 2006), ACC
 * @tparam NX Plant order: 1 or 2 (ESO dimension NX+1)
 * @tparam T  Scalar type (default: float)
 */
template<size_t NX, typename T = float, typename Store = detail::HeldEso<NX, T>>
class ESO {
    static_assert(NX == 1 || NX == 2, "ESO supports 1st- and 2nd-order plants only");

    template<size_t, typename, typename>
    friend class ESO;

    [[no_unique_address]] typename Store::Phi   Phi_{};                  ///< ZOH state map
    [[no_unique_address]] typename Store::Gamma Gamma_{};                ///< ZOH input map
    [[no_unique_address]] typename Store::L     L_{};                    ///< Current-estimator gain
    ColVec<NX + 1, T>                           z_{};                    ///< ESO state: [ŷ, ŷ̇, …, f̂]
    T                                           u_prev_{T{0}};           ///< Last applied command (post-saturation)
    bool                                        valid_{Store::constant}; ///< A successful discrete (Ts > 0) design

public:
    constexpr ESO() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!std::is_same_v<U, T> && !Store::constant)
    constexpr ESO(const design::ESOResult<NX, U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : ESO(result.template as<T>()) {}

    /**
     * @brief From the design result, which carries the sample time.
     *
     * A failed or continuous (Ts = 0) design is ill-formed under constant evaluation; at run
     * time it leaves the observer invalid (valid() false, state held).
     */
    constexpr ESO(const design::ESOResult<NX, T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        requires(!Store::constant)
        : Phi_(result.Phi),
          Gamma_(result.Gamma),
          L_(result.L),
          valid_(result.success && result.Ts > T{0}) {
        if ((!result.success || !(result.Ts > T{0})) && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_sample_time : result.status); // constant-evaluated (constinit / constexpr): a failed or continuous design does not compile
        }
    }

    template<typename U>
        requires(!Store::constant)
    constexpr ESO(const ESO<NX, U>& other)
        : Phi_(other.Phi_.template as<T>()),
          Gamma_(other.Gamma_.template as<T>()),
          L_(other.L_.template as<T>()),
          z_(other.z_.template as<T>()),
          u_prev_(static_cast<T>(other.u_prev_)),
          valid_(other.valid_) {}

    /** @brief Time update with the input the plant holds over the period: ẑ ← Φ ẑ + Γ u. */
    constexpr void predict(const ColVec<1, T>& u) {
        if (valid_) {
            u_prev_ = u[0];
            predict();
        }
    }

    /** @brief Time update with the last applied input (back_calculate or the last predict's u). */
    constexpr void predict() {
        if (valid_) {
            z_ = ColVec<NX + 1, T>((Phi_ * z_) + (Gamma_ * ColVec<1, T>(u_prev_)));
        }
    }

    /** @brief Measurement update: ẑ ← ẑ + L (y − ẑ₀). Does not advance time. */
    constexpr void update(const ColVec<1, T>& y) {
        if (valid_) {
            const T innovation = y[0] - z_[0];
            z_ = ColVec<NX + 1, T>(z_ + (L_ * ColVec<1, T>(innovation)));
        }
    }

    /**
     * @brief Fused tick: predict with @p u, correct with @p y, return [ŷ, (ẏ,) f̂]
     *
     * @param y Measurement taken now
     * @param u Input applied over the previous interval
     */
    constexpr const ColVec<NX + 1, T>& step(const ColVec<1, T>& y, const ColVec<1, T>& u) {
        predict(u);
        update(y);
        return z_;
    }

    /** @brief Fused tick with the last applied input (back_calculate or the last step's u). */
    constexpr const ColVec<NX + 1, T>& step(T y) {
        predict();
        update(y);
        return z_;
    }

    /**
     * @brief Record the input the plant received, for the next predict() / step(y)
     *
     * Same signature as the laws' back_calculate. A feedforward added downstream is already
     * part of @p u_sat, so @p u_ff does not enter.
     */
    constexpr void back_calculate(T u_sat, [[maybe_unused]] T u_ff = T{0}) { u_prev_ = u_sat; }

    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /** @brief ESO state [ŷ, (ẏ,) f̂]. */
    [[nodiscard]] constexpr const ColVec<NX + 1, T>& state() const { return z_; }

    /** @brief Estimated output ŷ = z₀. */
    [[nodiscard]] constexpr T output() const { return z_[0]; }

    /** @brief Lumped-disturbance estimate f̂ = z_{NX}. */
    [[nodiscard]] constexpr T disturbance() const { return z_[NX]; }

    /** @brief Command the ESO used last tick (post-back_calculate if you called it). */
    [[nodiscard]] constexpr T last_command() const { return u_prev_; }

    /** @brief Zero the state and the held input. */
    constexpr void reset() { reset(ColVec<NX + 1, T>{}); }

    /** @brief Start from @p z0 with no held input: e.g. [y₀, 0, 0] avoids the startup transient from ẑ = 0. */
    constexpr void reset(const ColVec<NX + 1, T>& z0) {
        z_ = z0;
        u_prev_ = T{0};
    }
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<size_t NX, typename T>
ESO(const design::ESOResult<NX, T>&) -> ESO<NX, T>;

namespace detail {

/** @brief The runtime an ESO design runs as: none for anything else. */
template<const auto& Res, typename T, typename R = std::remove_cvref_t<decltype(Res)>>
struct eso_runtime {};

/** @brief design::eso: with members, or on its constants. */
template<const auto& Res, typename T, size_t NX, typename U>
struct eso_runtime<Res, T, design::ESOResult<NX, U>> {
    using type = ESO<NX, T>;                              ///< With members
    using on_constants = ESO<NX, T, ConstantEso<Res, T>>; ///< On the design's constants
};

/** @brief An ESO design. */
template<const auto& Res, typename T>
concept eso_design = requires { typename eso_runtime<Res, T>::type; };

} // namespace detail

/** @brief Law<kRes, T> of a design::eso initialized at start-up: the ESO it runs as. A runtime of the same shape assigned replaces it. */
template<const auto& Res, typename T>
    requires detail::eso_design<Res, T>
struct Law<Res, T> : detail::eso_runtime<Res, T>::type {
    using Runtime = detail::eso_runtime<Res, T>::type; ///< The runtime the design runs as

    /** @brief The design's observer. */
    constexpr Law() : Runtime{Res} {}

    using Runtime::operator=;
};

/** @brief Law<kRes, T> of a `constexpr` design::eso: the same ESO on the design's matrices as constants of the type. */
template<const auto& Res, typename T>
    requires(detail::eso_design<Res, T> && detail::constant_design<Res>)
struct Law<Res, T> : detail::eso_runtime<Res, T>::on_constants {
    using Runtime = detail::eso_runtime<Res, T>::type; ///< The runtime with members: the type of a member a redesign assigns

    static_assert(detail::constant_design_runs(Res.success, Res.status, Res.Ts));
};

} // namespace damp
