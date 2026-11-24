// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private: a public header that exports it provides its names

/**
 * @file lq_law.hpp
 * @brief Law<kRes, T> of an LQ design: the runtime the design names, on members or on the
 *        design's matrices as constants of the type
 *
 * Each runtime's header (lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp) specializes
 * detail::lq_runtime for the designs it runs.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/state_feedback_core.hpp"
#include "damp/controllers/law.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/const_matrix.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {

namespace detail {

/** @brief The runtime an LQ design runs as, in the scalar @p T (specialized by lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp). */
template<typename R, typename T>
struct lq_runtime {};

/** @brief An LQ design: a result that has a runtime in lq_runtime. */
template<const auto& Res, typename T>
concept lq_design = requires {
    typename lq_runtime<damp::remove_cvref_t<decltype(Res)>, T>::type;
    Res.success;
};

} // namespace detail

/** @brief Law<kRes, T> of an LQ design initialized at start-up: the runtime it runs as (LQR, LQI, LQG, LQGI, OFLQG), default-constructed from the design. A runtime of the same shape assigned replaces it. */
template<const auto& Res, typename T>
    requires detail::lq_design<Res, T>
struct Law<Res, T> : detail::lq_runtime<damp::remove_cvref_t<decltype(Res)>, T>::type {
    using Runtime = detail::lq_runtime<damp::remove_cvref_t<decltype(Res)>, T>::type; ///< The runtime the design runs as

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res.template as<T>()} {}

    using Runtime::operator=;
};

namespace detail {

/** @brief A design a Law can run on constants: it succeeded and is discrete. Anything else does not compile, and the error names its status. */
template<const auto& Res>
inline constexpr bool constant_design_usable = [] {
    const bool ok = Res.success && discrete_gain(Res);
    if (!ok) {
        if constexpr (requires { Res.status; }) {
            refuse_failed_design(Res.success ? design::DesignStatus::invalid_sample_time : Res.status);
        } else {
            refuse_failed_design();
        }
    }
    return ok;
}();

/** @brief The law's matrices as constants of the type: the store of a `constexpr` design's Law. */
/** @brief A design's demand limit as constants of the type: what a `constexpr` design's Law holds, in no storage. */
template<const auto& D>
struct ConstDemand {
    static constexpr size_t state = D.state; ///< The state asked for (none: off)

    /** @brief Whether term @p j asks. */
    [[nodiscard]] static constexpr bool outer(size_t j) {
        return D.outer(j);
    }

    static constexpr auto lim = D.lim; ///< The most they may ask

    /** @brief The limit as a runtime that holds its matrices holds it: what limits() reports. */
    [[nodiscard]] static constexpr auto terms() {
        return D;
    }
};

template<const auto& Res, typename T>
struct ConstantGains {
    static constexpr bool constant = true;               ///< The matrices are constants of the type
    static constexpr auto design = Res.template as<T>(); ///< The design, in the runtime's scalar
    static constexpr auto kK = design.K;                 ///< K
    static constexpr auto kCi = [] {
        if constexpr (requires { design.integral.Ci; }) {
            if constexpr (damp::remove_cvref_t<decltype(design)>::kEstimated) {
                return design.integral.Ci; // on the measurement y
            } else {
                return design.integral.Ci_x; // on the measured state: Ci·C
            }
        } else {
            return Matrix<0, 0, T>{};
        }
    }(); ///< The integrated outputs
    static constexpr auto kCiu = [] {
        if constexpr (requires { design.integral.Ci_u; }) {
            return design.integral.Ci_u;
        } else {
            return Matrix<0, 0, T>{};
        }
    }(); ///< Their feedthrough
    static constexpr auto kTargetD = [] {
        if constexpr (requires { design.target.target_d; }) {
            return design.target.target_d;
        } else {
            return Matrix<0, 0, T>{};
        }
    }(); ///< The target per unit d̂
    static constexpr auto kTargetR = [] {
        if constexpr (requires { design.target.target_r; }) {
            return design.target.target_r;
        } else {
            return Matrix<0, 0, T>{};
        }
    }(); ///< The target per unit r
    static constexpr auto kDemand = [] {
        if constexpr (requires { design.limits; }) {
            return design.limits.demand;
        } else {
            return ::damp::detail::DemandTerms<damp::remove_cvref_t<decltype(design.K)>::cols(), T>{};
        }
    }(); ///< The limit on what the outer terms ask of a state (off without one)
    static constexpr T kFadeDecay = [] {
        if constexpr (requires { requires damp::remove_cvref_t<decltype(design)>::kHandover; }) {
            return damp::exp(-design.Ts / design.fade_time);
        } else {
            return T{0};
        }
    }(); ///< The handover fade's factor per tick, exp(−Ts/fade_time) (0 without a Handover)

    template<size_t, size_t, typename>
    using K = ConstMatrix<kK>; ///< The gain
    template<size_t, size_t, typename>
    using Ci = ConstMatrix<kCi>; ///< The integrated outputs
    template<size_t, size_t, typename>
    using Ci_u = ConstMatrix<kCiu>; ///< Their feedthrough
    template<size_t, size_t, typename>
    using TargetD = ConstMatrix<kTargetD>; ///< The target per unit d̂
    template<size_t, size_t, typename>
    using TargetR = ConstMatrix<kTargetR>; ///< The target per unit r
    template<size_t, typename>
    using Demand = ConstDemand<kDemand>; ///< The limit on what the outer terms ask of a state
};

/** @brief The estimator of a `constexpr` design's Law: none for a measured state. */
template<const auto& Res, typename T, typename R = damp::remove_cvref_t<decltype(Res)>>
struct constant_estimator {
    using type = NoEstimator; ///< NoEstimator
};

/** @brief A design with a filter runs its FixedGainObserver on constants of the type: NUK inputs, u and then a measured disturbance's. */
template<const auto& Res, typename T, size_t NX, size_t NU, typename U, size_t NI, size_t NY, size_t NXK, size_t NUK, size_t NW, size_t NV, size_t ND, size_t NZ, size_t NRU, size_t NM, size_t NDU, bool HO>
struct constant_estimator<Res, T, design::StateFeedbackResult<NX, NU, U, NI, NY, design::KalmanResult<NXK, NUK, NY, U, NW, NV, ND>, NZ, NRU, NM, NDU, HO>> {
    using type = FixedGainObserver<NXK, NUK, NY, T, ND, ConstantObserver<Res, T, NXK - ND>>; ///< The observer and the disturbance
};

} // namespace detail

/** @brief Law<kRes, T> of a `constexpr` LQ design: the same runtime on the design's matrices as constants of the type. */
template<const auto& Res, typename T>
    requires(detail::lq_design<Res, T> && detail::constant_design<Res>)
struct Law<Res, T> : detail::lq_runtime<damp::remove_cvref_t<decltype(Res)>, T>::template on_constants<typename detail::constant_estimator<Res, T>::type, detail::ConstantGains<Res, T>> {
    using Runtime = detail::lq_runtime<damp::remove_cvref_t<decltype(Res)>, T>::type; ///< The runtime with members: the type of a member a redesign assigns

    static_assert(detail::constant_design_usable<Res>);

    /**
     * @brief The design's law: its matrices and its demand limit are the type's, its integrator step
     *        and its other limits members (a limit set_limits refuses does not compile when
     *        constant-evaluated).
     */
    constexpr Law() {
        auto lim = this->limits();
        if constexpr (requires { this->integral.step; }) {
            this->integral.step = static_cast<T>(Res.integral.step);
            this->integral.unwind = Res.integral.unwind.template as<T>();
            this->tracking_time = static_cast<T>(Res.tracking_time);
            lim.integral = Res.limits.integral.template as<T>();
        }
        if constexpr (requires { Res.limits; }) {
            lim.input = Res.limits.input.template as<T>();
            lim.input_rate = Res.limits.input_rate.template as<T>();
        }
        if (!this->set_limits(lim) && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // a result whose limits were edited by hand
        }
    }

    /** @brief There is nothing to replace: the design is in the type. */
    template<typename Other>
        requires(!damp::is_same_v<damp::remove_cvref_t<Other>, Law>)
    constexpr Law& operator=(const Other& other) {
        static_assert(sizeof(Other) == 0, "This Law is a constexpr design's: its matrices are constants of the type, so there is nothing to replace. A member a redesign assigns is Law<kRes, T>::Runtime law{kRes.as<T>()};");
        static_cast<void>(other);
        return *this;
    }
};

} // namespace damp
