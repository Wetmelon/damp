// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lq_law.hpp
 * @brief Law<kRes, T> of an LQ or gain design: the runtime the design names, on members or on the
 *        design's matrices as constants of the type
 *
 * Each runtime's header (lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp) specializes
 * detail::lq_runtime for the designs it runs.
 */

#include "damp/backend.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/controllers/state_feedback_core.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/matrix/const_matrix.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp {

namespace detail {

/** @brief The runtime an LQ design runs as, in the scalar @p T (specialized by lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp). */
template<typename R, typename T>
struct lq_runtime {};

/** @brief An LQ or gain design: a result that has a runtime in lq_runtime. */
template<const auto& Res, typename T>
concept lq_design = requires {
    typename lq_runtime<damp::remove_cvref_t<decltype(Res)>, T>::type;
    Res.success;
};

} // namespace detail

/** @brief Law<kRes, T> of an LQ or gain design initialized at start-up: the runtime it runs as (LQR, LQI, LQG, LQGI, OFLQG), default-constructed from the design. A runtime of the same shape assigned replaces it. */
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
};

/** @brief The estimator of a `constexpr` design's Law: none for a measured state. */
template<const auto& Res, typename T, typename R = damp::remove_cvref_t<decltype(Res)>>
struct constant_estimator {
    using type = NoEstimator; ///< NoEstimator
};

/** @brief A design with a filter runs its FixedGainObserver on constants of the type. */
template<const auto& Res, typename T, size_t NX, size_t NU, typename U, size_t NI, size_t NY, size_t NXK, size_t NW, size_t NV, size_t ND, size_t NZ, size_t NRU>
struct constant_estimator<Res, T, design::StateFeedbackResult<NX, NU, U, NI, NY, design::KalmanResult<NXK, NU, NY, U, NW, NV, ND>, NZ, NRU>> {
    using type = FixedGainObserver<NXK, NU, NY, T, ND, ConstantObserver<Res, T, NXK - ND>>; ///< The observer and the disturbance
};

} // namespace detail

/** @brief Law<kRes, T> of a `constexpr` LQ or gain design: the same runtime on the design's matrices as constants of the type. */
template<const auto& Res, typename T>
    requires(detail::lq_design<Res, T> && detail::constant_design<Res>)
struct Law<Res, T> : detail::lq_runtime<damp::remove_cvref_t<decltype(Res)>, T>::template on_constants<typename detail::constant_estimator<Res, T>::type, detail::ConstantGains<Res, T>> {
    using Runtime = detail::lq_runtime<damp::remove_cvref_t<decltype(Res)>, T>::type; ///< The runtime with members: the type of a member a redesign assigns

    static_assert(detail::constant_design_usable<Res>);

    /** @brief The design's law. */
    constexpr Law() {
        if constexpr (requires { this->integral.step; }) {
            this->integral.step = static_cast<T>(Res.integral.step);
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
