// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file status.hpp
 * @brief Why a design failed: the status every design result carries
 *
 * A design result keeps `bool success` (so `static_assert(res.success)` reads plainly) and
 * adds `DesignStatus status`, which says what to change when it is false. Both are set in
 * one place, detail::set_status, so they cannot disagree: `success == (status == ok)`.
 *
 * ok means this design's own checks passed: finite inputs, weights, solver convergence, and,
 * where the design runs one, a closed-loop pole test (Hurwitz or Schur). It is not a
 * robustness certificate. The Riccati residual lives on RiccatiResidual, not as a separate
 * status: a small residual is necessary, not sufficient, for a usable gain.
 */

#include <cstdint>

#include "damp/backend.hpp"

namespace damp {
namespace design {

/**
 * @brief Why a design failed, one value per cause class (not per code branch)
 *
 * Estimator designs solve the dual Riccati equation, so a dual that cannot be stabilized is
 * reported as @c not_detectable: (A, C) is what to change.
 */
enum class DesignStatus : std::uint8_t {
    ok = 0,              ///< This design's checks passed (not a robustness certificate)
    not_stabilizable,    ///< (A, B) fails PBH, or the closed-loop poles fail this design's Hurwitz/Schur test
    not_detectable,      ///< (A, C) cannot be detected (estimator designs, through the dual)
    not_converged,       ///< Iteration cap reached, the Schur form did not isolate the stable subspace, or the result overflowed the range of T (a discretization whose unstable mode grows past it within one period)
    singular,            ///< Plant or numerics singular: rank-deficient B, R + BᵀSB, innovation covariance
    bad_weights,         ///< Design inputs inconsistent: Q / R not symmetric, Q not PSD, R not PD, cross term not PSD, a mode on the stability boundary Q does not weight, Q_integral not PD with integral action, an output index outside [0, NY), a model or noise entry that is not finite
    not_placeable,       ///< Pole placement cannot assign the requested spectrum
    invalid_sample_time, ///< Ts ≤ 0 or non-finite, Ts given for an already-discrete plant, or no design for this time domain
    invalid_parameter,   ///< A scalar design input outside its domain: a cutoff outside (0, Nyquist), a non-positive gain, time constant, damping or bandwidth
};

} // namespace design

namespace detail {

/** @brief Set a design result's status and success together: success == (status == ok). */
template<typename Result>
constexpr Result& set_status(Result& result, design::DesignStatus status) {
    result.status = status;
    result.success = (status == design::DesignStatus::ok);
    return result;
}

/** @brief A result of type @p Result refused with @p status. */
template<typename Result>
[[nodiscard]] constexpr Result refused(design::DesignStatus status) {
    Result result{};
    set_status(result, status);
    return result;
}

/** @brief The estimator's dual Riccati status: a dual that cannot be stabilized is an undetectable (A, C). */
[[nodiscard]] constexpr design::DesignStatus dual_status(design::DesignStatus status) {
    return status == design::DesignStatus::not_stabilizable ? design::DesignStatus::not_detectable : status;
}

/**
 * @brief One refusal per DesignStatus, so a refused constant-evaluated law names its cause
 *
 * None is constexpr: reached under constant evaluation (`static constinit LQR law{res};` on a
 * failed design), the call is the compile error, and the compiler's "call to non-constexpr
 * function" line carries the function's name, e.g. design_refused::not_stabilizable(). At run
 * time each traps, as refuse_failed_design() does.
 */
namespace design_refused {

/** @brief Refused: (A, B) cannot be stabilized, or the closed-loop certificate failed. */
[[noreturn]] inline void not_stabilizable() noexcept {
    refuse_failed_design();
}

/** @brief Refused: (A, C) cannot be detected. */
[[noreturn]] inline void not_detectable() noexcept {
    refuse_failed_design();
}

/** @brief Refused: an iteration cap was reached, or the Schur form did not isolate the stable subspace. */
[[noreturn]] inline void not_converged() noexcept {
    refuse_failed_design();
}

/** @brief Refused: the plant or the numerics are singular. */
[[noreturn]] inline void singular() noexcept {
    refuse_failed_design();
}

/** @brief Refused: the weights are inconsistent (Q / R not symmetric, Q not PSD, R not PD, …), or the model has a non-finite entry. */
[[noreturn]] inline void bad_weights() noexcept {
    refuse_failed_design();
}

/** @brief Refused: pole placement cannot assign the requested spectrum. */
[[noreturn]] inline void not_placeable() noexcept {
    refuse_failed_design();
}

/** @brief Refused: Ts ≤ 0 or non-finite, a Ts for a discrete plant, or a continuous design for a discrete runtime. */
[[noreturn]] inline void invalid_sample_time() noexcept {
    refuse_failed_design();
}

/** @brief Refused: a scalar design input outside its domain. */
[[noreturn]] inline void invalid_parameter() noexcept {
    refuse_failed_design();
}

} // namespace design_refused

/**
 * @brief Refuse a failed design under constant evaluation, naming @p status in the compile error
 *
 * Dispatches to the design_refused function of @p status. A status of ok (a refusal for a
 * reason the status does not carry) takes the generic refuse_failed_design().
 */
[[noreturn]] constexpr void refuse_failed_design(design::DesignStatus status) noexcept {
    switch (status) {
        case design::DesignStatus::not_stabilizable:
            design_refused::not_stabilizable();
        case design::DesignStatus::not_detectable:
            design_refused::not_detectable();
        case design::DesignStatus::not_converged:
            design_refused::not_converged();
        case design::DesignStatus::singular:
            design_refused::singular();
        case design::DesignStatus::bad_weights:
            design_refused::bad_weights();
        case design::DesignStatus::not_placeable:
            design_refused::not_placeable();
        case design::DesignStatus::invalid_sample_time:
            design_refused::invalid_sample_time();
        case design::DesignStatus::invalid_parameter:
            design_refused::invalid_parameter();
        case design::DesignStatus::ok:
        default:
            refuse_failed_design(); // inside the switch: a trailing call would make every path non-constant
    }
}

} // namespace detail
} // namespace damp
