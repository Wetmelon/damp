// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file model_result.hpp
 * @brief Results of designs whose product is a model: a transfer function or a state-space system
 *
 * Filter prototypes (design::lowpass_1st(fc), design::butterworth_lowpass, …) and model
 * builders such as design::sogi return the model with the same success / status pair the
 * gain designs carry, so an invalid input is refused instead of returning a placeholder model.
 */

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

namespace damp {
namespace design {

/**
 * @struct TransferFunctionResult
 * @brief A designed transfer function with its success / status
 *
 * @tparam Nnum Numerator coefficient count
 * @tparam Nden Denominator coefficient count
 * @tparam T    Scalar type
 */
template<size_t Nnum, size_t Nden, typename T = double>
struct [[nodiscard]] TransferFunctionResult {
    TransferFunction<Nnum, Nden, T> tf{};                                ///< The designed transfer function (zero when refused)
    bool                            success{false};                      ///< status == DesignStatus::ok
    DesignStatus                    status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr TransferFunctionResult<Nnum, Nden, U> as() const {
        return TransferFunctionResult<Nnum, Nden, U>{tf.template as<U>(), success, status};
    }
};

/**
 * @struct StateSpaceResult
 * @brief A designed state-space model with its success / status
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam NY Number of outputs
 * @tparam T  Scalar type
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
struct [[nodiscard]] StateSpaceResult {
    StateSpace<NX, NU, NY, T> sys{};                               ///< The designed model (zero when refused); sys.Ts is its domain
    bool                      success{false};                      ///< status == DesignStatus::ok
    DesignStatus              status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr StateSpaceResult<NX, NU, NY, U> as() const {
        return StateSpaceResult<NX, NU, NY, U>{sys.template as<U>(), success, status};
    }
};

namespace detail {

/** @brief A successful TransferFunctionResult holding @p tf. */
template<size_t Nnum, size_t Nden, typename T>
[[nodiscard]] constexpr TransferFunctionResult<Nnum, Nden, T> designed(const TransferFunction<Nnum, Nden, T>& tf) {
    return TransferFunctionResult<Nnum, Nden, T>{tf, true, DesignStatus::ok};
}

/** @brief A successful StateSpaceResult holding @p sys. */
template<size_t NX, size_t NU, size_t NY, typename T>
[[nodiscard]] constexpr StateSpaceResult<NX, NU, NY, T> designed(const StateSpace<NX, NU, NY, T>& sys) {
    return StateSpaceResult<NX, NU, NY, T>{sys, true, DesignStatus::ok};
}

} // namespace detail

} // namespace design
} // namespace damp
