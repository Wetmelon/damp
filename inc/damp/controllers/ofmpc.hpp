// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file ofmpc.hpp
 * @brief The offset-free MPC: design::ofmpc, its MPCResult, the OFMPC runtime
 *
 * design::ofmpc<NP, NC>(plant, [Ts,] weights, constraints, KalmanNoise{.process, .measurement, .disturbance},
 * disturbance) → MPCResult → OFMPC. step(r, y) is predict, update(y), control(r). The
 * low-latency tick is update(y), control(r), back_calculate, predict. The QP, the
 * limits and the measured-state MPC are in mpc.hpp.
 */

#include "damp/backend.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @brief The offset-free MPC on a discrete plant: a Kalman filter on [x; d] (@p noise, @p disturbance)
 *        and the QP posed about the steady-state target under d̂: runs as OFMPC.
 *
 * As design::mpc, with the measured outputs C and D and the controlled outputs the rows the
 * disturbance spec lists (all of them without a list): see "Offset-free control" there.
 */
template<size_t NP, size_t NC, size_t NT_max = 0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Disturbance, size_t ND = detail::disturbance_traits<Disturbance, NU, NY, NW>::states, size_t NZ = detail::mpc_outputs<NU, NY, NW, Disturbance>>
    requires(detail::is_disturbance_model<Disturbance>)
[[nodiscard]] constexpr auto ofmpc(
    const StateSpace<NX, NU, NY, T, NW, NV>&                 sys,
    const MPCWeights<NU, NZ, T>&                             weights,
    const MPCConstraints<NU, NZ, T>&                         constraints,
    const damp::type_identity_t<KalmanNoise<NW, NV, T, ND>>& noise,
    const Disturbance&                                       disturbance
) {
    return detail::mpc_design<true, NP, NC, NT_max>(sys, T{0}, weights, constraints, noise, disturbance);
}

/** @brief The sampled-data offset-free MPC: the continuous plant and cost, the QP at @p Ts (the noise densities): runs as OFMPC. */
template<size_t NP, size_t NC, size_t NT_max = 0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Disturbance, size_t ND = detail::disturbance_traits<Disturbance, NU, NY, NW>::states, size_t NZ = detail::mpc_outputs<NU, NY, NW, Disturbance>>
    requires(detail::is_disturbance_model<Disturbance>)
[[nodiscard]] constexpr auto ofmpc(
    const StateSpace<NX, NU, NY, T, NW, NV>&                 plant,
    damp::type_identity_t<T>                                 Ts,
    const MPCWeights<NU, NZ, T>&                             weights,
    const MPCConstraints<NU, NZ, T>&                         constraints,
    const damp::type_identity_t<KalmanNoise<NW, NV, T, ND>>& noise,
    const Disturbance&                                       disturbance
) {
    return detail::mpc_design<false, NP, NC, NT_max>(plant, Ts, weights, constraints, noise, disturbance);
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Offset-free constrained linear MPC: a Kalman filter on [x; d] and the QP about the
 *        steady-state target under d̂ (design::ofmpc)
 *
 * The estimator is a FixedGainObserver with design::kalman's fixed gain, holding deviations
 * from the target. The tick is update(y) then control(r), or step(r, y) for both; @p r is the
 * controlled outputs' reference (held, or a Preview). back_calculate(u_sat) records what went
 * out, predict() moves the observer with it. A failed design does not compile when
 * constant-evaluated.
 *
 * @note Compare with MATLAB®'s mpcmove(mpcobj, xc, ym, r) with its default disturbance model.
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, typename T, typename Solver, size_t NT_max, typename Estimator>
class OFMPC : private detail::MPCCore<NX, NU, NY, NP, NC, 0, T, Solver, NT_max, Estimator> {
    using Core = detail::MPCCore<NX, NU, NY, NP, NC, 0, T, Solver, NT_max, Estimator>;

public:
    using Reference = Core::Reference;       ///< r: held, or a Preview
    using Preview = Core::Preview;           ///< The reference at steps 1 … NP
    using InputPreview = Core::InputPreview; ///< The input target per move

    static constexpr size_t NYM = Core::NYM; ///< Measured outputs
    static constexpr size_t NDE = Core::NDE; ///< Estimated disturbance states
    static constexpr size_t NZ1 = Core::NZ1; ///< QP decision size (moves + slack)
    static constexpr size_t NI1 = Core::NI1; ///< QP constraint rows
    static_assert(NDE > 0, "OFMPC runs design::ofmpc: the filter estimates a disturbance.");

    constexpr OFMPC() = default;

    /** @brief From design::ofmpc (its filter with any noise inputs). A failed design is a compile error when constant-evaluated. */
    template<typename Kalman>
    constexpr explicit OFMPC(const design::MPCResult<typename Core::Artifacts, Kalman, typename Core::TargetResult>& result) : Core(result) {}

    /** @brief From design::ofmpc with a specific solver instance (for stateful solver policies). */
    template<typename Kalman>
    constexpr OFMPC(const design::MPCResult<typename Core::Artifacts, Kalman, typename Core::TargetResult>& result, const Solver& solver) : Core(result, solver) {}

    /** @brief The filter's measurement update with @p y. */
    constexpr void update(const ColVec<NYM, T>& y) { Core::update(y); }

    /** @brief The move about the target for @p r under d̂: u = u_prev + Δu₀, clamped to the u and Δu boxes. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const Reference& r) { return Core::control(r, ColVec<0, T>{}); }

    /** @brief One tick: predict(), update(@p y), then control(@p r). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const Reference& r, const ColVec<NYM, T>& y) { return Core::step(r, ColVec<0, T>{}, y); }

    using Core::applied;
    using Core::back_calculate;
    using Core::disturbance_estimate;
    using Core::estimator;
    using Core::input_target;
    using Core::input_target_preview;
    using Core::last_command;
    using Core::last_cost;
    using Core::last_iterations;
    using Core::last_slack;
    using Core::last_status;
    using Core::limits;
    using Core::predict;
    using Core::qp_solver;
    using Core::reset;
    using Core::set_input_target;
    using Core::set_limits;
    using Core::state_estimate;
};

/** @brief OFMPC law{res}: the sizes, scalar, filter and target of design::ofmpc's result. */
template<size_t NXM, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, size_t NT, size_t NXK, size_t NYM, size_t NW, size_t NV, size_t NXE>
OFMPC(const design::MPCResult<design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT>, design::KalmanResult<NXK, NU, NYM, T, NW, NV, ND>, design::SteadyStateTargetResult<NXE, NU, ND, NY, T>>&)
    -> OFMPC<NXM, NU, NY, NP, NC, T, design::WarmStartActiveSetSolver<(NC * NU) + 1, (4 * NC * NU) + (2 * NP * NY) + NT + 1, T>, NT, FixedGainObserver<NXK, NU, NYM, T, ND>>;

/** @brief OFMPC law{res, solver}: as OFMPC law{res}, with that solver. */
template<size_t NXM, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, size_t NT, size_t NXK, size_t NYM, size_t NW, size_t NV, size_t NXE, typename Solver>
OFMPC(const design::MPCResult<design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT>, design::KalmanResult<NXK, NU, NYM, T, NW, NV, ND>, design::SteadyStateTargetResult<NXE, NU, ND, NY, T>>&, const Solver&) -> OFMPC<NXM, NU, NY, NP, NC, T, Solver, NT, FixedGainObserver<NXK, NU, NYM, T, ND>>;

namespace detail {

/** @brief design::ofmpc runs as OFMPC. */
template<size_t NXM, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename U, size_t NT, size_t NXK, size_t NYM, size_t NW, size_t NV, size_t NXE, typename T>
struct mpc_runtime<design::MPCResult<design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, U, NT>, design::KalmanResult<NXK, NU, NYM, U, NW, NV, ND>, design::SteadyStateTargetResult<NXE, NU, ND, NY, U>>, T> {
    using type = OFMPC<NXM, NU, NY, NP, NC, T, design::WarmStartActiveSetSolver<(NC * NU) + 1, (4 * NC * NU) + (2 * NP * NY) + NT + 1, T>, NT, FixedGainObserver<NXK, NU, NYM, T, ND>>; ///< OFMPC
    static constexpr bool filtered = true;                                                                                                                                              ///< Its filter can run on the design's constants
    template<const auto& Res>
    using on_constants = OFMPC<NXM, NU, NY, NP, NC, T, design::WarmStartActiveSetSolver<(NC * NU) + 1, (4 * NC * NU) + (2 * NP * NY) + NT + 1, T>, NT, FixedGainObserver<NXK, NU, NYM, T, ND, ConstantObserver<Res, T, NXK - ND>>>; ///< OFMPC with its filter on the design's constants
};

} // namespace detail

} // namespace damp
