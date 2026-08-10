// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file simulate.hpp
 * @brief Closed-loop simulation helpers (host-only)
 *
 * @defgroup simulate Closed-Loop Simulation
 * @brief Host closed-loop simulation of plants with continuous or sampled controllers
 *
 * Host-only. Result storage uses @c std::vector; this header is not part of
 * the embeddable @c control.hpp surface. Pull via @c workbench.hpp on the
 * workstation — never into target firmware.
 *
 * Provides functions to simulate closed-loop systems using the ODE solvers
 * from solver.hpp and integrators from integrator.hpp.
 *
 * Ts consistency.
 * - Continuous LTI plants: @c simulate_lti / @c simulate_lti_siso expect
 *   @c sys.Ts == 0; the plant is integrated in continuous time while the
 *   controller ZOH period is the separate @p Ts argument of @c simulate_sampled.
 * - Discrete plants: @c simulate_discrete advances with @c sys.Ts (time labels
 *   only — the map is one step per call). Pass a discrete model whose @c Ts
 *   matches the control period you intend to validate.
 * - Multi-rate: @c Ts is the controller sample period; @p fine_solver's step
 *   is the plant substep and must be ≤ @c Ts.
 *
 * Usage:
 * @code
 *   // Nonlinear plant: dx/dt = f(t, x, u)
 *   auto plant = [](double t, auto x, auto u) { return ...; };
 *   auto output = [](auto x) { return ...; };
 *   auto controller = [&](auto y) { return -K * y; };
 *
 *   auto result = simulate(plant, output, controller, solver, x0, {0.0, 10.0});
 *
 *   // SISO PID (or any control(r,y) SisoController) with a constant reference:
 *   auto pid = PIController<double>{design::pid(2.0, 5.0, 0.0).discretize(Ts)};
 *   auto sil = simulate_lti_siso(plant_ss, pid, 1.0, fine_solver, Ts, x0, {0.0, 5.0});
 * @endcode
 */

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include "damp/backend.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/systems/state_space.hpp"

namespace damp::sim {

/**
 * @brief Adapt a SISO `control(r, y)` controller for callables that want
 *        `u = f(y)` or `u = f(t, y)` (the older simulate* contracts).
 *
 * Holds a constant reference @p r. Models the usual digital-control SIL setup:
 * the same fixed-rate object that runs on target (`control(r, y)`) is exercised
 * on the host. Sample time is already baked into discrete controllers; @p Ts is
 * retained only as the plant ZOH period for the multi-rate simulator.
 *
 * @tparam C Controller type with `control(T, T) -> T` (e.g. PIDController)
 * @tparam T Scalar type
 */
template<typename C, typename T = double>
struct SisoReferenceAdapter {
    C* c{};
    T  r{T{0}};
    T  Ts{T{0}}; ///< Controller sample period (ZOH interval); not passed to control()

    [[nodiscard]] ColVec<1, T> operator()(const ColVec<1, T>& y) const {
        return ColVec<1, T>{c->control(r, y[0])};
    }

    [[nodiscard]] ColVec<1, T> operator()(T /*t*/, const ColVec<1, T>& y) const {
        return (*this)(y);
    }
};

/**
 * @brief Build a @ref SisoReferenceAdapter for a SISO `control(r,y)` controller
 */
template<typename C, typename T>
[[nodiscard]] constexpr SisoReferenceAdapter<std::remove_reference_t<C>, T>
siso_ref(C& c, T r, T Ts) {
    return SisoReferenceAdapter<std::remove_reference_t<C>, T>{&c, r, Ts};
}

/**
 * @brief Result of a closed-loop simulation
 *
 * @tparam NX Number of plant states
 * @tparam NU Number of control inputs
 * @tparam NY Number of plant outputs
 * @tparam T  Scalar type
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
struct SimulationResult {
    std::vector<T>             t; ///< Time points
    std::vector<ColVec<NX, T>> x; ///< State history
    std::vector<ColVec<NY, T>> y; ///< Output history
    std::vector<ColVec<NU, T>> u; ///< Control input history
};

/**
 * @brief Simulate a nonlinear plant with a controller in closed loop
 *
 * The plant dynamics are dx/dt = plant(t, x, u), the output is y = output(x),
 * and the controller computes u = controller(y). At each solver step, the
 * controller output is computed from the current output, then held constant
 * for the duration of the integration step (zero-order hold).
 *
 * @tparam Plant      Callable: (T t, ColVec<NX,T> x, ColVec<NU,T> u) -> ColVec<NX,T>
 * @tparam Output     Callable: (ColVec<NX,T> x) -> ColVec<NY,T>
 * @tparam Controller Callable: (ColVec<NY,T> y) -> ColVec<NU,T>
 * @tparam Solver     FixedStepSolver or AdaptiveStepSolver
 *
 * @param plant      Plant dynamics: dx/dt = plant(t, x, u)
 * @param output     Output function: y = output(x)
 * @param controller Controller: u = controller(y)
 * @param solver     ODE solver instance
 * @param x0         Initial state
 * @param t_span     Simulation time span {t0, tf}
 * @return SimulationResult with time, state, output, and input history
 */
template<size_t NX, size_t NU, size_t NY, typename T, typename Plant, typename Output, typename Controller, typename Solver>
SimulationResult<NX, NU, NY, T> simulate(
    Plant&&                 plant,
    Output&&                output,
    Controller&&            controller,
    const Solver&           solver,
    const ColVec<NX, T>&    x0,
    const damp::pair<T, T>& t_span
) {
    // Capture controller output so the plant closure can use it
    ColVec<NU, T> u_current{};

    // Wrap plant + controller into f(t, x) for the solver
    auto f = [&](T t, const ColVec<NX, T>& x) -> ColVec<NX, T> {
        return plant(t, x, u_current);
    };

    // Solve, but intercept each step to record outputs and update controller
    SimulationResult<NX, NU, NY, T> sim;

    // Compute initial output and control
    ColVec<NY, T> y0 = output(x0);
    u_current = controller(y0);

    // Record initial conditions
    sim.t.push_back(t_span.first);
    sim.x.push_back(x0);
    sim.y.push_back(y0);
    sim.u.push_back(u_current);

    // Create a copy of the solver so we can add a step callback
    auto sim_solver = solver;
    sim_solver.set_on_step([&](T t, const ColVec<NX, T>& x) {
        ColVec<NY, T> y = output(x);
        u_current = controller(y);

        sim.t.push_back(t);
        sim.x.push_back(x);
        sim.y.push_back(y);
        sim.u.push_back(u_current);
    });

    sim_solver.solve(f, x0, t_span);

    return sim;
}

/**
 * @brief Simulate a continuous plant under a discrete (sampled) controller — multi-rate.
 *
 * The closed loop runs at two rates: the controller is evaluated once per control
 * period @p Ts and its output held (zero-order hold) while @p fine_solver
 * integrates the continuous plant across that period, recording every sub-step. So the
 * trajectory keeps full resolution between control updates (smooth, accurate continuous
 * response) even though the controller only ticks at @p Ts. This is the digital-control
 * pattern — a current loop or servo updating at @p Ts while the physical plant evolves
 * continuously — and it decouples the integration rate from the control rate, unlike
 * @ref simulate where the controller fires on every solver step.
 *
 * @tparam Plant      (T t, ColVec\<NX\> x, ColVec\<NU\> u) -> ColVec\<NX\>   (dx/dt)
 * @tparam Output     (ColVec\<NX\> x) -> ColVec\<NY\>
 * @tparam Controller (T t, ColVec\<NY\> y) -> ColVec\<NU\>   (held across each period)
 * @tparam Solver     a FixedStepSolver whose step is the (fine) integration step
 *
 * @param plant        Continuous plant f(t, x, u) -> dx/dt
 * @param output       Measurement map y = g(x)
 * @param controller   Sampled law u = c(t, y), held for each control period
 * @param fine_solver  Fixed-step continuous solver (sub-step within Ts)
 * @param Ts           Control period [s]
 * @param x0           Initial plant state
 * @param t_span       Integration interval {t0, tf}
 * @param log_stride   Record every N-th fine sample within a control period (default 1)
 */
template<size_t NX, size_t NU, size_t NY, typename T, typename Plant, typename Output, typename Controller, typename Solver>
SimulationResult<NX, NU, NY, T> simulate_sampled(
    Plant&&                 plant,
    Output&&                output,
    Controller&&            controller,
    const Solver&           fine_solver,
    T                       Ts,
    const ColVec<NX, T>&    x0,
    const damp::pair<T, T>& t_span,
    std::size_t             log_stride = 1
) {
    SimulationResult<NX, NU, NY, T> sim;
    ColVec<NX, T>                   x = x0;
    T                               t = t_span.first;
    const std::size_t               stride = (log_stride == 0) ? 1 : log_stride;

    while (t < t_span.second - T(1e-12)) {
        const T             t_next = (t + Ts < t_span.second) ? t + Ts : t_span.second;
        const ColVec<NY, T> y = output(x);
        const ColVec<NU, T> u = controller(t, y); // zero-order hold across this period

        if (sim.t.empty()) { // record the initial sample once
            sim.t.push_back(t);
            sim.x.push_back(x);
            sim.y.push_back(y);
            sim.u.push_back(u);
        }

        const auto f = [&](T tt, const ColVec<NX, T>& xx) { return plant(tt, xx, u); };
        const auto seg = fine_solver.solve(f, x, {t, t_next});
        for (size_t j = 1; j < seg.size(); ++j) { // skip seg[0] (== current sample)
            const bool keep = (j % stride == 0) || (j + 1 == seg.size());
            if (!keep) {
                continue;
            }
            sim.t.push_back(seg.t[j]);
            sim.x.push_back(seg.x[j]);
            sim.y.push_back(output(seg.x[j]));
            sim.u.push_back(u);
        }
        x = seg.x.back();
        t = t_next;
    }
    return sim;
}

/**
 * @brief Simulate a nonlinear plant with state-feedback controller
 *
 * Convenience overload where the controller receives the full state x
 * rather than the output y. Output is computed for recording only.
 *
 * @param plant      dx/dt = plant(t, x, u)
 * @param output     y = output(x) — for recording only
 * @param controller u = controller(x) — state feedback
 * @param solver     ODE solver instance
 * @param x0         Initial state
 * @param t_span     {t0, tf}
 */
template<size_t NX, size_t NU, size_t NY, typename T, typename Plant, typename Output, typename Controller, typename Solver>
SimulationResult<NX, NU, NY, T> simulate_state_feedback(
    Plant&&                 plant,
    Output&&                output,
    Controller&&            controller,
    const Solver&           solver,
    const ColVec<NX, T>&    x0,
    const damp::pair<T, T>& t_span
) {
    ColVec<NU, T> u_current{};

    auto f = [&](T t, const ColVec<NX, T>& x) -> ColVec<NX, T> {
        return plant(t, x, u_current);
    };

    SimulationResult<NX, NU, NY, T> sim;

    ColVec<NY, T> y0 = output(x0);
    u_current = controller(x0);

    sim.t.push_back(t_span.first);
    sim.x.push_back(x0);
    sim.y.push_back(y0);
    sim.u.push_back(u_current);

    auto sim_solver = solver;
    sim_solver.set_on_step([&](T t, const ColVec<NX, T>& x) {
        ColVec<NY, T> y = output(x);
        u_current = controller(x);

        sim.t.push_back(t);
        sim.x.push_back(x);
        sim.y.push_back(y);
        sim.u.push_back(u_current);
    });

    sim_solver.solve(f, x0, t_span);

    return sim;
}

/**
 * @brief Simulate a continuous LTI system with a controller
 *
 * Uses StateSpace matrices (A, B, C, D) as the plant.
 * Plant output is y = Cx + Du with the ZOH control that was applied over the
 * step just completed; the controller then updates u for the next interval
 * (one fixed-point update closes direct feedthrough when D ≠ 0).
 *
 * @param sys        Continuous-time state-space system (Ts == 0)
 * @param controller u = controller(y)
 * @param solver     ODE solver
 * @param x0         Initial state
 * @param t_span     {t0, tf}
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Controller, typename Solver>
SimulationResult<NX, NU, NY, T> simulate_lti(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    Controller&&                             controller,
    const Solver&                            solver,
    const ColVec<NX, T>&                     x0,
    const damp::pair<T, T>&                  t_span
) {
    // Specialized path so y can include D·u (the generic simulate() output map
    // is y = f(x) only).
    ColVec<NU, T> u_current{};

    auto f = [&](T /*t*/, const ColVec<NX, T>& x) -> ColVec<NX, T> {
        return sys.A * x + sys.B * u_current;
    };

    SimulationResult<NX, NU, NY, T> sim;

    ColVec<NY, T> y0 = sys.C * x0 + sys.D * u_current;
    u_current = controller(y0);
    y0 = sys.C * x0 + sys.D * u_current;

    sim.t.push_back(t_span.first);
    sim.x.push_back(x0);
    sim.y.push_back(y0);
    sim.u.push_back(u_current);

    auto sim_solver = solver;
    sim_solver.set_on_step([&](T t, const ColVec<NX, T>& x) {
        ColVec<NY, T> y = sys.C * x + sys.D * u_current;
        u_current = controller(y);
        y = sys.C * x + sys.D * u_current;

        sim.t.push_back(t);
        sim.x.push_back(x);
        sim.y.push_back(y);
        sim.u.push_back(u_current);
    });

    sim_solver.solve(f, x0, t_span);
    return sim;
}

/**
 * @brief SIL: continuous SISO LTI plant under a sampled SISO controller with constant reference
 *
 * Multi-rate digital-control pattern: controller ticks every @p Ts via
 * `control(r, y)` (ZOH on u); @p fine_solver integrates the plant between
 * ticks. Same controller object as on target — Design Is Deploy / gate-and-validate.
 *
 * @param sys          Continuous plant (SISO: NU = NY = 1)
 * @param controller   Object with `control(T r, T y) -> T` (e.g. PIDController after discretize)
 * @param r            Constant reference
 * @param fine_solver  Fixed-step solver for plant substeps
 * @param Ts           Controller sample period [s] (ZOH interval)
 * @param x0           Initial state
 * @param t_span       {t0, tf}
 */
template<size_t NX, size_t NW, size_t NV, typename T, typename Controller, typename Solver>
[[nodiscard]] SimulationResult<NX, 1, 1, T> simulate_lti_siso(
    const StateSpace<NX, 1, 1, T, NW, NV>& sys,
    Controller&                            controller,
    T                                      r,
    const Solver&                          fine_solver,
    T                                      Ts,
    const ColVec<NX, T>&                   x0,
    const damp::pair<T, T>&                t_span
)
    requires requires(Controller& c, T a, T b) {
        { c.control(a, b) } -> std::convertible_to<T>;
    }
{
    auto plant = [&](T /*t*/, const ColVec<NX, T>& x, const ColVec<1, T>& u) -> ColVec<NX, T> {
        return sys.A * x + sys.B * u;
    };
    auto output = [&](const ColVec<NX, T>& x) -> ColVec<1, T> {
        return ColVec<1, T>{sys.C * x};
    };
    auto adapted = siso_ref(controller, r, Ts);
    return simulate_sampled<NX, 1, 1, T>(plant, output, adapted, fine_solver, Ts, x0, t_span);
}

/**
 * @brief Simulate a discrete-time nonlinear plant with a controller
 *
 * Steps a user-supplied discrete map x[k+1] = f(k, x[k], u[k]) directly — no ODE
 * solver, no continuous-time integration. This is the path for plants that are
 * natively discrete and nonlinear (a switching-model converter, a sampled-data
 * map, a difference equation identified from data), where the continuous
 * `simulate()` overloads don't apply.
 *
 * Timing matches `simulate_discrete()`: at each step the output is computed from
 * the current state, the controller produces u from that output, then the plant
 * map advances the state (zero-order hold on u across the step).
 *
 * @tparam NX Number of plant states
 * @tparam NU Number of control inputs
 * @tparam NY Number of plant outputs
 * @tparam Plant      Callable: (size_t k, ColVec<NX,T> x, ColVec<NU,T> u) -> ColVec<NX,T>  (returns x[k+1])
 * @tparam Output     Callable: (ColVec<NX,T> x) -> ColVec<NY,T>
 * @tparam Controller Callable: (ColVec<NY,T> y) -> ColVec<NU,T>
 *
 * @param plant      Discrete dynamics: x[k+1] = plant(k, x, u)
 * @param output     Output map: y = output(x)
 * @param controller Controller: u = controller(y)
 * @param x0         Initial state
 * @param Ts         Sample time [s] — used only to populate the time vector
 * @param n_steps    Number of steps to simulate
 * @return SimulationResult with time, state, output, and input history (n_steps + 1 samples)
 */
template<size_t NX, size_t NU, size_t NY, typename T, typename Plant, typename Output, typename Controller>
SimulationResult<NX, NU, NY, T> simulate_discrete_nonlinear(
    Plant&&              plant,
    Output&&             output,
    Controller&&         controller,
    const ColVec<NX, T>& x0,
    T                    Ts,
    size_t               n_steps
) {
    SimulationResult<NX, NU, NY, T> sim;
    sim.t.reserve(n_steps + 1);
    sim.x.reserve(n_steps + 1);
    sim.y.reserve(n_steps + 1);
    sim.u.reserve(n_steps + 1);

    ColVec<NX, T> x = x0;
    T             t = T(0);

    for (size_t k = 0; k <= n_steps; ++k) {
        ColVec<NY, T> y = output(x);
        ColVec<NU, T> u = controller(y);

        sim.t.push_back(t);
        sim.x.push_back(x);
        sim.y.push_back(y);
        sim.u.push_back(u);

        if (k < n_steps) {
            x = plant(k, x, u); // x[k+1] = f(k, x[k], u[k])
            t += Ts;
        }
    }

    return sim;
}

/**
 * @brief Simulate a discrete-time system with a controller
 *
 * Steps the discrete system x_{k+1} = Ax_k + Bu_k with u_k = controller(y_k)
 * and plant output y_k = Cx_k + Du_k. When D ≠ 0 the controller is evaluated
 * once with the held previous u (zero at k = 0), then y is recomputed with the
 * new u so recorded outputs include direct feedthrough.
 *
 * @param sys        Discrete-time state-space system (Ts > 0)
 * @param controller u = controller(y)
 * @param x0         Initial state
 * @param n_steps    Number of time steps to simulate
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Controller>
SimulationResult<NX, NU, NY, T> simulate_discrete(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    Controller&&                             controller,
    const ColVec<NX, T>&                     x0,
    size_t                                   n_steps
) {
    SimulationResult<NX, NU, NY, T> sim;
    sim.t.reserve(n_steps + 1);
    sim.x.reserve(n_steps + 1);
    sim.y.reserve(n_steps + 1);
    sim.u.reserve(n_steps + 1);

    ColVec<NX, T> x = x0;
    ColVec<NU, T> u{};
    T             t = T(0);

    for (size_t k = 0; k <= n_steps; ++k) {
        ColVec<NY, T> y = sys.C * x + sys.D * u;
        u = controller(y);
        y = sys.C * x + sys.D * u;

        sim.t.push_back(t);
        sim.x.push_back(x);
        sim.y.push_back(y);
        sim.u.push_back(u);

        if (k < n_steps) {
            x = sys.A * x + sys.B * u;
            t += sys.Ts;
        }
    }

    return sim;
}

} // namespace damp::sim
