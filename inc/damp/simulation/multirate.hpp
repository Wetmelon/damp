// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file multirate.hpp
 * @brief Multi-rate closed-loop simulation harness (Phase 1 / roadmap #18)
 *
 * Host-only (`std::vector` result storage) — not embeddable.
 *
 * Real deployments are multi-rate: an ISR-level inner loop at kilohertz and
 * RTOS outer loops at hertz. A faithful SIL must reproduce the rate hierarchy
 * and its inter-rate effects — sample-and-hold across rate boundaries and the
 * one-sample transport delay an outer task sees on plant state. That is where
 * most "works in sim, oscillates on hardware" gaps come from.
 *
 * Model
 * -----
 * - @c Ts_base is the fastest discrete tick (base rate).
 * - Each discrete law fires every @c divisor base ticks (integer ≥ 1).
 * - The plant is continuous: a fine-step ODE solver advances it between base
 *   ticks with the held plant input (ZOH from the fastest law that drives @c u).
 * - All rate crossings are sample-and-hold: slower→faster and control→plant.
 *
 * Order at base tick @c k (time @c t = t0 + k·Ts_base):
 * 1. Sample plant output @c y = output(x).
 * 2. Run discrete laws that fire this tick (slow → fast for a cascade: outer
 *    updates the held reference, then inner updates the held plant input).
 * 3. Integrate the plant with held @c u over @c [t, t + Ts_base].
 *
 * Single-rate is the special case @c divisor = 1 for every law (and matches
 * simulate_sampled when the control law is the same).
 *
 * @code
 * // Two-rate cascade: outer every 10 base ticks, inner every tick
 * auto outer = [&](double, const ColVec<1>& y) {
 *     return ColVec<1>{outer_pi.control(r, y[0])};  // → i*
 * };
 * auto inner = [&](double, const ColVec<1>& r_i, const ColVec<1>& y) {
 *     return ColVec<1>{inner_pi.control(r_i[0], y[0])};  // → u
 * };
 * RK4<2> rk4;
 * FixedStepSolver plant_solver(rk4, 1e-5);
 * auto tr = simulate_two_rate<2, 1, 1, 1>(
 *     plant, output, outer, 10, inner, 1, plant_solver, 1e-4, x0, {0.0, 1.0});
 * @endcode
 *
 * @see simulate_sampled for the single-rate digital-control pattern
 * @see "Digital Control of Dynamic Systems" (Franklin, Powell & Workman, 3rd ed.),
 *      multirate sampling chapter
 */

#include <cstddef>
#include <utility>
#include <vector>

#include "damp/backend.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/simulation/simulate.hpp"
#include "damp/systems/state_space.hpp"

namespace damp::sim {

/**
 * @brief True when base tick @p tick is a firing instant for @p divisor
 *
 * Divisor 0 never fires (invalid schedule). Tick 0 fires for every positive
 * divisor so held signals are initialized on the first sample.
 */
[[nodiscard]] constexpr bool fires_at(std::size_t tick, std::size_t divisor) noexcept {
    return divisor != 0 && (tick % divisor) == 0;
}

/**
 * @brief Configuration for multi-rate closed-loop runs
 */
template<typename T = double>
struct MultiRateConfig {
    T           Ts_base{};         ///< Fastest discrete period [s] (must be > 0)
    std::size_t log_stride{1};     ///< Keep every N-th fine sample within a base period
    std::size_t max_base_ticks{0}; ///< 0 = unlimited; otherwise stop after this many base ticks
};

/**
 * @brief Result of a two-rate cascade simulation
 *
 * Extends the plant trajectory with the held outer reference @c r at each
 * logged sample (ZOH of the outer law between its firing instants).
 */
template<size_t NX, size_t NU, size_t NY, size_t NR, typename T = double>
struct TwoRateSimulationResult {
    std::vector<T>             t{};
    std::vector<ColVec<NX, T>> x{};
    std::vector<ColVec<NY, T>> y{};
    std::vector<ColVec<NU, T>> u{};
    std::vector<ColVec<NR, T>> r{}; ///< Held outer → inner reference
    bool                       truncated{false};
    std::size_t                base_tick_count{0};
};

/**
 * @brief Multi-rate closed-loop simulation with @ref MultiRateConfig
 *
 * At each base tick @c k the callable
 * @c u = control_step(k, t, y) runs once; @c u is held while @p fine_solver
 * integrates the continuous plant over the next @c Ts_base. Implement multi-rate
 * scheduling inside @p control_step with @ref fires_at (or free-form logic).
 *
 * @tparam Plant        (T t, ColVec\<NX\> x, ColVec\<NU\> u) → ColVec\<NX\>  (dx/dt)
 * @tparam Output       (ColVec\<NX\> x) → ColVec\<NY\>
 * @tparam ControlStep  (std::size_t k, T t, ColVec\<NY\> y) → ColVec\<NU\>
 * @tparam Solver       FixedStepSolver (plant sub-step)
 *
 * @see simulate_sampled — equivalent when control_step ignores @c k and runs every tick
 */
template<
    size_t NX,
    size_t NU,
    size_t NY,
    typename T,
    typename Plant,
    typename Output,
    typename ControlStep,
    typename Solver>
[[nodiscard]] SimulationResult<NX, NU, NY, T> simulate_multirate(
    Plant&&                 plant,
    Output&&                output,
    ControlStep&&           control_step,
    const Solver&           fine_solver,
    const ColVec<NX, T>&    x0,
    const damp::pair<T, T>& t_span,
    MultiRateConfig<T>      config
) {
    SimulationResult<NX, NU, NY, T> sim{};
    if (!(config.Ts_base > T{0}) || !(t_span.second > t_span.first)) {
        return sim;
    }

    ColVec<NX, T>     x = x0;
    T                 t = t_span.first;
    const std::size_t stride = (config.log_stride == 0) ? 1 : config.log_stride;
    std::size_t       k = 0;

    while (t < t_span.second - static_cast<T>(1e-12)) {
        if (config.max_base_ticks != 0 && k >= config.max_base_ticks) {
            break;
        }

        const T             t_next = (t + config.Ts_base < t_span.second) ? t + config.Ts_base : t_span.second;
        const ColVec<NY, T> y = output(x);
        const ColVec<NU, T> u = control_step(k, t, y);

        if (sim.t.empty()) {
            sim.t.push_back(t);
            sim.x.push_back(x);
            sim.y.push_back(y);
            sim.u.push_back(u);
        }

        const auto f = [&](T tt, const ColVec<NX, T>& xx) { return plant(tt, xx, u); };
        const auto seg = fine_solver.solve(f, x, {t, t_next});
        for (size_t j = 1; j < seg.size(); ++j) {
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
        ++k;
    }
    return sim;
}

/**
 * @brief Multi-rate closed-loop simulation (Ts_base + log_stride convenience)
 */
template<
    size_t NX,
    size_t NU,
    size_t NY,
    typename T,
    typename Plant,
    typename Output,
    typename ControlStep,
    typename Solver>
[[nodiscard]] SimulationResult<NX, NU, NY, T> simulate_multirate(
    Plant&&                 plant,
    Output&&                output,
    ControlStep&&           control_step,
    const Solver&           fine_solver,
    T                       Ts_base,
    const ColVec<NX, T>&    x0,
    const damp::pair<T, T>& t_span,
    std::size_t             log_stride = 1
) {
    MultiRateConfig<T> cfg{.Ts_base = Ts_base, .log_stride = log_stride};
    return simulate_multirate<NX, NU, NY, T>(
        damp::forward<Plant>(plant),
        damp::forward<Output>(output),
        damp::forward<ControlStep>(control_step),
        fine_solver,
        x0,
        t_span,
        cfg
    );
}

/**
 * @brief Two-rate cascade with @ref MultiRateConfig
 *
 * Outer fires every @p outer_div base ticks and produces a held reference
 * @c r for the inner law. Inner fires every @p inner_div base ticks and
 * produces the held plant input @c u. Both laws see the plant sample at the
 * start of the base tick (same ADC/snapshot semantics as a shared sample
 * instant on hardware).
 *
 * Typical cascade: @p inner_div = 1 (ISR), @p outer_div ≫ 1 (RTOS task).
 * The outer task only observes plant state on its own firing grid — the
 * inter-sample evolution is invisible to it (transport delay / hold).
 *
 * @tparam NR     Dimension of the outer → inner reference
 * @tparam Outer  (T t, ColVec\<NY\> y) → ColVec\<NR\>
 * @tparam Inner  (T t, ColVec\<NR\> r, ColVec\<NY\> y) → ColVec\<NU\>
 */
template<
    size_t NX,
    size_t NU,
    size_t NY,
    size_t NR,
    typename T,
    typename Plant,
    typename Output,
    typename Outer,
    typename Inner,
    typename Solver>
[[nodiscard]] TwoRateSimulationResult<NX, NU, NY, NR, T> simulate_two_rate(
    Plant&&                 plant,
    Output&&                output,
    Outer&&                 outer,
    std::size_t             outer_div,
    Inner&&                 inner,
    std::size_t             inner_div,
    const Solver&           fine_solver,
    const ColVec<NX, T>&    x0,
    const damp::pair<T, T>& t_span,
    MultiRateConfig<T>      config
) {
    TwoRateSimulationResult<NX, NU, NY, NR, T> sim{};
    if (!(config.Ts_base > T{0}) || !(t_span.second > t_span.first)) {
        return sim;
    }
    if (outer_div == 0 || inner_div == 0) {
        return sim;
    }

    ColVec<NX, T>     x = x0;
    ColVec<NR, T>     r{};
    ColVec<NU, T>     u{};
    T                 t = t_span.first;
    const std::size_t stride = (config.log_stride == 0) ? 1 : config.log_stride;
    std::size_t       k = 0;

    while (t < t_span.second - static_cast<T>(1e-12)) {
        if (config.max_base_ticks != 0 && k >= config.max_base_ticks) {
            sim.truncated = true;
            break;
        }

        const T             t_next = (t + config.Ts_base < t_span.second) ? t + config.Ts_base : t_span.second;
        const ColVec<NY, T> y = output(x);

        // Slow → fast: outer reference first, then inner plant command.
        if (fires_at(k, outer_div)) {
            r = outer(t, y);
        }
        if (fires_at(k, inner_div)) {
            u = inner(t, r, y);
        }

        if (sim.t.empty()) {
            sim.t.push_back(t);
            sim.x.push_back(x);
            sim.y.push_back(y);
            sim.u.push_back(u);
            sim.r.push_back(r);
        }

        const auto f = [&](T tt, const ColVec<NX, T>& xx) { return plant(tt, xx, u); };
        const auto seg = fine_solver.solve(f, x, {t, t_next});
        for (size_t j = 1; j < seg.size(); ++j) {
            const bool keep = (j % stride == 0) || (j + 1 == seg.size());
            if (!keep) {
                continue;
            }
            sim.t.push_back(seg.t[j]);
            sim.x.push_back(seg.x[j]);
            sim.y.push_back(output(seg.x[j]));
            sim.u.push_back(u);
            sim.r.push_back(r);
        }
        x = seg.x.back();
        t = t_next;
        ++k;
    }
    sim.base_tick_count = k;
    return sim;
}

/**
 * @brief Two-rate cascade: outer reference law + inner plant-input law
 *
 * @param plant        Continuous plant f(t, x, u) -> dx/dt
 * @param output       Measurement map y = g(x)
 * @param outer        Outer reference law (fires every outer_div base ticks)
 * @param outer_div    Outer period as multiple of Ts_base
 * @param inner        Inner plant-input law (fires every inner_div base ticks)
 * @param inner_div    Inner period as multiple of Ts_base
 * @param fine_solver  Fixed-step continuous solver
 * @param Ts_base      Base (fastest) discrete period [s]
 * @param x0           Initial plant state
 * @param t_span       Integration interval {t0, tf}
 * @param log_stride   Record every N-th fine sample within a base period (default 1)
 */
template<
    size_t NX,
    size_t NU,
    size_t NY,
    size_t NR,
    typename T,
    typename Plant,
    typename Output,
    typename Outer,
    typename Inner,
    typename Solver>
[[nodiscard]] TwoRateSimulationResult<NX, NU, NY, NR, T> simulate_two_rate(
    Plant&&                 plant,
    Output&&                output,
    Outer&&                 outer,
    std::size_t             outer_div,
    Inner&&                 inner,
    std::size_t             inner_div,
    const Solver&           fine_solver,
    T                       Ts_base,
    const ColVec<NX, T>&    x0,
    const damp::pair<T, T>& t_span,
    std::size_t             log_stride = 1
) {
    MultiRateConfig<T> cfg{.Ts_base = Ts_base, .log_stride = log_stride};
    return simulate_two_rate<NX, NU, NY, NR, T>(
        damp::forward<Plant>(plant),
        damp::forward<Output>(output),
        damp::forward<Outer>(outer),
        outer_div,
        damp::forward<Inner>(inner),
        inner_div,
        fine_solver,
        x0,
        t_span,
        cfg
    );
}

/**
 * @brief Two-rate cascade on a continuous LTI plant (A, B, C)
 *
 * Plant RHS is @c Ax + Bu; control sampling uses @c y = Cx (same as
 * simulate_sampled — no D feedthrough of the not-yet-computed @c u).
 */
template<
    size_t NX,
    size_t NU,
    size_t NY,
    size_t NR,
    size_t NW,
    size_t NV,
    typename T,
    typename Outer,
    typename Inner,
    typename Solver>
[[nodiscard]] TwoRateSimulationResult<NX, NU, NY, NR, T> simulate_two_rate_lti(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    Outer&&                                  outer,
    std::size_t                              outer_div,
    Inner&&                                  inner,
    std::size_t                              inner_div,
    const Solver&                            fine_solver,
    T                                        Ts_base,
    const ColVec<NX, T>&                     x0,
    const damp::pair<T, T>&                  t_span,
    std::size_t                              log_stride = 1
) {
    auto plant = [&](T /*t*/, const ColVec<NX, T>& x, const ColVec<NU, T>& u) -> ColVec<NX, T> {
        return sys.A * x + sys.B * u;
    };
    auto output = [&](const ColVec<NX, T>& x) -> ColVec<NY, T> { return sys.C * x; };

    return simulate_two_rate<NX, NU, NY, NR, T>(
        plant,
        output,
        damp::forward<Outer>(outer),
        outer_div,
        damp::forward<Inner>(inner),
        inner_div,
        fine_solver,
        Ts_base,
        x0,
        t_span,
        log_stride
    );
}

} // namespace damp::sim
