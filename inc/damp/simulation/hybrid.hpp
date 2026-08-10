// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file hybrid.hpp
 * @brief Event-driven piecewise-LTI simulation (Exact ZOH between events)
 *
 * Host-only — not embeddable (`std::vector` result storage). Pull via
 * @c workbench.hpp on the workstation; never into target firmware.
 *
 * Piecewise-linear plants (switched converters, etc.): each mode is a continuous
 * StateSpace with @c Ts == 0; between events the plant is advanced with
 * Exact (matrix exponential, input held over the event interval — that
 * interval is the effective ZOH period, independent of any discrete @c Ts on
 * the mode tables). Nonlinear plants (drone, flexible mechanisms) stay on
 * @ref simulate / simulate_sampled — do not force them into LTI modes.
 *
 * Log density:
 * - LogPolicy::EventsOnly — sample at each event (PWM edges, control ticks)
 * - LogPolicy::FixedStride — also sample every @c h_fine inside segments
 *   (dense tank waveforms); prefer NPY export for large runs, not Plotly HTML
 *
 * @code
 * PiecewiseLTI<2, 2, 1, 1> plant{.modes = {mode_on, mode_off}};
 * auto schedule = [&](double t, std::size_t mode, auto const& x, auto const& u) {
 *     return next_pwm_edge(t); // absolute time
 * };
 * auto on_event = [&](double t, std::size_t mode, auto const& x, auto const& y, auto const& u) {
 *     return HybridEventAction<2, 1, 1>{.mode = 1 - mode, .u = u, .x = x};
 * };
 * auto tr = simulate_piecewise_lti(plant, 0, x0, u0, {0.0, 1e-3}, schedule, on_event);
 * @endcode
 *
 * @see Exact
 */

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

#include "damp/backend.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/simulation/cached_zoh.hpp"
#include "damp/simulation/integrator.hpp"
#include "damp/systems/state_space.hpp"

namespace damp::sim {

/**
 * @brief How often to record samples along a hybrid trajectory
 */
enum class LogPolicy : std::uint8_t {
    EventsOnly,  ///< At each event (and endpoints)
    FixedStride, ///< Every h_fine inside segments, plus events
};

/**
 * @brief Configuration for event-driven hybrid / multi-rate style runs
 */
template<typename T = double>
struct EventSimConfig {
    LogPolicy   log = LogPolicy::EventsOnly;
    T           h_fine = T{0};   ///< Used when log == FixedStride; must be > 0
    T           min_step = T{0}; ///< Skip / clamp steps smaller than this (0 = none)
    std::size_t max_events = 1'000'000;
};

/**
 * @brief Result of a piecewise-LTI hybrid simulation
 *
 * Same shape as @ref SimulationResult plus the active mode index at each sample.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
struct HybridSimulationResult {
    std::vector<T>             t{};
    std::vector<ColVec<NX, T>> x{};
    std::vector<ColVec<NY, T>> y{};
    std::vector<ColVec<NU, T>> u{};
    std::vector<std::uint32_t> mode{};
    bool                       truncated{false}; ///< Hit max_events before tf
    std::size_t                event_count{0};
};

/**
 * @brief Plant: NModes continuous LTI pieces (same state/input/output dimensions)
 */
template<size_t NModes, size_t NX, size_t NU, size_t NY, typename T = double>
struct PiecewiseLTI {
    damp::array<StateSpace<NX, NU, NY, T>, NModes> modes{};
};

/**
 * @brief Action applied when an event fires (after Exact advance to the event time)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
struct HybridEventAction {
    std::size_t   mode{}; ///< Mode for the *next* segment
    ColVec<NU, T> u{};    ///< Input held over the next segment (ZOH)
    ColVec<NX, T> x{};    ///< State after optional projection / reset
};

namespace hybrid_detail {

template<size_t NX, size_t NU, size_t NY, typename T>
[[nodiscard]] constexpr ColVec<NY, T> output_of(
    const StateSpace<NX, NU, NY, T>& sys,
    const ColVec<NX, T>&             x,
    const ColVec<NU, T>&             u
) {
    return sys.C * x + sys.D * u;
}

template<size_t NX, size_t NU, size_t NY, typename T>
void push_sample(
    HybridSimulationResult<NX, NU, NY, T>& tr,
    T                                      t,
    const ColVec<NX, T>&                   x,
    const ColVec<NY, T>&                   y,
    const ColVec<NU, T>&                   u,
    std::uint32_t                          mode
) {
    tr.t.push_back(t);
    tr.x.push_back(x);
    tr.y.push_back(y);
    tr.u.push_back(u);
    tr.mode.push_back(mode);
}

template<size_t NX, size_t NU, typename T>
[[nodiscard]] ColVec<NX, T> exact_step(
    const Exact<NX, T>&      exact,
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const ColVec<NX, T>&     x,
    const ColVec<NU, T>&     u,
    T                        h
) {
    if (!(h > T{0})) {
        return x;
    }
    return exact.evolve(A, B, x, u, h).x;
}

} // namespace hybrid_detail

/**
 * @brief Simulate a piecewise-LTI plant with Exact integration between events
 *
 * @tparam NModes  Number of LTI modes
 * @tparam NX,NU,NY Dimensions
 * @tparam T       Scalar
 * @tparam Schedule  (T t, size_t mode, ColVec x, ColVec u) -> T t_next
 *                   Absolute time of the next event; must satisfy t_next > t
 *                   (or >= tf to finish). Clamped to tf.
 * @tparam OnEvent   (T t, size_t mode, ColVec x, ColVec y, ColVec u)
 *                   -> HybridEventAction  applied at each event at time t
 *                   (after the Exact advance onto t).
 *
 * @param plant     Mode table (continuous StateSpace, Ts == 0)
 * @param mode0     Initial mode index in [0, NModes)
 * @param x0        Initial state
 * @param u0        Initial input (held until first OnEvent change)
 * @param t_span    {t0, tf}
 * @param schedule  Next event time
 * @param on_event  Mode/input/state update at events
 * @param config    Logging and limits
 */
template<
    size_t NModes,
    size_t NX,
    size_t NU,
    size_t NY,
    typename T,
    typename Schedule,
    typename OnEvent>
[[nodiscard]] HybridSimulationResult<NX, NU, NY, T> simulate_piecewise_lti(
    const PiecewiseLTI<NModes, NX, NU, NY, T>& plant,
    std::size_t                                mode0,
    const ColVec<NX, T>&                       x0,
    const ColVec<NU, T>&                       u0,
    const damp::pair<T, T>&                    t_span,
    Schedule&&                                 schedule,
    OnEvent&&                                  on_event,
    EventSimConfig<T>                          config = {}
) {
    HybridSimulationResult<NX, NU, NY, T> tr{};
    if (!(t_span.second > t_span.first) || mode0 >= NModes) {
        return tr;
    }

    const Exact<NX, T> exact{};
    const T            t0 = t_span.first;
    const T            tf = t_span.second;

    std::size_t   mode = mode0;
    ColVec<NX, T> x = x0;
    ColVec<NU, T> u = u0;
    T             t = t0;

    const auto& sys0 = plant.modes[mode];
    hybrid_detail::push_sample(tr, t, x, hybrid_detail::output_of(sys0, x, u), u, static_cast<std::uint32_t>(mode));

    std::size_t events = 0;
    while (t < tf) {
        if (events >= config.max_events) {
            tr.truncated = true;
            break;
        }

        T t_next = schedule(t, mode, x, u);
        if (t_next > tf) {
            t_next = tf;
        }
        if (!(t_next > t)) {
            // Schedule returned a non-advancing time — stop cleanly at tf if we can
            if (t < tf) {
                // Last open segment to tf with current mode
                const auto& sys = plant.modes[mode];
                x = hybrid_detail::exact_step(exact, sys.A, sys.B, x, u, tf - t);
                t = tf;
                hybrid_detail::push_sample(
                    tr, t, x, hybrid_detail::output_of(sys, x, u), u, static_cast<std::uint32_t>(mode)
                );
            }
            break;
        }

        const auto& sys = plant.modes[mode];
        const T     h_seg = t_next - t;

        if (config.log == LogPolicy::FixedStride && config.h_fine > T{0}) {
            // Full fine steps share one CachedZoh (expm once per segment stride);
            // any remainder shorter than h_fine still uses Exact.
            const auto    zoh = CachedZoh<NX, NU, T>::from(sys.A, sys.B, config.h_fine);
            T             t_cursor = t;
            ColVec<NX, T> x_cursor = x;
            while (t_cursor < t_next) {
                T h = config.h_fine;
                if (t_cursor + h > t_next) {
                    h = t_next - t_cursor;
                }
                if (config.min_step > T{0} && h < config.min_step && t_cursor + h < t_next) {
                    h = (config.min_step < (t_next - t_cursor)) ? config.min_step : (t_next - t_cursor);
                }
                if (h == config.h_fine) {
                    x_cursor = zoh.step(x_cursor, u);
                } else {
                    x_cursor = hybrid_detail::exact_step(exact, sys.A, sys.B, x_cursor, u, h);
                }
                t_cursor += h;
                if (t_cursor < t_next) {
                    hybrid_detail::push_sample(
                        tr,
                        t_cursor,
                        x_cursor,
                        hybrid_detail::output_of(sys, x_cursor, u),
                        u,
                        static_cast<std::uint32_t>(mode)
                    );
                }
            }
            x = x_cursor;
            t = t_next;
        } else {
            x = hybrid_detail::exact_step(exact, sys.A, sys.B, x, u, h_seg);
            t = t_next;
        }

        const ColVec<NY, T> y = hybrid_detail::output_of(plant.modes[mode], x, u);
        hybrid_detail::push_sample(tr, t, x, y, u, static_cast<std::uint32_t>(mode));

        if (t >= tf) {
            break;
        }

        // Event at t: update mode / u / x for the next segment
        HybridEventAction<NX, NU, NY, T> act = on_event(t, mode, x, y, u);
        if (act.mode >= NModes) {
            tr.truncated = true;
            break;
        }
        mode = act.mode;
        u = act.u;
        x = act.x;
        ++events;
        tr.event_count = events;

        // Re-log after projection / mode change (same t, new discrete state)
        const ColVec<NY, T> y2 = hybrid_detail::output_of(plant.modes[mode], x, u);
        if (!tr.t.empty() && tr.t.back() == t) {
            tr.x.back() = x;
            tr.y.back() = y2;
            tr.u.back() = u;
            tr.mode.back() = static_cast<std::uint32_t>(mode);
        } else {
            hybrid_detail::push_sample(tr, t, x, y2, u, static_cast<std::uint32_t>(mode));
        }
    }

    return tr;
}

/**
 * @brief Exact open-loop step of a single continuous LTI system
 *
 * Host helper for piecewise-Exact SILs (also used by unit tests).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] ColVec<NX, T> exact_lti_step(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const ColVec<NX, T>&                     x,
    const ColVec<NU, T>&                     u,
    T                                        h
) {
    Exact<NX, T> exact{};
    return hybrid_detail::exact_step(exact, sys.A, sys.B, x, u, h);
}

/**
 * @brief Locate a state zero-crossing of g(x) on an Exact LTI segment
 *
 * Advances from @p x0 under constant @p u with the given continuous LTI system.
 * Assumes g(x0) and g(x(t1)) have opposite signs (or a root in between). Returns
 * the time in (t0, t1] where g crosses zero (bisection on Exact evaluations).
 *
 * @tparam G Callable (ColVec<NX,T> const&) -> T
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename G>
[[nodiscard]] T locate_zero_crossing_exact(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const ColVec<NX, T>&                     x0,
    const ColVec<NU, T>&                     u,
    T                                        t0,
    T                                        t1,
    G&&                                      g,
    T                                        tol = static_cast<T>(1e-12),
    int                                      max_iter = 60
) {
    Exact<NX, T>  exact{};
    T             ta = t0;
    T             tb = t1;
    ColVec<NX, T> xa = x0;
    T             ga = g(xa);
    ColVec<NX, T> xb = hybrid_detail::exact_step(exact, sys.A, sys.B, x0, u, t1 - t0);
    T             gb = g(xb);
    if (ga * gb > T{0}) {
        return t1; // no sign change — return endpoint
    }
    for (int i = 0; i < max_iter; ++i) {
        const T tm = static_cast<T>(0.5) * (ta + tb);
        if (tb - ta <= tol) {
            return tm;
        }
        const ColVec<NX, T> xm = hybrid_detail::exact_step(exact, sys.A, sys.B, x0, u, tm - t0);
        const T             gm = g(xm);
        if (ga * gm <= T{0}) {
            tb = tm;
            gb = gm;
        } else {
            ta = tm;
            ga = gm;
        }
    }
    return static_cast<T>(0.5) * (ta + tb);
}

} // namespace damp::sim
