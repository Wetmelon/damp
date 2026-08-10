// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file time_response.hpp
 * @brief Step / impulse / lsim and stepinfo / lsiminfo (host-only)
 *
 * Host tooling: result storage uses @c std::vector. Not part of the embeddable
 * @c control.hpp surface — pull via @c workbench.hpp or this header directly.
 *
 * Ts consistency. Continuous plants are ZOH-discretized onto the uniform
 * grid spacing @c dt = time[1]-time[0]. Discrete plants (@c sys.Ts > 0) are
 * iterated with their own sample time; the time vector must match that period
 * (within default_tol) or the response is empty. Do not mix a discrete
 * model with an unrelated host grid.
 */

#include <cstddef>
#include <vector>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"


namespace damp {
namespace analysis {

// ============================================================================
// Time-Domain Responses
// ============================================================================

/**
 * @brief Multi-channel time-domain response sampled on a time grid
 *
 * For a system with NU inputs and NY outputs, `y[k](i, j)` is output i at time
 * step k in response to a canonical input (unit step / unit impulse) applied to
 * input channel j alone. This mirrors MATLAB®'s (Nt × Ny × Nu) response array.
 *
 * @tparam NY Number of outputs
 * @tparam NU Number of inputs
 */
template<size_t NY, size_t NU, typename T = double>
struct TimeResponse {
    std::vector<T>                 t; ///< Time points (s)
    std::vector<Matrix<NY, NU, T>> y; ///< y[k](i,j): output i from a canonical input on channel j
};

/**
 * @brief Result of a single-trajectory simulation: time, output, and state history
 *
 * Used by `lsim` (forced response to a given input) and `initial` (free
 * response). `y[k]` is the NY-vector output and `x[k]` the full state at step k.
 *
 * @tparam NX Number of states
 * @tparam NY Number of outputs
 */
template<size_t NX, size_t NY, typename T = double>
struct LsimResult {
    std::vector<T>             t; ///< Time points (s)
    std::vector<ColVec<NY, T>> y; ///< Output vector at each time point
    std::vector<ColVec<NX, T>> x; ///< State vector at each time point
};

namespace detail {

/**
 * @brief ZOH-discretize a continuous system onto a uniform time grid.
 *
 * Pass-through if already discrete and the grid spacing matches @c sys.Ts.
 * Returns @c nullopt when the grid is empty, @c dt ≤ 0, discrete @c Ts
 * disagrees with the grid, or (Tustin-only path) discretization fails — ZOH
 * always succeeds for finite matrices.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] damp::optional<StateSpace<NX, NU, NY, T, NW, NV>> discretize_on_grid(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const std::vector<T>&                    time
) {
    if (time.empty()) {
        return damp::nullopt;
    }
    T dt = T{1};
    if (time.size() > 1) {
        dt = time[1] - time[0];
    } else if (sys.is_discrete()) {
        dt = sys.Ts;
    }
    if (!(dt > T{0})) {
        return damp::nullopt;
    }
    if (sys.is_discrete()) {
        const T tol = default_tol<T>() * (T{1} + damp::abs(sys.Ts));
        if (damp::abs(dt - sys.Ts) > tol) {
            return damp::nullopt; // host grid ≠ plant Ts
        }
        return sys;
    }
    return damp::discretize(sys, dt, DiscretizationMethod::ZOH);
}

//! Iterate a discrete system from x0 with constant input u, collecting y[k]=Cx+Du for n steps.
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] std::vector<ColVec<NY, T>> const_input_response(
    const StateSpace<NX, NU, NY, T, NW, NV>& dsys,
    const ColVec<NX, T>&                     x0,
    const ColVec<NU, T>&                     u,
    size_t                                   n
) {
    std::vector<ColVec<NY, T>> y;
    y.reserve(n);
    ColVec<NX, T> x = x0;
    for (size_t k = 0; k < n; ++k) {
        y.push_back(ColVec<NY, T>(dsys.C * x + dsys.D * u));
        x = ColVec<NX, T>(dsys.A * x + dsys.B * u);
    }
    return y;
}

//! Scatter a single-channel response (column j) into a multi-channel TimeResponse.
template<size_t NY, size_t NU, typename T>
void set_column(TimeResponse<NY, NU, T>& r, size_t j, const std::vector<ColVec<NY, T>>& yj) {
    for (size_t k = 0; k < yj.size(); ++k) {
        for (size_t i = 0; i < NY; ++i) {
            r.y[k](i, j) = yj[k][i];
        }
    }
}

} // namespace detail

/**
 * @brief Step response of a (MIMO) state-space system
 *
 * Applies a unit step to each input channel in turn (others held at zero) from
 * zero initial state and records every output, so `y[k](i, j)` is output i at
 * step k due to a step on input j. Continuous systems are ZOH-discretized on the
 * (uniform) time grid; discrete systems are iterated directly. MATLAB®
 * equivalent: `step(sys, t)`.
 *
 * @param sys  State-space system (continuous or discrete)
 * @param time Uniformly spaced time vector (e.g. analysis::linspace(0, tf, n))
 * @return TimeResponse with per-channel output history; empty if grid/Ts invalid
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] TimeResponse<NY, NU, T> step(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const std::vector<T>&                    time
) {
    TimeResponse<NY, NU, T> r;
    const auto              dsys_opt = detail::discretize_on_grid(sys, time);
    if (!dsys_opt) {
        return r;
    }
    const auto& dsys = *dsys_opt;
    r.t = time;
    r.y.assign(time.size(), Matrix<NY, NU, T>{});
    for (size_t j = 0; j < NU; ++j) {
        ColVec<NU, T> uj{};
        uj[j] = T{1};
        detail::set_column(r, j, detail::const_input_response(dsys, ColVec<NX, T>{}, uj, time.size()));
    }
    return r;
}

/**
 * @brief Impulse response of a (MIMO) state-space system
 *
 * For each input channel j, computes the free response from initial state
 * x₀ = B(:,j), giving y(t) = C·e^{At}·B(:,j) (the strictly-proper part; the
 * D·δ(t) feedthrough term is not plotted). MATLAB® equivalent: `impulse(sys, t)`.
 *
 * @param sys  State-space system
 * @param time Uniformly spaced time vector
 * @return TimeResponse with per-channel output history; empty if grid/Ts invalid
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] TimeResponse<NY, NU, T> impulse(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const std::vector<T>&                    time
) {
    TimeResponse<NY, NU, T> r;
    const auto              dsys_opt = detail::discretize_on_grid(sys, time);
    if (!dsys_opt) {
        return r;
    }
    const auto& dsys = *dsys_opt;
    r.t = time;
    r.y.assign(time.size(), Matrix<NY, NU, T>{});
    for (size_t j = 0; j < NU; ++j) {
        ColVec<NX, T> x0j{};
        for (size_t i = 0; i < NX; ++i) {
            x0j[i] = sys.B(i, j);
        }
        detail::set_column(r, j, detail::const_input_response(dsys, x0j, ColVec<NU, T>{}, time.size()));
    }
    return r;
}

/**
 * @brief Initial-condition (free) response of a (MIMO) state-space system
 *
 * Computes y(t) = C·e^{At}·x₀ for the unforced system (u=0). MATLAB®
 * equivalent: `initial(sys, x0, t)`.
 *
 * @param sys  State-space system
 * @param x0   Initial state
 * @param time Uniformly spaced time vector
 * @return LsimResult with time, output, and state history; empty if grid/Ts invalid
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] LsimResult<NX, NY, T> initial(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const ColVec<NX, T>&                     x0,
    const std::vector<T>&                    time
) {
    LsimResult<NX, NY, T> r;
    const auto            dsys_opt = detail::discretize_on_grid(sys, time);
    if (!dsys_opt) {
        return r;
    }
    const auto& dsys = *dsys_opt;
    r.t = time;
    r.y.reserve(time.size());
    r.x.reserve(time.size());
    ColVec<NX, T> x = x0;
    for (size_t k = 0; k < time.size(); ++k) {
        r.x.push_back(x);
        r.y.push_back(ColVec<NY, T>(dsys.C * x));
        x = ColVec<NX, T>(dsys.A * x);
    }
    return r;
}

/**
 * @brief Forced time response of a (MIMO) state-space system to an input signal
 *
 * Simulates ẋ = Ax + Bu, y = Cx + Du driven by the supplied input samples,
 * from initial state x₀. Continuous systems are ZOH-discretized on the (uniform)
 * time grid — exact at the sample points for piecewise-constant input — and
 * discrete systems are iterated directly. MATLAB® equivalent: `lsim(sys, u, t, x0)`.
 *
 * @param sys  State-space system (continuous or discrete)
 * @param u    Input samples, one ColVec\<NU\> per time point (length == time.size())
 * @param time Uniformly spaced time vector
 * @param x0   Initial state (defaults to zero)
 * @return LsimResult with time, output, and state history; empty on length/Ts mismatch
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] LsimResult<NX, NY, T> lsim(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const std::vector<ColVec<NU, T>>&        u,
    const std::vector<T>&                    time,
    const ColVec<NX, T>&                     x0 = {}
) {
    LsimResult<NX, NY, T> r;
    // Guard length mismatch: never read past the shorter of u and time.
    if (u.size() != time.size() || time.empty()) {
        return r; // empty result (caller can check r.t.empty())
    }

    const auto dsys_opt = detail::discretize_on_grid(sys, time);
    if (!dsys_opt) {
        return r;
    }
    const auto& dsys = *dsys_opt;

    r.t = time;
    r.y.reserve(time.size());
    r.x.reserve(time.size());

    ColVec<NX, T> x = x0;
    for (size_t k = 0; k < time.size(); ++k) {
        r.x.push_back(x);
        r.y.push_back(ColVec<NY, T>(dsys.C * x + dsys.D * u[k]));
        x = ColVec<NX, T>(dsys.A * x + dsys.B * u[k]);
    }
    return r;
}

/**
 * @brief Single-input convenience overload of lsim taking a scalar input signal
 *
 * MATLAB® allows `lsim(sys, u, t)` with a plain vector u for single-input systems.
 */
template<size_t NX, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] LsimResult<NX, NY, T> lsim(
    const StateSpace<NX, 1, NY, T, NW, NV>& sys,
    const std::vector<T>&                   u,
    const std::vector<T>&                   time,
    const ColVec<NX, T>&                    x0 = {}
) {
    std::vector<ColVec<1, T>> uv;
    uv.reserve(u.size());
    for (const T& ui : u) {
        uv.push_back(ColVec<1, T>{ui});
    }
    return lsim(sys, uv, time, x0);
}

/**
 * @brief Step response of a SISO transfer function (MATLAB® `step(tf, t)`)
 *
 * Empty if companion realization fails or the time grid is invalid.
 */
template<size_t Nnum, size_t Nden, typename T>
[[nodiscard]] TimeResponse<1, 1, T> step(
    const TransferFunction<Nnum, Nden, T>& tf,
    const std::vector<T>&                  time
) {
    const auto ss = tf.to_state_space();
    if (!ss) {
        return {};
    }
    return step(*ss, time);
}

/**
 * @brief Impulse response of a SISO transfer function (MATLAB® `impulse(tf, t)`)
 *
 * Empty if companion realization fails or the time grid is invalid.
 */
template<size_t Nnum, size_t Nden, typename T>
[[nodiscard]] TimeResponse<1, 1, T> impulse(
    const TransferFunction<Nnum, Nden, T>& tf,
    const std::vector<T>&                  time
) {
    const auto ss = tf.to_state_space();
    if (!ss) {
        return {};
    }
    return impulse(*ss, time);
}

// ============================================================================
// Response Characteristics
// ============================================================================

/**
 * @brief Step-response characteristics of a single output signal
 *
 * Sample-resolution metrics (no sub-sample interpolation) measured against the
 * initial value y(0) and the supplied steady-state value. MATLAB® equivalent:
 * `stepinfo`.
 */
template<typename T = double>
struct StepInfo {
    T rise_time{};     ///< Time to go from 10% to 90% of the total change
    T settling_time{}; ///< Last time the response stays within the settling band of yfinal
    T settling_min{};  ///< Minimum value once the response first enters the settling band
    T settling_max{};  ///< Maximum value once the response first enters the settling band
    T overshoot{};     ///< Percent overshoot beyond yfinal (0 if none)
    T undershoot{};    ///< Percent undershoot below the initial value (0 if none)
    T peak{};          ///< Peak absolute value of the response
    T peak_time{};     ///< Time of the peak absolute value
};

/**
 * @brief Compute step-response characteristics from an output/time signal
 *
 * @param y           Output samples
 * @param t           Matching time samples (same length as y)
 * @param yfinal      Steady-state value the response settles to
 * @param settle_frac Settling band as a fraction of |yfinal - y0| (default 0.02 → ±2%)
 */
template<typename T = double>
[[nodiscard]] StepInfo<T> stepinfo(
    const std::vector<T>& y,
    const std::vector<T>& t,
    T                     yfinal,
    T                     settle_frac = T(0.02)
) {
    StepInfo<T> info;
    if (y.empty()) {
        return info;
    }
    const T y0 = y.front();
    const T span = yfinal - y0;
    const T aspan = damp::abs(span);

    // Rise time: first crossing of 10% then 90% of the total change.
    const T lo = y0 + T(0.1) * span;
    const T hi = y0 + T(0.9) * span;
    T       t_lo = t.front();
    T       t_hi = t.front();
    bool    got_lo = false;
    bool    got_hi = false;
    for (size_t k = 0; k < y.size(); ++k) {
        const T reach_lo = (span >= T(0)) ? (y[k] >= lo) : (y[k] <= lo);
        const T reach_hi = (span >= T(0)) ? (y[k] >= hi) : (y[k] <= hi);
        if (!got_lo && reach_lo) {
            t_lo = t[k];
            got_lo = true;
        }
        if (!got_hi && reach_hi) {
            t_hi = t[k];
            got_hi = true;
            break;
        }
    }
    info.rise_time = t_hi - t_lo;

    // Settling time: last instant the response is outside the ±settle_frac band.
    const T band = settle_frac * aspan;
    info.settling_time = t.front();
    for (size_t k = 0; k < y.size(); ++k) {
        if (damp::abs(y[k] - yfinal) > band) {
            info.settling_time = (k + 1 < t.size()) ? t[k + 1] : t[k];
        }
    }
    // settling_min/max: extremes after first entering the band.
    info.settling_min = yfinal;
    info.settling_max = yfinal;
    bool entered = false;
    for (size_t k = 0; k < y.size(); ++k) {
        if (!entered && damp::abs(y[k] - yfinal) <= band) {
            entered = true;
            info.settling_min = y[k];
            info.settling_max = y[k];
        }
        if (entered) {
            info.settling_min = damp::min(info.settling_min, y[k]);
            info.settling_max = damp::max(info.settling_max, y[k]);
        }
    }

    // Peak (absolute), overshoot, undershoot.
    info.peak = damp::abs(y.front());
    info.peak_time = t.front();
    T ymax = y.front();
    T ymin = y.front();
    for (size_t k = 0; k < y.size(); ++k) {
        if (damp::abs(y[k]) > info.peak) {
            info.peak = damp::abs(y[k]);
            info.peak_time = t[k];
        }
        ymax = damp::max(ymax, y[k]);
        ymin = damp::min(ymin, y[k]);
    }
    if (aspan > T(0)) {
        const T over = (span >= T(0)) ? (ymax - yfinal) : (yfinal - ymin);
        const T under = (span >= T(0)) ? (y0 - ymin) : (ymax - y0);
        info.overshoot = damp::max(T(0), over / aspan * T(100));
        info.undershoot = damp::max(T(0), under / aspan * T(100));
    }
    return info;
}

/**
 * @brief Step-response characteristics of a SISO system (MATLAB® `stepinfo(sys)`)
 *
 * Runs `step(sys, time)` and summarizes the single output. The steady-state
 * value is taken as the last sample. For MIMO systems, extract the desired
 * `y[k](i, j)` channel and call the signal overload.
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] StepInfo<T> stepinfo(
    const StateSpace<NX, 1, 1, T, NW, NV>& sys,
    const std::vector<T>&                  time
) {
    const auto     resp = step(sys, time);
    std::vector<T> y;
    y.reserve(resp.y.size());
    for (const auto& yk : resp.y) {
        y.push_back(yk(0, 0));
    }
    return stepinfo(y, resp.t, y.empty() ? T(0) : y.back());
}

/**
 * @brief Transient characteristics of an arbitrary response signal
 *
 * Like `stepinfo` but makes no step assumption: reports settling relative to a
 * given final value plus the signal's extremes. MATLAB® equivalent: `lsiminfo`.
 */
template<typename T = double>
struct LsimInfo {
    T settling_time{}; ///< Last time the response leaves the settling band of yfinal
    T settling_min{};  ///< Minimum after first entering the band
    T settling_max{};  ///< Maximum after first entering the band
    T min{};           ///< Global minimum of the signal
    T max{};           ///< Global maximum of the signal
    T min_time{};      ///< Time of the global minimum
    T max_time{};      ///< Time of the global maximum
};

/**
 * @brief Compute transient characteristics from an output/time signal
 *
 * @param y           Output samples
 * @param t           Matching time samples
 * @param yfinal      Reference value for the settling band
 * @param settle_frac Settling band as a fraction of |yfinal| (default 0.02)
 */
template<typename T = double>
[[nodiscard]] LsimInfo<T> lsiminfo(
    const std::vector<T>& y,
    const std::vector<T>& t,
    T                     yfinal,
    T                     settle_frac = T(0.02)
) {
    LsimInfo<T> info;
    if (y.empty()) {
        return info;
    }
    const T band = settle_frac * damp::abs(yfinal);

    info.min = y.front();
    info.max = y.front();
    info.min_time = t.front();
    info.max_time = t.front();
    info.settling_time = t.front();
    for (size_t k = 0; k < y.size(); ++k) {
        if (y[k] < info.min) {
            info.min = y[k];
            info.min_time = t[k];
        }
        if (y[k] > info.max) {
            info.max = y[k];
            info.max_time = t[k];
        }
        if (damp::abs(y[k] - yfinal) > band) {
            info.settling_time = (k + 1 < t.size()) ? t[k + 1] : t[k];
        }
    }

    info.settling_min = yfinal;
    info.settling_max = yfinal;
    bool entered = false;
    for (size_t k = 0; k < y.size(); ++k) {
        if (!entered && damp::abs(y[k] - yfinal) <= band) {
            entered = true;
            info.settling_min = y[k];
            info.settling_max = y[k];
        }
        if (entered) {
            info.settling_min = damp::min(info.settling_min, y[k]);
            info.settling_max = damp::max(info.settling_max, y[k]);
        }
    }
    return info;
}

} // namespace analysis
} // namespace damp
