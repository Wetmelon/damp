// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file solver.hpp
 * @brief ODE solver wrappers (host)
 *
 * @defgroup solver ODE Solver Wrappers
 * @brief Host fixed- and adaptive-step ODE solvers (RK4, DP45, events) over integrators
 *
 * High-level solvers that manage time-stepping, result storage, event
 * detection, and zero-crossing location on top of the integrators in
 * integrator.hpp. Host-only (allocating result storage).
 *
 * Fixed step — pick any integrator, hold h constant:
 * @code
 *   auto result = fixed_solve\<RK4\>(f, x0, {0.0, 10.0}, 0.01);
 *   // or: FixedStepSolver solver(RK4<2>{}, 0.01); solver.solve(...);
 * @endcode
 *
 * Adaptive step — only pairs modelling @ref AdaptiveStepIntegrator
 * (DP45 non-stiff, RK23 cheap non-stiff, TRBDF2 stiff). Error control matches
 * SciPy `solve_ivp` (atol=1e-6, rtol=1e-3, weighted RMS ≤ 1), with Hairer
 * initial-step selection, Gustafsson PI step control, and scalar-or-vector tols:
 * @code
 *   auto result = adaptive_solve\<DP45\>(f, x0, {0.0, 10.0});
 *   auto stiff  = adaptive_solve\<TRBDF2\>(f, x0, {0.0, 10.0}, {.atol=1e-8, .rtol=1e-6});
 *   // or: AdaptiveStepSolver solver(DP45<2>{}, {.first_step=0.01}); ...
 *   // per-component: AdaptiveOptions<2>{.atol = ColVec<2>{1e-12, 1e-3}, ...}
 * @endcode
 */

#include <cstddef>
#include <functional>
#include <iterator>
#include <limits>
#include <vector>

#include "damp/backend.hpp"
#include "damp/matrix/colvec.hpp"
#include "integrator.hpp"

namespace damp::sim {

namespace detail {

/**
 * @brief Number of steps of size @p h that cover [t0, tf]
 *
 * A trailing partial step counts once; a remainder within 64·eps·n of a whole
 * step count (the rounding of (tf − t0)/h, e.g. 0.3/0.1 = 2.9999999999999996)
 * does not add a sliver step. Callers place step k at t0 + k·h (not a running
 * sum) and land the last one exactly on tf, so time labels do not drift.
 */
template<typename T>
[[nodiscard]] constexpr size_t step_count(T t0, T tf, T h) {
    const T n = (tf - t0) / h;
    if (!(n > T{0})) {
        return 0;
    }
    return static_cast<size_t>(damp::ceil(n - (T{64} * eps<T>() * n)));
}

/// Time of step @p k of @p n on the uniform grid t0 + k·h, with the last step exactly at tf.
template<typename T>
[[nodiscard]] constexpr T grid_time(T t0, T tf, T h, size_t k, size_t n) {
    return (k >= n) ? tf : t0 + (static_cast<T>(k) * h);
}

} // namespace detail

/**
 * @brief Result of an ODE solve operation
 *
 * Stores the full time and state history from a simulation.
 * Supports range-based for loops over (time, state) pairs.
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct SolveResult {
    std::vector<T>             t; ///< Time points
    std::vector<ColVec<NX, T>> x; ///< State vectors at each time point

    bool   success = true; ///< false on invalid input (h ≤ 0, tf < t0) or adaptive-control failure
    size_t nfev = 0;       ///< Number of integrator evolve() calls (accepted + rejected + event location)
    size_t n_accepted = 0; ///< Accepted steps (adaptive solvers only)
    size_t n_rejected = 0; ///< Rejected steps (adaptive solvers only)

    /// Iterator for range-based for: for (const auto& [t, x] : result)
    struct Iterator {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = damp::pair<T, const ColVec<NX, T>&>;
        using reference = value_type;

        Iterator(const std::vector<T>* t_vec, const std::vector<ColVec<NX, T>>* x_vec, size_t idx)
            : t_vec_(t_vec), x_vec_(x_vec), idx_(idx) {}

        reference operator*() const { return {(*t_vec_)[idx_], (*x_vec_)[idx_]}; }

        Iterator& operator++() {
            ++idx_;
            return *this;
        }

        Iterator operator++(int) {
            Iterator tmp = *this;
            ++idx_;
            return tmp;
        }

        bool operator==(const Iterator& other) const { return idx_ == other.idx_; }
        bool operator!=(const Iterator& other) const { return idx_ != other.idx_; }

    private:
        const std::vector<T>*             t_vec_;
        const std::vector<ColVec<NX, T>>* x_vec_;
        size_t                            idx_;
    };

    Iterator begin() const { return Iterator(&t, &x, 0); }
    Iterator end() const { return Iterator(&t, &x, t.size()); }

    /// Number of recorded time steps
    size_t size() const { return t.size(); }
};

/**
 * @brief Fixed-step ODE solver
 *
 * Wraps any integrator with a fixed time step. Supports optional callbacks
 * for step notification, early stopping, and event detection.
 *
 * @tparam NX             Number of states
 * @tparam T              Scalar type
 * @tparam IntegratorType Integrator struct (e.g., RK4<NX,T>)
 */
template<size_t NX, typename T, typename IntegratorType>
class FixedStepSolver {
public:
    explicit FixedStepSolver(IntegratorType integrator, T step_size = T(0.01))
        : integrator_(integrator), h_(step_size) {}

    /**
     * @brief Solve dx/dt = f(t, x) over [t0, tf]
     *
     * @param f      Right-hand side callable: f(t, x) -> ColVec<NX,T>
     * @param x0     Initial state
     * @param t_span Pair {t0, tf}
     * @return SolveResult with full time/state history
     */
    template<typename F>
    SolveResult<NX, T> solve(const F& f, const ColVec<NX, T>& x0, const damp::pair<T, T>& t_span) const {
        // Reset multi-step integrators if they have a reset() method
        if constexpr (requires { integrator_.reset(); }) {
            integrator_.reset();
        }

        SolveResult<NX, T> result;
        const T            t0 = t_span.first;
        const T            tf = t_span.second;
        T                  t = t0;
        ColVec<NX, T>      x = x0;

        result.t.push_back(t);
        result.x.push_back(x);
        if (!(h_ > T{0}) || !(tf >= t0)) {
            result.success = false; // forward integration with a positive step only
            return result;
        }

        // Step k lands on t0 + k·h (no running-sum drift); the last lands exactly on tf.
        const size_t n_steps = detail::step_count(t0, tf, h_);
        result.t.reserve(n_steps + 1);
        result.x.reserve(n_steps + 1);

        for (size_t k = 1; k <= n_steps; ++k) {
            const T                  t_next = detail::grid_time(t0, tf, h_, k, n_steps);
            IntegrationResult<NX, T> step_result = integrator_.evolve(f, x, t, t_next - t);
            ++result.nfev;

            t = t_next;
            x = step_result.x;
            result.t.push_back(t);
            result.x.push_back(x);

            // Step callback
            if (on_step_) {
                on_step_(t, x);
            }

            // Event detection
            if (event_detector_) {
                auto [triggered, new_state, stop] = event_detector_(t, x);
                if (triggered) {
                    x = new_state;
                    result.x.back() = x;
                }
                if (stop) {
                    break;
                }
            }

            // Stop condition
            if (stop_condition_ && stop_condition_(t, x)) {
                break;
            }
        }

        return result;
    }

    /// Set callback invoked after each accepted step
    template<typename Callback>
    void set_on_step(Callback&& cb) {
        on_step_ = damp::forward<Callback>(cb);
    }

    /// Set condition to stop simulation early: returns true to stop
    template<typename Condition>
    void set_stop_condition(Condition&& cond) {
        stop_condition_ = damp::forward<Condition>(cond);
    }

    /**
     * @brief Set event detector
     *
     * Callable returning a struct/tuple of (bool triggered, ColVec<NX,T> new_state, bool stop).
     * If triggered, state is replaced. If stop, simulation ends.
     */
    template<typename Detector>
    void set_event_detector(Detector&& det) {
        event_detector_ = damp::forward<Detector>(det);
    }

    void set_step_size(T h) { h_ = h; }
    T    get_step_size() const { return h_; }

private:
    mutable IntegratorType integrator_;
    T                      h_;

    std::function<void(T, const ColVec<NX, T>&)>                                   on_step_;
    std::function<bool(T, const ColVec<NX, T>&)>                                   stop_condition_;
    std::function<damp::tuple<bool, ColVec<NX, T>, bool>(T, const ColVec<NX, T>&)> event_detector_;
};
/// CTAD deduction guides for FixedStepSolver
template<template<size_t, typename> class Int, size_t NX, typename T>
FixedStepSolver(Int<NX, T>, T) -> FixedStepSolver<NX, T, Int<NX, T>>;

template<template<size_t, typename> class Int, size_t NX, typename T>
FixedStepSolver(Int<NX, T>) -> FixedStepSolver<NX, T, Int<NX, T>>;

/**
 * @brief Fixed-step solve in one call
 *
 * @tparam Method Integrator template (e.g. RK4, Heun, TRBDF2)
 * @param f       Right-hand side: f(t, x) → ColVec
 * @param x0      Initial state
 * @param t_span  {t0, tf}
 * @param h       [s] Step size
 */
template<template<size_t, typename> class Method, size_t NX, typename T, typename F>
[[nodiscard]] SolveResult<NX, T> fixed_solve(const F& f, const ColVec<NX, T>& x0, const damp::pair<T, T>& t_span, T h) {
    FixedStepSolver solver(Method<NX, T>{}, h);
    return solver.solve(f, x0, t_span);
}

/**
 * @brief Scalar-or-vector local-error tolerance for adaptive ODE control
 *
 * Construct from a scalar (broadcast to every state) or a ColVec for
 * per-component control. Matches SciPy `solve_ivp` accepting scalar or
 * array-valued `atol` / `rtol`.
 *
 * @tparam NX State dimension
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct ErrorTolerance {
    ColVec<NX, T> v{};

    /// Broadcast a scalar tolerance to all states
    constexpr ErrorTolerance(T s = T{0}) {
        for (size_t i = 0; i < NX; ++i) {
            v[i] = s;
        }
    }

    /// Per-component tolerances
    constexpr ErrorTolerance(const ColVec<NX, T>& c) : v(c) {}

    /// Brace list: one value broadcasts; NX values set components
    constexpr ErrorTolerance(std::initializer_list<T> il) {
        if (il.size() == 1) {
            const T s = *il.begin();
            for (size_t i = 0; i < NX; ++i) {
                v[i] = s;
            }
        } else {
            v = ColVec<NX, T>(il);
        }
    }

    [[nodiscard]] constexpr T  operator[](size_t i) const { return v[i]; }
    [[nodiscard]] constexpr T& operator[](size_t i) { return v[i]; }
};

/**
 * @brief Options for adaptive-step ODE integration
 *
 * Defaults match SciPy `solve_ivp` / `RK45`:
 * - atol = 1e-6, rtol = 1e-3 (scalar broadcast, or pass a ColVec per component)
 * - max_step unbounded (`numeric_limits::max()`)
 * - step growth factors safety=0.9, min_factor=0.2, max_factor=10
 * - first_step = 0 → Hairer / SciPy `select_initial_step` from f(t0, x0)
 * - Gustafsson PI step controller after the first accepted step
 *
 * @note Compare with scipy.integrate.solve_ivp(…, atol, rtol, first_step, max_step).
 * @see "Solving Ordinary Differential Equations I" (Hairer, Nørsett & Wanner), §II.4
 * @see Gustafsson, "Control-theoretic techniques for stepsize selection…," ACM TOMS 20(4), 1994
 *
 * @tparam NX State dimension
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct AdaptiveOptions {
    ErrorTolerance<NX, T> atol{static_cast<T>(1e-6)};               ///< Absolute local-error tolerance (scalar or ColVec)
    ErrorTolerance<NX, T> rtol{static_cast<T>(1e-3)};               ///< Relative local-error tolerance (scalar or ColVec)
    T                     first_step = T{0};                        ///< Initial h; 0 → Hairer select_initial_step
    T                     min_step = T{0};                          ///< Floor on h; 0 → 10·ε·max(|t0|, |tf|)
    T                     max_step = std::numeric_limits<T>::max(); ///< Ceiling on h (SciPy: inf)
    T                     safety = static_cast<T>(0.9);             ///< Step-size safety factor
    T                     min_factor = static_cast<T>(0.2);         ///< Min h shrink/grow factor per attempt
    T                     max_factor = T{10};                       ///< Max h grow factor (SciPy RK45)
    size_t                max_nfev = 1'000'000;                     ///< Cap on evolve() calls
    bool                  fail_on_min_step = true;                  ///< Fail if error still exceeds tol at min_step
};

/**
 * @brief Adaptive-step ODE solver
 *
 * Wraps an @ref AdaptiveStepIntegrator and adjusts h so the weighted RMS
 * local error stays ≤ 1 (SciPy / Hairer style):
 *
 *     scale_i = atol_i + rtol_i · max(|x_i|, |x⁺_i|)
 *     error_norm = rms(error_i / scale_i)
 *
 * Step-size control:
 * - First h from Hairer §II.4 / SciPy `select_initial_step` when
 *   `first_step == 0` (order-aware, uses f(t0,x0) and a trial f).
 * - After the first accepted step: Gustafsson PI controller
 *   `factor ∝ err^{-0.7/k} · err_prev^{0.4/k}` (k = error_order).
 * - On reject / first step: pure I-control `factor ∝ err^{-1/k}`.
 *
 * Zero crossings are located by bisection on the integrated trajectory (each
 * probe is one integrator step from the start of the accepted step); the
 * earliest crossing among all registered functions ends the step, at the first
 * bisection point past the sign change. Fixed-step integrators
 * (RK4, …) are rejected at compile time — their `error` is always zero and
 * would freeze h at the initial guess.
 *
 * Prefer adaptive_solve for the common case; keep this class when you
 * need events, zero-crossings, or per-step callbacks.
 *
 * | Method  | Use when                         |
 * |---------|----------------------------------|
 * | DP45    | Non-stiff default (MATLAB® ode45) |
 * | RK23    | Cheap non-stiff (ode23)          |
 * | TRBDF2  | Stiff / switched plants (ode23tb)|
 *
 * @note Compare with scipy.integrate.solve_ivp atol/rtol error control.
 * @see "Solving Ordinary Differential Equations I" (Hairer, Nørsett & Wanner), §II.4
 * @see Gustafsson, "Control-theoretic techniques for stepsize selection…," ACM TOMS 20(4), 1994
 *
 * @tparam NX             Number of states
 * @tparam T              Scalar type
 * @tparam IntegratorType Integrator modelling @ref AdaptiveStepIntegrator
 */
template<size_t NX, typename T, typename IntegratorType>
    requires AdaptiveStepIntegrator<IntegratorType, NX, T>
class AdaptiveStepSolver {
public:
    /**
     * @brief Construct with SciPy-style options (preferred)
     */
    explicit AdaptiveStepSolver(IntegratorType integrator, AdaptiveOptions<NX, T> opts = {})
        : integrator_(damp::move(integrator)), opts_(opts) {}

    /**
     * @brief Construct with first step and absolute/relative tolerances (scalar broadcast)
     */
    explicit AdaptiveStepSolver(IntegratorType integrator, T first_step, T atol, T rtol)
        : integrator_(damp::move(integrator)), opts_{.atol = atol, .rtol = rtol, .first_step = first_step} {}

    /**
     * @brief Solve dx/dt = f(t, x) over [t0, tf] with adaptive stepping
     *
     * @param f      Right-hand side callable: f(t, x) -> ColVec<NX,T>
     * @param x0     Initial state
     * @param t_span Pair {t0, tf}
     * @return SolveResult with full time/state history
     */
    template<typename F>
    SolveResult<NX, T> solve(const F& f, const ColVec<NX, T>& x0, const damp::pair<T, T>& t_span) const {
        if constexpr (requires { integrator_.reset(); }) {
            integrator_.reset();
        }

        SolveResult<NX, T> result;
        const T            t0 = t_span.first;
        const T            tf = t_span.second;
        T                  t = t0;
        ColVec<NX, T>      x = x0;

        result.t.push_back(t);
        result.x.push_back(x);
        if (!(tf >= t0)) {
            result.success = false; // forward integration only
            return result;
        }
        const T span = tf - t0;

        // Auto min-step: a few ULPs of the largest |t| in the span (relative, so the
        // floor scales with the time units).
        T h_min = (opts_.min_step > T{0}) ? opts_.min_step : T{10} * eps<T>() * damp::max(damp::abs(t0), damp::abs(tf));
        T h_max = opts_.max_step;
        // std::clamp is UB when lo > hi — swap rather than poison the controller.
        if (h_min > h_max) {
            damp::swap(h_min, h_max);
        }

        T h;
        if (opts_.first_step > T{0}) {
            h = opts_.first_step;
        } else {
            const ColVec<NX, T> f0 = f(t0, x0);
            h = select_initial_step(f, t0, x0, f0, IntegratorType::error_order, opts_.atol, opts_.rtol, span, h_max);
        }
        h = damp::clamp(h, h_min, h_max);

        result.t.reserve(1000);
        result.x.reserve(1000);

        // Record an accepted point and run the callbacks; true ⇒ stop requested.
        const auto accept = [&](T t_new, const ColVec<NX, T>& x_new) {
            t = t_new;
            x = x_new;
            result.t.push_back(t);
            result.x.push_back(x);
            ++result.n_accepted;
            if (on_step_) {
                on_step_(t, x);
            }
            if (event_detector_) {
                auto [triggered, new_state, stop] = event_detector_(t, x);
                if (triggered) {
                    x = new_state;
                    result.x.back() = x;
                }
                if (stop) {
                    return true;
                }
            }
            return stop_condition_ && stop_condition_(t, x);
        };

        constexpr size_t max_consecutive_rejections = 100;
        // Local error ~ O(h^{error_order}); error_norm is tol-scaled (target 1).
        // Pure I: h_new/h ≈ error_norm^{-1/error_order}.
        // Gustafsson PI after an accepted step: err^{-0.7/k} · err_prev^{0.4/k}.
        constexpr T inv_k = T{1} / T{IntegratorType::error_order};
        size_t      consecutive_rejections = 0;
        T           error_old = T{-1}; // < 0 ⇒ no previous accepted error (use pure I)
        bool        last_rejected = false;

        while (t < tf) {
            if (result.nfev >= opts_.max_nfev) {
                result.success = false;
                return result;
            }

            const T remaining = tf - t;
            // Propose h, never smaller than h_min while there is a full min-step
            // of room left; always finish exactly at tf on the last partial step.
            T target_step = damp::min(h, remaining);
            if (remaining > h_min && target_step < h_min) {
                target_step = h_min;
            }
            const bool last_step = !(target_step < remaining);
            if (last_step) {
                target_step = remaining;
            }
            if (!(target_step > T{0})) {
                result.success = false;
                return result;
            }
            const T t_end = last_step ? tf : t + target_step;

            IntegrationResult<NX, T> step_result = integrator_.evolve(f, x, t, target_step);
            ++result.nfev;

            const T error_norm = weighted_rms_error(x, step_result.x, step_result.error, opts_.atol, opts_.rtol);
            // A non-finite estimate (stage overflow, f undefined) is a rejection that must shrink h.
            const bool finite = damp::isfinite(error_norm);
            const bool within_tol = finite && error_norm <= T{1};
            // "At min step" only for a full h_min attempt.
            const bool at_min_step = target_step <= h_min * (T{1} + T{10} * eps<T>()) && remaining > h_min;

            if (!within_tol && at_min_step && (opts_.fail_on_min_step || !finite)) {
                result.success = false;
                return result;
            }

            const T factor = step_factor(error_norm, error_old, within_tol, last_rejected, inv_k);

            if (!(within_tol || at_min_step)) {
                // Reject: retry with a smaller step (keep error_old for the next accept).
                h = damp::clamp(target_step * factor, h_min, h_max);
                ++consecutive_rejections;
                ++result.n_rejected;
                last_rejected = true;
                if (consecutive_rejections >= max_consecutive_rejections) {
                    result.success = false;
                    return result;
                }
                continue;
            }

            // Accept (within tol, or forced at min step when fail_on_min_step is false).
            consecutive_rejections = 0;
            last_rejected = false;
            if (error_norm > T{0}) {
                error_old = error_norm;
            }
            h = damp::clamp(target_step * factor, h_min, h_max);

            const CrossingResult crossing = earliest_crossing(f, x, t, step_result.x, t_end, result.nfev);
            const bool           stop = crossing.found ? accept(crossing.time, crossing.x) : accept(t_end, step_result.x);
            if (stop) {
                return result;
            }
        }

        return result;
    }

    /// Add a zero-crossing function to monitor
    template<typename ZCF>
    void add_zero_crossing(ZCF&& zcf) {
        zero_crossings_.push_back(damp::forward<ZCF>(zcf));
    }

    template<typename Callback>
    void set_on_step(Callback&& cb) {
        on_step_ = damp::forward<Callback>(cb);
    }

    template<typename Condition>
    void set_stop_condition(Condition&& cond) {
        stop_condition_ = damp::forward<Condition>(cond);
    }

    template<typename Detector>
    void set_event_detector(Detector&& det) {
        event_detector_ = damp::forward<Detector>(det);
    }

    void                                 set_options(AdaptiveOptions<NX, T> opts) { opts_ = opts; }
    [[nodiscard]] AdaptiveOptions<NX, T> get_options() const { return opts_; }

    void set_atol(T atol) { opts_.atol = ErrorTolerance<NX, T>{atol}; }
    void set_atol(const ColVec<NX, T>& atol) { opts_.atol = ErrorTolerance<NX, T>{atol}; }
    void set_rtol(T rtol) { opts_.rtol = ErrorTolerance<NX, T>{rtol}; }
    void set_rtol(const ColVec<NX, T>& rtol) { opts_.rtol = ErrorTolerance<NX, T>{rtol}; }

    [[nodiscard]] ErrorTolerance<NX, T> get_atol() const { return opts_.atol; }
    [[nodiscard]] ErrorTolerance<NX, T> get_rtol() const { return opts_.rtol; }

private:
    struct CrossingResult {
        bool          found = false;
        T             time = T(0);
        ColVec<NX, T> x{};
    };

    /// Hairer §II.4 / SciPy select_initial_step — order-aware h0 from f(t0,y0).
    template<typename F>
    [[nodiscard]] static T select_initial_step(
        const F&                     f,
        T                            t0,
        const ColVec<NX, T>&         y0,
        const ColVec<NX, T>&         f0,
        int                          order,
        const ErrorTolerance<NX, T>& atol,
        const ErrorTolerance<NX, T>& rtol,
        T                            span,
        T                            h_max
    ) {
        if (!(span > T{0})) {
            return T{0};
        }

        T d0_sq = T{0};
        T d1_sq = T{0};
        for (size_t i = 0; i < NX; ++i) {
            const T scale = atol[i] + damp::abs(y0[i]) * rtol[i];
            const T denom = (scale > T{0}) ? scale : std::numeric_limits<T>::min();
            const T ys = y0[i] / denom;
            const T fs = f0[i] / denom;
            d0_sq += ys * ys;
            d1_sq += fs * fs;
        }
        const T d0 = damp::sqrt(d0_sq / T{NX});
        const T d1 = damp::sqrt(d1_sq / T{NX});

        T h0 = (d0 < static_cast<T>(1e-5) || d1 < static_cast<T>(1e-5)) ? static_cast<T>(1e-6) : static_cast<T>(0.01) * d0 / d1;
        h0 = damp::min({h0, span, h_max});
        if (!(h0 > T{0})) {
            return damp::min(span, h_max);
        }

        const ColVec<NX, T> y1 = y0 + h0 * f0;
        const ColVec<NX, T> f1 = f(t0 + h0, y1);

        T d2_sq = T{0};
        for (size_t i = 0; i < NX; ++i) {
            const T scale = atol[i] + damp::abs(y0[i]) * rtol[i];
            const T denom = (scale > T{0}) ? scale : std::numeric_limits<T>::min();
            const T dfs = (f1[i] - f0[i]) / denom;
            d2_sq += dfs * dfs;
        }
        const T d2 = damp::sqrt(d2_sq / T{NX}) / h0;

        T h1;
        if (d1 <= static_cast<T>(1e-15) && d2 <= static_cast<T>(1e-15)) {
            h1 = damp::max(static_cast<T>(1e-6), h0 * static_cast<T>(1e-3));
        } else {
            // Hairer: h1^p · max(d1, d2) = 0.01 with p the local-error order. SciPy
            // writes 1/(error_estimator_order + 1), which is the same p (5 for DP45).
            h1 = damp::pow(static_cast<T>(0.01) / damp::max(d1, d2), T{1} / static_cast<T>(order));
        }

        return damp::min({T{100} * h0, h1, span, h_max});
    }

    /// Gustafsson PI on accept (with prior error); pure I on reject / first step;
    /// the minimum factor after a non-finite error estimate.
    [[nodiscard]] T step_factor(T error_norm, T error_old, bool within_tol, bool last_rejected, T inv_k) const {
        if (!damp::isfinite(error_norm)) {
            return opts_.min_factor;
        }
        T factor = opts_.max_factor;
        if (error_norm > T{0}) {
            if (within_tol && error_old > T{0} && damp::isfinite(error_old)) {
                // PI: safety · err^{-0.7/k} · err_old^{0.4/k}
                factor = opts_.safety * damp::pow(error_norm, -static_cast<T>(0.7) * inv_k) * damp::pow(error_old, static_cast<T>(0.4) * inv_k);
            } else {
                // Pure I: safety · err^{-1/k}
                factor = opts_.safety * damp::pow(error_norm, -inv_k);
            }
            factor = damp::clamp(factor, opts_.min_factor, opts_.max_factor);
        }
        // After a rejection, do not grow on the next accepted step (SciPy RK).
        if (within_tol && last_rejected) {
            factor = damp::min(factor, T{1});
        }
        return factor;
    }

    /// SciPy-style weighted RMS: rms(err_i / (atol_i + rtol_i * max(|y_i|, |ynew_i|)))
    /// Scale is floored away from zero so atol=rtol=0 with a zero state cannot 0/0 → NaN.
    [[nodiscard]] static constexpr T weighted_rms_error(
        const ColVec<NX, T>&         x,
        const ColVec<NX, T>&         x_new,
        const ColVec<NX, T>&         err,
        const ErrorTolerance<NX, T>& atol,
        const ErrorTolerance<NX, T>& rtol
    ) {
        T sum_sq = T{0};
        for (size_t i = 0; i < NX; ++i) {
            const T scale = atol[i] + rtol[i] * damp::max(damp::abs(x[i]), damp::abs(x_new[i]));
            const T denom = (scale > T{0}) ? scale : std::numeric_limits<T>::min();
            const T r = err[i] / denom;
            sum_sq += r * r;
        }
        return damp::sqrt(sum_sq / T{NX});
    }

    [[nodiscard]] static constexpr bool opposite_signs(T a, T b) {
        return (a > T{0} && b < T{0}) || (a < T{0} && b > T{0});
    }

    /// Earliest sign change of any registered zero-crossing function over the step [t_a, t_b].
    template<typename F>
    [[nodiscard]] CrossingResult earliest_crossing(
        const F&             f,
        const ColVec<NX, T>& x_a,
        T                    t_a,
        const ColVec<NX, T>& x_b,
        T                    t_b,
        size_t&              nfev
    ) const {
        CrossingResult best{};
        for (const auto& zcf : zero_crossings_) {
            const T z_a = zcf(t_a, x_a);
            if (!opposite_signs(z_a, zcf(t_b, x_b))) {
                continue;
            }
            const CrossingResult c = locate_crossing(f, zcf, x_a, t_a, z_a, x_b, t_b, nfev);
            if (!best.found || c.time < best.time) {
                best = c;
            }
        }
        return best;
    }

    /// Bisection on the integrated trajectory x(τ) = evolve(f, x_a, t_a, τ − t_a).
    /// Keeps z(lo) on z_a's side and returns hi — the first point at or past the
    /// crossing — so the recorded state has crossed and is not re-detected.
    template<typename F>
    [[nodiscard]] CrossingResult locate_crossing(
        const F&                                         f,
        const std::function<T(T, const ColVec<NX, T>&)>& zcf,
        const ColVec<NX, T>&                             x_a,
        T                                                t_a,
        T                                                z_a,
        const ColVec<NX, T>&                             x_b,
        T                                                t_b,
        size_t&                                          nfev
    ) const {
        T             lo = t_a;
        T             hi = t_b;
        ColVec<NX, T> x_hi = x_b;
        // Resolve the time to a few ULPs of the step's end-point magnitude.
        const T t_tol = T{4} * eps<T>() * damp::max(damp::abs(t_a), damp::abs(t_b));
        for (int iter = 0; iter < 200 && (hi - lo) > t_tol; ++iter) {
            const T mid = lo + ((hi - lo) / T{2});
            if (!(mid > lo && mid < hi)) {
                break;
            }
            const ColVec<NX, T> x_mid = integrator_.evolve(f, x_a, t_a, mid - t_a).x;
            ++nfev;
            const T z_mid = zcf(mid, x_mid);
            if (z_mid != T{0} && !opposite_signs(z_a, z_mid)) {
                lo = mid;
            } else {
                hi = mid;
                x_hi = x_mid;
            }
        }
        return {true, hi, x_hi};
    }

    mutable IntegratorType integrator_;
    AdaptiveOptions<NX, T> opts_;

    std::function<void(T, const ColVec<NX, T>&)>                                   on_step_;
    std::function<bool(T, const ColVec<NX, T>&)>                                   stop_condition_;
    std::function<damp::tuple<bool, ColVec<NX, T>, bool>(T, const ColVec<NX, T>&)> event_detector_;
    std::vector<std::function<T(T, const ColVec<NX, T>&)>>                         zero_crossings_;
};

/// CTAD deduction guides for AdaptiveStepSolver
template<template<size_t, typename> class Int, size_t NX, typename T>
AdaptiveStepSolver(Int<NX, T>, AdaptiveOptions<NX, T>) -> AdaptiveStepSolver<NX, T, Int<NX, T>>;

template<template<size_t, typename> class Int, size_t NX, typename T>
AdaptiveStepSolver(Int<NX, T>, T, T, T) -> AdaptiveStepSolver<NX, T, Int<NX, T>>;

template<template<size_t, typename> class Int, size_t NX, typename T>
AdaptiveStepSolver(Int<NX, T>) -> AdaptiveStepSolver<NX, T, Int<NX, T>>;

/**
 * @brief Adaptive-step solve in one call
 *
 * Instantiates `Method<NX,T>`, wraps it in AdaptiveStepSolver, and
 * integrates. Method must model @ref AdaptiveStepIntegrator (DP45, RK23,
 * TRBDF2). Defaults match SciPy `solve_ivp`.
 *
 * @code
 *   auto y = adaptive_solve\<DP45\>(f, x0, {0.0, 1.0});                         // SciPy defaults
 *   auto z = adaptive_solve\<TRBDF2\>(f, x0, {0.0, 1.0}, {.atol=1e-8, .rtol=1e-6});
 *   auto w = adaptive_solve\<DP45\>(f, x0, {0.0, 1.0},
 *       AdaptiveOptions<2>{.atol = ColVec<2>{1e-12, 1e-3}, .rtol = 1e-6});
 * @endcode
 *
 * @tparam Method Integrator template (DP45, RK23, or TRBDF2)
 * @param f       Right-hand side: f(t, x) → ColVec
 * @param x0      Initial state
 * @param t_span  {t0, tf}
 * @param opts    Adaptive tolerances and step limits
 */
template<template<size_t, typename> class Method, size_t NX, typename T, typename F>
    requires AdaptiveStepIntegrator<Method<NX, T>, NX, T>
[[nodiscard]] SolveResult<NX, T> adaptive_solve(
    const F&                f,
    const ColVec<NX, T>&    x0,
    const damp::pair<T, T>& t_span,
    AdaptiveOptions<NX, T>  opts = {}
) {
    AdaptiveStepSolver solver(Method<NX, T>{}, opts);
    return solver.solve(f, x0, t_span);
}

} // namespace damp::sim
