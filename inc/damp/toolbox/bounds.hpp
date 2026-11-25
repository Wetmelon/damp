// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file bounds.hpp
 * @brief Per-channel interval limits — the one primitive for box constraints on
 *        a state, output, or input vector.
 *
 * Every "limit" in the library is the same shape: a closed interval
 * @f$[\,\text{lower}_i,\ \text{upper}_i\,]@f$ per channel. A state-limit set is a
 * `Bounds\<Nx\>`; an output/input-limit set is a `Bounds\<Nu\>` — same type, different
 * vector. The ends are column vectors, so they enter matrix expressions directly.
 * A scalar saturation pair is a `Bounds<1>`, written `{lo, hi}`: the PID output and
 * integrator limits, the ESC and relay-autotune clamps, the SOGI-FLL frequency range.
 * The vector form carries the MPC input, move and output limits and the MHE state
 * limits. The symmetric @f$|x_i|\le \text{max}_i@f$ convention is Bounds::symmetric.
 *
 * Scope is just the hard constraint. A soft/weighted constraint (MPC: a `Q`/`R`
 * cost coefficient or a slack penalty per channel) composes a `Bounds` with a
 * weight; it is not a member here.
 *
 * @see controllers/pid.hpp (Bounds<1> on output + integrator),
 *      design/input_limits.hpp (InputLimits), design/state_limits.hpp (StateLimits),
 *      estimation/mhe.hpp (MHELimits).
 */


#include "damp/backend.hpp"
#include "damp/design/options.hpp" // clamp, numeric_limits
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"

namespace damp {

/**
 * @brief A per-channel closed-interval box constraint.
 *
 * An aggregate: `{.lower = -12.0, .upper = 12.0}` for one channel (a ColVec<1> takes its
 * scalar), `{.lower = {-1.0, -2.0}, .upper = {1.0, 2.0}}` for two. An end left out stays
 * unbounded. The default is unbounded: every end at the widest finite value, ∓max(). The
 * library builds under -ffinite-math-only, where producing ±inf is undefined, so
 * the ends are never ±inf; a zero-initialised box would clamp every channel to 0.
 *
 * @tparam N Number of channels
 * @tparam T Scalar type (floating point)
 */
template<size_t N, typename T = double>
struct Bounds {
    ColVec<N, T> lower = Matrix<N, 1, T>::constant(-damp::numeric_limits<T>::max()); ///< lower(i) ≤ x(i)
    ColVec<N, T> upper = Matrix<N, 1, T>::constant(damp::numeric_limits<T>::max());  ///< x(i) ≤ upper(i)

    /** @brief Symmetric box [−mag, +mag] per channel (the |xᵢ| ≤ magᵢ convention). */
    [[nodiscard]] static constexpr Bounds symmetric(const ColVec<N, T>& mag) {
        return Bounds{
            .lower = -mag,
            .upper = mag,
        };
    }

    /**
     * @brief Per-channel clamp into the box.
     *
     * damp::clamp per channel, which has no precondition: an interval whose lower end is above its
     * upper one (which valid() refuses at design, but a member set at run time can hold) returns its
     * upper end on every backend.
     */
    [[nodiscard]] constexpr ColVec<N, T> saturate(const ColVec<N, T>& x) const {
        ColVec<N, T> out{};
        for (size_t i = 0; i < N; ++i) {
            out(i) = damp::clamp(x(i), lower(i), upper(i));
        }
        return out;
    }

    /** @brief True iff every channel is within its interval. */
    [[nodiscard]] constexpr bool contains(const ColVec<N, T>& x) const {
        for (size_t i = 0; i < N; ++i) {
            if (x(i) < lower(i) || x(i) > upper(i)) {
                return false;
            }
        }
        return true;
    }

    /** @brief True iff the box is usable: every end finite (∓max() for none) and no lower end above its upper one. */
    [[nodiscard]] constexpr bool valid() const {
        for (size_t i = 0; i < N; ++i) {
            if (!damp::isfinite(lower(i)) || !damp::isfinite(upper(i)) || !(lower(i) <= upper(i))) {
                return false;
            }
        }
        return true;
    }

    /** @brief Clamp of a scalar into the SISO interval. */
    [[nodiscard]] constexpr T saturate(T x) const
        requires(N == 1)
    {
        return damp::clamp(x, lower(0), upper(0));
    }

    /** @brief True iff the scalar is within the SISO interval. */
    [[nodiscard]] constexpr bool contains(T x) const
        requires(N == 1)
    {
        return (x >= lower(0)) && (x <= upper(0));
    }

    /**
     * @brief The same box in scalar type U. A finite end at or beyond max()/4 in T is unbounded and
     * becomes ±max() of U, and so does a finite end past the range of a narrower U: a plain cast
     * of either to float would overflow. A NaN or ±Inf end stays one, for valid() to refuse.
     */
    template<typename U>
    [[nodiscard]] constexpr Bounds<N, U> as() const {
        Bounds<N, U> out{};
        for (size_t i = 0; i < N; ++i) {
            out.lower(i) = convert_end<U>(lower(i));
            out.upper(i) = convert_end<U>(upper(i));
        }
        return out;
    }

private:
    template<typename U>
    [[nodiscard]] static constexpr U convert_end(T v) {
        if (!damp::isfinite(v)) {
            return static_cast<U>(v); // a NaN or ±Inf end stays one, for valid() to refuse, read before a comparison -ffast-math may fold
        }
        constexpr T unbounded = damp::numeric_limits<T>::max() / T{4};
        if (v >= unbounded) {
            return damp::numeric_limits<U>::max();
        }
        if (v <= -unbounded) {
            return -damp::numeric_limits<U>::max();
        }
        if constexpr (damp::numeric_limits<U>::max_exponent < damp::numeric_limits<T>::max_exponent) {
            constexpr T widest = static_cast<T>(damp::numeric_limits<U>::max()); // U is narrower: its range fits in T
            return static_cast<U>(damp::clamp(v, -widest, widest));
        }
        return static_cast<U>(v);
    }
};

/**
 * @brief A limit on what a state-feedback law's outer terms ask of one of its states: a cascade's
 *        outer limit, inside one gain
 *
 * Read the law u = u* − K·z, z = [x − x*; ξ], as a cascade on the state s: u = −K_s·(z_s − z_c) − (the
 * other terms), where the outer terms (the plant states and integrators listed in @ref outer) ask for
 * z_c = −(Σ K_j·z_j)/K_s. The law clamps z_c into [@ref lower, @ref upper], per input that acts on s (K_s ≠ 0). A
 * balancing robot's position, speed and integral ask for a lean; this bounds the lean they may ask
 * for, not the lean itself, which no law without a prediction can bound. The integrators hold
 * against the clamp as against the InputLimits', or unwind with a tracking_time. At a steady state the outer terms carry the input that holds it (none
 * for a robot balanced upright; a spring's force through its position): a limit below that holds the
 * plant short of its reference.
 *
 * On the tick the law is computed as ever, u = u* − K·z, a plain product with the gain, and the clamp
 * corrects it afterwards. For each input i that acts on s, the outer terms' share of uᵢ,
 * oᵢ = −Σ K_ij·z_j over @ref outer, is formed again, and where oᵢ/K_is falls outside its box it is
 * replaced by K_is·clamp(oᵢ/K_is). K itself never changes. The law is piecewise linear: with K split by
 * columns into K_out, the outer terms, and K_in = K − K_out, and k_s the column on s,
 * @f[
 *   u = u^* - K_{\mathrm{in}}\, z + k_s \circ \operatorname{sat}_{\mathrm{lim}}\!\left(\frac{-K_{\mathrm{out}}\, z}{k_s}\right),
 * @f]
 * entry by entry over the inputs, which inside the box is u* − K·z. Outside it only the outer terms
 * are bounded: the gain on s and the other terms keep their full stiffness. An input whose gain on s
 * is rounding (|K_is| ≤ 4ε of its row's largest) is left as it is. Off, the limit costs one
 * comparison; on, a second pass over each row of K, NU·(NX + NI) multiply-adds, and NU divisions.
 *
 * What it does not promise. The split is a reading of the gain, not part of its design: the gain is
 * optimal while the clamp does not bind, and while it binds the loop holds a saturation the Riccati
 * equation never saw (a Lur'e system), which nothing certifies. The outer terms then act through a
 * bounded demand; where they must stabilize a mode that is unstable on its own, that may not bring it
 * back. A disturbance can still carry the state itself past the box. A bound on the state with a
 * guarantee is design::mpc's StateLimits, or an ActionGovernor around the law.
 *
 * Published work has each part. A pole-placed state feedback with integral action on an unstable
 * plant has been read as the equivalent cascade of loops, there to protect the actuator with an
 * override. The pitch setpoint an outer position loop asks of a multirotor has been saturated, with an
 * anti-windup compensator synthesized on generalized sector conditions. Nested saturations clamp the
 * intermediate control of a chain of integrators, with gains chosen so that the saturated loop is
 * provably stable. The setpoint read out of one finished optimal gain and clamped per input, with the
 * integrators held against it, is not in that work. docs/demand_lim.md derives the construction, the
 * loop it makes, and a regional certificate Damp does not compute yet.
 *
 * A design takes it among its trailing options (design::state_feedback), the terms that ask listed
 * by index in @ref outer: the plant states 0 … NX − 1, then the integrators NX … NX + NI − 1 (the
 * delay states of a plant whose input arrives late are never outer).
 * `DemandLimit{.state = 2, .outer = {0, 1, 4}, .lower = -0.26, .upper = 0.26}`: the position, the
 * speed and the integrator ask for the lean, x₂, at most 0.26 either way. @ref outer is an array,
 * so the limit takes M, how many terms ask, from its list; the ends are scalars, one state's, and
 * one left out is unbounded. A state that is not the plant's, an index past the terms, one listed
 * twice, the state among its own outer terms, or a box Bounds::valid refuses (an end NaN or ±Inf,
 * or the lower above the upper) is refused with DesignStatus::invalid_parameter.
 *
 * A runtime that holds its matrices holds it with its other limits, read with `limits()` and set at
 * run time with `set_limits()` (`auto l = law.limits(); l.demand = DemandLimit{…};
 * law.set_limits(l);`, and `l.demand = {}` turns it off); the Law of a `constexpr` design holds
 * it as a constant of the type, so it costs that law nothing, and its set_limits refuses another.
 *
 * @see "Antiwindup and override control for exponentially unstable systems with actuator
 *      constraints" (Glattfelder & Schaufelberger, 2005), 16th IFAC World Congress, §4: a state
 *      feedback read as a cascade
 * @see "Fixed-dynamics antiwindup design: application to pitch-limited position control of multirotor
 *      unmanned aerial vehicles" (Ghignoni, Invernizzi & Lovera, 2021), IEEE Trans. Control Systems
 *      Technology: a saturated setpoint from an outer loop
 * @see "Global stabilization and restricted tracking for multiple integrators with bounded controls"
 *      (Teel, 1992), Systems & Control Letters 18(3): nested saturations
 *
 * @tparam M The terms that ask, listed in @ref outer
 */
template<size_t M, typename T = double>
struct DemandLimit {
    static_assert(M > 0, "DemandLimit lists the terms that ask by index: .outer = {i, ...}, at least one.");
    using option_kind = design::kind::DemandLimit; ///< Its kind of design option

    size_t state{static_cast<size_t>(-1)};         ///< s: the state asked for, an index of the plant state x (required)
    size_t outer[M]{};                             ///< The terms that ask, by index: plant states 0 … NX − 1, then the integrators (an array, so M is its list's length)
    T      lower{-damp::numeric_limits<T>::max()}; ///< The lowest they may ask of x_s, as a deviation from the reference point (the target with one; default unbounded)
    T      upper{damp::numeric_limits<T>::max()};  ///< The highest they may ask of x_s, as that deviation (default unbounded)
};

namespace detail {

/**
 * @brief A DemandLimit as a law holds it: one flag per term of the gain, the plant states then the
 *        integrators
 *
 * Made from a DemandLimit, whose list it reads; default-constructed it is off. An index past the N
 * terms, one listed twice, a state that is none of them or a box that is not one (an end NaN or
 * ±Inf, or the lower above the upper) leaves it unusable: valid() refuses it at design and in
 * set_limits, and the law reads it as off.
 */
template<size_t N, typename T = double>
class DemandTerms {
public:
    static constexpr size_t none = static_cast<size_t>(-1); ///< No state: the limit is off

    size_t       state{none}; ///< s: the state asked for (none: off)
    Bounds<1, T> lim{};       ///< The most the outer terms may ask of x_s

    /** @brief Off. */
    constexpr DemandTerms() = default;

    /** @brief The limit @p d on N terms, its outer terms flagged. Its box is tested in U, before as<T>() reads an infinite end as unbounded. */
    template<size_t M, typename U>
    constexpr DemandTerms(const DemandLimit<M, U>& d) // NOLINT(google-explicit-constructor): limits.demand = DemandLimit{…}
        : state{d.state}, lim{box(d).template as<T>()}, usable_{d.state < N && box(d).valid()} {
        for (size_t j = 0; j < M; ++j) {
            const size_t i = d.outer[j];
            if (i >= N || outer_[i]) {
                usable_ = false; // past the terms, or listed twice
            } else {
                outer_[i] = true;
            }
        }
        if (!usable_) {
            state = none; // the law reads it as off
        }
    }

    /** @brief Whether term @p j asks. */
    [[nodiscard]] constexpr bool outer(size_t j) const {
        return outer_[j];
    }

    /** @brief The same limit in scalar type U. */
    template<typename U>
    [[nodiscard]] constexpr DemandTerms<N, U> as() const {
        DemandTerms<N, U> out{};
        out.state = state;
        out.lim = lim.template as<U>();
        out.outer_ = outer_;
        out.usable_ = usable_;
        return out;
    }

    /** @brief Usable on a plant of @p NX states: off, or its list read whole, a plant state that is not one of its own outer terms, and a box (Bounds::valid). */
    template<size_t NX>
    [[nodiscard]] constexpr bool valid() const {
        return usable_ && ((state == none) || (state < NX && !outer_[state] && lim.valid()));
    }

    /** @brief Whether @p o is the same limit: both off, or the same state, outer terms and box (both valid, so their ends are numbers). */
    [[nodiscard]] constexpr bool same(const DemandTerms& o) const {
        if (state == none || o.state == none) {
            return state == o.state;
        }
        if (state != o.state || lim.lower(0) != o.lim.lower(0) || lim.upper(0) != o.lim.upper(0)) {
            return false;
        }
        for (size_t j = 0; j < N; ++j) {
            if (outer_[j] != o.outer_[j]) {
                return false;
            }
        }
        return true;
    }

private:
    template<size_t, typename>
    friend class DemandTerms;

    /** @brief The box of a DemandLimit's ends, in its own scalar. */
    template<size_t M, typename U>
    [[nodiscard]] static constexpr Bounds<1, U> box(const DemandLimit<M, U>& d) {
        return Bounds<1, U>{.lower = ColVec<1, U>{d.lower}, .upper = ColVec<1, U>{d.upper}};
    }

    damp::array<bool, N> outer_{};      ///< One flag per term
    bool                 usable_{true}; ///< The list was read whole
};

} // namespace detail

} // namespace damp
