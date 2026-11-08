// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file transforms.hpp
 * @brief Three-phase reference-frame transforms (Clarke/Park)
 */


#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"

namespace damp {

/**
 * @brief Clarke/Park and related three-phase frame transforms (abc, αβ, dq)
 *
 * @defgroup transforms Reference-Frame Transforms
 * @brief Clarke/Park and related three-phase frame transforms (abc, αβ, dq)
 *
 * Clarke, Park, and their inverses map between the reference frames used in
 * field-oriented control (FOC) and grid-tie synchronization:
 *
 *   - abc  three-phase stationary frame (@c ColVec\<3, T\>)
 *   - αβ   two-phase stationary frame (@ref AlphaBeta)
 *   - dq   rotor-synchronous rotating frame (@ref DirectQuadrature)
 *   - zero scalar common-mode / zero-sequence (@ref zero_sequence) — not a third
 *     frame type; optional on the inverse path for four-wire reconstruction
 *
 * Phase quantities live in a @c ColVec\<3, T\> so they compose with the rest of
 * the linear-algebra library; the αβ and dq pairs use named structs so the two
 * orthogonal components never get silently swapped.
 *
 * The scaling convention is part of the αβ/dq type (a template parameter), so
 * the compiler refuses to mix amplitude-invariant and power-invariant quantities
 * — e.g. forming power from a voltage and current of different conventions is a
 * compile error, and SVM accepts only amplitude-invariant volts. Transforms
 * starting from raw abc take the convention as a (defaulted) template argument;
 * everything downstream deduces it from the typed operand.
 *
 * Also provides the phasor-domain symmetrical-component (Fortescue) transform for
 * unbalanced / fault analysis.
 *
 * @see https://en.wikipedia.org/wiki/Alpha%E2%80%93beta_transformation
 * @see https://en.wikipedia.org/wiki/Direct-quadrature-zero_transformation
 */

/**
 * @brief Scaling convention for the Clarke/Park family
 * @ingroup transforms
 *
 * Selects the magnitude normalisation of the αβ/dq transforms:
 *   - AmplitudeInvariant (2/3): αβ/dq magnitude equals the peak phase
 *     amplitude. The default; three-phase power is @f$ \frac{3}{2}(v_\alpha i_\alpha
 *     + v_\beta i_\beta) @f$.
 *   - PowerInvariant (√(2/3), "Concordia"): the transform is orthonormal on the
 *     αβ plane, so power is preserved directly — @f$ p = v_\alpha i_\alpha +
 *     v_\beta i_\beta @f$ with no 3/2 factor.
 *
 * @see https://en.wikipedia.org/wiki/Alpha%E2%80%93beta_transformation
 */
enum class Convention : uint8_t {
    AmplitudeInvariant,
    PowerInvariant,
};

/**
 * @brief Direct-quadrature (rotor-frame) component pair
 * @ingroup transforms
 *
 * Behaves as the complex number @f$ d + jq @f$: abs() / arg() give its polar form,
 * conj() / the complex operator* treat it as a phasor, and the additive operators
 * act component-wise. The @p C scaling convention is carried in the type so dq
 * quantities of different conventions cannot be combined.
 */
template<typename T = float, Convention C = Convention::AmplitudeInvariant>
struct DirectQuadrature {
    T d{}; ///< d-axis (direct) component
    T q{}; ///< q-axis (quadrature) component

    /** @brief The scaling convention this quantity was produced with. */
    static constexpr Convention convention = C;

    /** @brief Phase angle atan2(q, d) [rad]. */
    [[nodiscard]] constexpr T arg() const { return damp::atan2(q, d); }

    /** @brief Magnitude √(d² + q²). */
    [[nodiscard]] constexpr T abs() const { return damp::hypot(d, q); }

    /** @brief Component-wise add. */
    constexpr DirectQuadrature& operator+=(const DirectQuadrature& o) {
        d += o.d;
        q += o.q;
        return *this;
    }

    /** @brief Component-wise subtract. */
    constexpr DirectQuadrature& operator-=(const DirectQuadrature& o) {
        d -= o.d;
        q -= o.q;
        return *this;
    }

    /** @brief Scale both components by @p s. */
    constexpr DirectQuadrature& operator*=(T s) {
        d *= s;
        q *= s;
        return *this;
    }

    /** @brief Divide both components by @p s. */
    constexpr DirectQuadrature& operator/=(T s) {
        d /= s;
        q /= s;
        return *this;
    }

    /** @brief Negation (d → −d, q → −q). */
    [[nodiscard]] constexpr DirectQuadrature operator-() const { return {-d, -q}; }

    /** @brief Complex conjugate d − jq. */
    [[nodiscard]] constexpr DirectQuadrature conj() const { return {d, -q}; }

    /** @brief Component-wise sum. */
    [[nodiscard]] friend constexpr DirectQuadrature operator+(DirectQuadrature a, const DirectQuadrature& b) { return a += b; }
    /** @brief Component-wise difference. */
    [[nodiscard]] friend constexpr DirectQuadrature operator-(DirectQuadrature a, const DirectQuadrature& b) { return a -= b; }
    /** @brief Scale by @p s. */
    [[nodiscard]] friend constexpr DirectQuadrature operator*(DirectQuadrature v, T s) { return v *= s; }
    /** @brief Scale by @p s. */
    [[nodiscard]] friend constexpr DirectQuadrature operator*(T s, DirectQuadrature v) { return v *= s; }
    /** @brief Divide by @p s. */
    [[nodiscard]] friend constexpr DirectQuadrature operator/(DirectQuadrature v, T s) { return v /= s; }

    /** @brief Complex product, treating both operands as d + jq. */
    [[nodiscard]] friend constexpr DirectQuadrature operator*(const DirectQuadrature& x, const DirectQuadrature& y) {
        return {
            (x.d * y.d) - (x.q * y.q),
            (x.d * y.q) + (x.q * y.d),
        };
    }

    /** @brief Exact component-wise equality. */
    [[nodiscard]] friend constexpr bool operator==(const DirectQuadrature&, const DirectQuadrature&) = default;
};

/**
 * @brief Alpha-beta (stationary-frame) component pair
 * @ingroup transforms
 *
 * Behaves as the complex number @f$ \alpha + j\beta @f$: abs() / arg() give its
 * polar form, conj() / the complex operator* treat it as a phasor, and the
 * additive operators act component-wise. The @p C scaling convention is carried in
 * the type so αβ quantities of different conventions cannot be combined.
 */
template<typename T = float, Convention C = Convention::AmplitudeInvariant>
struct AlphaBeta {
    T alpha{}; ///< α-axis component (aligned with phase a)
    T beta{};  ///< β-axis component (leads α by 90°)

    /** @brief The scaling convention this quantity was produced with. */
    static constexpr Convention convention = C;

    /** @brief Phase angle atan2(β, α) [rad]. */
    [[nodiscard]] constexpr T arg() const { return damp::atan2(beta, alpha); }

    /** @brief Magnitude √(α² + β²). */
    [[nodiscard]] constexpr T abs() const { return damp::hypot(alpha, beta); }

    /** @brief Component-wise add. */
    constexpr AlphaBeta& operator+=(const AlphaBeta& o) {
        alpha += o.alpha;
        beta += o.beta;
        return *this;
    }
    /** @brief Component-wise subtract. */
    constexpr AlphaBeta& operator-=(const AlphaBeta& o) {
        alpha -= o.alpha;
        beta -= o.beta;
        return *this;
    }
    /** @brief Scale both components by @p s. */
    constexpr AlphaBeta& operator*=(T s) {
        alpha *= s;
        beta *= s;
        return *this;
    }
    /** @brief Divide both components by @p s. */
    constexpr AlphaBeta& operator/=(T s) {
        alpha /= s;
        beta /= s;
        return *this;
    }

    /** @brief Negation (α → −α, β → −β). */
    [[nodiscard]] constexpr AlphaBeta operator-() const { return {-alpha, -beta}; }

    /** @brief Complex conjugate α − jβ. */
    [[nodiscard]] constexpr AlphaBeta conj() const { return {alpha, -beta}; }

    /** @brief Component-wise sum. */
    [[nodiscard]] friend constexpr AlphaBeta operator+(AlphaBeta a, const AlphaBeta& b) { return a += b; }
    /** @brief Component-wise difference. */
    [[nodiscard]] friend constexpr AlphaBeta operator-(AlphaBeta a, const AlphaBeta& b) { return a -= b; }
    /** @brief Scale by @p s. */
    [[nodiscard]] friend constexpr AlphaBeta operator*(AlphaBeta v, T s) { return v *= s; }
    /** @brief Scale by @p s. */
    [[nodiscard]] friend constexpr AlphaBeta operator*(T s, AlphaBeta v) { return v *= s; }
    /** @brief Divide by @p s. */
    [[nodiscard]] friend constexpr AlphaBeta operator/(AlphaBeta v, T s) { return v /= s; }

    /** @brief Complex product, treating both operands as α + jβ. */
    [[nodiscard]] friend constexpr AlphaBeta operator*(const AlphaBeta& x, const AlphaBeta& y) {
        return {
            (x.alpha * y.alpha) - (x.beta * y.beta),
            (x.alpha * y.beta) + (x.beta * y.alpha),
        };
    }

    /** @brief Exact component-wise equality. */
    [[nodiscard]] friend constexpr bool operator==(const AlphaBeta&, const AlphaBeta&) = default;
};

/**
 * @brief Zero-sequence (common-mode) scalar of a three-phase set
 * @ingroup transforms
 *
 * Projection onto the all-ones axis. Amplitude-invariant (default):
 * @f$ v_0 = (a + b + c) / 3 @f$. Power-invariant (Concordia):
 * @f$ v_0 = (a + b + c) / \sqrt{3} @f$.
 *
 * Equal bias on all three phases lands entirely here and cancels from αβ;
 * independent per-phase offsets still leak into αβ. For four-wire current
 * control, regulate this scalar next to @ref DirectQuadrature (Park does not
 * rotate it). Modulation zero-sequence injection is a separate scalar on the
 * pole commands — not this measurement.
 *
 * @tparam C  Scaling convention (default amplitude-invariant)
 * @param abc Phase quantities {a, b, c}
 * @return Zero-sequence component under convention @p C
 */
template<typename T = float, Convention C = Convention::AmplitudeInvariant>
[[nodiscard]] constexpr T zero_sequence(const ColVec<3, T>& abc) {
    if constexpr (C == Convention::PowerInvariant) {
        return (abc[0] + abc[1] + abc[2]) * damp::numbers::inv_sqrt3_v<T>;
    } else {
        return (abc[0] + abc[1] + abc[2]) / T{3};
    }
}

/**
 * @brief Clarke transform (abc → αβ)
 * @ingroup transforms
 *
 * Rank-2 projection orthogonal to common mode. Amplitude-invariant by default:
 * @f[
 *   \alpha = \frac{2a - b - c}{3}, \qquad \beta = \frac{b - c}{\sqrt{3}}
 * @f]
 * Power-invariant (Concordia) uses orthonormal √(2/3) scaling on the αβ plane.
 * Common-mode content cancels from αβ; use @ref zero_sequence when the residual
 * is needed (four-wire, diagnostics).
 *
 * @note Compare with MATLAB®'s Clarke Transform block (Motor Control Blockset).
 *
 * @tparam C  Scaling convention (default amplitude-invariant)
 * @param abc Phase quantities {a, b, c}
 * @return αβ components tagged with convention @p C
 */
template<typename T = float, Convention C = Convention::AmplitudeInvariant>
[[nodiscard]] constexpr AlphaBeta<T, C> clarke_transform(const ColVec<3, T>& abc) {
    if constexpr (C == Convention::PowerInvariant) {
        constexpr T sqrt_2_3 = damp::numbers::sqrt2_v<T> * damp::numbers::inv_sqrt3_v<T>;
        return {
            .alpha = sqrt_2_3 * (abc[0] - ((abc[1] + abc[2]) / T{2})),
            .beta = (abc[1] - abc[2]) * damp::numbers::inv_sqrt2_v<T>,
        };
    } else {
        return {
            .alpha = ((T{2} * abc[0]) - abc[1] - abc[2]) / T{3},
            .beta = (abc[1] - abc[2]) * damp::numbers::inv_sqrt3_v<T>,
        };
    }
}

/**
 * @brief Inverse Clarke transform (αβ → abc)
 * @ingroup transforms
 *
 * Reconstructs three phases from αβ and an optional zero-sequence scalar
 * (default 0 — three-wire). Convention is taken from the αβ type.
 * Amplitude-invariant adds @p zero to every phase; power-invariant uses the
 * transpose of the orthonormal forward matrix.
 *
 * @note Compare with MATLAB®'s Inverse Clarke Transform block (Motor Control Blockset).
 *
 * @param ab   αβ components
 * @param zero Zero-sequence under the same convention as @p ab (default 0)
 * @return Phase quantities {a, b, c}
 */
template<typename T = float, Convention C = Convention::AmplitudeInvariant>
[[nodiscard]] constexpr ColVec<3, T> inverse_clarke_transform(const AlphaBeta<T, C>& ab, T zero = T{0}) {
    if constexpr (C == Convention::PowerInvariant) {
        constexpr T sqrt_2_3 = damp::numbers::sqrt2_v<T> * damp::numbers::inv_sqrt3_v<T>;
        constexpr T inv_sqrt2 = damp::numbers::inv_sqrt2_v<T>;
        constexpr T inv_sqrt6 = inv_sqrt2 * damp::numbers::inv_sqrt3_v<T>;
        const T     z = zero * damp::numbers::inv_sqrt3_v<T>;
        return {
            (sqrt_2_3 * ab.alpha) + z,
            (-inv_sqrt6 * ab.alpha) + (inv_sqrt2 * ab.beta) + z,
            (-inv_sqrt6 * ab.alpha) - (inv_sqrt2 * ab.beta) + z,
        };
    } else {
        const T half_sqrt3_beta = damp::numbers::sqrt3_v<T> * ab.beta / T{2};
        return {
            ab.alpha + zero,
            (-ab.alpha / T{2}) + half_sqrt3_beta + zero,
            (-ab.alpha / T{2}) - half_sqrt3_beta + zero,
        };
    }
}

/**
 * @brief Park transform (αβ → dq)
 * @ingroup transforms
 *
 * Rotates the stationary αβ frame into the rotor-synchronous dq frame (the
 * rotation is independent of the amplitude/power scaling, so the scaling
 * convention passes through on the type):
 * @f[
 *   d =  \alpha\cos\theta + \beta\sin\theta, \qquad
 *   q = -\alpha\sin\theta + \beta\cos\theta
 * @f]
 *
 * @p theta is the electrical angle of the d-axis in the αβ plane: a pure
 * αβ space vector at angle @p theta maps to `(d, q) = (|v|, 0)`. The q-axis
 * leads d by +90° electrical (standard FOC / IEEE d–q). With
 * amplitude-invariant scaling, a constant positive @f$ q @f$ (and @f$ d = 0 @f$)
 * produces a positive-sequence abc set that leads the d-axis by 90° (for
 * θ = ωt, phase a = q·cos(ωt + 90°)).
 *
 * @note Compare with MATLAB®'s Park Transform block (Motor Control Blockset).
 *
 * @param ab    αβ components
 * @param theta Rotor / d-axis electrical angle [rad]
 * @return dq components (same convention as @p ab)
 */
template<typename T = float, Convention C = Convention::AmplitudeInvariant>
[[nodiscard]] constexpr DirectQuadrature<T, C> park_transform(const AlphaBeta<T, C>& ab, T theta) {
    const auto [s, c] = damp::sincos(theta);

    return {
        .d = (ab.alpha * c) + (ab.beta * s),
        .q = (-ab.alpha * s) + (ab.beta * c),
    };
}

/**
 * @brief Inverse Park transform (dq → αβ)
 * @ingroup transforms
 *
 * Exact inverse of @ref park_transform under the same d-axis / +90° q-lead
 * convention:
 * @f[
 *   \alpha = d\cos\theta - q\sin\theta, \qquad
 *   \beta  = d\sin\theta + q\cos\theta
 * @f]
 *
 * @note Compare with MATLAB®'s Inverse Park Transform block (Motor Control Blockset).
 *
 * @param dq    dq components
 * @param theta Rotor / d-axis electrical angle [rad]
 * @return αβ components (same convention as @p dq)
 */
template<typename T = float, Convention C = Convention::AmplitudeInvariant>
[[nodiscard]] constexpr AlphaBeta<T, C> inverse_park_transform(const DirectQuadrature<T, C>& dq, T theta) {
    const auto [s, c] = damp::sincos(theta);

    return {
        .alpha = (dq.d * c) - (dq.q * s),
        .beta = (dq.d * s) + (dq.q * c),
    };
}

/**
 * @brief Fused Clarke-Park transform (abc → dq)
 * @ingroup transforms
 *
 * Clarke then Park: maps three-phase stationary quantities to the rotor frame;
 * the usual measurement-side step of an FOC loop.
 *
 * @tparam C  Scaling convention (default amplitude-invariant)
 * @param abc   Phase quantities {a, b, c}
 * @param theta Rotor electrical angle [rad]
 * @return dq components tagged with convention @p C
 */
template<typename T = float, Convention C = Convention::AmplitudeInvariant>
[[nodiscard]] constexpr DirectQuadrature<T, C> clarke_park_transform(const ColVec<3, T>& abc, T theta) {
    return park_transform(clarke_transform<T, C>(abc), theta);
}

/**
 * @brief Fused inverse Park-Clarke transform (dq → abc)
 * @ingroup transforms
 *
 * Inverse Park then inverse Clarke: maps rotor-frame quantities to three
 * phases; the usual command-side step of an FOC loop. Convention taken
 * from the input type. Optional @p zero is the same scalar as
 * @ref inverse_clarke_transform (default 0).
 *
 * @param dq    dq components
 * @param theta Rotor electrical angle [rad]
 * @param zero  Zero-sequence under the same convention as @p dq (default 0)
 * @return Phase quantities {a, b, c}
 */
template<typename T = float, Convention C = Convention::AmplitudeInvariant>
[[nodiscard]] constexpr ColVec<3, T> inverse_park_clarke_transform(
    const DirectQuadrature<T, C>& dq, T theta, T zero = T{0}
) {
    return inverse_clarke_transform(inverse_park_transform(dq, theta), zero);
}

/**
 * @brief Instantaneous active and reactive power
 * @ingroup transforms
 *
 * Three-phase instantaneous power by the Akagi p–q theory. From the αβ voltage
 * and current vectors:
 * @f[
 *   p = k\,(v_\alpha i_\alpha + v_\beta i_\beta), \qquad
 *   q = k\,(v_\beta i_\alpha - v_\alpha i_\beta)
 * @f]
 * @f$ p @f$ is the real power transferred; @f$ q @f$ is the (Akagi) imaginary
 * power that circulates between phases without net transfer. The scale factor
 * @f$ k @f$ is fixed by the operands' convention — @f$ 3/2 @f$ for
 * amplitude-invariant (so the result is true three-phase watts/VAr), @f$ 1 @f$
 * for power-invariant — and voltage and current must share that convention (the
 * type system enforces it).
 *
 * @f$ q @f$ has the complex-power sign @f$ \operatorname{Im}(V I^*) @f$, positive
 * for an inductive load. The 1984 paper defines @f$ q @f$ with the opposite sign;
 * the sign here is the one used by Akagi, Watanabe and Aredes (2007).
 *
 * @see "Instantaneous Power Theory and Applications to Power Conditioning"
 *      (Akagi, Watanabe, Aredes, 2007)
 * @see H. Akagi, Y. Kanazawa, A. Nabae, "Instantaneous reactive power
 *      compensators comprising switching devices without energy storage
 *      components," IEEE Trans. Ind. Appl., vol. IA-20, no. 3, pp. 625-630, 1984.
 */
template<typename T = float>
struct InstantaneousPower {
    T p = {}; ///< [W]   instantaneous active (real) power
    T q = {}; ///< [VAr] instantaneous reactive (Akagi imaginary) power

    /** @brief Apparent power S = √(p² + q²) [VA] — the complex-power magnitude. */
    [[nodiscard]] constexpr T abs() const { return damp::hypot(p, q); }

    /** @brief Power-factor angle φ = atan2(q, p) [rad]. */
    [[nodiscard]] constexpr T arg() const { return damp::atan2(q, p); }

    /** @brief Power factor cos φ = p / S in [-1, 1] (0 when S = 0). */
    [[nodiscard]] constexpr T power_factor() const {
        const T s = abs();
        return s > T{0} ? p / s : T{0};
    }
};

/** @copydoc InstantaneousPower */
template<typename T = float, Convention C = Convention::AmplitudeInvariant>
[[nodiscard]] constexpr InstantaneousPower<T> instantaneous_power(const AlphaBeta<T, C>& v, const AlphaBeta<T, C>& i) {
    constexpr T k = (C == Convention::PowerInvariant) ? T{1} : static_cast<T>(1.5);
    // Complex power S = k · V · conj(I); real part is P, imaginary part is Q.
    const AlphaBeta<T, C> s = v * i.conj();
    return {.p = k * s.alpha, .q = k * s.beta};
}

/**
 * @brief Instantaneous active and reactive power from dq quantities
 * @ingroup transforms
 *
 * The rotor/synchronous-frame form of instantaneous_power():
 * @f[
 *   p = k\,(v_d i_d + v_q i_q), \qquad q = k\,(v_q i_d - v_d i_q)
 * @f]
 * Equivalent to the αβ form (both are frame-invariant scalars); use whichever
 * frame the signals are already in. Convention (hence @f$ k @f$) comes from the
 * operands' type.
 */
template<typename T = float, Convention C = Convention::AmplitudeInvariant>
[[nodiscard]] constexpr InstantaneousPower<T> instantaneous_power(const DirectQuadrature<T, C>& v, const DirectQuadrature<T, C>& i) {
    constexpr T k = (C == Convention::PowerInvariant) ? T{1} : static_cast<T>(1.5);
    // Complex power S = k · V · conj(I); real part is P, imaginary part is Q.
    const DirectQuadrature<T, C> s = v * i.conj();
    return {.p = k * s.d, .q = k * s.q};
}

/**
 * @brief Symmetrical (sequence) components of a three-phase phasor set
 * @ingroup transforms
 *
 * The Fortescue decomposition resolves an unbalanced set of three complex phasors
 * into three balanced sets: zero, positive, and negative sequence. Used for
 * unbalanced-grid / fault analysis and for sequence-domain current control.
 *
 * @tparam T Scalar type
 */
template<typename T = float>
struct SequenceComponents {
    complex<T> zero = {};     ///< Zero-sequence phasor V₀ (co-phasal)
    complex<T> positive = {}; ///< Positive-sequence phasor V₁ (abc rotation)
    complex<T> negative = {}; ///< Negative-sequence phasor V₂ (acb rotation)
};

/**
 * @brief Forward symmetrical-component (Fortescue) transform (abc → 012)
 * @ingroup transforms
 *
 * Decomposes three complex phase phasors into zero/positive/negative sequence
 * with the operator @f$ a = e^{j2\pi/3} @f$:
 * @f[
 *   \begin{bmatrix} V_0 \\ V_1 \\ V_2 \end{bmatrix}
 *   = \frac{1}{3}
 *   \begin{bmatrix} 1 & 1 & 1 \\ 1 & a & a^2 \\ 1 & a^2 & a \end{bmatrix}
 *   \begin{bmatrix} V_a \\ V_b \\ V_c \end{bmatrix}
 * @f]
 *
 * @param abc Three-phase phasors {Va, Vb, Vc}
 * @return Sequence components {zero, positive, negative}
 *
 * @see https://en.wikipedia.org/wiki/Symmetrical_components
 * @see C. L. Fortescue, "Method of symmetrical co-ordinates applied to the
 *      solution of polyphase networks," Trans. AIEE, vol. 37, no. 2, pp.
 *      1027-1140, 1918.
 */
template<typename T = float>
[[nodiscard]] constexpr SequenceComponents<T> symmetrical_components(const ColVec<3, complex<T>>& abc) {
    // Rotation operator a = e^{j120°} and its square a² = e^{j240°} = e^{-j120°}.
    constexpr complex<T> a = {static_cast<T>(-0.5), damp::numbers::sqrt3_v<T> / T{2}};
    constexpr complex<T> a2 = {static_cast<T>(-0.5), -damp::numbers::sqrt3_v<T> / T{2}};

    return {
        .zero = (abc[0] + abc[1] + abc[2]) / T{3},
        .positive = (abc[0] + (a * abc[1]) + (a2 * abc[2])) / T{3},
        .negative = (abc[0] + (a2 * abc[1]) + (a * abc[2])) / T{3},
    };
}

/**
 * @brief Inverse symmetrical-component transform (012 → abc)
 * @ingroup transforms
 *
 * Recombines sequence components back into phase phasors:
 * @f[
 *   \begin{bmatrix} V_a \\ V_b \\ V_c \end{bmatrix}
 *   = \begin{bmatrix} 1 & 1 & 1 \\ 1 & a^2 & a \\ 1 & a & a^2 \end{bmatrix}
 *   \begin{bmatrix} V_0 \\ V_1 \\ V_2 \end{bmatrix}
 * @f]
 *
 * @param s Sequence components {zero, positive, negative}
 * @return Three-phase phasors {Va, Vb, Vc}
 *
 * @see https://en.wikipedia.org/wiki/Symmetrical_components
 */
template<typename T = float>
[[nodiscard]] constexpr ColVec<3, complex<T>> inverse_symmetrical_components(const SequenceComponents<T>& s) {
    constexpr complex<T> a = {static_cast<T>(-0.5), damp::numbers::sqrt3_v<T> / T{2}};
    constexpr complex<T> a2 = {static_cast<T>(-0.5), -damp::numbers::sqrt3_v<T> / T{2}};

    return {
        s.zero + s.positive + s.negative,
        s.zero + (a2 * s.positive) + (a * s.negative),
        s.zero + (a * s.positive) + (a2 * s.negative),
    };
}

} // namespace damp
