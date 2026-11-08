// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file biquad.hpp
 * @brief Second-order sections: design::second_order_sections, lowpass_sections and
 *        highpass_sections → SecondOrderSectionsResult → BiquadCascade, and the Biquad section
 */


#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/filters/iir_design.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"

namespace damp {
namespace design {

/**
 * @brief An IIR filter as a cascade of second-order sections: the result BiquadCascade runs
 *
 * @f[
 *   H(z) = \prod_{k=1}^{NS} \frac{b_{0k} + b_{1k} z^{-1} + b_{2k} z^{-2}}{1 + a_{1k} z^{-1} + a_{2k} z^{-2}}
 * @f]
 * with the sections in cascade order and the overall gain in sections[0]. A first-order section
 * has b₂ = a₂ = 0; a section past the filter's order is the identity {1, 0, 0, 0, 0}.
 *
 * @tparam NS Number of sections
 * @tparam T  Scalar type
 */
template<size_t NS, typename T = double>
struct [[nodiscard]] SecondOrderSectionsResult {
    damp::array<SecondOrderCoeffs<T>, NS> sections{};                          ///< In cascade order, the gain in sections[0]
    bool                                  success{false};                      ///< status == DesignStatus::ok
    DesignStatus                          status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same sections in scalar U. */
    template<typename U>
    [[nodiscard]] constexpr SecondOrderSectionsResult<NS, U> as() const {
        SecondOrderSectionsResult<NS, U> out{};
        for (size_t i = 0; i < NS; ++i) {
            out.sections[i] = sections[i].template as<U>();
        }
        out.success = success;
        out.status = status;
        return out;
    }
};

namespace detail {

/** @brief One section's roots: none, one real root, or a pair (complex conjugate, or two reals). */
template<typename T>
struct SectionRoots {
    damp::array<damp::complex<T>, 2> r{};  ///< The roots, r[0 .. n)
    size_t                           n{0}; ///< 0, 1 or 2
};

/** @brief Up to M finite roots of a polynomial and its first nonzero coefficient. */
template<size_t M, typename T>
struct FiniteRoots {
    damp::array<damp::complex<T>, M> roots{};       ///< roots[0 .. count)
    size_t                           count{0};      ///< The degree once the leading zeros are dropped
    T                                lead{0};       ///< First nonzero coefficient; 0 for the zero polynomial
    bool                             success{true}; ///< The eigenvalue solve converged
};

/**
 * @brief Roots of @p c (highest power first), its leading exact zeros dropped: in z they are roots
 *        at infinity, which a section realizes as a delay.
 */
template<size_t N, typename T>
    requires(N >= 1)
[[nodiscard]] constexpr FiniteRoots<N - 1, T> finite_roots(const damp::array<T, N>& c) {
    FiniteRoots<N - 1, T> out{};
    if (c[0] != T{0}) {
        const auto r = poly_roots(c);
        out.roots = r.roots;
        out.count = N - 1;
        out.lead = c[0];
        out.success = r.success;
    } else if constexpr (N > 1) {
        damp::array<T, N - 1> tail{};
        for (size_t i = 0; i + 1 < N; ++i) {
            tail[i] = c[i + 1];
        }
        const auto inner = finite_roots(tail);
        for (size_t i = 0; i < inner.count; ++i) {
            out.roots[i] = inner.roots[i];
        }
        out.count = inner.count;
        out.lead = inner.lead;
        out.success = inner.success;
    }
    return out;
}

/** @brief Roots grouped into sections' worth: conjugate pairs, then reals two by two. */
template<size_t M, typename T>
struct RootGroups {
    damp::array<SectionRoots<T>, M> g{};           ///< g[0 .. count)
    size_t                          count{0};      ///< Number of groups
    bool                            success{true}; ///< Every complex root found its conjugate
};

/**
 * @brief Group @p n roots as zp2sos does: each complex root with its conjugate, and the real roots
 *        two by two, each with the one closest to it in absolute value.
 *
 * A root is real when |Im r| ≤ 100·ε·|r| (MATLAB®'s cplxpair tolerance). An odd real root is
 * left alone, the smallest in absolute value.
 */
template<size_t M, typename T>
[[nodiscard]] constexpr RootGroups<M, T> group_roots(const damp::array<damp::complex<T>, M>& r, size_t n) {
    RootGroups<M, T>     out{};
    damp::array<bool, M> used{};
    damp::array<T, M>    reals{};
    size_t               nr = 0;
    for (size_t i = 0; i < n; ++i) {
        if (damp::abs(r[i].imag()) <= T{100} * eps<T>() * damp::abs(r[i])) {
            reals[nr++] = r[i].real();
            used[i] = true;
        }
    }
    for (size_t i = 0; i < n; ++i) {
        if (used[i] || r[i].imag() < T{0}) {
            continue;
        }
        size_t best = M;
        T      best_d{0};
        for (size_t j = 0; j < n; ++j) {
            if (used[j] || !(r[j].imag() < T{0})) {
                continue;
            }
            const T d = damp::abs(r[j] - damp::conj(r[i]));
            if (best == M || d < best_d) {
                best = j;
                best_d = d;
            }
        }
        if (best == M) {
            out.success = false;
            return out;
        }
        used[i] = true;
        used[best] = true;
        const T re = (r[i].real() + r[best].real()) / T{2};
        const T im = (r[i].imag() - r[best].imag()) / T{2};
        out.g[out.count++] = SectionRoots<T>{{damp::complex<T>{re, im}, damp::complex<T>{re, -im}}, 2};
    }
    for (size_t i = 0; i < n; ++i) {
        if (!used[i]) {
            out.success = false;
            return out;
        }
    }
    if constexpr (M >= 2) {               // one root has nothing to sort or pair
        for (size_t i = 1; i < nr; ++i) { // by |r|, largest first
            const T v = reals[i];
            size_t  j = i;
            while (j > 0 && damp::abs(reals[j - 1]) < damp::abs(v)) {
                reals[j] = reals[j - 1];
                --j;
            }
            reals[j] = v;
        }
        for (size_t i = 0; i + 1 < nr; i += 2) {
            out.g[out.count++] = SectionRoots<T>{{damp::complex<T>{reals[i], T{0}}, damp::complex<T>{reals[i + 1], T{0}}}, 2};
        }
    }
    if (nr % 2 == 1) {
        out.g[out.count++] = SectionRoots<T>{{damp::complex<T>{reals[nr - 1], T{0}}, damp::complex<T>{}}, 1};
    }
    return out;
}

/** @brief Distance of a group's roots from the unit circle: min |1 − |r||. */
template<typename T>
[[nodiscard]] constexpr T circle_distance(const SectionRoots<T>& s) {
    T d = damp::abs(T{1} - damp::abs(s.r[0]));
    if (s.n == 2) {
        const T d1 = damp::abs(T{1} - damp::abs(s.r[1]));
        d = (d1 < d) ? d1 : d;
    }
    return d;
}

/** @brief Closest distance between a root of @p a and a root of @p b. */
template<typename T>
[[nodiscard]] constexpr T group_distance(const SectionRoots<T>& a, const SectionRoots<T>& b) {
    T    d{0};
    bool first = true;
    for (size_t i = 0; i < a.n; ++i) {
        for (size_t j = 0; j < b.n; ++j) {
            const T dij = damp::abs(a.r[i] - b.r[j]);
            if (first || dij < d) {
                d = dij;
                first = false;
            }
        }
    }
    return d;
}

/** @brief The monic polynomial in z of a group's roots, highest power first: {1, c₁, c₂}. */
template<typename T>
[[nodiscard]] constexpr damp::array<T, 3> monic(const SectionRoots<T>& s) {
    damp::array<T, 3> c{T{1}, T{0}, T{0}};
    if (s.n == 1) {
        c[1] = -s.r[0].real();
    } else if (s.n == 2) {
        c[1] = -(s.r[0] + s.r[1]).real();
        c[2] = (s.r[0] * s.r[1]).real();
    }
    return c;
}

/**
 * @brief The section Π(z − qᵢ)/Π(z − pⱼ) in z⁻¹: a section with fewer zeros than poles carries the
 *        difference as a delay, b = {0, …}.
 */
template<typename T>
[[nodiscard]] constexpr SecondOrderCoeffs<T> section_of(const SectionRoots<T>& poles, const SectionRoots<T>& zeros) {
    const auto        a = monic(poles);
    const auto        q = monic(zeros);
    damp::array<T, 3> b{};
    const size_t      shift = poles.n - zeros.n;
    for (size_t i = 0; i <= zeros.n; ++i) {
        b[i + shift] = q[i];
    }
    return SecondOrderCoeffs<T>{b[0], b[1], b[2], a[1], a[2]};
}

/**
 * @brief Pair zeros with poles and order the sections as zp2sos does by default ('up', scale 'none').
 *
 * The poles closest to the unit circle take the zeros closest to them first; the sections run
 * from the poles farthest from the unit circle to the closest, and the gain is in the first.
 */
template<size_t NS, size_t MZ, size_t MP, typename T>
[[nodiscard]] constexpr SecondOrderSectionsResult<NS, T> sections_from_roots(
    const damp::array<damp::complex<T>, MZ>& zeros,
    size_t                                   nz,
    const damp::array<damp::complex<T>, MP>& poles,
    size_t                                   np,
    T                                        gain
) {
    using Result = SecondOrderSectionsResult<NS, T>;
    if (nz > np || !damp::isfinite(gain)) { // more zeros than poles: not causal in z
        return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter);
    }
    const auto zg = group_roots(zeros, nz);
    const auto pg = group_roots(poles, np);
    if (!zg.success || !pg.success) {
        return ::damp::detail::refused<Result>(DesignStatus::not_converged);
    }

    const size_t n_zg = (zg.count < MZ) ? zg.count : MZ;
    const size_t n_pg = (pg.count < MP) ? pg.count : MP;

    damp::array<size_t, MP> order{}; // pole groups, closest to the unit circle first
    if constexpr (MP > 1) {          // one group has nothing to sort
        for (size_t i = 0; i < n_pg; ++i) {
            size_t j = i;
            while (j > 0 && circle_distance(pg.g[i]) < circle_distance(pg.g[order[j - 1]])) {
                order[j] = order[j - 1];
                --j;
            }
            order[j] = i;
        }
    }

    size_t z_pairs = 0; // zero pairs and pole pairs not yet placed: a zero pair needs a pole pair
    size_t p_pairs = 0;
    for (size_t i = 0; i < n_zg; ++i) {
        if (zg.g[i].n == 2) {
            ++z_pairs;
        }
    }
    for (size_t i = 0; i < n_pg; ++i) {
        if (pg.g[i].n == 2) {
            ++p_pairs;
        }
    }

    Result                out{};
    damp::array<bool, MZ> taken{};
    for (size_t k = 0; k < n_pg; ++k) {
        const auto& p = pg.g[order[k]];
        const bool  needs_pair = (p.n == 2) && (z_pairs > 0) && (z_pairs >= p_pairs);
        size_t      best = MZ;
        T           best_d{0};
        for (size_t j = 0; j < n_zg; ++j) {
            const auto& q = zg.g[j];
            if (taken[j] || q.n > p.n || (needs_pair && q.n != 2)) {
                continue;
            }
            const T d = group_distance(q, p);
            if (best == MZ || d < best_d) {
                best = j;
                best_d = d;
            }
        }
        SectionRoots<T> q{};
        if (best != MZ) {
            taken[best] = true;
            q = zg.g[best];
            z_pairs -= (q.n == 2) ? 1 : 0;
        }
        p_pairs -= (p.n == 2) ? 1 : 0;
        out.sections[n_pg - 1 - k] = section_of(p, q);
    }
    for (size_t j = 0; j < n_zg; ++j) {
        if (!taken[j]) {
            return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter);
        }
    }
    for (size_t k = n_pg; k < NS; ++k) {
        out.sections[k] = SecondOrderCoeffs<T>{T{1}, T{0}, T{0}, T{0}, T{0}};
    }
    out.sections[0].b0 *= gain;
    out.sections[0].b1 *= gain;
    out.sections[0].b2 *= gain;
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

} // namespace detail

/**
 * @brief Factor a discrete transfer function in z into second-order sections, for BiquadCascade
 *
 * A single direct form of order N ≳ 3 loses accuracy as its poles crowd near z = 1: the
 * coefficients' rounding moves the poles by much more than the rounding itself. The same filter
 * as a cascade of second-order sections keeps each pole pair in its own two coefficients.
 *
 * Pairing and order follow MATLAB®'s zp2sos defaults: the poles closest to the unit circle take
 * the zeros closest to them, the sections run from the poles farthest from the unit circle to the
 * closest ('up'), and the gain is in the first section (scale 'none').
 *
 * @p tf is in z, highest power first and aligned on the constant term, as TransferFunction
 * stores it. A numerator of lower degree than the denominator is a delay, kept: the sections
 * without zeros carry it (b = {0, 0, 1} for two poles, none of them zeros). Leading exact zeros
 * lower a polynomial's degree.
 *
 * @code
 * constexpr auto delay = fractional_delay<5>(4.6e-3, 1e-3);           // Thiran, 4.6 samples
 * constexpr auto sos   = design::second_order_sections(*delay);       // 3 sections
 * static constinit BiquadCascade<3, float> line{sos};
 * @endcode
 *
 * @note Compare with MATLAB®'s tf2sos (matlab::tf2sos takes MATLAB®'s z⁻¹ vectors).
 * @see "Discrete-Time Signal Processing" (Oppenheim & Schafer, 3rd ed., 2010), §6.3.2 (cascade form)
 *
 * @param tf Transfer function in z, order Nden − 1 ≥ 1, no more numerator than denominator coefficients
 * @return Nden / 2 sections; invalid_parameter if a coefficient is not finite, the denominator is
 *         zero, or its leading zeros leave it of lower degree than the numerator (not causal);
 *         not_converged if a root solve fails
 */
template<size_t Nnum, size_t Nden, typename T>
    requires(Nnum >= 1 && Nnum <= Nden && Nden >= 2)
[[nodiscard]] constexpr SecondOrderSectionsResult<Nden / 2, T> second_order_sections(const TransferFunction<Nnum, Nden, T>& tf) {
    constexpr size_t NS = Nden / 2;
    using Result = SecondOrderSectionsResult<NS, T>;
    for (size_t i = 0; i < Nnum; ++i) {
        if (!damp::isfinite(tf.num[i])) {
            return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter);
        }
    }
    for (size_t i = 0; i < Nden; ++i) {
        if (!damp::isfinite(tf.den[i])) {
            return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter);
        }
    }
    const auto zr = detail::finite_roots(tf.num);
    const auto pr = detail::finite_roots(tf.den);
    if (pr.lead == T{0}) {
        return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter);
    }
    if (!zr.success || !pr.success) {
        return ::damp::detail::refused<Result>(DesignStatus::not_converged);
    }
    return detail::sections_from_roots<NS>(zr.roots, zr.count, pr.roots, pr.count, zr.lead / pr.lead);
}

/**
 * @brief Second-order sections of a discrete ZPK, H(z) = k·Π(z − zᵢ)/Π(z − pⱼ), for BiquadCascade
 *
 * As the transfer-function form: fewer zeros than poles is a delay, kept in the sections.
 *
 * @note Compare with MATLAB®'s zp2sos(z, p, k).
 *
 * @param sys Zeros, poles (at least one, and no fewer than zeros) and gain in z
 * @return (Np + 1) / 2 sections; invalid_parameter if a root or the gain is not finite;
 *         not_converged if a complex root has no conjugate
 */
template<size_t Nz, size_t Np, typename T>
    requires(Nz <= Np && Np >= 1)
[[nodiscard]] constexpr SecondOrderSectionsResult<(Np + 1) / 2, T> second_order_sections(const ZPK<Nz, Np, T>& sys) {
    constexpr size_t NS = (Np + 1) / 2;
    using Result = SecondOrderSectionsResult<NS, T>;
    for (size_t i = 0; i < Nz; ++i) {
        if (!damp::isfinite(sys.zeros[i].real()) || !damp::isfinite(sys.zeros[i].imag())) {
            return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter);
        }
    }
    for (size_t j = 0; j < Np; ++j) {
        if (!damp::isfinite(sys.poles[j].real()) || !damp::isfinite(sys.poles[j].imag())) {
            return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter);
        }
    }
    return detail::sections_from_roots<NS>(sys.zeros, Nz, sys.poles, Np, sys.gain);
}

namespace detail {

/** @brief Qₖ of the Butterworth pair k (1 ≤ k ≤ N/2): 1 / (2 sin((2k − 1)π / 2N)). */
template<typename T>
[[nodiscard]] constexpr T butterworth_Q(size_t k, size_t N) {
    return T{1} / (T{2} * damp::sin(static_cast<T>((2 * k) - 1) * damp::numbers::pi_v<T> / static_cast<T>(2 * N)));
}

/** @brief A first-order section as a SecondOrderCoeffs, b₂ = a₂ = 0. */
template<typename T>
[[nodiscard]] constexpr SecondOrderCoeffs<T> first_order_section(const FirstOrderCoeffs<T>& c) {
    return SecondOrderCoeffs<T>{c.b0, c.b1, T{0}, c.a1, T{0}, c.success, c.status};
}

} // namespace detail

/**
 * @brief Nth-order Butterworth low-pass at cutoff @p fc as second-order sections, any order
 *
 * Each section is the prewarped Tustin @ref lowpass_2nd at the same corner, Qₖ = 1/(2 sin((2k − 1)π/2N)),
 * with a @ref lowpass_1st section for odd N, so the cascade is MATLAB®'s butter(N, 2·fc·Ts) and the
 * −3 dB corner lands on @p fc. The sections run from the lowest Q to the highest ('up'), each with
 * unity DC gain.
 *
 * @code
 * static constinit BiquadCascade<3, float> lpf{design::lowpass_sections<6>(200.0f, 1.0f / 10000.0f)};
 * @endcode
 *
 * @note Compare with MATLAB®'s [z, p, k] = butter(N, 2*fc*Ts); sos = zp2sos(z, p, k).
 *
 * @tparam N Filter order (≥ 1)
 * @param fc [Hz] Cutoff frequency, 0 < fc < 1/(2·Ts)
 * @param Ts [s] Sample time, > 0
 * @return ⌈N/2⌉ sections; invalid_sample_time if Ts ≤ 0, invalid_parameter if fc is outside (0, Nyquist)
 */
template<size_t N, typename T = double>
    requires(N >= 1)
[[nodiscard]] constexpr SecondOrderSectionsResult<(N + 1) / 2, T> lowpass_sections(T fc, T Ts) {
    using Result = SecondOrderSectionsResult<(N + 1) / 2, T>;
    if (!detail::valid_cutoff(fc, Ts)) {
        return ::damp::detail::refused<Result>(detail::cutoff_status(Ts));
    }
    Result out{};
    size_t s = 0;
    if constexpr (N % 2 == 1) {
        out.sections[s++] = detail::first_order_section(lowpass_1st(fc, Ts));
    }
    for (size_t k = N / 2; k >= 1; --k) {
        out.sections[s++] = lowpass_2nd(fc, detail::butterworth_Q<T>(k, N), Ts);
    }
    for (const auto& c : out.sections) {
        if (!c.success) {
            return ::damp::detail::refused<Result>(c.status);
        }
    }
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

/**
 * @brief Nth-order Butterworth high-pass at cutoff @p fc as second-order sections, any order
 *
 * The high-pass counterpart of @ref lowpass_sections: @ref highpass_2nd sections at
 * the Butterworth Qₖ and a @ref highpass_1st section for odd N, each with unity gain at Nyquist.
 *
 * @note Compare with MATLAB®'s [z, p, k] = butter(N, 2*fc*Ts, 'high'); sos = zp2sos(z, p, k).
 *
 * @tparam N Filter order (≥ 1)
 * @param fc [Hz] Cutoff frequency, 0 < fc < 1/(2·Ts)
 * @param Ts [s] Sample time, > 0
 * @return ⌈N/2⌉ sections; invalid_sample_time if Ts ≤ 0, invalid_parameter if fc is outside (0, Nyquist)
 */
template<size_t N, typename T = double>
    requires(N >= 1)
[[nodiscard]] constexpr SecondOrderSectionsResult<(N + 1) / 2, T> highpass_sections(T fc, T Ts) {
    using Result = SecondOrderSectionsResult<(N + 1) / 2, T>;
    if (!detail::valid_cutoff(fc, Ts)) {
        return ::damp::detail::refused<Result>(detail::cutoff_status(Ts));
    }
    Result out{};
    size_t s = 0;
    if constexpr (N % 2 == 1) {
        out.sections[s++] = detail::first_order_section(highpass_1st(fc, Ts));
    }
    for (size_t k = N / 2; k >= 1; --k) {
        out.sections[s++] = highpass_2nd(fc, detail::butterworth_Q<T>(k, N), Ts);
    }
    for (const auto& c : out.sections) {
        if (!c.success) {
            return ::damp::detail::refused<Result>(c.status);
        }
    }
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

} // namespace design

/**
 * @brief Second-order IIR (biquad) section runtime.
 *
 * Runs a SecondOrderCoeffs section in Direct Form I:
 *   y[n] = b0·x[n] + b1·x[n-1] + b2·x[n-2] − a1·y[n-1] − a2·y[n-2]
 *
 * Pair with any design:: biquad designer (notch, bandpass, highpass_2nd,
 * peaking, lowshelf, highshelf, lowpass_2nd).
 *
 * @code
 * Biquad<float> notch{design::notch(50.0f, 5.0f, 1.0f / 1000.0f)};
 * float clean = notch(sample);
 * @endcode
 */
template<typename T = float>
class Biquad {
public:
    constexpr Biquad() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr Biquad(const design::SecondOrderCoeffs<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : Biquad(result.template as<T>()) {}

    constexpr Biquad(const design::SecondOrderCoeffs<T>& c) // NOLINT(google-explicit-constructor): assign-from-design
        : b0_(c.b0), b1_(c.b1), b2_(c.b2), a1_(c.a1), a2_(c.a2) {
        if (!c.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(c.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /** @brief Convert across scalar precision, preserving coefficients and delay line. */
    template<typename U>
    constexpr explicit Biquad(const Biquad<U>& o)
        : b0_(static_cast<T>(o.b0_)),
          b1_(static_cast<T>(o.b1_)),
          b2_(static_cast<T>(o.b2_)),
          a1_(static_cast<T>(o.a1_)),
          a2_(static_cast<T>(o.a2_)),
          x1_(static_cast<T>(o.x1_)),
          x2_(static_cast<T>(o.x2_)),
          y1_(static_cast<T>(o.y1_)),
          y2_(static_cast<T>(o.y2_)) {}

    /** @brief Process one sample. */
    constexpr T operator()(T x) {
        const T y = (b0_ * x) + (b1_ * x1_) + (b2_ * x2_) - (a1_ * y1_) - (a2_ * y2_);
        x2_ = x1_;
        x1_ = x;
        y2_ = y1_;
        y1_ = y;
        return y;
    }

    /** @brief Reset the internal delay line. */
    constexpr void reset() {
        x1_ = x2_ = y1_ = y2_ = T{0};
    }

    /**
     * @brief Replace the coefficients in place, keeping the delay line — for adaptive
     * filters that retune (e.g. grid-frequency tracking) without losing state.
     */
    constexpr void set_coefficients(const design::SecondOrderCoeffs<T>& c) {
        b0_ = c.b0;
        b1_ = c.b1;
        b2_ = c.b2;
        a1_ = c.a1;
        a2_ = c.a2;
    }

    /** @brief Most recent output y[n-1]; read/adjust to seed anti-windup unwinding. */
    [[nodiscard]] constexpr T last_output() const { return y1_; }
    constexpr void            set_last_output(T y1) { y1_ = y1; }

private:
    template<typename>
    friend class Biquad;

    T b0_{1}, b1_{0}, b2_{0}, a1_{0}, a2_{0};
    T x1_{0}, x2_{0}, y1_{0}, y2_{0};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
Biquad(const design::SecondOrderCoeffs<T>&) -> Biquad<T>;

/**
 * @brief Cascade of second-order sections (SOS) for higher-order IIR filters.
 *
 * Chains NSections biquads in series. Cascading is the numerically preferred
 * realization for higher-order IIR filters (vs. a single high-order section).
 * Build it from design::lowpass_sections / highpass_sections (Butterworth, any
 * order) or design::second_order_sections (any transfer function in z).
 *
 * @code
 * static constinit BiquadCascade<2, float> lpf{design::lowpass_sections<4>(100.0f, 1.0f / 1000.0f)};
 * float smooth = lpf(sample);
 * @endcode
 *
 * @tparam NSections Number of biquad sections
 * @tparam T         Scalar type
 */
template<size_t NSections, typename T = float>
class BiquadCascade {
public:
    constexpr BiquadCascade() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr BiquadCascade(const design::SecondOrderSectionsResult<NSections, U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : BiquadCascade(result.template as<T>()) {}

    /** @brief From a sections design. A failed design does not compile when constant-evaluated; at run time it outputs 0. */
    constexpr BiquadCascade(const design::SecondOrderSectionsResult<NSections, T>& result) { // NOLINT(google-explicit-constructor): assign-from-design
        if (!result.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        for (size_t i = 0; i < NSections; ++i) {
            sections_[i] = Biquad<T>(result.sections[i]);
        }
    }

    /** @brief From hand-built sections in cascade order. */
    constexpr explicit BiquadCascade(const damp::array<design::SecondOrderCoeffs<T>, NSections>& sections) {
        for (size_t i = 0; i < NSections; ++i) {
            sections_[i] = Biquad<T>(sections[i]);
        }
    }

    /** @brief Process one sample through all sections in series. */
    constexpr T operator()(T x) {
        for (size_t i = 0; i < NSections; ++i) {
            x = sections_[i](x);
        }
        return x;
    }

    /** @brief Process only the first @p n_active sections (bypass the rest as identity). */
    constexpr T operator()(T x, size_t n_active) {
        const size_t n = (n_active < NSections) ? n_active : NSections;
        for (size_t i = 0; i < n; ++i) {
            x = sections_[i](x);
        }
        return x;
    }

    /** @brief Replace coefficients of section @p i (keeps that section's delay line). */
    constexpr void set_section(size_t i, const design::SecondOrderCoeffs<T>& c) {
        if (i < NSections) {
            sections_[i].set_coefficients(c);
        }
    }

    /** @brief Access one section (e.g. reset a single biquad after retune). */
    [[nodiscard]] constexpr Biquad<T>& section(size_t i) { return sections_[i]; }

    [[nodiscard]] constexpr const Biquad<T>& section(size_t i) const { return sections_[i]; }

    /** @brief Reset every section. */
    constexpr void reset() {
        for (size_t i = 0; i < NSections; ++i) {
            sections_[i].reset();
        }
    }

    /** @brief Reset and load identity (bypass) coefficients into every section. */
    constexpr void reset_identity() {
        constexpr design::SecondOrderCoeffs<T> id{T{1}, T{0}, T{0}, T{0}, T{0}};
        for (size_t i = 0; i < NSections; ++i) {
            sections_[i] = Biquad<T>(id);
        }
    }

private:
    damp::array<Biquad<T>, NSections> sections_{};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<size_t NSections, typename T>
BiquadCascade(const design::SecondOrderSectionsResult<NSections, T>&) -> BiquadCascade<NSections, T>;
} // namespace damp
