// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file fir.hpp
 * @brief Direct-form FIR runtime and window-method coefficient design (`fir1`)
 *
 * Example: low-pass FIR for a 1 kHz sample rate, 100 Hz cutoff (fc = 0.1 cycles/sample)
 * @code
 * #include "damp/filters/fir.hpp"
 * using namespace damp;
 *
 * constexpr auto des = design::fir1<64, double>(
 *     31, 0.1, design::FirType::Lowpass, design::FirWindow::Hamming);
 * static_assert(des.success);
 *
 * FirFilter<64, float> fir(des.as<float>());
 * // In ISR: float y = fir(x);
 * @endcode
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp {
namespace design {

/**
 * @brief Ideal frequency-selective response for @ref fir1 / @ref fir_window
 */
enum class FirType {
    Lowpass,  ///< Pass [0, fc], stop (fc, 0.5]
    Highpass, ///< Stop [0, fc], pass (fc, 0.5]
    Bandpass, ///< Pass [f_lo, f_hi]
    Bandstop  ///< Stop [f_lo, f_hi]
};

/**
 * @brief Window applied to the ideal truncated impulse response
 *
 * @see Oppenheim & Schafer, "Discrete-Time Signal Processing", window method
 */
enum class FirWindow {
    Rectangular, ///< w[n] = 1 (narrowest main lobe, highest sidelobes)
    Hann,        ///< raised-cosine, 0 endpoints (a.k.a. Hanning)
    Hamming,     ///< 0.54 − 0.46 cos (MATLAB® fir1 default)
    Blackman     ///< three-term cosine window
};

/**
 * @brief Window-method FIR design result
 *
 * Coefficients are stored b[0]…b[n_taps−1] with b[0] on the newest sample
 * (same convention as @ref FirFilter). Passband gain is normalized to ≈ 1.
 *
 * @tparam MaxTaps Storage capacity for coefficients
 * @tparam T       Scalar type
 */
template<size_t MaxTaps, typename T = double>
struct FirDesignResult {
    damp::array<T, MaxTaps> b{};            ///< FIR coefficients (only first n_taps used)
    size_t                  n_taps{0};      ///< Active length (1…MaxTaps when success)
    bool                    success{false}; ///< true if design completed

    /**
     * @brief Convert coefficients to a different scalar type
     * @tparam U Target scalar type
     */
    template<typename U>
    [[nodiscard]] constexpr FirDesignResult<MaxTaps, std::remove_const_t<U>> as() const {
        FirDesignResult<MaxTaps, std::remove_const_t<U>> out{};
        using O = std::remove_const_t<U>;
        for (size_t i = 0; i < MaxTaps; ++i) {
            out.b[i] = static_cast<O>(b[i]);
        }
        out.n_taps = n_taps;
        out.success = success;
        return out;
    }
};

namespace detail {

/**
 * @brief Window weight w[n] for length @p N (order M = N−1)
 *
 * Endpoints use the standard discrete formulas with M = N−1 in the cosine
 * argument @f$ 2\pi n / M @f$. For N = 1 the window is identically 1.
 */
template<typename T>
[[nodiscard]] constexpr T fir_window_weight(FirWindow window, size_t n, size_t N) {
    if (N <= 1) {
        return T{1};
    }
    const T M = static_cast<T>(N - 1);
    const T tn = static_cast<T>(n);
    const T two_pi = T{2} * damp::numbers::pi_v<T>;
    switch (window) {
        case FirWindow::Rectangular:
            return T{1};
        case FirWindow::Hann:
            return static_cast<T>(0.5) * (T{1} - damp::cos(two_pi * tn / M));
        case FirWindow::Hamming:
            return static_cast<T>(0.54) - (static_cast<T>(0.46) * damp::cos(two_pi * tn / M));
        case FirWindow::Blackman:
            return static_cast<T>(0.42) - (static_cast<T>(0.5) * damp::cos(two_pi * tn / M))
                 + (static_cast<T>(0.08) * damp::cos(T{2} * two_pi * tn / M));
    }
    return T{1};
}

/**
 * @brief Ideal low-pass impulse response sample (unit sample rate, fc in cycles/sample)
 *
 * @f[
 *   h_d[n] = 2 f_c \,\mathrm{sinc}\bigl(2 f_c (n-\alpha)\bigr),\quad
 *   \alpha = (N-1)/2,\quad
 *   \mathrm{sinc}(x)=\sin(\pi x)/(\pi x)
 * @f]
 * with the removable singularity @f$ h_d[\alpha] = 2 f_c @f$ when @f$ n = \alpha @f$.
 *
 * @see Oppenheim & Schafer, Discrete-Time Signal Processing, window method for FIR design
 */
template<typename T>
[[nodiscard]] constexpr T ideal_lowpass_tap(size_t n, size_t N, T fc) {
    const T alpha = static_cast<T>(N - 1) / T{2};
    const T x = static_cast<T>(n) - alpha;
    if (x == T{0}) {
        return T{2} * fc;
    }
    const T arg = T{2} * damp::numbers::pi_v<T> * fc * x;
    return damp::sin(arg) / (damp::numbers::pi_v<T> * x);
}

/**
 * @brief Real (zero-phase) gain of a linear-phase FIR at digital frequency @p f
 *
 * @f$ G(f) = \sum_n b[n]\cos\bigl(2\pi f (n-\alpha)\bigr) @f$ with
 * @f$ \alpha=(N-1)/2 @f$. Used to normalize passband gain to 1.
 */
template<size_t MaxTaps, typename T>
[[nodiscard]] constexpr T fir_passband_gain(
    const damp::array<T, MaxTaps>& b, size_t N, T f_cycles
) {
    const T alpha = static_cast<T>(N - 1) / T{2};
    const T w = T{2} * damp::numbers::pi_v<T> * f_cycles;
    T       g = T{0};
    for (size_t n = 0; n < N; ++n) {
        g += b[n] * damp::cos(w * (static_cast<T>(n) - alpha));
    }
    return g;
}

template<typename T>
[[nodiscard]] constexpr bool valid_fc(T fc) {
    return (fc > T{0}) && (fc < static_cast<T>(0.5));
}

template<typename T>
[[nodiscard]] constexpr bool valid_band(T f_lo, T f_hi) {
    return valid_fc(f_lo) && valid_fc(f_hi) && (f_lo < f_hi);
}

/**
 * @brief Highpass / bandstop need an integer group delay sample (odd length / Type I)
 */
[[nodiscard]] constexpr size_t force_odd_length(size_t n_taps, size_t max_taps) {
    if ((n_taps % 2) == 1) {
        return n_taps;
    }
    if (n_taps + 1 <= max_taps) {
        return n_taps + 1;
    }
    if (n_taps >= 2) {
        return n_taps - 1;
    }
    return n_taps;
}

} // namespace detail

/**
 * @brief Window-method FIR design (normalized frequency, single cutoff)
 *
 * Designs an equiripple-free FIR by truncating the ideal impulse response and
 * multiplying by a window. Cutoff @p fc is in cycles per sample on
 * @f$ (0,\,0.5) @f$ (Nyquist = 0.5). Equivalent Hz form:
 * @f$ f_c = f_{\mathrm{Hz}} / f_s @f$.
 *
 * For @ref FirType::Highpass, an even @p n_taps is adjusted to the nearest odd
 * length ≤ @p MaxTaps (Type I linear phase needs a centre tap for spectral
 * inversion). Bandpass/bandstop use the two-cutoff overload.
 *
 * Coefficients are scaled so the real passband-centre gain is 1.
 *
 * @note Compare with MATLAB®'s fir1(n, Wn) (MATLAB uses order n = n_taps−1 and
 *       Wn on (0,1] with 1 = Nyquist; Damp uses n_taps and fc on (0, 0.5)).
 *
 * @see "Discrete-Time Signal Processing" (Oppenheim & Schafer), FIR window method
 * @see fir_window() alias
 *
 * @tparam MaxTaps Coefficient storage capacity
 * @tparam T       Scalar type (default double for design)
 * @param n_taps   Desired number of taps (1…MaxTaps)
 * @param fc       Cutoff frequency in cycles/sample, strictly in (0, 0.5)
 * @param type     Lowpass or Highpass (Bandpass/Bandstop require two cutoffs)
 * @param window   Window family (default Hamming)
 * @return FirDesignResult; success=false on invalid arguments
 */
template<size_t MaxTaps, typename T = double>
[[nodiscard]] constexpr FirDesignResult<MaxTaps, T> fir1(
    size_t    n_taps,
    T         fc,
    FirType   type = FirType::Lowpass,
    FirWindow window = FirWindow::Hamming
) {
    FirDesignResult<MaxTaps, T> result{};
    if (n_taps < 1 || n_taps > MaxTaps || !detail::valid_fc(fc)) {
        return result;
    }
    if (type == FirType::Bandpass || type == FirType::Bandstop) {
        return result; // need two cutoffs
    }

    size_t N = n_taps;
    if (type == FirType::Highpass) {
        N = detail::force_odd_length(n_taps, MaxTaps);
        if ((N % 2) == 0) {
            return result;
        }
    }

    // Ideal low-pass × window
    for (size_t n = 0; n < N; ++n) {
        const T hd = detail::ideal_lowpass_tap<T>(n, N, fc);
        const T w = detail::fir_window_weight<T>(window, n, N);
        result.b[n] = hd * w;
    }

    // Unit DC gain on the low-pass prototype (exact spectral inversion needs this)
    {
        T s = T{0};
        for (size_t n = 0; n < N; ++n) {
            s += result.b[n];
        }
        if (damp::abs(s) < static_cast<T>(1e-30)) {
            return result;
        }
        const T inv_s = T{1} / s;
        for (size_t n = 0; n < N; ++n) {
            result.b[n] *= inv_s;
        }
    }

    if (type == FirType::Highpass) {
        // Spectral inversion: h_hp = δ[n−α] − h_lp  → exact DC zero when sum(h_lp)=1
        const size_t alpha_i = (N - 1) / 2;
        for (size_t n = 0; n < N; ++n) {
            result.b[n] = -result.b[n];
        }
        result.b[alpha_i] += T{1};

        // Unit gain at Nyquist (preserves the DC null)
        const T g = detail::fir_passband_gain<MaxTaps, T>(result.b, N, static_cast<T>(0.5));
        if (damp::abs(g) < static_cast<T>(1e-30)) {
            return result;
        }
        const T inv_g = T{1} / g;
        for (size_t n = 0; n < N; ++n) {
            result.b[n] *= inv_g;
        }
    }

    result.n_taps = N;
    result.success = true;
    return result;
}

/**
 * @brief Window-method FIR design with two cutoffs (bandpass / bandstop)
 *
 * @p f_lo and @p f_hi are in cycles/sample on (0, 0.5) with @p f_lo < @p f_hi.
 * Bandstop lengths are forced odd (Type I). Passband-centre gain is normalized
 * to 1 (bandpass midband; bandstop at DC).
 *
 * @tparam MaxTaps Coefficient storage capacity
 * @tparam T       Scalar type
 * @param n_taps   Desired taps
 * @param f_lo     Lower cutoff (cycles/sample)
 * @param f_hi     Upper cutoff (cycles/sample)
 * @param type     Bandpass or Bandstop
 * @param window   Window family
 * @return FirDesignResult; success=false on invalid arguments
 */
template<size_t MaxTaps, typename T = double>
[[nodiscard]] constexpr FirDesignResult<MaxTaps, T> fir1(
    size_t    n_taps,
    T         f_lo,
    T         f_hi,
    FirType   type = FirType::Bandpass,
    FirWindow window = FirWindow::Hamming
) {
    FirDesignResult<MaxTaps, T> result{};
    if (n_taps < 1 || n_taps > MaxTaps || !detail::valid_band(f_lo, f_hi)) {
        return result;
    }
    if (type != FirType::Bandpass && type != FirType::Bandstop) {
        // Single-cutoff types: use the first cutoff only
        return fir1<MaxTaps, T>(n_taps, f_lo, type, window);
    }

    size_t N = n_taps;
    if (type == FirType::Bandstop) {
        N = detail::force_odd_length(n_taps, MaxTaps);
        if ((N % 2) == 0) {
            return result;
        }
    }

    for (size_t n = 0; n < N; ++n) {
        const T h_hi = detail::ideal_lowpass_tap<T>(n, N, f_hi);
        const T h_lo = detail::ideal_lowpass_tap<T>(n, N, f_lo);
        const T w = detail::fir_window_weight<T>(window, n, N);
        if (type == FirType::Bandpass) {
            // Ideal bandpass = LPF(f_hi) − LPF(f_lo)
            result.b[n] = (h_hi - h_lo) * w;
        } else {
            // Ideal bandstop = δ − bandpass = LPF(f_lo) + (δ − LPF(f_hi))
            result.b[n] = (h_lo - h_hi) * w;
        }
    }

    if (type == FirType::Bandstop) {
        const size_t alpha_i = (N - 1) / 2;
        result.b[alpha_i] += T{1};
    }

    const T f_pb = (type == FirType::Bandpass) ? ((f_lo + f_hi) / T{2}) : T{0};
    const T g = detail::fir_passband_gain<MaxTaps, T>(result.b, N, f_pb);
    if (damp::abs(g) < static_cast<T>(1e-30)) {
        return result;
    }
    const T inv_g = T{1} / g;
    for (size_t n = 0; n < N; ++n) {
        result.b[n] *= inv_g;
    }

    result.n_taps = N;
    result.success = true;
    return result;
}

/**
 * @brief Window-method FIR design from frequencies in Hz
 *
 * Converts @f$ f_c = f_{\mathrm{Hz}}/f_s @f$ then calls @ref fir1.
 *
 * @param n_taps  Number of taps
 * @param fc_hz   Cutoff [Hz]
 * @param fs_hz   Sample rate [Hz] (> 0); Nyquist = fs/2
 * @param type    Lowpass or Highpass
 * @param window  Window family
 */
template<size_t MaxTaps, typename T = double>
[[nodiscard]] constexpr FirDesignResult<MaxTaps, T> fir1_hz(
    size_t    n_taps,
    T         fc_hz,
    T         fs_hz,
    FirType   type = FirType::Lowpass,
    FirWindow window = FirWindow::Hamming
) {
    if (fs_hz <= T{0}) {
        return FirDesignResult<MaxTaps, T>{};
    }
    return fir1<MaxTaps, T>(n_taps, fc_hz / fs_hz, type, window);
}

/**
 * @brief Two-cutoff Hz overload for bandpass / bandstop
 */
template<size_t MaxTaps, typename T = double>
[[nodiscard]] constexpr FirDesignResult<MaxTaps, T> fir1_hz(
    size_t    n_taps,
    T         f_lo_hz,
    T         f_hi_hz,
    T         fs_hz,
    FirType   type = FirType::Bandpass,
    FirWindow window = FirWindow::Hamming
) {
    if (fs_hz <= T{0}) {
        return FirDesignResult<MaxTaps, T>{};
    }
    return fir1<MaxTaps, T>(n_taps, f_lo_hz / fs_hz, f_hi_hz / fs_hz, type, window);
}

/**
 * @brief Alias for @ref fir1 (descriptive name)
 * @see fir1()
 */
template<size_t MaxTaps, typename T = double>
[[nodiscard]] constexpr FirDesignResult<MaxTaps, T> fir_window(
    size_t    n_taps,
    T         fc,
    FirType   type = FirType::Lowpass,
    FirWindow window = FirWindow::Hamming
) {
    return fir1<MaxTaps, T>(n_taps, fc, type, window);
}

/**
 * @brief Two-cutoff alias for @ref fir1
 */
template<size_t MaxTaps, typename T = double>
[[nodiscard]] constexpr FirDesignResult<MaxTaps, T> fir_window(
    size_t    n_taps,
    T         f_lo,
    T         f_hi,
    FirType   type = FirType::Bandpass,
    FirWindow window = FirWindow::Hamming
) {
    return fir1<MaxTaps, T>(n_taps, f_lo, f_hi, type, window);
}

} // namespace design

/**
 * @brief Direct-form FIR filter runtime (tapped delay line + dot product)
 *
 * Implements
 * @f[
 *   y[n] = \sum_{k=0}^{N-1} b[k]\,x[n-k]
 * @f]
 * with a fixed-capacity circular buffer. @p MaxTaps bounds storage at compile
 * time; the active length @p N may be set at runtime with @p N ≤ MaxTaps.
 *
 * Pair with @ref design::fir1 / @ref design::fir_window for coefficients, or
 * load an arbitrary coefficient array.
 *
 * @tparam MaxTaps Maximum number of taps (buffer capacity)
 * @tparam T       Scalar type (default float for embedded)
 */
template<size_t MaxTaps, typename T = float>
class FirFilter {
public:
    static_assert(MaxTaps >= 1, "FirFilter needs MaxTaps >= 1");

    /// Identity (single tap b0 = 1) until @ref init
    constexpr FirFilter() {
        b_[0] = T{1};
    }

    /// Construct from a design result (no-op if !success — leaves identity)
    constexpr explicit FirFilter(const design::FirDesignResult<MaxTaps, T>& des) { init(des); }

    /**
     * @brief Load coefficients from a design result
     * @return true if des.success and length is valid
     */
    constexpr bool init(const design::FirDesignResult<MaxTaps, T>& des) {
        if (!des.success || des.n_taps < 1 || des.n_taps > MaxTaps) {
            return false;
        }
        return init(des.b, des.n_taps);
    }

    /**
     * @brief Load @p n coefficients from an array (clamped to [1, MaxTaps])
     * @param coeffs Coefficient array (at least @p n valid entries)
     * @param n      Active length
     * @return true if n is in range
     */
    constexpr bool init(const damp::array<T, MaxTaps>& coeffs, size_t n) {
        if (n < 1 || n > MaxTaps) {
            return false;
        }
        n_taps_ = n;
        for (size_t i = 0; i < n; ++i) {
            b_[i] = coeffs[i];
        }
        for (size_t i = n; i < MaxTaps; ++i) {
            b_[i] = T{0};
        }
        reset();
        return true;
    }

    /// Process one sample
    constexpr T operator()(T x) {
        x_[write_] = x;
        T      y = T{0};
        size_t j = write_;
        for (size_t k = 0; k < n_taps_; ++k) {
            y += b_[k] * x_[j];
            j = (j == 0) ? (n_taps_ - 1) : (j - 1);
        }
        write_ = (write_ + 1) % n_taps_;
        return y;
    }

    /// Clear the delay line (coefficients unchanged)
    constexpr void reset() {
        x_ = {};
        write_ = 0;
    }

    [[nodiscard]] constexpr size_t n_taps() const { return n_taps_; }

    /// Read coefficient b[k] (0 if k ≥ n_taps)
    [[nodiscard]] constexpr T coefficient(size_t k) const {
        return (k < n_taps_) ? b_[k] : T{0};
    }

private:
    damp::array<T, MaxTaps> b_{};
    damp::array<T, MaxTaps> x_{};
    size_t                  n_taps_{1};
    size_t                  write_{0};
};

/**
 * @brief Normalized LMS adaptive FIR filter
 *
 * Updates coefficients each sample with the normalized LMS rule
 * @f[
 *   \mathbf{b} \leftarrow \mathbf{b} + \frac{\mu\,e[n]}
 *     {\varepsilon + \|\mathbf{x}\|^2}\,\mathbf{x},\quad
 *   e[n] = d[n] - \mathbf{b}^\top\mathbf{x}
 * @f]
 * where @p d is the desired response and @p x is the tapped input vector.
 * Useful for system identification, adaptive noise cancellation, and line
 * enhancement.
 *
 * @tparam MaxTaps Maximum taps
 * @tparam T       Scalar type
 *
 * @see "Adaptive Filter Theory" (Haykin), NLMS algorithm
 */
template<size_t MaxTaps, typename T = float>
class NlmsFilter {
public:
    static_assert(MaxTaps >= 1, "NlmsFilter needs MaxTaps >= 1");

    constexpr NlmsFilter() = default;

    /**
     * @brief Configure length and step size
     * @param n_taps Active taps (clamped to [1, MaxTaps])
     * @param mu     Step size, typically in (0, 2) for NLMS stability
     * @param eps    Regularization @f$ \varepsilon > 0 @f$ (avoids /0)
     */
    constexpr void init(size_t n_taps, T mu = static_cast<T>(0.5), T eps = static_cast<T>(1e-6)) {
        n_taps_ = damp::clamp(n_taps, size_t{1}, MaxTaps);
        mu_ = mu;
        eps_ = (eps > T{0}) ? eps : static_cast<T>(1e-6);
        reset();
    }

    /**
     * @brief Filter @p x toward desired @p d; update weights; return output y
     *
     * @param x Input sample (pushed into the tapped delay line)
     * @param d Desired response for this sample
     * @return Filter output y = bᵀx before the weight update
     */
    constexpr T operator()(T x, T d) {
        x_[write_] = x;

        // y = b · x  (same circular indexing as FirFilter)
        T      y = T{0};
        T      power = T{0};
        size_t j = write_;
        for (size_t k = 0; k < n_taps_; ++k) {
            const T xk = x_[j];
            y += b_[k] * xk;
            power += xk * xk;
            j = (j == 0) ? (n_taps_ - 1) : (j - 1);
        }

        const T e = d - y;
        const T step = (mu_ * e) / (eps_ + power);
        j = write_;
        for (size_t k = 0; k < n_taps_; ++k) {
            b_[k] += step * x_[j];
            j = (j == 0) ? (n_taps_ - 1) : (j - 1);
        }

        write_ = (write_ + 1) % n_taps_;
        last_error_ = e;
        return y;
    }

    /// Zero coefficients and delay line
    constexpr void reset() {
        b_ = {};
        x_ = {};
        write_ = 0;
        last_error_ = T{0};
    }

    [[nodiscard]] constexpr size_t n_taps() const { return n_taps_; }
    [[nodiscard]] constexpr T      last_error() const { return last_error_; }
    [[nodiscard]] constexpr T      coefficient(size_t k) const {
        return (k < n_taps_) ? b_[k] : T{0};
    }

    /// Copy coefficients into a fixed FIR (freeze adaptation)
    [[nodiscard]] constexpr design::FirDesignResult<MaxTaps, T> coefficients() const {
        design::FirDesignResult<MaxTaps, T> out{};
        for (size_t i = 0; i < n_taps_; ++i) {
            out.b[i] = b_[i];
        }
        out.n_taps = n_taps_;
        out.success = true;
        return out;
    }

private:
    damp::array<T, MaxTaps> b_{};
    damp::array<T, MaxTaps> x_{};
    size_t                  n_taps_{MaxTaps};
    size_t                  write_{0};
    T                       mu_{static_cast<T>(0.5)};
    T                       eps_{static_cast<T>(1e-6)};
    T                       last_error_{T{0}};
};

} // namespace damp
