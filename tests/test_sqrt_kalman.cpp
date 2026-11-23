// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// The recursive filters carry the triangular factor S of their covariance, P = S Sᵀ
// (estimation/sqrt_covariance.hpp). These tests hold that form to what the covariance form could
// not do (a prior far wider than a reading's noise, in float), to what it did where it was sound
// (a Joseph-form reference kept here as the oracle), and to the statistics of the model itself.
// Float is the first case of every template.

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/estimation/ekf.hpp"
#include "damp/estimation/fixed_gain_observer.hpp" // IWYU pragma: keep (defines Law<kRes, T>::Recursive)
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/sqrt_covariance.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/matrix/pattern_matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

/// Standard normals that are the same on every host: a 64-bit LCG (Knuth's MMIX constants) through Box–Muller.
struct Normal {
    std::uint64_t state;

    double uniform() {
        state = (state * 6364136223846793005ULL) + 1442695040888963407ULL;
        return (static_cast<double>(state >> 11U) + 0.5) * (1.0 / 9007199254740992.0);
    }
    double operator()() {
        const double u1 = uniform();
        const double u2 = uniform();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    }
};

/// Largest entry of |a − b|, in double.
template<typename A, typename B>
double gap(const A& a, const B& b) {
    double g = 0.0;
    for (size_t r = 0; r < A::rows(); ++r) {
        for (size_t c = 0; c < A::cols(); ++c) {
            g = std::fmax(g, std::fabs(static_cast<double>(a(r, c)) - static_cast<double>(b(r, c))));
        }
    }
    return g;
}

/// Largest |entry|, in double.
template<typename A>
double largest(const A& a) {
    double g = 0.0;
    for (size_t r = 0; r < A::rows(); ++r) {
        for (size_t c = 0; c < A::cols(); ++c) {
            g = std::fmax(g, std::fabs(static_cast<double>(a(r, c))));
        }
    }
    return g;
}

template<size_t R, size_t C, typename T>
bool all_finite(const Matrix<R, C, T>& m) {
    return mat::all_of(m, [](T v) { return damp::isfinite(v); });
}

/// A factor as the filter keeps it: zeros above the diagonal, nothing negative on it.
template<size_t N, typename T>
bool is_factor(const Matrix<N, N, T>& S) {
    for (size_t i = 0; i < N; ++i) {
        if (S(i, i) < T{0}) {
            return false;
        }
        for (size_t j = i + 1; j < N; ++j) {
            if (S(i, j) != T{0}) {
                return false;
            }
        }
    }
    return true;
}

/// The bits of a float or double.
template<typename T>
auto bits_of(T x) {
    if constexpr (sizeof(T) == 4) {
        return std::bit_cast<std::uint32_t>(x);
    } else {
        return std::bit_cast<std::uint64_t>(x);
    }
}

/// Every entry of @p a and @p b has the same bits.
template<size_t R, size_t C, typename T>
bool same_bits(const Matrix<R, C, T>& a, const Matrix<R, C, T>& b) {
    for (size_t r = 0; r < R; ++r) {
        for (size_t c = 0; c < C; ++c) {
            if (bits_of(a(r, c)) != bits_of(b(r, c))) {
                return false;
            }
        }
    }
    return true;
}

/// A quiet NaN and +Inf made at run time: a constant can be folded away under -ffast-math.
template<typename T>
T runtime_nan() {
    volatile T zero = T{0};
    return zero / zero;
}
template<typename T>
T runtime_inf() {
    volatile T zero = T{0};
    return T{1} / zero;
}

/**
 * The covariance form, in double: the oracle. P ← A P Aᵀ + G Q Gᵀ, and the Joseph update
 * P ← (I − K C) P (I − K C)ᵀ + K R Kᵀ with K = P Cᵀ (C P Cᵀ + R)⁻¹.
 */
template<size_t NX, size_t NU>
struct JosephReference {
    Matrix<NX, NX> A{};
    Matrix<NX, NU> B{};
    Matrix<NX, NX> GQG{};
    ColVec<NX>     x{};
    Matrix<NX, NX> P{};

    void predict(const ColVec<NU>& u) {
        x = ColVec<NX>(A * x) + ColVec<NX>(B * u);
        P = Matrix<NX, NX>(A * P * A.t()) + GQG;
        P = Matrix<NX, NX>(P + P.t()) * 0.5;
    }

    template<size_t NY>
    bool update(const ColVec<NY>& y, const Matrix<NY, NX>& C, const Matrix<NY, NY>& R) {
        const Matrix<NY, NY> S = Matrix<NY, NY>(C * P * C.t()) + R;
        const auto           Kt = mat::cholesky_solve(S, Matrix<NY, NX>(C * P));
        if (!Kt) {
            return false;
        }
        const Matrix<NX, NY> K = Kt.value().t();
        x += K * ColVec<NY>(y - (C * x));
        const Matrix<NX, NX> IKC = Matrix<NX, NX>::identity() - Matrix<NX, NX>(K * C);
        P = Matrix<NX, NX>(IKC * P * IKC.t()) + Matrix<NX, NX>(K * R * K.t());
        P = Matrix<NX, NX>(P + P.t()) * 0.5;
        return true;
    }
};

/// The reference on a filter's own model, in the filter's rounded numbers.
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
JosephReference<NX, NU> reference_of(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Matrix<NW, NW, T>& Q, const ColVec<NX, T>& x0, const Matrix<NX, NX, T>& P0) {
    JosephReference<NX, NU> ref{};
    ref.A = sys.A.template as<double>();
    ref.B = sys.B.template as<double>();
    const Matrix<NX, NW> G = sys.G.template as<double>();
    ref.GQG = G * Q.template as<double>() * G.t();
    ref.x = x0.template as<double>();
    ref.P = P0.template as<double>();
    return ref;
}

// ------------------------------------------------------------------------------------------------
// A state and the bias of its sensor: the model whose covariance the float filter lost
// ------------------------------------------------------------------------------------------------

constexpr double kReadingVariance = 1e-6;

/// x₀ decays slowly and is driven; x₁ is the sensor's bias; the one reading is their sum.
template<typename T>
StateSpace<2, 1, 1, T, 2, 1> bias_model() {
    StateSpace<2, 1, 1, T, 2, 1> sys{};
    sys.A = {{static_cast<T>(0.9999), T{0}}, {T{0}, T{1}}};
    sys.B = {{static_cast<T>(1e-4)}, {T{0}}};
    sys.C = {{T{1}, T{1}}};
    sys.G = Matrix<2, 2, T>::identity();
    sys.H = Matrix<1, 1, T>::identity();
    sys.Ts = T{1};
    return sys;
}

template<typename T>
Matrix<2, 2, T> bias_process_noise() {
    return Matrix<2, 2, T>::diagonal({static_cast<T>(1e-8), static_cast<T>(1e-12)});
}

struct BiasRun {
    int    refused{0};      ///< Updates that returned false
    bool   factor_ok{true}; ///< The factor finite and triangular after every tick
    double output_rms{0.0}; ///< RMS of C (x̂ − x) over the run
};

/// 20 000 ticks of predict(u) then update(y) from P₀ = ratio·R·I, on readings every scalar sees alike.
template<typename T>
BiasRun run_bias_filter(double ratio) {
    const auto                     sys = bias_model<T>();
    KalmanFilter<2, 1, 1, T, 2, 1> kf{
        sys,
        {.Q = bias_process_noise<T>(), .R = static_cast<T>(kReadingVariance), .P0 = Matrix<2, 2, T>::identity() * static_cast<T>(ratio * kReadingVariance)},
    };
    Normal       noise{20261114ULL};
    double       x = 0.3;
    const double bias = 0.05;
    BiasRun      out{};
    double       sum = 0.0;
    const int    ticks = 20000;
    for (int k = 0; k < ticks; ++k) {
        const double u = std::sin(1e-3 * k);
        x = (0.9999 * x) + (1e-4 * u) + (1e-4 * noise());
        const double y = x + bias + (std::sqrt(kReadingVariance) * noise());
        kf.predict(ColVec<1, T>{static_cast<T>(u)});
        if (!kf.update(ColVec<1, T>{static_cast<T>(y)})) {
            ++out.refused;
        }
        const auto& S = kf.covariance_factor();
        out.factor_ok = out.factor_ok && all_finite(S) && is_factor(S);
        const double e = (static_cast<double>(kf.state()[0]) + static_cast<double>(kf.state()[1])) - (x + bias);
        sum += e * e;
    }
    out.output_rms = std::sqrt(sum / ticks);
    return out;
}

// ------------------------------------------------------------------------------------------------
// A three-state model for the comparisons with the covariance form
// ------------------------------------------------------------------------------------------------

/// A double integrator and a lag behind it, the first state read.
template<typename T>
StateSpace<3, 1, 1, T, 3, 1> chain_model() {
    StateSpace<3, 1, 1, T, 3, 1> sys{};
    sys.A = {{T{1}, static_cast<T>(0.1), static_cast<T>(0.005)}, {T{0}, T{1}, static_cast<T>(0.1)}, {T{0}, T{0}, static_cast<T>(0.95)}};
    sys.B = {{T{0}}, {static_cast<T>(0.005)}, {static_cast<T>(0.1)}};
    sys.C = {{T{1}, T{0}, T{0}}};
    sys.G = Matrix<3, 3, T>::identity();
    sys.H = Matrix<1, 1, T>::identity();
    sys.Ts = static_cast<T>(0.1);
    return sys;
}

template<typename T>
Matrix<3, 3, T> chain_process_noise() {
    // Not diagonal: the factor of Q has an entry below its diagonal.
    Matrix<3, 3, T> Q = Matrix<3, 3, T>::diagonal({static_cast<T>(1e-4), static_cast<T>(1e-3), static_cast<T>(1e-2)});
    Q(2, 1) = static_cast<T>(1e-3);
    Q(1, 2) = static_cast<T>(1e-3);
    return Q;
}

/// A mass on a spring, continuous.
constexpr StateSpace<2, 1, 1, double, 2, 1> mass_spring() {
    StateSpace<2, 1, 1, double, 2, 1> plant{};
    plant.A = {{0.0, 1.0}, {-4.0, -0.4}};
    plant.B = {{0.0}, {1.0}};
    plant.C = {{1.0, 0.0}};
    return plant;
}

/// The Kalman design on that plant with its input 2.2 samples late: [x; v(k−1); v(k−2); v(k−3)].
inline constexpr auto kDelayedDesign = design::kalman(mass_spring(), 0.01, design::KalmanNoise{.process = Matrix<2, 2>{{1e-4, 0.0}, {0.0, 1e-2}}, .measurement = 1e-4}, InputDelay<3>{0.022});

inline constexpr auto kSecondSensorPattern = SparsityPattern<2, 3>{}.unit_diagonal(0, 0, 1).values(1, 1, 1, 2); // [1 0 0; 0 · ·]

/// Which reading a comparison run takes each tick, after the model's own.
enum class Reading {
    own,            ///< update(y) alone: the model's dense C, a scalar
    vector_diag,    ///< update(y, C, R): two rows, a diagonal R
    vector_full,    ///< update(y, C, R): two rows, a full R
    pattern_scalar, ///< update(y, C, R): one row of a C on its pattern
    pattern_vector, ///< update(y, C, R): two rows of a C on its pattern, a full R
    feedthrough,    ///< update(y, C, D, R): two rows and the held input
};

/// The worst gaps to the reference over a run, in units of eps<T>() times the scale of the data.
struct Agreement {
    double covariance{0.0}; ///< max |P − P_ref| / (eps · max |P_ref|)
    double state{0.0};      ///< max |x̂ − x̂_ref| / (eps · (max |x̂_ref| + √max P_ref))
    bool   every_update{true};
};

template<typename T>
Agreement run_against_joseph(Reading reading) {
    const auto                     sys = chain_model<T>();
    const Matrix<3, 3, T>          Q = chain_process_noise<T>();
    const Matrix<3, 3, T>          P0 = Matrix<3, 3, T>::diagonal({T{1}, T{2}, static_cast<T>(0.5)});
    const ColVec<3, T>             x0{static_cast<T>(0.2), static_cast<T>(-0.1), T{0}};
    const T                        r_own = static_cast<T>(0.04);
    const T                        r_row = static_cast<T>(0.09);
    KalmanFilter<3, 1, 1, T, 3, 1> kf{sys, {.Q = Q, .R = r_own, .x0 = x0, .P0 = P0}};
    auto                           ref = reference_of(sys, Q, x0, P0);

    const Matrix<2, 3, T>                  C2{{T{1}, static_cast<T>(0.5), T{0}}, {T{0}, T{1}, static_cast<T>(-0.3)}};
    const Matrix<2, 1, T>                  D2{{static_cast<T>(0.25)}, {T{0}}};
    const Matrix<2, 2, T>                  R_diag = Matrix<2, 2, T>::diagonal({static_cast<T>(0.09), static_cast<T>(0.25)});
    const Matrix<2, 2, T>                  R_full{{static_cast<T>(0.09), static_cast<T>(0.06)}, {static_cast<T>(0.06), static_cast<T>(0.25)}};
    PatternMatrix<kSecondSensorPattern, T> Cp{};
    Cp.template at<1, 1>() = static_cast<T>(0.8);
    Cp.template at<1, 2>() = static_cast<T>(-0.4);
    const Matrix<2, 3, T>         Cp_dense = Cp.to_matrix();
    const Matrix<1, 3, T>         Cp_row{{T{0}, static_cast<T>(0.8), static_cast<T>(-0.4)}};
    static constexpr auto         kRowPattern = SparsityPattern<1, 3>{}.values(0, 1, 1, 2);
    PatternMatrix<kRowPattern, T> Cp1{};
    Cp1.template at<0, 1>() = static_cast<T>(0.8);
    Cp1.template at<0, 2>() = static_cast<T>(-0.4);

    Normal     noise{77ULL};
    ColVec<3>  truth = x0.template as<double>();
    Agreement  out{};
    const auto note = [&] {
        const double scale = largest(ref.P);
        out.covariance = std::fmax(out.covariance, gap(kf.covariance(), ref.P) / (static_cast<double>(eps<T>()) * scale));
        out.state = std::fmax(out.state, gap(kf.state(), ref.x) / (static_cast<double>(eps<T>()) * (largest(ref.x) + std::sqrt(scale))));
    };
    for (int k = 0; k < 300; ++k) {
        const T         u = static_cast<T>(std::sin(0.05 * k));
        const ColVec<1> ud{static_cast<double>(u)};
        truth = ColVec<3>(ref.A * truth) + ColVec<3>(ref.B * ud) + ColVec<3>{0.01 * noise(), 0.03 * noise(), 0.1 * noise()};
        kf.predict(ColVec<1, T>{u});
        ref.predict(ud);
        note();
        const auto two = [&](const Matrix<2, 3, T>& C, double feed) {
            const ColVec<2> y = ColVec<2>(C.template as<double>() * truth) + ColVec<2>{(0.3 * noise()) + feed, 0.5 * noise()};
            return y.template as<T>();
        };
        const ColVec<1, T> y_own{static_cast<T>(truth[0] + (0.2 * noise()))};
        bool               ok = kf.update(y_own) && ref.update(y_own.template as<double>(), sys.C.template as<double>(), Matrix<1, 1>{static_cast<double>(r_own)});
        note();
        switch (reading) {
            case Reading::own:
                break;
            case Reading::vector_diag: {
                const ColVec<2, T> y = two(C2, 0.0);
                ok = ok && kf.update(y, C2, R_diag) && ref.update(y.template as<double>(), C2.template as<double>(), R_diag.template as<double>());
                break;
            }
            case Reading::vector_full: {
                const ColVec<2, T> y = two(C2, 0.0);
                ok = ok && kf.update(y, C2, R_full) && ref.update(y.template as<double>(), C2.template as<double>(), R_full.template as<double>());
                break;
            }
            case Reading::pattern_scalar: {
                const ColVec<1, T> y{static_cast<T>((0.8 * truth[1]) - (0.4 * truth[2]) + (0.3 * noise()))};
                ok = ok && kf.update(y, Cp1, r_row) && ref.update(y.template as<double>(), Cp_row.template as<double>(), Matrix<1, 1>{static_cast<double>(r_row)});
                break;
            }
            case Reading::pattern_vector: {
                const ColVec<2, T> y = two(Cp_dense, 0.0);
                ok = ok && kf.update(y, Cp, R_full) && ref.update(y.template as<double>(), Cp_dense.template as<double>(), R_full.template as<double>());
                break;
            }
            case Reading::feedthrough: {
                const ColVec<2, T> y = two(C2, 0.25 * static_cast<double>(u));
                const ColVec<2>    y_less_feed = y.template as<double>() - ColVec<2>(D2.template as<double>() * ud);
                ok = ok && kf.update(y, C2, D2, R_full) && ref.update(y_less_feed, C2.template as<double>(), R_full.template as<double>());
                break;
            }
        }
        out.every_update = out.every_update && ok;
        note();
    }
    return out;
}

// ------------------------------------------------------------------------------------------------
// The same state and bias through the ExtendedKalmanFilter
// ------------------------------------------------------------------------------------------------

/// ẋ₀ = −x₀ + u and a constant bias x₁, the same model in either scalar.
template<typename T>
struct DecayAndBias {
    [[nodiscard]] ColVec<2, T> operator()(const ColVec<2, T>& x, const ColVec<1, T>& u) const { return ColVec<2, T>{u[0] - x[0], T{0}}; }
};

/// The reading: the state plus the sensor's bias.
template<typename T>
struct BiasedReading {
    [[nodiscard]] ColVec<1, T>                    operator()(const ColVec<2, T>& x, const ColVec<1, T>& /*u*/) const { return ColVec<1, T>{x[0] + x[1]}; }
    [[nodiscard]] MeasurementJacobian<2, 1, 1, T> jacobian(const ColVec<2, T>& /*x*/, const ColVec<1, T>& /*u*/) const {
        return {.C = Matrix<1, 2, T>{{T{1}, T{1}}}};
    }
};

struct EkfRun {
    int    refused{0};
    bool   factor_ok{true};
    double output_rms{0.0};
};

template<typename T>
EkfRun run_ekf(double ratio) {
    const double                                         r = 1e-6;
    const ExtendedKalmanFilterConfig<DecayAndBias<T>, T> config{
        .Ts = static_cast<T>(0.01),
        .Q = Matrix<2, 2, T>::diagonal({static_cast<T>(1e-8), static_cast<T>(1e-12)}),
        .P0 = Matrix<2, 2, T>::identity() * static_cast<T>(ratio * r),
    };
    ExtendedKalmanFilter<DecayAndBias<T>, T> ekf{DecayAndBias<T>{}, config};
    Normal                                   noise{555ULL};
    double                                   x = 0.3;
    const double                             bias = 0.05;
    const double                             decay = std::exp(-0.01);
    EkfRun                                   out{};
    double                                   sum = 0.0;
    const int                                ticks = 5000;
    for (int k = 0; k < ticks; ++k) {
        const double u = std::sin(2e-3 * k);
        x = (decay * x) + ((1.0 - decay) * u) + (1e-4 * noise());
        const double y = x + bias + (std::sqrt(r) * noise());
        ekf.predict(ColVec<1, T>{static_cast<T>(u)});
        if (!ekf.update(BiasedReading<T>{}, ColVec<1, T>{static_cast<T>(y)}, static_cast<T>(r))) {
            ++out.refused;
        }
        out.factor_ok = out.factor_ok && all_finite(ekf.covariance_factor()) && is_factor(ekf.covariance_factor());
        const double e = (static_cast<double>(ekf.state()[0]) + static_cast<double>(ekf.state()[1])) - (x + bias);
        sum += e * e;
    }
    out.output_rms = std::sqrt(sum / ticks);
    return out;
}

// ------------------------------------------------------------------------------------------------
// Perfect measurements: a reading, or a row of one, with no noise
// ------------------------------------------------------------------------------------------------

/// A chain filter five ticks in: a full factor, every state uncertain.
template<typename T>
KalmanFilter<3, 1, 1, T, 3, 1> settled_chain_filter() {
    KalmanFilter<3, 1, 1, T, 3, 1> kf{chain_model<T>(), {.Q = chain_process_noise<T>(), .R = static_cast<T>(0.04), .P0 = Matrix<3, 3, T>::identity()}};
    for (int k = 0; k < 5; ++k) {
        kf.predict(ColVec<1, T>{static_cast<T>(0.3)});
        REQUIRE(kf.update(ColVec<1, T>{static_cast<T>(0.1 * k)}));
    }
    return kf;
}

/// The Joseph form in double from where a filter stands now: for one update.
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
JosephReference<NX, NU> reference_at(const KalmanFilter<NX, NU, NY, T, NW, NV>& kf) {
    JosephReference<NX, NU> ref{};
    ref.x = kf.state().template as<double>();
    ref.P = detail::covariance_of_factor(kf.covariance_factor().template as<double>());
    return ref;
}

/// |C x̂ − y| of row M, in eps<T>() of the row's terms |y| + Σ|C x̂|.
template<size_t M, size_t NX, size_t NY, typename T>
double pin_error(const Matrix<NY, NX, T>& C, const ColVec<NX, T>& x, const ColVec<NY, T>& y) {
    static_assert(M < NY);
    double sum = 0.0;
    double terms = std::fabs(static_cast<double>(y[M]));
    for (size_t i = 0; i < NX; ++i) {
        const double t = static_cast<double>(C(M, i)) * static_cast<double>(x[i]);
        sum += t;
        terms += std::fabs(t);
    }
    return std::fabs(sum - static_cast<double>(y[M])) / (static_cast<double>(eps<T>()) * terms);
}

/// The variance of a one-row reading C x in a filter's factor, ‖Sᵀcᵀ‖², in double.
template<size_t NX, typename T>
double reading_variance(const Matrix<1, NX, T>& C, const Matrix<NX, NX, T>& S) {
    double sum = 0.0;
    for (size_t k = 0; k < NX; ++k) {
        double f = 0.0;
        for (size_t i = k; i < NX; ++i) {
            f += static_cast<double>(C(0, i)) * static_cast<double>(S(i, k));
        }
        sum += f * f;
    }
    return sum;
}

/// The gaps of a filter to the Joseph form after the same update, in eps<T>() of the data.
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
double gap_to(const KalmanFilter<NX, NU, NY, T, NW, NV>& kf, const JosephReference<NX, NU>& ref) {
    const auto   e = static_cast<double>(eps<T>());
    const double scale = largest(ref.P);
    return std::fmax(gap(kf.covariance(), ref.P) / (e * scale), gap(kf.state(), ref.x) / (e * (largest(ref.x) + std::sqrt(scale))));
}

} // namespace

TEST_SUITE("KalmanFilter on the factor of its covariance") {
    TEST_CASE_TEMPLATE("a prior up to 10¹² times the reading's variance: no update refused, the double filter's error", T, float, double) {
        // The covariance form in float lost P on this model from P₀/R ≈ 10⁵: indefinite within
        // three ticks, and every reading after that refused (double from 10⁹).
        for (int decade = 2; decade <= 12; ++decade) {
            const double  ratio = std::pow(10.0, decade);
            const BiasRun run = run_bias_filter<T>(ratio);
            const BiasRun in_double = run_bias_filter<double>(ratio);
            CAPTURE(decade);
            CHECK(run.refused == 0);
            CHECK(run.factor_ok);
            // The reading's noise is 1e-3 and the filter averages it down; the float filter's
            // output error is the double one's to within 5 %.
            CHECK(in_double.output_rms < 5e-4);
            CHECK(run.output_rms <= 1.05 * in_double.output_rms);
        }
    }

    TEST_CASE_TEMPLATE("where the covariance form was sound the two agree, tick for tick", T, float, double) {
        // Three states, the first read every tick, then a second reading: dense and patterned C,
        // one row and two, diagonal and full R, with and without feedthrough. The reference is
        // the Joseph form in double on the filter's own rounded model, so the gap is the filter's
        // arithmetic. Over 300 ticks it is at most 10.5 eps·‖P‖ on the covariance and
        // 5.1 eps·(‖x̂‖ + √‖P‖) on the estimate, in float and in double, with and without
        // -ffast-math and fused multiply-adds; the bound is 64 eps.
        const double bound = 64.0;
        for (const Reading reading :
             {Reading::own, Reading::vector_diag, Reading::vector_full, Reading::pattern_scalar, Reading::pattern_vector, Reading::feedthrough}) {
            const Agreement a = run_against_joseph<T>(reading);
            CAPTURE(static_cast<int>(reading));
            CAPTURE(a.covariance);
            CAPTURE(a.state);
            CHECK(a.every_update);
            CHECK(a.covariance <= bound);
            CHECK(a.state <= bound);
        }
    }
}

TEST_SUITE("KalmanFilter on the factor: the model's statistics") {
    TEST_CASE_TEMPLATE("from 100 P∞ the running covariance converges to the design's steady state", T, float, double) {
        const auto sys = chain_model<double>();
        const auto design = design::kalman(sys, design::KalmanNoise{.process = chain_process_noise<double>(), .measurement = 0.04});
        REQUIRE(design.success);
        KalmanFilter<3, 1, 1, T, 3, 1> kf{
            sys.template as<T>(),
            {.Q = chain_process_noise<double>().template as<T>(), .R = static_cast<T>(0.04), .P0 = Matrix<3, 3>(design.P * 100.0).template as<T>()},
        };
        CHECK(gap(kf.covariance(), design.P * 100.0) <= 64.0 * static_cast<double>(eps<T>()) * largest(design.P) * 100.0);
        // P does not depend on the readings: P∞ is the a priori covariance, the one after predict.
        for (int k = 0; k < 400; ++k) {
            REQUIRE(kf.update(ColVec<1, T>{T{0}}));
            kf.predict(ColVec<1, T>{T{0}});
        }
        // The filter's model and noise are the design's rounded to T, so its fixed point is P∞ to
        // that rounding through the Riccati map: at most 6.9 eps·‖P∞‖ measured, in float and in
        // double, against a bound of 64 eps.
        const double relative = gap(kf.covariance(), design.P) / (static_cast<double>(eps<T>()) * largest(design.P));
        CAPTURE(relative);
        CHECK(relative <= 64.0);
    }

    TEST_CASE_TEMPLATE("with the model's own noise the error is the size its covariance says (NEES)", T, float, double) {
        // 200 independent runs of 60 ticks from x₀ ~ N(x̂₀, P₀), process and reading noise the
        // model's own, seed 4242. At the last update e = x − x̂ and NEES = eᵀP⁻¹e = ‖S⁻¹e‖², a
        // χ² of 3 degrees of freedom when the filter is consistent. The mean over the runs is
        // then χ²(600)/200: mean 3, standard deviation √(2·3/200) = 0.173, and the band is three
        // of those, [2.48, 3.52].
        const auto            sys = chain_model<T>();
        const Matrix<3, 3, T> Q = chain_process_noise<T>();
        const T               r = static_cast<T>(0.04);
        const Matrix<3, 3, T> P0 = Matrix<3, 3, T>::diagonal({T{1}, T{2}, static_cast<T>(0.5)});
        const ColVec<3, T>    x0{static_cast<T>(0.2), static_cast<T>(-0.1), T{0}};
        const Matrix<3, 3>    A = sys.A.template as<double>();
        const Matrix<3, 1>    B = sys.B.template as<double>();
        Matrix<3, 3>          Lq{};
        REQUIRE(detail::semidefinite_factor(Q.template as<double>(), Lq));
        Normal     noise{4242ULL};
        const auto draw = [&] { return ColVec<3>{noise(), noise(), noise()}; };
        const int  runs = 200;
        double     sum = 0.0;
        for (int run = 0; run < runs; ++run) {
            KalmanFilter<3, 1, 1, T, 3, 1> kf{sys, {.Q = Q, .R = r, .x0 = x0, .P0 = P0}};
            const ColVec<3>                n0 = draw();
            ColVec<3>                      truth = x0.template as<double>() + ColVec<3>{n0[0], std::numbers::sqrt2 * n0[1], (1.0 / std::numbers::sqrt2) * n0[2]};
            for (int k = 0; k < 60; ++k) {
                const double u = std::sin(0.05 * k);
                truth = ColVec<3>(A * truth) + ColVec<3>(B * ColVec<1>{u}) + ColVec<3>(Lq * draw());
                const double y = truth[0] + (std::sqrt(static_cast<double>(r)) * noise());
                kf.predict(ColVec<1, T>{static_cast<T>(u)});
                REQUIRE(kf.update(ColVec<1, T>{static_cast<T>(y)}));
            }
            // z = S⁻¹ e by forward substitution on the factor
            const Matrix<3, 3> S = kf.covariance_factor().template as<double>();
            ColVec<3>          z = truth - kf.state().template as<double>();
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < i; ++j) {
                    z[i] -= S(i, j) * z[j];
                }
                z[i] /= S(i, i);
                sum += z[i] * z[i];
            }
        }
        const double nees = sum / runs;
        CAPTURE(nees);
        CHECK(nees >= 2.48);
        CHECK(nees <= 3.52);
    }
}

TEST_SUITE("KalmanFilter on the factor: singular covariances and refused readings") {
    TEST_CASE_TEMPLATE("P₀ of all zeros: the filter starts certain and the process noise opens it", T, float, double) {
        const auto                     sys = chain_model<T>();
        const Matrix<3, 3, T>          Q = chain_process_noise<T>();
        const T                        r = static_cast<T>(0.04);
        KalmanFilter<3, 1, 1, T, 3, 1> kf{sys, {.Q = Q, .R = r, .P0 = Matrix<3, 3, T>::zeros()}};
        auto                           ref = reference_of(sys, Q, ColVec<3, T>{}, Matrix<3, 3, T>::zeros());
        CHECK(largest(kf.covariance_factor()) == 0.0);
        // Certain of x̂ = 0: a reading moves nothing, and is not refused (R is positive definite).
        REQUIRE(kf.update(ColVec<1, T>{T{1}}));
        CHECK(largest(kf.state()) == 0.0);
        CHECK(largest(kf.covariance_factor()) == 0.0);
        double worst = 0.0;
        for (int k = 0; k < 50; ++k) {
            const T y = static_cast<T>(0.1 * std::cos(0.3 * k));
            kf.predict(ColVec<1, T>{static_cast<T>(0.5)});
            ref.predict(ColVec<1>{0.5});
            REQUIRE(kf.update(ColVec<1, T>{y}));
            REQUIRE(ref.update(ColVec<1>{static_cast<double>(y)}, sys.C.template as<double>(), Matrix<1, 1>{static_cast<double>(r)}));
            worst = std::fmax(worst, gap(kf.covariance(), ref.P) / (static_cast<double>(eps<T>()) * largest(ref.P)));
            CHECK(is_factor(kf.covariance_factor()));
        }
        CAPTURE(worst);
        CHECK(worst <= 64.0); // eps·‖P‖ against the Joseph form in double: 8.1 measured
    }

    TEST_CASE_TEMPLATE("a state known exactly, with no noise on it: its row of the factor stays zero and no reading moves it", T, float, double) {
        // The bias model with the bias known: P₀ = diag(p, 0) and Q = diag(q, 0).
        const auto                     sys = bias_model<T>();
        const Matrix<2, 2, T>          Q = Matrix<2, 2, T>::diagonal({static_cast<T>(1e-8), T{0}});
        const Matrix<2, 2, T>          P0 = Matrix<2, 2, T>::diagonal({static_cast<T>(1e-2), T{0}});
        const ColVec<2, T>             x0{T{0}, static_cast<T>(0.05)};
        const T                        r = static_cast<T>(kReadingVariance);
        KalmanFilter<2, 1, 1, T, 2, 1> kf{sys, {.Q = Q, .R = r, .x0 = x0, .P0 = P0}};
        auto                           ref = reference_of(sys, Q, x0, P0);
        Normal                         noise{9ULL};
        double                         worst = 0.0;
        for (int k = 0; k < 500; ++k) {
            const T u = static_cast<T>(std::sin(1e-2 * k));
            const T y = static_cast<T>(0.35 + (1e-3 * noise()));
            kf.predict(ColVec<1, T>{u});
            ref.predict(ColVec<1>{static_cast<double>(u)});
            REQUIRE(kf.update(ColVec<1, T>{y}));
            REQUIRE(ref.update(ColVec<1>{static_cast<double>(y)}, sys.C.template as<double>(), Matrix<1, 1>{static_cast<double>(r)}));
            const auto& S = kf.covariance_factor();
            REQUIRE(S(1, 0) == T{0});
            REQUIRE(S(1, 1) == T{0});
            REQUIRE(kf.state()[1] == x0[1]); // exactly
            worst = std::fmax(worst, std::fabs(static_cast<double>(kf.covariance()(0, 0)) - ref.P(0, 0)) / (static_cast<double>(eps<T>()) * ref.P(0, 0)));
        }
        CAPTURE(worst);
        CHECK(worst <= 64.0); // eps·P₀₀ against the Joseph form in double: 12.7 measured
    }

    TEST_CASE_TEMPLATE("a delay-register model, its Q singular: the registers hold the inputs exactly, with zero rows in the factor", T, float, double) {
        // A mass on a spring whose input arrives 1.5 samples late: [x; v(k−1); v(k−2)], the
        // process noise on the plant only (G has no row for a register).
        const auto model = *discretize(mass_spring(), 0.01, InputDelay<2>{0.015});
        const auto sys = model.template as<T>();
        static_assert(std::is_same_v<decltype(sys), const StateSpace<4, 1, 1, T, 2, 1>>);
        const Matrix<2, 2, T> Q = Matrix<2, 2, T>::diagonal({static_cast<T>(1e-4), static_cast<T>(1e-2)});
        Matrix<4, 4, T>       P0{};
        P0(0, 0) = static_cast<T>(1e-2);
        P0(1, 1) = static_cast<T>(1e-1);
        const T                        r = static_cast<T>(1e-4);
        KalmanFilter<4, 1, 1, T, 2, 1> kf{sys, {.Q = Q, .R = r, .P0 = P0}};
        auto                           ref = reference_of(sys, Q, ColVec<4, T>{}, P0);
        Normal                         noise{31ULL};
        double                         worst = 0.0;
        T                              u_before = T{0};
        for (int k = 0; k < 300; ++k) {
            const T u = static_cast<T>(std::sin(0.07 * k));
            const T y = static_cast<T>(0.1 * std::sin(0.02 * k) + (1e-2 * noise()));
            kf.predict(ColVec<1, T>{u});
            ref.predict(ColVec<1>{static_cast<double>(u)});
            REQUIRE(kf.update(ColVec<1, T>{y}));
            REQUIRE(ref.update(ColVec<1>{static_cast<double>(y)}, sys.C.template as<double>(), Matrix<1, 1>{static_cast<double>(r)}));
            const auto& S = kf.covariance_factor();
            for (size_t j = 0; j < 4; ++j) {
                REQUIRE(S(2, j) == T{0});
                REQUIRE(S(3, j) == T{0});
            }
            REQUIRE(kf.state()[2] == u);        // v(k−1): the input just applied, exactly
            REQUIRE(kf.state()[3] == u_before); // v(k−2)
            u_before = u;
            worst = std::fmax(worst, gap(kf.covariance(), ref.P) / (static_cast<double>(eps<T>()) * largest(ref.P)));
        }
        CAPTURE(worst);
        CHECK(worst <= 128.0); // eps·‖P‖ against the Joseph form in double, four states: 25.9 measured
    }

    TEST_CASE_TEMPLATE("Law<kRes, T>::Recursive, the filter on a design's constants: the same factor, the registers' zeros known at compile time", T, float, double) {
        // The design's steady state seeds both filters; the registers have no covariance in it.
        static_assert(kDelayedDesign.success);
        typename Law<kDelayedDesign, T>::Recursive on_constants;
        KalmanFilter<5, 1, 1, T, 2, 1>             held{kDelayedDesign};
        static_assert(sizeof(on_constants) < sizeof(held));
        const auto sys = kDelayedDesign.sys.template as<T>();
        auto       ref = reference_of(sys, kDelayedDesign.Q.template as<T>(), ColVec<5, T>{}, kDelayedDesign.P.template as<T>());
        const T    r = static_cast<T>(kDelayedDesign.R(0, 0));
        // One factor at compile time, the other at run time, from the same P∞ in double: the same to a few ulps.
        CHECK(gap(on_constants.covariance_factor(), held.covariance_factor()) <= 4.0 * static_cast<double>(eps<T>()) * largest(held.covariance_factor()));
        Normal noise{63ULL};
        double worst = 0.0;
        double apart = 0.0;
        for (int k = 0; k < 300; ++k) {
            const T u = static_cast<T>(std::sin(0.07 * k));
            const T y = static_cast<T>((0.1 * std::sin(0.02 * k)) + (1e-2 * noise()));
            on_constants.predict(ColVec<1, T>{u});
            held.predict(ColVec<1, T>{u});
            ref.predict(ColVec<1>{static_cast<double>(u)});
            REQUIRE(on_constants.update(ColVec<1, T>{y}));
            REQUIRE(held.update(ColVec<1, T>{y}));
            REQUIRE(ref.update(ColVec<1>{static_cast<double>(y)}, sys.C.template as<double>(), Matrix<1, 1>{static_cast<double>(r)}));
            const auto& S = on_constants.covariance_factor();
            REQUIRE(is_factor(S));
            for (size_t i = 2; i < 5; ++i) {
                for (size_t j = 0; j < 5; ++j) {
                    REQUIRE(S(i, j) == T{0});
                }
            }
            REQUIRE(on_constants.state()[2] == u); // v(k−1), exactly
            worst = std::fmax(worst, gap(on_constants.covariance(), ref.P) / (static_cast<double>(eps<T>()) * largest(ref.P)));
            apart = std::fmax(apart, gap(on_constants.covariance(), held.covariance()) / (static_cast<double>(eps<T>()) * largest(ref.P)));
            apart = std::fmax(apart, gap(on_constants.state(), held.state()) / (static_cast<double>(eps<T>()) * (largest(ref.x) + std::sqrt(largest(ref.P)))));
        }
        CAPTURE(worst);
        CAPTURE(apart);
        CHECK(worst <= 128.0); // eps·‖P‖ against the Joseph form in double, five states: 11.4 measured
        CHECK(apart <= 64.0);  // the two stores differ only in the order of their sums: 8.4 measured
    }

    TEST_CASE_TEMPLATE("a reading that cannot be used returns false and leaves the estimate and the factor bit for bit", T, float, double) {
        const auto                     sys = chain_model<T>();
        KalmanFilter<3, 1, 1, T, 3, 1> kf{sys, {.Q = chain_process_noise<T>(), .R = static_cast<T>(0.04), .P0 = Matrix<3, 3, T>::identity()}};
        for (int k = 0; k < 5; ++k) {
            kf.predict(ColVec<1, T>{static_cast<T>(0.3)});
            REQUIRE(kf.update(ColVec<1, T>{static_cast<T>(0.1 * k)}));
        }
        const ColVec<3, T>    x_before = kf.state();
        const Matrix<3, 3, T> S_before = kf.covariance_factor();
        const auto            untouched = [&] { return same_bits(kf.state(), x_before) && same_bits(kf.covariance_factor(), S_before); };

        const Matrix<1, 3, T> C1{{T{1}, T{0}, static_cast<T>(0.5)}};
        const ColVec<1, T>    y1{static_cast<T>(0.4)};
        const Matrix<2, 3, T> C2{{T{1}, static_cast<T>(0.5), T{0}}, {T{0}, T{1}, static_cast<T>(-0.3)}};
        const ColVec<2, T>    y2{static_cast<T>(0.4), static_cast<T>(-0.2)};
        const Matrix<2, 2, T> R2{{static_cast<T>(0.09), static_cast<T>(0.06)}, {static_cast<T>(0.06), static_cast<T>(0.25)}};

        SUBCASE("R that is not positive semidefinite: negative, indefinite, a zero variance with a covariance") {
            CHECK_FALSE(kf.update(y1, C1, static_cast<T>(-0.04)));
            CHECK(untouched());
            CHECK_FALSE(kf.update(y2, C2, Matrix<2, 2, T>{{T{1}, T{2}}, {T{2}, T{1}}}));
            CHECK(untouched());
            CHECK_FALSE(kf.update(y2, C2, Matrix<2, 2, T>{{static_cast<T>(0.09), static_cast<T>(0.06)}, {static_cast<T>(0.06), T{0}}}));
            CHECK(untouched());
            CHECK_FALSE(kf.update(y2, C2, Matrix<2, 2, T>{{T{0}, static_cast<T>(0.06)}, {static_cast<T>(0.06), static_cast<T>(0.25)}}));
            CHECK(untouched());
        }
        SUBCASE("R that is not finite") {
            for (const T bad : {runtime_nan<T>(), runtime_inf<T>()}) {
                CHECK_FALSE(kf.update(y1, C1, bad));
                CHECK(untouched());
                Matrix<2, 2, T> R = R2;
                R(1, 0) = bad;
                R(0, 1) = bad;
                CHECK_FALSE(kf.update(y2, C2, R));
                CHECK(untouched());
                R = R2;
                R(1, 1) = bad;
                CHECK_FALSE(kf.update(y2, C2, R));
                CHECK(untouched());
            }
        }
        SUBCASE("C that is not finite, in the second row of two: the first row is not applied either") {
            for (const T bad : {runtime_nan<T>(), runtime_inf<T>()}) {
                Matrix<2, 3, T> C = C2;
                C(1, 2) = bad;
                CHECK_FALSE(kf.update(y2, C, R2));
                CHECK(untouched());
                PatternMatrix<kSecondSensorPattern, T> Cp{};
                Cp.template at<1, 1>() = static_cast<T>(0.8);
                Cp.template at<1, 2>() = bad;
                CHECK_FALSE(kf.update(y2, Cp, R2));
                CHECK(untouched());
            }
        }
        SUBCASE("an innovation variance past the scalar's range, in the second row of two") {
            Matrix<2, 3, T> C = C2;
            C(1, 1) = damp::numeric_limits<T>::max() / T{4}; // finite, and its square is not
            CHECK_FALSE(kf.update(y2, C, R2));
            CHECK(untouched());
        }
        SUBCASE("the model's own reading with an R that is not a covariance") {
            KalmanFilter<3, 1, 1, T, 3, 1> blind{sys, {.Q = chain_process_noise<T>(), .R = static_cast<T>(-0.04), .P0 = Matrix<3, 3, T>::identity()}};
            blind.predict(ColVec<1, T>{static_cast<T>(0.3)});
            const ColVec<3, T>    x = blind.state();
            const Matrix<3, 3, T> S = blind.covariance_factor();
            CHECK_FALSE(blind.update(ColVec<1, T>{T{1}}));
            CHECK(same_bits(blind.state(), x));
            CHECK(same_bits(blind.covariance_factor(), S));
            CHECK(blind.step(ColVec<1, T>{T{1}}, ColVec<1, T>{T{0}})[0] == blind.state()[0]); // the prediction stands
        }
        // The filter is still the one it was: the readings it refused are usable once valid.
        CHECK(kf.update(y1, C1, static_cast<T>(0.04)));
        CHECK(kf.update(y2, C2, R2));
        CHECK_FALSE(untouched());
    }

    TEST_CASE_TEMPLATE("a P₀ or a Q that is not a covariance: the identity and no process noise at run time", T, float, double) {
        const auto            sys = chain_model<T>();
        const Matrix<3, 3, T> indefinite{{T{1}, T{2}, T{0}}, {T{2}, T{1}, T{0}}, {T{0}, T{0}, T{1}}};
        Matrix<3, 3, T>       not_finite = Matrix<3, 3, T>::identity();
        not_finite(2, 1) = runtime_nan<T>();
        not_finite(1, 2) = runtime_nan<T>();
        for (const auto& bad : {indefinite, not_finite}) {
            KalmanFilter<3, 1, 1, T, 3, 1> kf{sys, {.Q = chain_process_noise<T>(), .R = static_cast<T>(0.04), .P0 = bad}};
            CHECK(gap(kf.covariance_factor(), Matrix<3, 3, T>::identity()) == 0.0);
            kf.reset(ColVec<3, T>{}, Matrix<3, 3, T>::identity() * T{4});
            CHECK(gap(kf.covariance_factor(), Matrix<3, 3, T>::identity() * T{2}) == 0.0);
            kf.reset(ColVec<3, T>{}, bad);
            CHECK(gap(kf.covariance_factor(), Matrix<3, 3, T>::identity()) == 0.0);

            KalmanFilter<3, 1, 1, T, 3, 1> quiet{sys, {.Q = bad, .R = static_cast<T>(0.04), .P0 = Matrix<3, 3, T>::zeros()}};
            quiet.predict(ColVec<1, T>{T{1}});
            CHECK(largest(quiet.covariance_factor()) == 0.0); // nothing added
            CHECK(all_finite(quiet.state()));
        }
    }

    TEST_CASE("the filter in another scalar carries the factor over, and reset takes a singular P₀") {
        const auto                          sys = chain_model<double>();
        KalmanFilter<3, 1, 1, double, 3, 1> kf{sys, {.Q = chain_process_noise<double>(), .R = 0.04, .P0 = Matrix<3, 3>::diagonal({1.0, 2.0, 0.5})}};
        kf.predict(ColVec<1>{0.3});
        REQUIRE(kf.update(ColVec<1>{0.2}));
        KalmanFilter<3, 1, 1, float, 3, 1> kf_f{kf};
        CHECK(gap(kf_f.covariance_factor(), kf.covariance_factor()) <= static_cast<double>(eps<float>()) * largest(kf.covariance_factor()));
        CHECK(gap(kf_f.state(), kf.state()) <= static_cast<double>(eps<float>()) * largest(kf.state()));
        kf_f.predict(ColVec<1, float>{0.3F});
        kf.predict(ColVec<1>{0.3});
        REQUIRE(kf_f.update(ColVec<1, float>{0.25F}));
        REQUIRE(kf.update(ColVec<1>{0.25}));
        CHECK(gap(kf_f.covariance(), kf.covariance()) <= 64.0 * static_cast<double>(eps<float>()) * largest(kf.covariance()));

        Matrix<3, 3> rank_one{}; // x₁ = 2 x₀ exactly, x₂ known
        rank_one(0, 0) = 1.0;
        rank_one(0, 1) = 2.0;
        rank_one(1, 0) = 2.0;
        rank_one(1, 1) = 4.0;
        kf.reset(ColVec<3>{1.0, 2.0, 3.0}, rank_one);
        CHECK(gap(kf.covariance(), rank_one) <= 4.0 * eps<double>());
        CHECK(kf.covariance_factor()(1, 1) == 0.0);
        CHECK(kf.covariance_factor()(2, 2) == 0.0);
        CHECK(kf.state()[2] == 3.0);
    }

    TEST_CASE("the filter ticks in constant evaluation") {
        static constexpr auto kTicked = [] {
            StateSpace<2, 1, 1, double, 2, 1> sys{};
            sys.A = {{1.0, 0.1}, {0.0, 1.0}};
            sys.B = {{0.005}, {0.1}};
            sys.C = {{1.0, 0.0}};
            sys.G = Matrix<2, 2>::identity();
            sys.H = Matrix<1, 1>::identity();
            sys.Ts = 0.1;
            KalmanFilter kf{sys, {.Q = Matrix<2, 2>::identity() * 0.01, .R = 0.1}};
            kf.predict(ColVec<1>{0.0});
            const bool ok = kf.update(ColVec<1>{1.0});
            return damp::make_tuple(ok, kf.covariance(), kf.covariance_factor());
        }();
        static_assert(damp::get<0>(kTicked));
        static_assert(damp::get<2>(kTicked)(0, 1) == 0.0 && damp::get<2>(kTicked)(0, 0) > 0.0 && damp::get<2>(kTicked)(1, 1) > 0.0);
        // P = A I Aᵀ + 0.01 I = [1.02 0.1; 0.1 1.01], then the reading of x₀ with R = 0.1.
        const Matrix<2, 2>& P = damp::get<1>(kTicked);
        CHECK(P(0, 0) == doctest::Approx(1.02 - ((1.02 * 1.02) / 1.12)).epsilon(1e-13));
        CHECK(P(1, 0) == doctest::Approx(0.1 - ((0.1 * 1.02) / 1.12)).epsilon(1e-13));
        CHECK(P(1, 1) == doctest::Approx(1.01 - ((0.1 * 0.1) / 1.12)).epsilon(1e-13));
    }
}

TEST_SUITE("KalmanFilter on the factor: perfect measurements") {
    // A reading with no noise (R = 0, or a zero-variance row of R) is taken at the rounding of
    // its own prior variance: the estimate reads it to eps, and the factor keeps eps times the
    // prior variance in what it reads. The covariance form took these while C P Cᵀ + R was
    // positive definite, and so does the factor.
    //
    // How closely a perfect row is read, |C x̂ − y| in eps of its terms: at most 6.7 measured,
    // in float and in double, with and without -ffast-math and fused multiply-adds.
    constexpr double kPinned = 16.0;

    TEST_CASE_TEMPLATE("R = 0, one reading: the estimate reads it to a few eps, and the factor keeps eps of the prior variance there", T, float, double) {
        auto                  kf = settled_chain_filter<T>();
        const Matrix<1, 3, T> C{{T{1}, T{0}, static_cast<T>(0.5)}};
        const ColVec<1, T>    y{static_cast<T>(0.4)};
        const double          prior = reading_variance(C, kf.covariance_factor());
        auto                  ref = reference_at(kf);
        REQUIRE(kf.update(y, C, T{0}));
        REQUIRE(ref.update(y.template as<double>(), C.template as<double>(), Matrix<1, 1>{0.0}));
        CHECK(all_finite(kf.state()));
        CHECK(all_finite(kf.covariance_factor()));
        CHECK(is_factor(kf.covariance_factor()));
        CHECK(pin_error<0>(C, kf.state(), y) <= kPinned);
        // eps·prior/(1 + eps) in exact arithmetic; the factor holds it through a cancellation.
        const double left = reading_variance(C, kf.covariance_factor()) / (static_cast<double>(eps<T>()) * prior);
        CAPTURE(left);
        CHECK(left >= 0.5);
        CHECK(left <= 2.0);
        // The Joseph form with R = 0 leaves no variance there: the two differ by that eps·prior.
        // In eps·‖P‖ and eps·(‖x̂‖ + √‖P‖), here and for one perfect row of two: 6.3 measured.
        CHECK(gap_to(kf, ref) <= 64.0);
    }

    TEST_CASE_TEMPLATE("two readings and a singular R: a zero-variance row with no cross term, either row, both rows, and one noise shared by both", T, float, double) {
        const Matrix<2, 3, T> C{{T{1}, static_cast<T>(0.5), T{0}}, {T{0}, T{1}, static_cast<T>(-0.3)}};
        const ColVec<2, T>    y{static_cast<T>(0.4), static_cast<T>(-0.2)};
        const T               r0 = static_cast<T>(0.09);
        const T               r1 = static_cast<T>(0.25);
        SUBCASE("the second row perfect") {
            auto                  kf = settled_chain_filter<T>();
            auto                  ref = reference_at(kf);
            const Matrix<2, 2, T> R = Matrix<2, 2, T>::diagonal({r0, T{0}});
            REQUIRE(kf.update(y, C, R));
            REQUIRE(ref.update(y.template as<double>(), C.template as<double>(), R.template as<double>()));
            CHECK(pin_error<1>(C, kf.state(), y) <= kPinned);
            CHECK(gap_to(kf, ref) <= 64.0);
        }
        SUBCASE("the first row perfect") {
            auto                  kf = settled_chain_filter<T>();
            auto                  ref = reference_at(kf);
            const Matrix<2, 2, T> R = Matrix<2, 2, T>::diagonal({T{0}, r1});
            REQUIRE(kf.update(y, C, R));
            REQUIRE(ref.update(y.template as<double>(), C.template as<double>(), R.template as<double>()));
            CHECK(pin_error<0>(C, kf.state(), y) <= kPinned);
            CHECK(gap_to(kf, ref) <= 64.0);
        }
        SUBCASE("both rows perfect") {
            auto kf = settled_chain_filter<T>();
            auto ref = reference_at(kf);
            REQUIRE(kf.update(y, C, Matrix<2, 2, T>::zeros()));
            REQUIRE(ref.update(y.template as<double>(), C.template as<double>(), Matrix<2, 2>::zeros()));
            CHECK(pin_error<0>(C, kf.state(), y) <= kPinned);
            CHECK(pin_error<1>(C, kf.state(), y) <= kPinned);
            CHECK(gap_to(kf, ref) <= 128.0); // each row leaves eps of its prior variance, and the rows are correlated: 30.5 measured
            CHECK(is_factor(kf.covariance_factor()));
        }
        SUBCASE("one noise shared by both rows: their difference is a perfect reading") {
            auto                  kf = settled_chain_filter<T>();
            auto                  ref = reference_at(kf);
            const Matrix<2, 2, T> R{{r0, r0}, {r0, r0}};
            REQUIRE(kf.update(y, C, R));
            REQUIRE(ref.update(y.template as<double>(), C.template as<double>(), R.template as<double>()));
            const Matrix<1, 3, T> difference{{C(1, 0) - C(0, 0), C(1, 1) - C(0, 1), C(1, 2) - C(0, 2)}};
            CHECK(pin_error<0>(difference, kf.state(), ColVec<1, T>{y[1] - y[0]}) <= kPinned);
            CHECK(gap_to(kf, ref) <= 64.0);
        }
    }

    TEST_CASE_TEMPLATE("design::kalman with no measurement noise runs on KalmanFilter: every update is taken, through update and through step", T, float, double) {
        StateSpace<2, 1, 2, double, 2, 2> sys{};
        sys.A = {{0.9, 0.1}, {0.0, 0.8}};
        sys.B = {{0.0}, {1.0}};
        sys.C = Matrix<2, 2>::identity();
        sys.G = Matrix<2, 2>::identity();
        sys.H = Matrix<2, 2>::identity();
        sys.Ts = 0.01;
        const auto design = design::kalman(sys, design::KalmanNoise{.process = Matrix<2, 2>::identity() * 0.01, .measurement = Matrix<2, 2>::zeros()});
        REQUIRE(design.success);
        KalmanFilter<2, 1, 2, T, 2, 2> kf{design};
        const Matrix<2, 2, T>          C = Matrix<2, 2, T>::identity();
        for (int k = 0; k < 50; ++k) {
            const ColVec<1, T> u{static_cast<T>(std::sin(0.2 * k))};
            const ColVec<2, T> y{static_cast<T>(0.5 + (0.1 * std::cos(0.3 * k))), static_cast<T>(-0.2 + (0.05 * std::sin(0.7 * k)))};
            if (k % 2 == 0) {
                kf.predict(u);
                REQUIRE(kf.update(y));
            } else {
                const ColVec<2, T> before = kf.state();
                (void)kf.step(y, u);
                CHECK_FALSE(same_bits(kf.state(), before)); // the update was not skipped
            }
            CAPTURE(k);
            REQUIRE(pin_error<0>(C, kf.state(), y) <= kPinned);
            REQUIRE(pin_error<1>(C, kf.state(), y) <= kPinned);
            REQUIRE(is_factor(kf.covariance_factor()));
            REQUIRE(all_finite(kf.covariance_factor()));
        }
    }

    TEST_CASE_TEMPLATE("a perfect reading of what is known exactly is refused: the estimate and the factor bit for bit", T, float, double) {
        // The bias known exactly: P₀ = diag(p, 0).
        const auto                     sys = bias_model<T>();
        const Matrix<2, 2, T>          P0 = Matrix<2, 2, T>::diagonal({static_cast<T>(1e-2), T{0}});
        KalmanFilter<2, 1, 1, T, 2, 1> kf{sys, {.Q = bias_process_noise<T>(), .R = static_cast<T>(kReadingVariance), .x0 = ColVec<2, T>{T{0}, static_cast<T>(0.05)}, .P0 = P0}};
        const ColVec<2, T>             x_before = kf.state();
        const Matrix<2, 2, T>          S_before = kf.covariance_factor();
        const auto                     untouched = [&] { return same_bits(kf.state(), x_before) && same_bits(kf.covariance_factor(), S_before); };
        // One row, the known state alone.
        CHECK_FALSE(kf.update(ColVec<1, T>{static_cast<T>(0.07)}, Matrix<1, 2, T>{{T{0}, T{1}}}, T{0}));
        CHECK(untouched());
        // Two rows, the second the known state with no noise: the first is not applied either.
        CHECK_FALSE(kf.update(ColVec<2, T>{static_cast<T>(0.3), static_cast<T>(0.07)}, Matrix<2, 2, T>::identity(), Matrix<2, 2, T>::diagonal({static_cast<T>(1e-4), T{0}})));
        CHECK(untouched());
        // With noise the same reading is taken, and leaves the known state where it is.
        CHECK(kf.update(ColVec<2, T>{static_cast<T>(0.3), static_cast<T>(0.07)}, Matrix<2, 2, T>::identity(), Matrix<2, 2, T>::diagonal({static_cast<T>(1e-4), static_cast<T>(1e-4)})));
        CHECK(kf.state()[1] == x_before[1]);
    }

    TEST_CASE_TEMPLATE("a perfect reading repeated with nothing in between: once what it reads is down to rounding it is refused, and no update moves the estimate", T, float, double) {
        // x₀ + x₁ read with no noise. The first reading leaves eps of the prior variance in the
        // sum; a later one finds only the rounding of S's own entries there, which is not a
        // variance to divide by. The test for that is relative to the terms of C S.
        KalmanFilter<2, 1, 1, T, 2, 1> kf{
            bias_model<T>(),
            {.Q = bias_process_noise<T>(), .R = static_cast<T>(kReadingVariance), .P0 = Matrix<2, 2, T>::diagonal({static_cast<T>(1e-2), static_cast<T>(4e-2)})},
        };
        const Matrix<1, 2, T> C{{T{1}, T{1}}};
        const ColVec<1, T>    y{static_cast<T>(0.35)};
        REQUIRE(kf.update(y, C, T{0}));
        const ColVec<2, T> x_first = kf.state();
        CHECK(pin_error<0>(C, kf.state(), y) <= kPinned);
        bool refused = false;
        for (int repeat = 0; repeat < 3 && !refused; ++repeat) {
            const ColVec<2, T>    x_before = kf.state();
            const Matrix<2, 2, T> S_before = kf.covariance_factor();
            if (!kf.update(y, C, T{0})) {
                refused = true;
                CHECK(same_bits(kf.state(), x_before));
                CHECK(same_bits(kf.covariance_factor(), S_before));
            }
        }
        CHECK(refused);
        // No reading in the sequence moved the states apart: the innovation was at rounding.
        CHECK(gap(kf.state(), x_first) <= 64.0 * static_cast<double>(eps<T>()) * 0.35);
        CHECK(pin_error<0>(C, kf.state(), y) <= kPinned);
        CHECK(is_factor(kf.covariance_factor()));
        // Process noise opens the sum again, and the reading is taken.
        kf.predict(ColVec<1, T>{T{0}});
        CHECK(kf.update(y, C, T{0}));
    }
}

TEST_SUITE("ExtendedKalmanFilter on the factor of its covariance") {
    TEST_CASE_TEMPLATE("a prior 10¹⁰ times the reading's variance: no update refused, the double filter's error", T, float, double) {
        for (const double ratio : {1e4, 1e10}) {
            const EkfRun run = run_ekf<T>(ratio);
            const EkfRun in_double = run_ekf<double>(ratio);
            CAPTURE(ratio);
            CHECK(run.refused == 0);
            CHECK(run.factor_ok);
            CHECK(in_double.output_rms < 5e-4);
            CHECK(run.output_rms <= 1.05 * in_double.output_rms);
        }
    }

    TEST_CASE("reset takes a singular P₀, and a P₀ that is not a covariance becomes the identity") {
        ExtendedKalmanFilter<DecayAndBias<double>, double> ekf{DecayAndBias<double>{}, {.Ts = 0.01, .Q = Matrix<2, 2>::zeros(), .P0 = Matrix<2, 2>::zeros()}};
        CHECK(largest(ekf.covariance()) == 0.0);
        ekf.predict(ColVec<1>{1.0});
        CHECK(largest(ekf.covariance()) == 0.0);
        REQUIRE(ekf.update(BiasedReading<double>{}, ColVec<1>{2.0}, 0.5)); // certain: the reading moves nothing
        CHECK(ekf.state()[1] == 0.0);
        ekf.reset(ColVec<2>{}, Matrix<2, 2>{{1.0, 2.0}, {2.0, 1.0}});
        CHECK(gap(ekf.covariance(), Matrix<2, 2>::identity()) == 0.0);
        CHECK_FALSE(ekf.update(BiasedReading<double>{}, ColVec<1>{2.0}, -0.5)); // R is not a covariance
        CHECK(gap(ekf.covariance(), Matrix<2, 2>::identity()) == 0.0);
        REQUIRE(ekf.update(BiasedReading<double>{}, ColVec<1>{2.0}, 0.0)); // a perfect reading of x₀ + x₁
        CHECK(ekf.state()[0] + ekf.state()[1] == doctest::Approx(2.0).epsilon(1e-14));
        CHECK(ekf.state()[0] == doctest::Approx(1.0).epsilon(1e-14));
    }
}

TEST_SUITE("the square-root covariance kernel") {
    TEST_CASE_TEMPLATE("semidefinite_factor: S Sᵀ = P for definite, singular and zero P; an indefinite or non-finite P is refused", T, float, double) {
        const auto e = static_cast<double>(eps<T>());
        const auto factors = [&](const Matrix<3, 3, T>& P) {
            Matrix<3, 3, T> S = Matrix<3, 3, T>::identity() * T{7};
            const bool      ok = detail::semidefinite_factor(P, S);
            return ok && is_factor(S) && gap(detail::covariance_of_factor(S), P) <= 8.0 * e * largest(P);
        };
        CHECK(factors(Matrix<3, 3, T>{{T{4}, T{1}, static_cast<T>(0.5)}, {T{1}, T{3}, static_cast<T>(0.25)}, {static_cast<T>(0.5), static_cast<T>(0.25), T{2}}}));
        CHECK(factors(Matrix<3, 3, T>::zeros()));
        CHECK(factors(Matrix<3, 3, T>::diagonal({T{1}, T{0}, T{2}})));                               // a state known exactly, in the middle
        CHECK(factors(Matrix<3, 3, T>{{T{1}, T{2}, T{3}}, {T{2}, T{4}, T{6}}, {T{3}, T{6}, T{9}}})); // rank one
        CHECK(factors(Matrix<3, 3, T>{{T{1}, T{1}, T{0}}, {T{1}, T{1}, T{0}}, {T{0}, T{0}, T{5}}})); // two states perfectly correlated

        // A variance 10⁻²⁴ of the largest is a column like any other: the pivot's test is its own.
        Matrix<3, 3, T> units = Matrix<3, 3, T>::diagonal({static_cast<T>(1e12), T{1}, static_cast<T>(1e-12)});
        units(2, 1) = static_cast<T>(5e-7);
        units(1, 2) = static_cast<T>(5e-7);
        Matrix<3, 3, T> S{};
        REQUIRE(detail::semidefinite_factor(units, S));
        CHECK(static_cast<double>(S(2, 2)) == doctest::Approx(std::sqrt(0.75e-12)).epsilon(8.0 * e));
        CHECK(static_cast<double>(S(0, 0)) == doctest::Approx(1e6).epsilon(8.0 * e));

        // A negative pivot of roundoff size is a zero: P = v vᵀ rounded, less a few eps on one entry.
        Matrix<3, 3, T> rounded{{T{1}, T{2}, T{3}}, {T{2}, T{4}, T{6}}, {T{3}, T{6}, T{9}}};
        rounded(1, 1) = T{4} * (T{1} - (T{4} * eps<T>()));
        CHECK(detail::semidefinite_factor(rounded, S));
        CHECK(S(1, 1) == T{0});

        const Matrix<3, 3, T> before = S;
        const auto            refused = [&](const Matrix<3, 3, T>& P) { return !detail::semidefinite_factor(P, S) && same_bits(S, before); };
        CHECK(refused(Matrix<3, 3, T>{{T{1}, T{2}, T{0}}, {T{2}, T{1}, T{0}}, {T{0}, T{0}, T{1}}})); // a negative pivot
        CHECK(refused(Matrix<3, 3, T>{{T{0}, T{1}, T{0}}, {T{1}, T{0}, T{0}}, {T{0}, T{0}, T{1}}})); // a zero variance with a covariance
        CHECK(refused(Matrix<3, 3, T>::diagonal({T{1}, static_cast<T>(-1e-3), T{1}})));              // a negative variance
        for (const T bad : {runtime_nan<T>(), runtime_inf<T>()}) {
            Matrix<3, 3, T> P = Matrix<3, 3, T>::identity();
            P(1, 1) = bad;
            CHECK(refused(P));
            P = Matrix<3, 3, T>::identity();
            P(2, 0) = bad;
            CHECK(refused(P));
            P = Matrix<3, 3, T>::diagonal({T{0}, T{1}, T{1}}); // below a zero pivot
            P(2, 0) = bad;
            CHECK(refused(P));
        }
    }

    TEST_CASE_TEMPLATE("triangularize: S Sᵀ = W Wᵀ + E Eᵀ, S lower triangular with a non-negative diagonal", T, float, double) {
        const auto e = static_cast<double>(eps<T>());
        Normal     noise{3ULL};
        const auto random = [&]<size_t R, size_t C>(Matrix<R, C, T>& m, double scale) {
            for (size_t i = 0; i < R; ++i) {
                for (size_t j = 0; j < C; ++j) {
                    m(i, j) = static_cast<T>(scale * noise());
                }
            }
        };
        const auto check = [&]<bool Trapezoid, size_t N, size_t M>(Matrix<N, N, T> W, Matrix<N, M, T> E) {
            const Matrix<N, N> Wd = W.template as<double>();
            const Matrix<N, M> Ed = E.template as<double>();
            const Matrix<N, N> expect = Matrix<N, N>(Wd * Wd.t()) + Matrix<N, N>(Ed * Ed.t());
            detail::triangularize<Trapezoid>(W, E);
            CHECK(is_factor(W));
            // Rows of very different sizes: each entry of P to a few eps of its own rows' sizes.
            const Matrix<N, N> got = detail::covariance_of_factor(W.template as<double>());
            double             worst = 0.0;
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    worst = std::fmax(worst, std::fabs(got(i, j) - expect(i, j)) / (e * std::sqrt(expect(i, i) * expect(j, j)) + 1e-300));
                }
            }
            CAPTURE(worst);
            CHECK(worst <= 16.0 * static_cast<double>(N + M));
        };
        Matrix<4, 4, T> W{};
        Matrix<4, 3, T> E{};
        random(W, 1.0);
        random(E, 0.5);
        check.template operator()<false>(W, E);
        // Rows from 10⁶ to 10⁻⁶: states in different units.
        for (size_t j = 0; j < 4; ++j) {
            W(0, j) *= static_cast<T>(1e6);
            W(3, j) *= static_cast<T>(1e-6);
        }
        for (size_t j = 0; j < 3; ++j) {
            E(0, j) *= static_cast<T>(1e6);
            E(3, j) *= static_cast<T>(1e-6);
        }
        check.template operator()<false>(W, E);
        // A zero row (a register), and a row that is already its diagonal entry, negative.
        Matrix<4, 4, T> Wz{};
        Matrix<4, 3, T> Ez{};
        random(Wz, 1.0);
        random(Ez, 1.0);
        for (size_t j = 0; j < 4; ++j) {
            Wz(1, j) = T{0};
            Wz(0, j) = (j == 0) ? T{-2} : T{0};
        }
        for (size_t j = 0; j < 3; ++j) {
            Ez(1, j) = T{0};
            Ez(0, j) = T{0};
        }
        check.template  operator()<false>(Wz, Ez);
        Matrix<4, 4, T> Wout = Wz;
        Matrix<4, 3, T> Eout = Ez;
        detail::triangularize<false>(Wout, Eout);
        CHECK(Wout(0, 0) == T{2});
        for (size_t j = 0; j < 4; ++j) {
            CHECK(Wout(1, j) == T{0}); // exactly
        }
        // A lower-trapezoidal E: what is above its diagonal is not read.
        Matrix<4, 4, T> Et{};
        random(Et, 1.0);
        Matrix<4, 4, T> Et_clean = Et;
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = i + 1; j < 4; ++j) {
                Et_clean(i, j) = T{0};
                Et(i, j) = runtime_nan<T>();
            }
        }
        Matrix<4, 4, T> Wt = W;
        Matrix<4, 4, T> Wt_clean = W;
        detail::triangularize<true>(Wt, Et);
        detail::triangularize<false>(Wt_clean, Et_clean);
        CHECK(all_finite(Wt));
        CHECK(gap(Wt, Wt_clean) <= 64.0 * e * largest(Wt_clean));
    }

    TEST_CASE_TEMPLATE("process_noise_factor: Γ Γᵀ = G Q Gᵀ for fewer, as many and more noise inputs than states", T, float, double) {
        const auto e = static_cast<double>(eps<T>());
        const auto check = [&]<size_t NX, size_t NW>(const Matrix<NX, NW, T>& G, const Matrix<NW, NW, T>& Q) {
            auto                 Gamma = detail::process_noise_factor(G, Q);
            const Matrix<NX, NW> Gd = G.template as<double>();
            const Matrix<NX, NX> expect = Gd * Q.template as<double>() * Gd.t();
            const auto           Gm = Gamma.template as<double>();
            CHECK(gap(Matrix<NX, NX>(Gm * Gm.t()), expect) <= 64.0 * e * largest(expect));
            return Gamma;
        };
        const Matrix<3, 3, T> Q3{{T{4}, T{1}, T{0}}, {T{1}, T{2}, T{0}}, {T{0}, T{0}, T{0}}}; // one input with no variance
        const Matrix<2, 2, T> Q2{{T{4}, T{1}}, {T{1}, T{2}}};
        const Matrix<3, 2, T> G_few{{T{1}, T{0}}, {static_cast<T>(0.5), T{1}}, {T{0}, T{0}}}; // no row for the third state
        const auto            few = check(G_few, Q2);
        static_assert(std::is_same_v<decltype(few), const Matrix<3, 2, T>>);
        CHECK(few(2, 0) == T{0});
        CHECK(few(2, 1) == T{0});
        const auto square = check(Matrix<3, 3, T>{{T{1}, T{0}, T{2}}, {T{0}, T{1}, T{0}}, {T{1}, T{1}, T{1}}}, Q3);
        CHECK(is_factor(square));
        const auto more = check(Matrix<2, 3, T>{{T{1}, T{0}, T{2}}, {T{0}, T{1}, T{1}}}, Q3);
        static_assert(std::is_same_v<decltype(more), const Matrix<2, 2, T>>);
        CHECK(is_factor(more));
        // A Q that is not a covariance: no noise.
        CHECK(largest(detail::process_noise_factor(G_few, Matrix<2, 2, T>{{T{1}, T{2}}, {T{2}, T{1}}})) == 0.0);
        Matrix<3, 2, T> G_bad = G_few;
        G_bad(0, 1) = runtime_inf<T>();
        CHECK(largest(detail::process_noise_factor(G_bad, Q2)) == 0.0);
    }
}
