// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <limits>

#include "damp/backend.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/design/model_inverse.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

/**
 * Two coupled channels, minimum phase: output 0 has relative degree 1 and output 1 relative
 * degree 2; the zeros are 0.5 and −0.3, the poles 0.2 … 0.9. The inputs are mixed by a
 * nonsingular M, and the second channel reads the first one's state.
 */
constexpr StateSpace<5, 2, 2> mimo_minimum_phase() {
    const Matrix<5, 5> A{
        {0.0, 1.0, 0.0, 0.0, 0.0},
        {-0.48, 1.4, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1.0},
        {0.1, 0.0, 0.126, -0.95, 1.8},
    };
    const Matrix<5, 2> B{{0.0, 0.0}, {1.0, 0.5}, {0.0, 0.0}, {0.0, 0.0}, {-0.3, 1.0}};
    const Matrix<2, 5> C{{-0.5, 1.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.3, 1.0, 0.0}};
    return StateSpace<5, 2, 2>{.A = A, .B = B, .C = C, .Ts = 0.1};
}

/** (z − 2)/((z − 0.5)(z − 0.6)(z + 0.2)): relative degree 2, one zero outside the unit circle. */
constexpr StateSpace<3, 1, 1> siso_nonminimum_phase() {
    return StateSpace<3, 1, 1>{
        .A = Matrix<3, 3>{{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {-0.06, -0.08, 0.9}},
        .B = Matrix<3, 1>{{0.0}, {0.0}, {1.0}},
        .C = Matrix<1, 3>{{-2.0, 1.0, 0.0}},
        .Ts = 0.1,
    };
}

/**
 * (z² − 2z + 2.5)/((z − 0.5)(z − 0.6)(z + 0.2)(z − 0.3)): relative degree 2, a complex pair of
 * zeros 1 ± 1.22i outside the unit circle, which the Schur form holds as one 2×2 block.
 */
constexpr StateSpace<4, 1, 1> siso_complex_zeros() {
    return StateSpace<4, 1, 1>{
        .A = Matrix<4, 4>{{0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}, {0.018, -0.036, -0.35, 1.2}},
        .B = Matrix<4, 1>{{0.0}, {0.0}, {0.0}, {1.0}},
        .C = Matrix<1, 4>{{2.5, -2.0, 1.0, 0.0}},
        .Ts = 0.1,
    };
}

/** 1/(s(s + 1)(s + 2)) under ZOH at 0.1 s: relative degree 1, sampling zeros at −3.47 and −0.25. */
StateSpace<3, 1, 1> sampled_third_order() {
    const StateSpace<3, 1, 1> plant{
        .A = Matrix<3, 3>{{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, -2.0, -3.0}},
        .B = Matrix<3, 1>{{0.0}, {0.0}, {1.0}},
        .C = Matrix<1, 3>{{1.0, 0.0, 0.0}},
    };
    return *discretize(plant, 0.1);
}

template<size_t N, size_t NX, size_t NU>
struct Trajectory {
    damp::array<ColVec<NX>, N> x{};
    damp::array<ColVec<NU>, N> u{};
    damp::array<ColVec<NU>, N> y{};
};

/** The model's response from rest to an input that runs between samples @p on and @p off. */
template<size_t N, size_t NX, size_t NU, size_t NY>
Trajectory<N, NX, NY> simulate(const StateSpace<NX, NU, NY>& sys, size_t on, size_t off) {
    Trajectory<N, NX, NY> t{};
    ColVec<NX>            x{};
    for (size_t k = 0; k < N; ++k) {
        ColVec<NU> u{};
        if (k > on && k < off) {
            for (size_t i = 0; i < NU; ++i) {
                u(i) = damp::sin((0.3 * static_cast<double>(k)) + static_cast<double>(i)) - (0.5 * static_cast<double>(i));
            }
        }
        t.x[k] = x;
        t.u[k] = u;
        t.y[k] = ColVec<NY>(sys.C * x);
        x = ColVec<NX>((sys.A * x) + (sys.B * u));
    }
    return t;
}

/** y*ₖ … y*ₖ₊NP₋₁, held at the last sample past the end. */
template<size_t NP, size_t N, size_t NY>
damp::array<ColVec<NY>, NP> preview_at(const damp::array<ColVec<NY>, N>& y, size_t k) {
    damp::array<ColVec<NY>, NP> p{};
    for (size_t j = 0; j < NP; ++j) {
        p[j] = y[damp::min(k + j, N - 1)];
    }
    return p;
}

template<size_t N>
double max_abs(const ColVec<N>& v) {
    double m = 0.0;
    for (size_t i = 0; i < N; ++i) {
        m = damp::max(m, damp::abs(v(i)));
    }
    return m;
}

template<size_t R, size_t C>
double max_abs(const Matrix<R, C>& M) {
    double m = 0.0;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            m = damp::max(m, damp::abs(M(i, j)));
        }
    }
    return m;
}

/** Run the inverse over @p t and return the largest x* and u* errors against it, up to sample @p last. */
template<size_t NP, size_t N, size_t NX, size_t NU>
damp::array<double, 2> inverse_error(OutputInverse<NX, NU, NP, double>& inv, const Trajectory<N, NX, NU>& t, size_t last) {
    damp::array<double, 2> worst{};
    for (size_t k = 0; k <= last; ++k) {
        const auto p = inv.step(preview_at<NP>(t.y, k));
        worst[0] = damp::max(worst[0], max_abs(ColVec<NX>(p.x - t.x[k])));
        worst[1] = damp::max(worst[1], max_abs(ColVec<NU>(p.u - t.u[k])));
    }
    return worst;
}

template<size_t N, size_t NX, size_t NU>
damp::array<double, 2> trajectory_size(const Trajectory<N, NX, NU>& t) {
    damp::array<double, 2> size{};
    for (size_t k = 0; k < N; ++k) {
        size[0] = damp::max(size[0], max_abs(t.x[k]));
        size[1] = damp::max(size[1], max_abs(t.u[k]));
    }
    return size;
}

constexpr double kEps = damp::eps<double>();

/**
 * Three columns of an orthogonal N × N reflector Q: b₀ = Q·e₀, b₂ = Q·e₁ and b₁ = Q·(e₀ + τ·e₁) =
 * b₀ + τ·b₂ up to rounding, so σ₃/σ₁ is roundoff. In the order [b₀, b₁, b₂] every diagonal of an
 * unpivoted QR clears its column's floor; in the order [b₀, b₂, b₁] the last one does not.
 */
template<size_t N>
Matrix<N, 3> hidden_dependency(double tau, bool b1_last) {
    const damp::array<double, 4> w{0.3, -0.5, 0.7, 0.4};
    double                       ww = 0.0;
    for (size_t i = 0; i < N; ++i) {
        ww += w[i] * w[i];
    }
    Matrix<N, 3> B{};
    for (size_t i = 0; i < N; ++i) {
        const double q0 = ((i == 0) ? 1.0 : 0.0) - (2.0 * w[i] * w[0] / ww); // Q(i, 0)
        const double q1 = ((i == 1) ? 1.0 : 0.0) - (2.0 * w[i] * w[1] / ww); // Q(i, 1)
        B(i, 0) = q0;
        B(i, b1_last ? 2 : 1) = q0 + (tau * q1);
        B(i, b1_last ? 1 : 2) = q1;
    }
    return B;
}

constexpr auto kMimoInverse = design::output_inverse<3>(mimo_minimum_phase());
static_assert(kMimoInverse.success);
static_assert(kMimoInverse.relative_degree[0] == 1 && kMimoInverse.relative_degree[1] == 2);
static_assert(kMimoInverse.unstable_zeros == 0);

constexpr auto kMimoFeedforward = design::state_feedforward(mimo_minimum_phase());
static_assert(kMimoFeedforward.success);

} // namespace

TEST_SUITE("Model inverse: plant-inversion feedforward and the stable output inverse") {
    TEST_CASE("state_feedforward recovers the input of a model trajectory") {
        const auto sys = mimo_minimum_phase();
        const auto t = simulate<60>(sys, 3, 40);
        const auto res = design::state_feedforward(sys);
        REQUIRE(res.success);
        const StateFeedforward ff{res};
        double                 worst = 0.0;
        double                 bound = 0.0;
        for (size_t k = 0; k + 1 < 60; ++k) {
            const ColVec<2> u = ff.control(t.x[k], t.x[k + 1]);
            worst = damp::max(worst, max_abs(u - t.u[k]));
            // Rounding of the simulated step and of the two products, in absolute values.
            const double step = (max_abs(sys.A) * max_abs(t.x[k]) * 5.0) + (max_abs(sys.B) * max_abs(t.u[k]) * 2.0);
            bound = damp::max(bound, (max_abs(res.B_pinv) * 5.0 * (step + max_abs(t.x[k]) + max_abs(t.x[k + 1]))) + (max_abs(res.drift) * 5.0 * max_abs(t.x[k])));
        }
        CHECK(worst <= 64.0 * kEps * bound);
        CHECK(worst > 0.0); // the products round: the check is not vacuous
    }

    TEST_CASE("state_feedforward is the least-squares input off the model's range") {
        // x*ₖ₊₁ − A·x*ₖ = B·u + w with w orthogonal to range(B): u is recovered, w ignored.
        const auto             sys = mimo_minimum_phase();
        const StateFeedforward ff{design::state_feedforward(sys)};
        const ColVec<5>        x{0.2, -0.1, 0.4, 0.3, -0.2};
        const ColVec<2>        u{0.7, -1.1};
        const ColVec<5>        w{1.0, 0.0, 0.5, -0.25, 0.0}; // Bᵀw = 0: B's nonzero rows are 1 and 4
        const ColVec<5>        x_next = (sys.A * x) + (sys.B * u) + w;
        const ColVec<2>        got = ff.control(x, x_next);
        CHECK(max_abs(got - u) <= 64.0 * kEps * 4.0);
    }

    TEST_CASE("state_feedforward refuses a rank-deficient B and a continuous model") {
        const StateSpace<3, 2, 1> alike{
            .A = Matrix<3, 3>::identity() * 0.9,
            .B = Matrix<3, 2>{{1.0, 2.0}, {2.0, 4.0}, {0.0, 0.0}},
            .C = Matrix<1, 3>{{1.0, 0.0, 0.0}},
            .Ts = 0.1,
        };
        CHECK(design::state_feedforward(alike).status == design::DesignStatus::singular);
        auto continuous = mimo_minimum_phase();
        continuous.Ts = 0.0;
        CHECK(design::state_feedforward(continuous).status == design::DesignStatus::invalid_sample_time);
    }

    TEST_CASE_TEMPLATE("state_feedforward and output_inverse refuse a B whose dependency hides in the column order", T, float, double) {
        // B = [b₀, b₀ + τ·b₂, b₂] has rank 2. The unpivoted QR kept every |R(j, j)| above its
        // column's floor, and state_feedforward returned success with |B⁺| up to 2·10¹⁶; the
        // same columns as [b₀, b₂, b₁] were refused. Both orders are refused by mat::rank now.
        for (const double tau : {1e-4, 1e-6}) {
            for (const bool b1_last : {false, true}) {
                INFO("tau ", tau, ", b1 last ", b1_last);
                const Matrix<4, 3, T>        B = hidden_dependency<4>(tau, b1_last).template as<T>();
                const StateSpace<4, 3, 1, T> sys{
                    .A = Matrix<4, 4, T>::identity() * T(0.9),
                    .B = B,
                    .C = Matrix<1, 4, T>{{T(1), T(0), T(0), T(0)}},
                    .Ts = T(0.1),
                };
                const auto ff = design::state_feedforward(sys);
                CHECK_FALSE(ff.success);
                CHECK(ff.status == design::DesignStatus::singular);

                // The decoupling matrix of y = x, three outputs of relative degree 1, is B itself.
                const StateSpace<3, 3, 3, T> square{
                    .A = Matrix<3, 3, T>::identity() * T(0.5),
                    .B = hidden_dependency<3>(tau, b1_last).template as<T>(),
                    .C = Matrix<3, 3, T>::identity(),
                    .Ts = T(0.1),
                };
                CHECK(design::output_inverse<3>(square).status == design::DesignStatus::singular);
            }
        }
        // Columns of different sizes are not a dependency.
        const StateSpace<4, 3, 1, T> scaled{
            .A = Matrix<4, 4, T>::identity() * T(0.9),
            .B = Matrix<4, 3, T>{{T(1), T(0), T(0)}, {T(0), T(1e-3), T(0)}, {T(0), T(0), T(1e3)}, {T(0), T(0), T(0)}},
            .C = Matrix<1, 4, T>{{T(1), T(0), T(0), T(0)}},
            .Ts = T(0.1),
        };
        CHECK(design::state_feedforward(scaled).success);
    }

    TEST_CASE("output_inverse recovers x and u of a minimum-phase MIMO model, vector relative degree (1, 2)") {
        const auto sys = mimo_minimum_phase();
        const auto t = simulate<80>(sys, 5, 45);
        const auto res = design::output_inverse<3>(sys, damp::array<size_t, 2>{0, 1});
        REQUIRE(res.success);
        CHECK(res.relative_degree[0] == 1);
        CHECK(res.relative_degree[1] == 2);
        CHECK(res.unstable_zeros == 0);
        CHECK(res.truncation == 0.0);
        OutputInverse inv{res};
        const auto    err = inverse_error(inv, t, 80 - 3);
        const auto    size = trajectory_size(t);
        CHECK(err[0] <= 1024.0 * kEps * size[0]);
        CHECK(err[1] <= 1024.0 * kEps * (size[1] + (max_abs(res.u_zero) * size[0])));
    }

    TEST_CASE("output_inverse solves a zero outside the unit circle over the preview") {
        const auto sys = siso_nonminimum_phase();
        const auto t = simulate<120>(sys, 40, 70);
        const auto res = design::output_inverse<44>(sys);
        REQUIRE(res.success);
        CHECK(res.relative_degree[0] == 2);
        CHECK(res.unstable_zeros == 1);
        CHECK(res.truncation == doctest::Approx(damp::pow(0.5, 42.0)).epsilon(1e-9)); // 1/|z|^(NP − r), z = 2
        OutputInverse inv{res};
        const auto    err = inverse_error(inv, t, 120 - 44);
        const auto    size = trajectory_size(t);
        // The preview holds the reference past its end: the bound is the truncation, not ε.
        CHECK(err[0] <= ((64.0 * res.truncation) + (1024.0 * kEps)) * size[0]);
        CHECK(err[1] <= ((64.0 * res.truncation) + (1024.0 * kEps)) * (size[1] + (max_abs(res.u_zero) * size[0])));
    }

    TEST_CASE("output_inverse solves a complex pair of zeros outside the unit circle") {
        const auto sys = siso_complex_zeros();
        const auto t = simulate<160>(sys, 60, 90);
        const auto res = design::output_inverse<60>(sys);
        REQUIRE(res.success);
        CHECK(res.relative_degree[0] == 2);
        CHECK(res.unstable_zeros == 2);
        CHECK(res.truncation == doctest::Approx(damp::pow(1.0 / damp::sqrt(2.5), 58.0)).epsilon(1e-9)); // |z| = √2.5
        OutputInverse inv{res};
        const auto    err = inverse_error(inv, t, 160 - 60);
        const auto    size = trajectory_size(t);
        CHECK(err[0] <= ((64.0 * res.truncation) + (1024.0 * kEps)) * size[0]);
        CHECK(err[1] <= ((64.0 * res.truncation) + (1024.0 * kEps)) * (size[1] + (max_abs(res.u_zero) * size[0])));
    }

    TEST_CASE("output_inverse inverts a sampled third-order plant through its sampling zero at −3.47") {
        const auto sys = sampled_third_order();
        const auto t = simulate<200>(sys, 20, 60);
        const auto res = design::output_inverse<30>(sys);
        REQUIRE(res.success);
        CHECK(res.relative_degree[0] == 1);
        CHECK(res.unstable_zeros == 1);
        OutputInverse inv{res};
        const auto    err = inverse_error(inv, t, 200 - 30);
        const auto    size = trajectory_size(t);
        CHECK(err[0] <= ((64.0 * res.truncation) + (1024.0 * kEps)) * size[0]);
        CHECK(err[1] <= ((64.0 * res.truncation) + (1024.0 * kEps)) * (size[1] + (max_abs(res.u_zero) * size[0])));

        // The two feedforwards agree: u* is the plant-inversion input between consecutive x*.
        OutputInverse          again{res};
        const StateFeedforward ff{design::state_feedforward(sys)};
        auto                   p = again.step(preview_at<30>(t.y, 0));
        double                 worst = 0.0;
        for (size_t k = 0; k + 30 < 200; ++k) {
            const auto next = again.step(preview_at<30>(t.y, k + 1));
            worst = damp::max(worst, damp::abs(ff.control(p.x, next.x)(0) - p.u(0)));
            p = next;
        }
        CHECK(worst <= ((64.0 * res.truncation) + (1024.0 * kEps)) * (size[1] + (max_abs(res.u_zero) * size[0])));
    }

    TEST_CASE("output_inverse holds the steady state of design::steady_state_target at rest") {
        const auto sys = sampled_third_order();
        const auto res = design::output_inverse<30>(sys);
        REQUIRE(res.success);
        const auto target = design::steady_state_target(sys);
        REQUIRE(target.success);
        // Per unit output: x_hold is O(1); u_hold carries x_hold's rounding through D*⁻¹ (u_zero).
        const double tol_x = 1024.0 * kEps;
        const double tol_u = 1024.0 * kEps * (1.0 + max_abs(res.u_zero));
        for (size_t i = 0; i < 3; ++i) {
            CHECK(damp::abs(res.x_hold(i, 0) - target.target_r(i, 0)) <= tol_x);
        }
        CHECK(damp::abs(res.u_hold(0, 0) - target.target_r(3, 0)) <= tol_u);

        // A held reference in float: x* = x_hold·y*, the zero dynamics stay at rest.
        OutputInverse<3, 1, 30, float>    inv{res};
        damp::array<ColVec<1, float>, 30> held{};
        for (auto& v : held) {
            v = ColVec<1, float>{0.25f};
        }
        for (int k = 0; k < 5; ++k) {
            const auto p = inv.step(held);
            CHECK(p.x(0) == doctest::Approx(0.25f));
            CHECK(damp::abs(p.u(0)) <= 1e-6f);
        }
        CHECK(max_abs(Matrix<3, 1>(inv.zero_dynamics_state.template as<double>())) == 0.0);
    }

    TEST_CASE("Law<kRes, float> runs an output inverse") {
        Law<kMimoInverse, float>         law{};
        damp::array<ColVec<2, float>, 3> held{};
        for (auto& v : held) {
            v = ColVec<2, float>{0.5f, -0.25f};
        }
        const auto p = law.step(held);
        const auto y = mimo_minimum_phase().C.template as<float>() * p.x;
        CHECK(y(0, 0) == doctest::Approx(0.5f).epsilon(1e-5));
        CHECK(y(1, 0) == doctest::Approx(-0.25f).epsilon(1e-5));
    }

    TEST_CASE("output_inverse refuses what it cannot invert") {
        // The zero at −1 of a sampled double integrator: no bounded inverse.
        const StateSpace<2, 1, 1> double_integrator{
            .A = Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}},
            .B = Matrix<2, 1>{{0.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
        };
        CHECK(design::output_inverse<8>(*discretize(double_integrator, 0.01)).status == design::DesignStatus::singular);

        // Two outputs that read the same combination: the decoupling matrix is singular.
        const StateSpace<2, 2, 2> twin{
            .A = Matrix<2, 2>{{0.5, 0.0}, {0.0, 0.6}},
            .B = Matrix<2, 2>::identity(),
            .C = Matrix<2, 2>{{1.0, 1.0}, {2.0, 2.0}},
            .Ts = 0.1,
        };
        CHECK(design::output_inverse<3>(twin).status == design::DesignStatus::singular);

        // An output no input reaches.
        const StateSpace<2, 1, 1> unreached{
            .A = Matrix<2, 2>{{0.5, 0.0}, {0.0, 0.6}},
            .B = Matrix<2, 1>{{1.0}, {0.0}},
            .C = Matrix<1, 2>{{0.0, 1.0}},
            .Ts = 0.1,
        };
        CHECK(design::output_inverse<3>(unreached).status == design::DesignStatus::singular);

        // A preview that ends before the input reaches the output (r = 2 needs NP ≥ 3).
        CHECK(design::output_inverse<2>(siso_nonminimum_phase()).status == design::DesignStatus::invalid_parameter);

        // A continuous model, and controlled outputs out of range or listed twice.
        auto continuous = mimo_minimum_phase();
        continuous.Ts = 0.0;
        CHECK(design::output_inverse<3>(continuous).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::output_inverse<3>(mimo_minimum_phase(), damp::array<size_t, 2>{0, 2}).status == design::DesignStatus::bad_weights);
        CHECK(design::output_inverse<3>(mimo_minimum_phase(), damp::array<size_t, 2>{1, 1}).status == design::DesignStatus::bad_weights);
    }
}

TEST_SUITE("model inverse: roundoff in A − B·D*⁻¹C* and non-finite models") {
    TEST_CASE("output_inverse in float holds the steady state, states rescaled or not") {
        // A − B·D*⁻¹C* cancels whole entries to roundoff; balancing that noise spread the exponents
        // over decades, so float was refused, and with rescaled states passed the certificate with
        // C·x_hold off the identity by 1.8 %.
        const auto   sys = mimo_minimum_phase();
        Matrix<5, 5> S{};
        Matrix<5, 5> S_inv{};
        const double d[5] = {3.0, 1.0 / 3.0, 7.0, 1.0 / 7.0, 5.0};
        for (size_t i = 0; i < 5; ++i) {
            S(i, i) = d[i];
            S_inv(i, i) = 1.0 / d[i];
        }
        const StateSpace<5, 2, 2> scaled{.A = S_inv * sys.A * S, .B = S_inv * sys.B, .C = sys.C * S, .Ts = 0.1};
        for (const auto& model : {sys.as<float>(), scaled.as<float>()}) {
            const auto inv = design::output_inverse<3>(model);
            REQUIRE(inv.success);
            const Matrix<2, 2, float> Cx = model.C * inv.x_hold;
            CHECK(Cx(0, 0) == doctest::Approx(1.0f).epsilon(1e-4));
            CHECK(Cx(1, 1) == doctest::Approx(1.0f).epsilon(1e-4));
            CHECK(damp::abs(Cx(0, 1)) < 1e-4f);
            CHECK(damp::abs(Cx(1, 0)) < 1e-4f);
        }
    }

    TEST_CASE("state_feedforward refuses a non-finite model with bad_weights") {
        auto sys = mimo_minimum_phase();
        sys.A(0, 0) = std::numeric_limits<double>::quiet_NaN();
        const auto ff = design::state_feedforward(sys);
        CHECK_FALSE(ff.success);
        CHECK(ff.status == design::DesignStatus::bad_weights);
    }
}
