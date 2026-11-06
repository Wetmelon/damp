// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/design/pole_placement.hpp"
#include "damp/design/stability.hpp"
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

using Cplx = damp::complex<double>;

// A continuous plant with the given A and B (outputs = states), for the placement tests.
template<size_t NX, size_t NU>
constexpr StateSpace<NX, NU, NX> plant(const Matrix<NX, NX>& A, const Matrix<NX, NU>& B, double Ts = 0.0) {
    return StateSpace<NX, NU, NX>{.A = A, .B = B, .C = Matrix<NX, NX>::identity(), .Ts = Ts};
}

// Real poles as the complex pole list place() takes.
template<size_t NX>
constexpr std::array<Cplx, NX> poles_of(const std::array<double, NX>& p) {
    std::array<Cplx, NX> out{};
    for (size_t i = 0; i < NX; ++i) {
        out[i] = Cplx{p[i], 0.0};
    }
    return out;
}

// Largest distance between the (sorted, real) spectrum of A − B·K and the
// requested real poles.
template<size_t NX, size_t NU>
double placement_error(const Matrix<NX, NX>& A, const Matrix<NX, NU>& B, const std::array<double, NX>& poles) {
    const auto res = design::place(plant(A, B), poles_of(poles));
    REQUIRE(res.success);
    Matrix<NX, NX> ABK = A - (B * res.K);
    auto           r = mat::compute_eigenvalues(ABK);

    std::array<double, NX> got{};
    for (size_t i = 0; i < NX; ++i) {
        got[i] = r.values[i].real();
    }
    std::array<double, NX> want = poles;
    std::sort(got.begin(), got.end());
    std::sort(want.begin(), want.end());
    double worst = 0.0;
    for (size_t i = 0; i < NX; ++i) {
        worst = damp::max(worst, std::abs(got[i] - want[i]));
    }
    return worst;
}

} // namespace

TEST_SUITE("pole_placement") {
    TEST_CASE("single-input placement assigns the spectrum exactly") {
        SUBCASE("double integrator, poles {-1, -2}") {
            Matrix<2, 2> A = {{0.0, 1.0}, {0.0, 0.0}};
            Matrix<2, 1> B = {{0.0}, {1.0}};
            CHECK(placement_error(A, B, {-1.0, -2.0}) < 1e-9);
        }
        SUBCASE("3x3 companion, poles {-2, -3, -4}") {
            Matrix<3, 3> A = {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {-1.0, -2.0, -3.0}};
            Matrix<3, 1> B = {{0.0}, {0.0}, {1.0}};
            CHECK(placement_error(A, B, {-2.0, -3.0, -4.0}) < 1e-9);
        }
    }

    TEST_CASE("multi-input placement assigns the spectrum exactly") {
        SUBCASE("3x2") {
            Matrix<3, 3> A = {{1.0, 2.0, 0.0}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}};
            Matrix<3, 2> B = {{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}};
            CHECK(placement_error(A, B, {-1.0, -2.0, -3.0}) < 1e-8);
        }
        SUBCASE("4x2") {
            Matrix<4, 4> A = {{0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}, {1.0, 2.0, 3.0, 4.0}};
            Matrix<4, 2> B = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 1.0}};
            CHECK(placement_error(A, B, {-1.0, -2.0, -3.0, -4.0}) < 1e-8);
        }
        SUBCASE("repeated pole with multiplicity == NU is assignable") {
            Matrix<4, 4> A = {{0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}, {1.0, 2.0, 3.0, 4.0}};
            Matrix<4, 2> B = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 1.0}};
            CHECK(placement_error(A, B, {-2.0, -2.0, -3.0, -3.0}) < 1e-7);
        }
    }

    TEST_CASE("square B (NU == NX) places via mat::solve on B") {
        Matrix<2, 2> A = {{1.0, 2.0}, {3.0, 4.0}};
        Matrix<2, 2> B = Matrix<2, 2>::identity();
        CHECK(placement_error(A, B, {-1.0, -5.0}) < 1e-10);
    }

    TEST_CASE("single-input place matches Ackermann (unique SI gain)") {
        Matrix<3, 3> A = {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {-1.0, -2.0, -3.0}};
        Matrix<3, 1> B = {{0.0}, {0.0}, {1.0}};
        const auto   Kp = design::place(plant(A, B), {-2.0, -3.0, -4.0});
        REQUIRE(Kp.success);

        // Ackermann by hand: s³ + 9s² + 26s + 24 on the companion form with a = (1, 2, 3).
        const std::array<double, 3> Kack{24.0 - 1.0, 26.0 - 2.0, 9.0 - 3.0};
        for (size_t j = 0; j < 3; ++j) {
            CHECK(Kp.K(0, j) == doctest::Approx(Kack[j]).epsilon(1e-9));
        }
    }

    TEST_CASE("place rejects non-assignable problems") {
        Matrix<3, 3> A = {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {-1.0, -2.0, -3.0}};

        SUBCASE("multiplicity exceeding the input count (two inputs, triple pole)") {
            Matrix<3, 2> B = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
            const auto   res = design::place(plant(A, B), {-2.0, -2.0, -2.0});
            CHECK_FALSE(res.success);
            CHECK(res.status == design::DesignStatus::not_placeable);
        }
        SUBCASE("rank-deficient B") {
            Matrix<3, 2> B = {{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}}; // rank 1
            const auto   res = design::place(plant(A, B), {-1.0, -2.0, -3.0});
            CHECK_FALSE(res.success);
            CHECK(res.status == design::DesignStatus::singular);
        }
    }

    // Largest distance between the spectrum of A − B·K and a complex pole set,
    // comparing as sorted (real, |imag|) pairs.
    static constexpr auto placement_error_cplx = []<size_t NX, size_t NU>(
                                                     const Matrix<NX, NX>&                        A,
                                                     const Matrix<NX, NU>&                        B,
                                                     const std::array<damp::complex<double>, NX>& poles
                                                 ) {
        const auto res = design::place(plant(A, B), poles);
        REQUIRE(res.success);
        Matrix<NX, NX> ABK = A - (B * res.K);
        auto           r = mat::compute_eigenvalues(ABK);

        std::array<damp::pair<double, double>, NX> got;
        std::array<damp::pair<double, double>, NX> want;
        for (size_t i = 0; i < NX; ++i) {
            got[i] = {r.values[i].real(), std::abs(r.values[i].imag())};
            want[i] = {poles[i].real(), std::abs(poles[i].imag())};
        }
        std::sort(got.begin(), got.end());
        std::sort(want.begin(), want.end());
        double worst = 0.0;
        for (size_t i = 0; i < NX; ++i) {
            worst = damp::max(worst, std::abs(got[i].first - want[i].first));
            worst = damp::max(worst, std::abs(got[i].second - want[i].second));
        }
        return worst;
    };

    TEST_CASE("complex-conjugate placement assigns the spectrum exactly") {
        SUBCASE("3x2: one complex pair + one real") {
            Matrix<3, 3> A = {{1.0, 2.0, 0.0}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}};
            Matrix<3, 2> B = {{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}};
            CHECK(placement_error_cplx(A, B, std::array<Cplx, 3>{Cplx(-1, 2), Cplx(-1, -2), Cplx(-3, 0)}) < 1e-6);
        }
        SUBCASE("4x2: two complex pairs") {
            Matrix<4, 4> A = {{0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}, {1.0, 2.0, 3.0, 4.0}};
            Matrix<4, 2> B = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 1.0}};
            CHECK(placement_error_cplx(A, B, std::array<Cplx, 4>{Cplx(-1, 1), Cplx(-1, -1), Cplx(-2, 3), Cplx(-2, -3)}) < 1e-6);
        }
        SUBCASE("4x2: mixed real + complex pair") {
            Matrix<4, 4> A = {{0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}, {1.0, 2.0, 3.0, 4.0}};
            Matrix<4, 2> B = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 1.0}};
            CHECK(placement_error_cplx(A, B, std::array<Cplx, 4>{Cplx(-1, 2), Cplx(-1, -2), Cplx(-3, 0), Cplx(-4, 0)}) < 1e-6);
        }
        SUBCASE("square 2x2 complex pair") {
            Matrix<2, 2> A = {{1.0, 2.0}, {3.0, 4.0}};
            Matrix<2, 2> B = Matrix<2, 2>::identity();
            CHECK(placement_error_cplx(A, B, std::array<Cplx, 2>{Cplx(-1, 2), Cplx(-1, -2)}) < 1e-9);
        }
        SUBCASE("all-real array forwards to the real (conditioned) path") {
            // Exercise the complex overload with an all-real spectrum.
            Matrix<3, 3> A3 = {{1.0, 2.0, 0.0}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}};
            Matrix<3, 2> B = {{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}};
            CHECK(placement_error_cplx(A3, B, std::array<Cplx, 3>{Cplx(-1, 0), Cplx(-2, 0), Cplx(-3, 0)}) < 1e-8);
        }
    }

    TEST_CASE("place rejects a dangling (unpaired) complex pole") {
        Matrix<3, 3> A = {{1.0, 2.0, 0.0}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}};
        Matrix<3, 2> B = {{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}};
        // -1+2j without its conjugate is not a real closed-loop spectrum.
        CHECK_FALSE(design::place(plant(A, B), {Cplx(-1, 2), Cplx(-3, 0), Cplx(-4, 0)}).success);
    }

    TEST_CASE("matlab::place forwards to design::place") {
        Matrix<3, 3>        A = {{1.0, 2.0, 0.0}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}};
        Matrix<3, 2>        B = {{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}};
        std::array<Cplx, 3> poles = {Cplx(-1, 2), Cplx(-1, -2), Cplx(-3, 0)};

        const design::PlaceResult<3, 2> Km = matlab::place(A, B, poles); // the design:: result type
        const auto                      Kd = design::place(plant(A, B), poles);
        REQUIRE(Km.success);
        REQUIRE(Kd.success);
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(Km.K(i, j) == doctest::Approx(Kd.K(i, j)));
            }
        }
    }

    TEST_CASE("place is constexpr-evaluable") {
        constexpr bool ok = []() consteval {
            Matrix<2, 2> A = {{0.0, 1.0}, {0.0, 0.0}};
            Matrix<2, 1> B = {{0.0}, {1.0}};
            const auto   res = design::place(plant(A, B), {-1.0, -2.0});
            if (!res.success) {
                return false;
            }
            // Double-integrator with poles {-1,-2}: char poly s²+3s+2 ⇒ K = [2, 3].
            return damp::abs(res.K(0, 0) - 2.0) < 1e-9 && damp::abs(res.K(0, 1) - 3.0) < 1e-9;
        }();
        static_assert(ok, "place must work at compile time");
        CHECK(ok);
    }

    TEST_CASE("place(sys, poles, Ts) matches ZOH + z-map + place on the discrete plant") {
        Matrix<2, 2>        A = {{0.0, 1.0}, {0.0, 0.0}};
        Matrix<2, 1>        B = {{0.0}, {1.0}};
        const double        Ts = 0.01;
        std::array<Cplx, 2> poles_s{-4.0, -6.0};
        const auto          sys_c = plant(A, B);

        const auto res = design::place(sys_c, poles_s, Ts);
        REQUIRE(res.success);
        CHECK(res.Ts == doctest::Approx(Ts));

        const auto          sys_d = *discretize(sys_c, Ts, DiscretizationMethod::ZOH);
        std::array<Cplx, 2> poles_z{damp::exp(-4.0 * Ts), damp::exp(-6.0 * Ts)};
        const auto          manual = design::place(sys_d, poles_z);
        REQUIRE(manual.success);
        for (size_t j = 0; j < 2; ++j) {
            CHECK(res.K(0, j) == doctest::Approx(manual.K(0, j)).epsilon(1e-9));
        }

        // e is the achieved discrete spectrum of Ad − Bd K, near e^{s Ts}.
        std::array<double, 2> got{res.e[0].real(), res.e[1].real()};
        std::sort(got.begin(), got.end());
        CHECK(std::abs(got[0] - poles_z[1].real()) < 1e-8);
        CHECK(std::abs(got[1] - poles_z[0].real()) < 1e-8);
        CHECK(res.is_stable());

        const auto bad_ts = design::place(sys_c, poles_s, 0.0);
        CHECK_FALSE(bad_ts.success);
        CHECK(bad_ts.status == design::DesignStatus::invalid_sample_time);
    }

    TEST_CASE("a discrete plant with a trailing Ts is refused, not discretized twice") {
        // The trap this API closes: the same (plant, poles, Ts) call on a model that is
        // already discrete must not ZOH it a second time.
        Matrix<2, 2> A = {{0.0, 1.0}, {0.0, 0.0}};
        Matrix<2, 1> B = {{0.0}, {1.0}};
        const double Ts = 0.01;
        const auto   sys_d = *discretize(plant(A, B), Ts, DiscretizationMethod::ZOH);
        REQUIRE(sys_d.Ts == doctest::Approx(Ts));

        const auto res = design::place(sys_d, {-5.0, -7.0}, Ts);
        CHECK_FALSE(res.success);
        CHECK(res.status == design::DesignStatus::invalid_sample_time);

        const auto jordan = design::place_jordan(sys_d, damp::array<design::JordanBlock<>, 1>{{{{-5.0, 0.0}, 2}}}, Ts);
        CHECK(jordan.status == design::DesignStatus::invalid_sample_time);

        // A discrete plant takes z-plane poles directly.
        const auto z = design::place(sys_d, {damp::exp(-5.0 * Ts), damp::exp(-7.0 * Ts)});
        CHECK(z.success);
        CHECK(z.Ts == doctest::Approx(Ts));
    }

    TEST_CASE("single-input repeated poles form one Jordan chain (the Ackermann gain)") {
        // Double integrator, double pole at -3: s² + 6s + 9 ⇒ K = [9, 6].
        Matrix<2, 2> A = {{0.0, 1.0}, {0.0, 0.0}};
        Matrix<2, 1> B = {{0.0}, {1.0}};
        const auto   res = design::place(plant(A, B), {-3.0, -3.0});
        REQUIRE(res.success);
        CHECK(res.K(0, 0) == doctest::Approx(9.0).epsilon(1e-9));
        CHECK(res.K(0, 1) == doctest::Approx(6.0).epsilon(1e-9));

        // Deadbeat on a discrete triple integrator: all three z-poles at 0.
        Matrix<3, 3> Ad = {{1.0, 1.0, 0.5}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}};
        Matrix<3, 1> Bd = {{1.0 / 6.0}, {0.5}, {1.0}};
        const auto   db = design::place(plant(Ad, Bd, 1.0), {0.0, 0.0, 0.0});
        REQUIRE(db.success);
        const Matrix<3, 3> Acl = Ad - (Bd * db.K);
        CHECK((Acl * Acl * Acl).norm() < 1e-9); // nilpotent
    }

    TEST_CASE("PlaceResult deploys through LQR") {
        Matrix<2, 2>   A = {{1.0, 0.01}, {0.0, 1.0}};
        Matrix<2, 1>   B = {{0.00005}, {0.01}};
        constexpr auto res = design::place(StateSpace<2, 1, 2>{.A = {{1.0, 0.01}, {0.0, 1.0}}, .B = {{0.00005}, {0.01}}, .C = Matrix<2, 2>::identity(), .Ts = 0.01}, {0.9, 0.95});
        static_assert(res.success);

        static constinit LQR law{res.as<float>()};
        static_assert(std::is_same_v<decltype(law), LQR<2, 1, float>>);
        LQR<2, 1, float> assigned = res; // assign-from-design
        CHECK(law.K(0, 0) == doctest::Approx(static_cast<float>(res.K(0, 0))));
        CHECK(assigned.K(0, 1) == doctest::Approx(static_cast<float>(res.K(0, 1))));

        const auto e = stability::closed_loop_poles(A, B, res.K);
        CHECK(std::min(e[0].real(), e[1].real()) == doctest::Approx(0.9).epsilon(1e-9));
    }

    TEST_CASE("a repeated pole within NU is placed when the spectrum also has a complex pair") {
        const Matrix<4, 4> A = {{0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}, {1.0, 2.0, 3.0, 4.0}};
        const Matrix<4, 2> B = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 1.0}};
        for (const auto& poles : {std::array<Cplx, 4>{Cplx{-2.0, 0.0}, Cplx{-2.0, 0.0}, Cplx{-1.0, 1.0}, Cplx{-1.0, -1.0}}, std::array<Cplx, 4>{Cplx{-1.0, 1.0}, Cplx{-1.0, -1.0}, Cplx{-1.0, 1.0}, Cplx{-1.0, -1.0}}}) {
            const auto res = design::place(plant(A, B), poles);
            REQUIRE(res.success);
            const auto ev = mat::compute_eigenvalues(A - (B * res.K));
            for (const auto& p : poles) {
                double best = 1e300;
                for (size_t i = 0; i < 4; ++i) {
                    best = std::min(best, abs(Cplx(ev.values[i]) - p));
                }
                CHECK(best < 1e-6);
            }
        }
    }

    TEST_CASE("distinct poles closer than an absolute 1e-7 stay distinct (fast sample, slow plant)") {
        // Two thermal nodes (τ ≈ 2000 s and 200 s) at Ts = 1e-4: the z-plane images of
        // s = −1e-3 and −2e-3 differ by 1e-7, and both were placed at the first.
        const StateSpace<2, 1, 1> th{.A = {{-1.0 / 2000.0, 1.0 / 2000.0}, {1.0 / 200.0, -2.0 / 200.0}}, .B = {{1.0 / 2000.0}, {0.0}}, .C = {{0.0, 1.0}}};
        const auto                res = design::place(th, {-1e-3, -2e-3}, 1e-4);
        REQUIRE(res.success);
        const double s0 = std::log(res.e[0].real()) / 1e-4;
        const double s1 = std::log(res.e[1].real()) / 1e-4;
        CHECK(std::min(s0, s1) == doctest::Approx(-2e-3).epsilon(1e-4));
        CHECK(std::max(s0, s1) == doctest::Approx(-1e-3).epsilon(1e-4));

        // Slow s-plane poles on a double integrator are distinct, not a repeated pole.
        const auto slow = design::place(plant(Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}}, Matrix<2, 1>{{0.0}, {1.0}}), {-1e-8, -2e-8});
        REQUIRE(slow.success);
        CHECK(slow.K(0, 0) == doctest::Approx(2e-16).epsilon(1e-6));
        CHECK(slow.K(0, 1) == doctest::Approx(3e-8).epsilon(1e-6));
    }

    TEST_CASE("single-input nearly repeated poles take the unique gain") {
        // A cluster just past the 64·ε equality is not a repeat, and KNV's eigenvectors for it are
        // nearly parallel; with one input the gain is unique and Ackermann's formula gives it.
        const Matrix<3, 3> A3 = {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, 0.0}};
        const Matrix<3, 1> B3 = {{0.0}, {0.0}, {1.0}};
        const auto         tri = design::place(plant(A3, B3), {-1.0, -1.0 - 1e-9, -1.0 - 2e-9});
        REQUIRE(tri.success);
        CHECK(tri.K(0, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(tri.K(0, 1) == doctest::Approx(3.0).epsilon(1e-6));
        CHECK(tri.K(0, 2) == doctest::Approx(3.0).epsilon(1e-6));

        // Four integrators, poles −1 − k·1e-12: the computed eig(A − BK) spread by ε^{1/4} like a
        // chain's, so the certificate counts the cluster's multiplicity, not the 64·ε repeats.
        const Matrix<4, 4> A4 = {{0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}, {0.0, 0.0, 0.0, 0.0}};
        const Matrix<4, 1> B4 = {{0.0}, {0.0}, {0.0}, {1.0}};
        const auto         quad = design::place(plant(A4, B4), {-1.0, -1.0 - 1e-12, -1.0 - 2e-12, -1.0 - 3e-12});
        REQUIRE(quad.success);
        CHECK(quad.K(0, 0) == doctest::Approx(1.0).epsilon(1e-9)); // (s + 1)⁴
        CHECK(quad.K(0, 1) == doctest::Approx(4.0).epsilon(1e-9));
        CHECK(quad.K(0, 2) == doctest::Approx(6.0).epsilon(1e-9));
        CHECK(quad.K(0, 3) == doctest::Approx(4.0).epsilon(1e-9));

        // A double root as an eigensolver returns it, −2 ± 1e-12j: K = [4, 4].
        const auto pair = design::place(plant(Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}}, Matrix<2, 1>{{0.0}, {1.0}}), {Cplx{-2.0, 1e-12}, Cplx{-2.0, -1e-12}});
        REQUIRE(pair.success);
        CHECK(pair.K(0, 0) == doctest::Approx(4.0).epsilon(1e-9));
        CHECK(pair.K(0, 1) == doctest::Approx(4.0).epsilon(1e-9));
    }

    TEST_CASE("float poles equal up to rounding are one repeated pole") {
        const StateSpace<2, 1, 1, float> di{.A = {{0.0f, 1.0f}, {0.0f, 0.0f}}, .B = {{0.0f}, {1.0f}}, .C = {{1.0f, 0.0f}}};
        const float                      p = -1005.30969f;
        const auto                       res = design::place(di, {damp::complex<float>{p, 0.0f}, damp::complex<float>{std::nextafter(p, 0.0f), 0.0f}});
        REQUIRE(res.success);
        CHECK(res.K(0, 0) == doctest::Approx(p * p).epsilon(1e-5));
        CHECK(res.K(0, 1) == doctest::Approx(-2.0f * p).epsilon(1e-5));
    }

    TEST_CASE("a non-finite pole is refused") {
        const Matrix<4, 4> A = {{0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}, {1.0, 2.0, 3.0, 4.0}};
        const Matrix<4, 2> B = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 1.0}};
        const auto         res = design::place(plant(A, B), {std::nan(""), -2.0, -3.0, -4.0});
        CHECK_FALSE(res.success);
        CHECK(res.status == design::DesignStatus::invalid_parameter);
        CHECK(res.K(0, 0) == 0.0);
    }

    TEST_CASE("a failed design is not stable") {
        const auto res = design::place(plant(Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}}, Matrix<2, 1>{{0.0}, {0.0}}, 0.01), {0.5, 0.6});
        REQUIRE_FALSE(res.success);
        CHECK_FALSE(res.is_stable()); // e is all zeros, inside the unit circle
    }

    TEST_CASE("a continuous PlaceResult is refused: zero gain at run time, no compile constant-evaluated") {
        constexpr auto continuous = design::place(StateSpace<2, 1, 2>{.A = {{0.0, 1.0}, {0.0, 0.0}}, .B = {{0.0}, {1.0}}, .C = Matrix<2, 2>::identity()}, {-2.0, -3.0});
        static_assert(continuous.success && continuous.Ts == 0.0);
        // static constinit LQR law{continuous}; // ill-formed: refuse_failed_design under constant evaluation
        const LQR law{continuous}; // a continuous gain on a discrete tick: loaded as a failed design
        CHECK(law.K(0, 0) == 0.0);
        CHECK(law.K(0, 1) == 0.0);
    }
}
