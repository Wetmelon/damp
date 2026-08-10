// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstdlib>

#include "damp/backend.hpp"
#include "damp/design/pole_placement.hpp"
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

// Largest distance between the (sorted, real) spectrum of A − B·K and the
// requested real poles.
template<size_t NX, size_t NU>
double placement_error(const Matrix<NX, NX>& A, const Matrix<NX, NU>& B, const std::array<double, NX>& poles) {
    auto Kopt = design::place(A, B, poles);
    REQUIRE(Kopt.has_value());
    Matrix<NX, NX> ABK = A - (B * Kopt.value());
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
        auto         Kp = design::place(A, B, std::array<double, 3>{-2.0, -3.0, -4.0});
        REQUIRE(Kp.has_value());

        std::array<damp::complex<double>, 3> pc = {
            damp::complex<double>(-2.0, 0.0), damp::complex<double>(-3.0, 0.0), damp::complex<double>(-4.0, 0.0)
        };
        auto Ka = matlab::acker(A, B, pc);
        REQUIRE(Ka.has_value());

        for (size_t j = 0; j < 3; ++j) {
            CHECK(Kp.value()(0, j) == doctest::Approx(Ka.value()(0, j)).epsilon(1e-9));
        }
    }

    TEST_CASE("place rejects non-assignable problems") {
        Matrix<3, 3> A = {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {-1.0, -2.0, -3.0}};

        SUBCASE("multiplicity exceeding the input count (single input, triple pole)") {
            Matrix<3, 1> B = {{0.0}, {0.0}, {1.0}};
            CHECK_FALSE(design::place(A, B, std::array<double, 3>{-2.0, -2.0, -2.0}).has_value());
        }
        SUBCASE("rank-deficient B") {
            Matrix<3, 2> B = {{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}}; // rank 1
            CHECK_FALSE(design::place(A, B, std::array<double, 3>{-1.0, -2.0, -3.0}).has_value());
        }
    }

    // Largest distance between the spectrum of A − B·K and a complex pole set,
    // comparing as sorted (real, |imag|) pairs.
    static constexpr auto placement_error_cplx = []<size_t NX, size_t NU>(
                                                     const Matrix<NX, NX>&                        A,
                                                     const Matrix<NX, NU>&                        B,
                                                     const std::array<damp::complex<double>, NX>& poles
                                                 ) {
        auto Kopt = design::place(A, B, poles);
        REQUIRE(Kopt.has_value());
        Matrix<NX, NX> ABK = A - (B * Kopt.value());
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
        using Cplx = damp::complex<double>;
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
        using Cplx = damp::complex<double>;
        Matrix<3, 3> A = {{1.0, 2.0, 0.0}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}};
        Matrix<3, 2> B = {{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}};
        // -1+2j without its conjugate is not a real closed-loop spectrum.
        CHECK_FALSE(design::place(A, B, std::array<Cplx, 3>{Cplx(-1, 2), Cplx(-3, 0), Cplx(-4, 0)}).has_value());
    }

    TEST_CASE("matlab::place forwards to design::place") {
        using Cplx = damp::complex<double>;
        Matrix<3, 3>        A = {{1.0, 2.0, 0.0}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}};
        Matrix<3, 2>        B = {{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}};
        std::array<Cplx, 3> poles = {Cplx(-1, 2), Cplx(-1, -2), Cplx(-3, 0)};

        auto Km = matlab::place(A, B, poles);
        auto Kd = design::place(A, B, poles);
        REQUIRE(Km.has_value());
        REQUIRE(Kd.has_value());
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(Km.value()(i, j) == doctest::Approx(Kd.value()(i, j)));
            }
        }
    }

    TEST_CASE("place is constexpr-evaluable") {
        constexpr bool ok = []() consteval {
            Matrix<2, 2> A = {{0.0, 1.0}, {0.0, 0.0}};
            Matrix<2, 1> B = {{0.0}, {1.0}};
            auto         K = design::place(A, B, std::array<double, 2>{-1.0, -2.0});
            if (!K) {
                return false;
            }
            // Double-integrator with poles {-1,-2}: char poly s²+3s+2 ⇒ K = [2, 3].
            return damp::abs(K.value()(0, 0) - 2.0) < 1e-9 && damp::abs(K.value()(0, 1) - 3.0) < 1e-9;
        }();
        static_assert(ok, "place must work at compile time");
        CHECK(ok);
    }

    TEST_CASE("place(A,B,poles,Ts) matches manual ZOH + z-map + place_discrete") {
        Matrix<2, 2>           A = {{0.0, 1.0}, {0.0, 0.0}};
        Matrix<2, 1>           B = {{0.0}, {1.0}};
        const double           Ts = 0.01;
        damp::array<double, 2> poles_s{-4.0, -6.0};

        auto K_auto = design::place(A, B, poles_s, Ts);
        REQUIRE(K_auto.has_value());

        StateSpace<2, 1, 2>    sys_c{.A = A, .B = B, .C = Matrix<2, 2>::identity(), .Ts = 0.0};
        const auto             sys_d = *discretize(sys_c, Ts, DiscretizationMethod::ZOH);
        damp::array<double, 2> poles_z{damp::exp(poles_s[0] * Ts), damp::exp(poles_s[1] * Ts)};
        auto                   K_manual = design::place_discrete(sys_d.A, sys_d.B, poles_z);
        REQUIRE(K_manual.has_value());

        for (size_t j = 0; j < 2; ++j) {
            CHECK(K_auto.value()(0, j) == doctest::Approx(K_manual.value()(0, j)).epsilon(1e-9));
        }

        // Discrete closed-loop poles of Ad − Bd K near e^{s Ts}
        Matrix<2, 2>          AdBK = sys_d.A - (sys_d.B * K_auto.value());
        auto                  ev = mat::compute_eigenvalues(AdBK);
        std::array<double, 2> got{ev.values[0].real(), ev.values[1].real()};
        std::array<double, 2> want{poles_z[0], poles_z[1]};
        std::sort(got.begin(), got.end());
        std::sort(want.begin(), want.end());
        CHECK(std::abs(got[0] - want[0]) < 1e-8);
        CHECK(std::abs(got[1] - want[1]) < 1e-8);

        CHECK_FALSE(design::place(A, B, poles_s, 0.0).has_value());
    }

    TEST_CASE("place_discrete(Ad,Bd,poles_s,Ts) maps s-poles without re-c2d") {
        Matrix<2, 2>        A = {{0.0, 1.0}, {0.0, 0.0}};
        Matrix<2, 1>        B = {{0.0}, {1.0}};
        const double        Ts = 0.01;
        StateSpace<2, 1, 2> sys_c{.A = A, .B = B, .C = Matrix<2, 2>::identity(), .Ts = 0.0};
        const auto          sys_d = *discretize(sys_c, Ts, DiscretizationMethod::ZOH);

        damp::array<double, 2> poles_s{-5.0, -7.0};
        auto                   K_map = design::place_discrete(sys_d.A, sys_d.B, poles_s, Ts);
        auto                   K_z = design::place_discrete(
            sys_d.A, sys_d.B, damp::array<double, 2>{damp::exp(-5.0 * Ts), damp::exp(-7.0 * Ts)}
        );
        REQUIRE(K_map.has_value());
        REQUIRE(K_z.has_value());
        for (size_t j = 0; j < 2; ++j) {
            CHECK(K_map.value()(0, j) == doctest::Approx(K_z.value()(0, j)).epsilon(1e-9));
        }
    }

    TEST_CASE("place(sys, poles) forwards A,B (domain follows sys.Ts)") {
        Matrix<2, 2>           A = {{0.0, 1.0}, {0.0, 0.0}};
        Matrix<2, 1>           B = {{0.0}, {1.0}};
        StateSpace             sys{.A = A, .B = B, .C = Matrix<1, 2>{{1.0, 0.0}}, .Ts = 0.0};
        damp::array<double, 2> poles{-1.0, -2.0};
        auto                   K_sys = design::place(sys, poles);
        auto                   K_ab = design::place(A, B, poles);
        REQUIRE(K_sys.has_value());
        REQUIRE(K_ab.has_value());
        CHECK(K_sys.value()(0, 0) == doctest::Approx(K_ab.value()(0, 0)));
        CHECK(K_sys.value()(0, 1) == doctest::Approx(K_ab.value()(0, 1)));
    }

    TEST_CASE("place_observer is dual of place") {
        Matrix<2, 2>           A = {{0.0, 1.0}, {0.0, 0.0}};
        Matrix<1, 2>           C = {{1.0, 0.0}};
        damp::array<double, 2> poles{-3.0, -4.0};

        auto L = design::place_observer(A, C, poles);
        auto Kt = design::place(A.transpose(), C.transpose(), poles);
        REQUIRE(L.has_value());
        REQUIRE(Kt.has_value());
        // L = place(Aᵀ, Cᵀ, p)ᵀ
        CHECK(L.value()(0, 0) == doctest::Approx(Kt.value()(0, 0)).epsilon(1e-9));
        CHECK(L.value()(1, 0) == doctest::Approx(Kt.value()(0, 1)).epsilon(1e-9));

        Matrix<2, 2>          ALC = A - (L.value() * C);
        auto                  ev = mat::compute_eigenvalues(ALC);
        std::array<double, 2> got{ev.values[0].real(), ev.values[1].real()};
        std::array<double, 2> want{-3.0, -4.0};
        std::sort(got.begin(), got.end());
        std::sort(want.begin(), want.end());
        CHECK(std::abs(got[0] - want[0]) < 1e-8);
        CHECK(std::abs(got[1] - want[1]) < 1e-8);
    }

    TEST_CASE("place_observer(A,C,poles,Ts) matches discrete dual after ZOH") {
        Matrix<2, 2>           A = {{0.0, 1.0}, {0.0, 0.0}};
        Matrix<1, 2>           C = {{1.0, 0.0}};
        const double           Ts = 0.01;
        damp::array<double, 2> poles_s{-8.0, -10.0};

        auto L_auto = design::place_observer(A, C, poles_s, Ts);
        REQUIRE(L_auto.has_value());

        StateSpace<2, 1, 1>    sys_c{.A = A, .B = Matrix<2, 1>{{0.0}, {1.0}}, .C = C, .Ts = 0.0};
        const auto             sys_d = *discretize(sys_c, Ts, DiscretizationMethod::ZOH);
        damp::array<double, 2> poles_z{damp::exp(poles_s[0] * Ts), damp::exp(poles_s[1] * Ts)};
        auto                   L_manual = design::place_observer_discrete(sys_d.A, sys_d.C, poles_z);
        REQUIRE(L_manual.has_value());
        for (size_t i = 0; i < 2; ++i) {
            CHECK(L_auto.value()(i, 0) == doctest::Approx(L_manual.value()(i, 0)).epsilon(1e-9));
        }
    }
}
