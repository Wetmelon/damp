// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Invariant zeros: hand-built models with known zeros. The scipy comparisons are in
// test_oracle_lti.cpp.

#include <cmath>
#include <limits>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/zeros.hpp"
#include "damp/systems/zpk.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

/// Every expected root has a distinct computed root within tol·max(1, |root|).
template<size_t N, size_t M>
bool contains_roots(const damp::array<damp::complex<double>, N>& got, size_t n_got, const damp::array<damp::complex<double>, M>& want, double tol) {
    if (n_got != M) {
        return false;
    }
    damp::array<bool, N> used{};
    for (size_t i = 0; i < M; ++i) {
        bool found = false;
        for (size_t j = 0; j < n_got && !found; ++j) {
            if (!used[j] && damp::abs(got[j] - want[i]) <= tol * damp::max(1.0, damp::abs(want[i]))) {
                used[j] = true;
                found = true;
            }
        }
        if (!found) {
            return false;
        }
    }
    return true;
}

/// G(s) = [1/(s+1) 2/(s+2); 1/(s+3) 1/(s+4)]: det G·∏(s + k) = −(s² + 5s + 2), zeros (−5 ± √17)/2.
constexpr StateSpace<4, 2, 2> kSquare{
    .A = Matrix<4, 4>{{-1.0, 0.0, 0.0, 0.0}, {0.0, -2.0, 0.0, 0.0}, {0.0, 0.0, -3.0, 0.0}, {0.0, 0.0, 0.0, -4.0}},
    .B = Matrix<4, 2>{{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}},
    .C = Matrix<2, 4>{{1.0, 2.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 1.0}},
};

const damp::array<damp::complex<double>, 2> kSquareZeros{
    damp::complex<double>{(-5.0 + std::sqrt(17.0)) / 2.0, 0.0},
    damp::complex<double>{(-5.0 - std::sqrt(17.0)) / 2.0, 0.0},
};

} // namespace

TEST_SUITE("Transmission zeros") {
    TEST_CASE("a square MIMO model, constant-evaluated") {
        constexpr auto z = transmission_zeros(kSquare);
        static_assert(z.success && z.n_zeros == 2);
        CHECK(contains_roots(z.zeros, z.n_zeros, kSquareZeros, 1e-13));
    }

    TEST_CASE("matlab::tzero forwards to transmission_zeros") {
        const auto z = matlab::tzero(kSquare);
        REQUIRE(z.success);
        CHECK(contains_roots(z.zeros, z.n_zeros, kSquareZeros, 1e-13));
    }

    TEST_CASE("the zeros do not depend on the units of time, inputs, outputs or states") {
        // Time in milliseconds (A, B × 10³), inputs × 10⁶, outputs × 10⁻⁶, a state in mm: zeros × 10³.
        StateSpace<4, 2, 2> scaled = kSquare;
        scaled.A = scaled.A * 1e3;
        scaled.B = scaled.B * 1e9;
        scaled.C = scaled.C * 1e-6;
        for (size_t j = 0; j < 4; ++j) {
            scaled.A(0, j) *= 1e3;
            scaled.A(j, 0) *= 1e-3;
        }
        scaled.B(0, 0) *= 1e3;
        scaled.C(0, 0) *= 1e-3;
        const auto z = transmission_zeros(scaled);
        REQUIRE(z.success);
        damp::array<damp::complex<double>, 2> want{kSquareZeros[0] * 1e3, kSquareZeros[1] * 1e3};
        CHECK(contains_roots(z.zeros, z.n_zeros, want, 1e-12));
    }

    TEST_CASE("SISO zeros agree with to_zpk") {
        // (s − 2)(s + 3)/((s + 1)(s + 2)(s + 4)), controllable companion form
        const StateSpace<3, 1, 1> sys{
            .A = Matrix<3, 3>{{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {-8.0, -14.0, -7.0}},
            .B = Matrix<3, 1>{{0.0}, {0.0}, {1.0}},
            .C = Matrix<1, 3>{{-6.0, 1.0, 1.0}},
        };
        const auto tz = transmission_zeros(sys);
        const auto zp = to_zpk(sys);
        REQUIRE(tz.success);
        REQUIRE(zp.success);
        REQUIRE(tz.n_zeros == zp.n_zeros);
        damp::array<damp::complex<double>, 2> want{zp.zeros[0], zp.zeros[1]};
        CHECK(contains_roots(tz.zeros, tz.n_zeros, want, 1e-12));
        const damp::array<damp::complex<double>, 2> exact{damp::complex<double>{2.0, 0.0}, damp::complex<double>{-3.0, 0.0}};
        CHECK(contains_roots(tz.zeros, tz.n_zeros, exact, 1e-13));
    }

    TEST_CASE("a strictly proper chain has its zeros at infinity") {
        const StateSpace<3, 1, 1> chain{
            .A = Matrix<3, 3>{{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {-6.0, -11.0, -6.0}},
            .B = Matrix<3, 1>{{0.0}, {0.0}, {1.0}},
            .C = Matrix<1, 3>{{1.0, 0.0, 0.0}},
        };
        const auto z = transmission_zeros(chain);
        CHECK(z.success);
        CHECK(z.n_zeros == 0);
    }

    TEST_CASE("a small D gives the large zeros it implies; a D at roundoff level gives none") {
        // (s + 1)/((s + 1)(s + 2)(s + 3)) + D: zeros −1 and the roots of D(s² + 5s + 6) + 1.
        StateSpace<3, 1, 1> sys{
            .A = Matrix<3, 3>{{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {-6.0, -11.0, -6.0}},
            .B = Matrix<3, 1>{{0.0}, {0.0}, {1.0}},
            .C = Matrix<1, 3>{{1.0, 1.0, 0.0}},
        };
        sys.D(0, 0) = 1e-6;
        const auto big = transmission_zeros(sys);
        REQUIRE(big.success);
        const double                                im = std::sqrt(1e6 - 0.25);
        const damp::array<damp::complex<double>, 3> want{
            damp::complex<double>{-1.0, 0.0},
            damp::complex<double>{-2.5, im},
            damp::complex<double>{-2.5, -im},
        };
        CHECK(contains_roots(big.zeros, big.n_zeros, want, 1e-10));
        sys.D(0, 0) = 1e-20;
        const auto none = transmission_zeros(sys);
        REQUIRE(none.success);
        const damp::array<damp::complex<double>, 1> finite{damp::complex<double>{-1.0, 0.0}};
        CHECK(contains_roots(none.zeros, none.n_zeros, finite, 1e-13));
    }

    TEST_CASE("a sampled double integrator has its sampling zero at −1") {
        const StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .Ts = 0.1,
        };
        const auto                                  z = transmission_zeros(sys);
        const damp::array<damp::complex<double>, 1> want{damp::complex<double>{-1.0, 0.0}};
        REQUIRE(z.success);
        CHECK(contains_roots(z.zeros, z.n_zeros, want, 1e-13));
    }

    TEST_CASE("a non-square model: a common zero of both outputs, and its dual") {
        // f = (s − 2)/(s + 1) drives two first-order outputs: y = [1/(s+2); 1/(s+2) + 1/(s+5)]·f.
        const StateSpace<3, 1, 2> tall{
            .A = Matrix<3, 3>{{-1.0, 0.0, 0.0}, {-3.0, -2.0, 0.0}, {-3.0, 0.0, -5.0}},
            .B = Matrix<3, 1>{{1.0}, {1.0}, {1.0}},
            .C = Matrix<2, 3>{{0.0, 1.0, 0.0}, {0.0, 1.0, 1.0}},
        };
        const damp::array<damp::complex<double>, 1> want{damp::complex<double>{2.0, 0.0}};
        const auto                                  zt = transmission_zeros(tall);
        REQUIRE(zt.success);
        CHECK(contains_roots(zt.zeros, zt.n_zeros, want, 1e-13));
        const StateSpace<3, 2, 1> wide{.A = tall.A.transpose(), .B = tall.C.transpose(), .C = tall.B.transpose()};
        const auto                zw = transmission_zeros(wide);
        REQUIRE(zw.success);
        CHECK(contains_roots(zw.zeros, zw.n_zeros, want, 1e-13));
    }

    TEST_CASE("an unobservable mode is a zero of a tall model; an uncontrollable one is not") {
        const Matrix<3, 3>        A{{-1.0, 1.0, 0.0}, {0.0, -2.0, 1.0}, {0.0, 0.0, -3.0}}; // x₁ is not seen by C
        const StateSpace<3, 1, 2> unobservable{
            .A = A,
            .B = Matrix<3, 1>{{0.0}, {0.0}, {1.0}},
            .C = Matrix<2, 3>{{0.0, 1.0, 0.0}, {0.0, 1.0, 1.0}},
        };
        const damp::array<damp::complex<double>, 1> want{damp::complex<double>{-1.0, 0.0}};
        const auto                                  z = transmission_zeros(unobservable);
        REQUIRE(z.success);
        CHECK(contains_roots(z.zeros, z.n_zeros, want, 1e-13));
        // x₁ = −x₁ alone: uncontrollable, so [A − λI, B] loses a row at −1, but P(−1) keeps full column rank.
        const StateSpace<3, 1, 2> uncontrollable{
            .A = Matrix<3, 3>{{-1.0, 0.0, 0.0}, {0.0, -2.0, 1.0}, {0.0, 0.0, -3.0}},
            .B = Matrix<3, 1>{{0.0}, {0.0}, {1.0}},
            .C = Matrix<2, 3>{{1.0, 1.0, 0.0}, {0.0, 1.0, 1.0}},
        };
        const auto none = transmission_zeros(uncontrollable);
        REQUIRE(none.success);
        CHECK(none.n_zeros == 0);
    }

    TEST_CASE("a zero model output has no zeros") {
        const StateSpace<2, 1, 1> sys{.A = Matrix<2, 2>{{0.0, 1.0}, {-2.0, -3.0}}, .B = Matrix<2, 1>{{0.0}, {1.0}}};
        const auto                z = transmission_zeros(sys);
        CHECK(z.success);
        CHECK(z.n_zeros == 0);
    }

    TEST_CASE("a non-finite entry is refused as bad_weights") {
        StateSpace<4, 2, 2> bad = kSquare;
        bad.C(1, 2) = std::numeric_limits<double>::quiet_NaN();
        const auto z = transmission_zeros(bad);
        CHECK_FALSE(z.success);
        CHECK(z.status == design::DesignStatus::bad_weights);
    }

    TEST_CASE("float") {
        const auto z = transmission_zeros(kSquare.as<float>());
        REQUIRE(z.success);
        const auto zd = z.as<double>();
        CHECK(contains_roots(zd.zeros, zd.n_zeros, kSquareZeros, 1e-5));
    }
}
