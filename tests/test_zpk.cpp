// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0

#include <cmath>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

/// Match complex roots unordered within tol (pairwise greedy).
template<size_t N>
bool roots_match(
    const damp::array<damp::complex<double>, N>& a,
    const damp::array<damp::complex<double>, N>& b,
    double                                       tol = 1e-8
) {
    damp::array<bool, N> used{};
    for (size_t i = 0; i < N; ++i) {
        bool found = false;
        for (size_t j = 0; j < N; ++j) {
            if (used[j]) {
                continue;
            }
            if (damp::abs(a[i] - b[j]) <= tol) {
                used[j] = true;
                found = true;
                break;
            }
        }
        if (!found) {
            return false;
        }
    }
    return true;
}

} // namespace

TEST_SUITE("ZPK first-class LTI") {

    TEST_CASE("poly_roots degree 1 and 2") {
        constexpr auto r1 = poly_roots(damp::array<double, 2>{2.0, 1.0}); // x + 2 = 0
        CHECK(r1.success);
        CHECK(r1.roots[0].real() == doctest::Approx(-2.0));
        CHECK(r1.roots[0].imag() == doctest::Approx(0.0));

        // (x-1)(x-3) = x² - 4x + 3 → ascending [3, -4, 1]
        constexpr auto r2 = poly_roots(damp::array<double, 3>{3.0, -4.0, 1.0});
        CHECK(r2.success);
        damp::array<damp::complex<double>, 2> expect{
            damp::complex<double>{1.0, 0.0},
            damp::complex<double>{3.0, 0.0}
        };
        CHECK(roots_match(r2.roots, expect));
    }

    TEST_CASE("TF → ZPK → TF first-order") {
        // G = 1/(s+1)  num=[1], den=[1,1]
        constexpr TransferFunction<1, 2> tf{.num = {1.0}, .den = {1.0, 1.0}};
        constexpr auto                   art = to_zpk(tf);
        static_assert(art.success);
        CHECK(art.success);
        CHECK(art.zpk.gain == doctest::Approx(1.0));
        CHECK(art.zpk.poles[0].real() == doctest::Approx(-1.0));

        constexpr auto tf2 = art.zpk.to_transfer_function();
        CHECK(tf2.num[0] == doctest::Approx(1.0));
        CHECK(tf2.den[0] == doctest::Approx(1.0));
        CHECK(tf2.den[1] == doctest::Approx(1.0));
    }

    TEST_CASE("TF → ZPK → TF second-order underdamped") {
        // ωn=2, ζ=0.5 → den = s² + 2s + 4, num = 4
        constexpr TransferFunction<1, 3> tf{.num = {4.0}, .den = {4.0, 2.0, 1.0}};
        constexpr auto                   art = to_zpk(tf);
        static_assert(art.success);
        CHECK(art.zpk.gain == doctest::Approx(4.0));

        const double                          sigma = -0.5 * 2.0;
        const double                          wd = 2.0 * std::sqrt(1.0 - 0.5 * 0.5);
        damp::array<damp::complex<double>, 2> expect{
            damp::complex<double>{sigma, wd},
            damp::complex<double>{sigma, -wd}
        };
        CHECK(roots_match(art.zpk.poles, expect, 1e-9));

        auto tf2 = art.zpk.to_transfer_function();
        CHECK(tf2.num[0] == doctest::Approx(4.0));
        CHECK(tf2.den[0] == doctest::Approx(4.0));
        CHECK(tf2.den[1] == doctest::Approx(2.0));
        CHECK(tf2.den[2] == doctest::Approx(1.0));
    }

    TEST_CASE("TF → ZPK → TF lead-lag") {
        // C(s) = (s+2)/(s+10)  num=[2,1], den=[10,1]
        constexpr TransferFunction<2, 2> tf{.num = {2.0, 1.0}, .den = {10.0, 1.0}};
        constexpr auto                   art = tf2zpk(tf);
        static_assert(art.success);
        CHECK(art.zpk.gain == doctest::Approx(1.0));
        CHECK(art.zpk.zeros[0].real() == doctest::Approx(-2.0));
        CHECK(art.zpk.poles[0].real() == doctest::Approx(-10.0));

        auto tf2 = zpk2tf(art.zpk);
        CHECK(tf2.num[0] == doctest::Approx(2.0));
        CHECK(tf2.num[1] == doctest::Approx(1.0));
        CHECK(tf2.den[0] == doctest::Approx(10.0));
        CHECK(tf2.den[1] == doctest::Approx(1.0));

        // DC gain H(0) = 2/10 = 0.2
        auto dc = art.zpk.dcgain();
        REQUIRE(dc.has_value());
        CHECK(*dc == doctest::Approx(0.2));
    }

    TEST_CASE("ZPK cancel_matching / minreal") {
        // (s+1)(s+2) / (s+1)(s+3) with explicit factors
        ZPK<2, 2, double> sys{
            .zeros = {damp::complex<double>{-1.0, 0.0}, damp::complex<double>{-2.0, 0.0}},
            .poles = {damp::complex<double>{-1.0, 0.0}, damp::complex<double>{-3.0, 0.0}},
            .gain = 1.0
        };

        auto cr = sys.cancel_matching();
        CHECK(cr.n_cancelled == 1);

        // Active front: zero -2, pole -3; cancelled pair -1 at end
        CHECK(cr.zpk.zeros[0].real() == doctest::Approx(-2.0));
        CHECK(cr.zpk.poles[0].real() == doctest::Approx(-3.0));
        CHECK(cr.zpk.zeros[1].real() == doctest::Approx(-1.0));
        CHECK(cr.zpk.poles[1].real() == doctest::Approx(-1.0));

        // After cancel common (s+1): expansion still matches DC of (s+2)/(s+3)
        auto dc = cr.zpk.dcgain();
        REQUIRE(dc.has_value());
        CHECK(*dc == doctest::Approx(2.0 / 3.0));

        auto mr = design::minreal_zpk(sys);
        CHECK(mr.n_cancelled == 1);
    }

    TEST_CASE("SS biproper → ZPK") {
        // G = (s+2)/(s+1) → A=-1, B=1, C=1, D=1 (from TF biproper realization)
        constexpr TransferFunction<2, 2> tf{.num = {2.0, 1.0}, .den = {1.0, 1.0}};
        constexpr auto                   ss = tf.to_state_space().value();
        constexpr auto                   art = to_zpk(ss);
        static_assert(art.success);
        CHECK(art.n_zeros == 1);
        CHECK(art.gain == doctest::Approx(1.0));
        CHECK(art.poles[0].real() == doctest::Approx(-1.0));
        CHECK(art.zeros[0].real() == doctest::Approx(-2.0));

        auto zpk_opt = art.as_zpk<1>();
        REQUIRE(zpk_opt.has_value());
        auto dc = zpk_opt->dcgain();
        REQUIRE(dc.has_value());
        CHECK(*dc == doctest::Approx(2.0));
    }

    TEST_CASE("SS strictly proper → ZPK (first order)") {
        constexpr TransferFunction<1, 2> tf{.num = {1.0}, .den = {1.0, 1.0}};
        constexpr auto                   ss = tf.to_state_space().value();
        auto                             art = ss2zpk(ss);
        CHECK(art.success);
        CHECK(art.n_zeros == 0);
        CHECK(art.poles[0].real() == doctest::Approx(-1.0));
        CHECK(art.gain == doctest::Approx(1.0));
    }

    TEST_CASE("SS → TF Leverrier matches known plant") {
        // Double integrator-like: A = [0 1; 0 -1], B = [0;1], C = [1 0], D = 0
        // G(s) = 1 / (s(s+1)) = 1/(s²+s)
        constexpr StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{0.0, 1.0}, {0.0, -1.0}},
            .B = Matrix<2, 1>{{0.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>{{0.0}},
        };
        auto tf = ss2tf(sys);
        CHECK(tf.den[2] == doctest::Approx(1.0));
        CHECK(tf.den[1] == doctest::Approx(1.0));
        CHECK(tf.den[0] == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(tf.num[2] == doctest::Approx(0.0).epsilon(1e-10));
        // num should be 1 (constant) → ascending [1, 0, 0]
        CHECK(tf.num[0] == doctest::Approx(1.0));
        CHECK(tf.num[1] == doctest::Approx(0.0).epsilon(1e-10));

        auto art = to_zpk(sys);
        CHECK(art.success);
        CHECK(art.n_zeros == 0);
        CHECK(art.gain == doctest::Approx(1.0));
        damp::array<damp::complex<double>, 2> expect{
            damp::complex<double>{0.0, 0.0},
            damp::complex<double>{-1.0, 0.0}
        };
        CHECK(roots_match(art.poles, expect, 1e-8));
    }

    TEST_CASE("dcgain_discrete at z=1") {
        // H(z) = 0.5 (z+1)/(z-0.5) → H(1) = 0.5 * 2 / 0.5 = 2
        ZPK<1, 1, double> hz{
            .zeros = {damp::complex<double>{-1.0, 0.0}},
            .poles = {damp::complex<double>{0.5, 0.0}},
            .gain = 0.5
        };
        auto dc = hz.dcgain_discrete();
        REQUIRE(dc.has_value());
        CHECK(*dc == doctest::Approx(2.0));
    }

    TEST_CASE("constexpr smoke ZPK round-trip") {
        constexpr TransferFunction<2, 2> tf{.num = {3.0, 1.0}, .den = {5.0, 1.0}};
        constexpr auto                   art = to_zpk(tf);
        static_assert(art.success);
        constexpr auto tf2 = art.zpk.to_transfer_function();
        static_assert(tf2.num[1] == 1.0);
        constexpr auto dc = art.zpk.dcgain();
        static_assert(dc.has_value());
        CHECK(*dc == doctest::Approx(3.0 / 5.0));
    }

    TEST_CASE("zero leading num fails TF→ZPK") {
        TransferFunction<2, 2, double> bad{.num = {1.0, 0.0}, .den = {1.0, 1.0}};
        auto                           art = to_zpk(bad);
        CHECK_FALSE(art.success);
    }
}
