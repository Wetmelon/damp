// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0

#include <cmath>
#include <limits>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
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
        constexpr auto r1 = poly_roots(damp::array<double, 2>{1.0, 2.0}); // x + 2 = 0
        CHECK(r1.success);
        CHECK(r1.roots[0].real() == doctest::Approx(-2.0));
        CHECK(r1.roots[0].imag() == doctest::Approx(0.0));

        // (x-1)(x-3) = x² - 4x + 3 → {1, -4, 3}
        constexpr auto r2 = poly_roots(damp::array<double, 3>{1.0, -4.0, 3.0});
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
        constexpr TransferFunction<1, 3> tf{.num = {4.0}, .den = {1.0, 2.0, 4.0}};
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
        CHECK(tf2.den[0] == doctest::Approx(1.0));
        CHECK(tf2.den[1] == doctest::Approx(2.0));
        CHECK(tf2.den[2] == doctest::Approx(4.0));
    }

    TEST_CASE("TF → ZPK → TF lead-lag") {
        // C(s) = (s+2)/(s+10)  num=[1,2], den=[1,10]
        constexpr TransferFunction<2, 2> tf{.num = {1.0, 2.0}, .den = {1.0, 10.0}};
        constexpr auto                   art = to_zpk(tf);
        static_assert(art.success);
        CHECK(art.zpk.gain == doctest::Approx(1.0));
        CHECK(art.zpk.zeros[0].real() == doctest::Approx(-2.0));
        CHECK(art.zpk.poles[0].real() == doctest::Approx(-10.0));

        auto tf2 = art.zpk.to_transfer_function();
        CHECK(tf2.num[0] == doctest::Approx(1.0));
        CHECK(tf2.num[1] == doctest::Approx(2.0));
        CHECK(tf2.den[0] == doctest::Approx(1.0));
        CHECK(tf2.den[1] == doctest::Approx(10.0));

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
        constexpr TransferFunction<2, 2> tf{.num = {1.0, 2.0}, .den = {1.0, 1.0}};
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
        auto                             art = to_zpk(ss);
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
        CHECK(tf.den[0] == doctest::Approx(1.0));
        CHECK(tf.den[1] == doctest::Approx(1.0));
        CHECK(tf.den[2] == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(tf.num[0] == doctest::Approx(0.0).epsilon(1e-10));
        // num should be 1 (constant) → {0, 0, 1}
        CHECK(tf.num[2] == doctest::Approx(1.0));
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
        constexpr TransferFunction<2, 2> tf{.num = {1.0, 3.0}, .den = {1.0, 5.0}};
        constexpr auto                   art = to_zpk(tf);
        static_assert(art.success);
        constexpr auto tf2 = art.zpk.to_transfer_function();
        static_assert(tf2.num[0] == 1.0);
        constexpr auto dc = art.zpk.dcgain();
        static_assert(dc.has_value());
        CHECK(*dc == doctest::Approx(3.0 / 5.0));
    }

    TEST_CASE("zero leading num fails TF→ZPK") {
        TransferFunction<2, 2, double> bad{.num = {0.0, 1.0}, .den = {1.0, 1.0}};
        auto                           art = to_zpk(bad);
        CHECK_FALSE(art.success);
    }

    TEST_CASE("SS→ZPK keeps a finite zero of a small-gain plant") {
        // G = 1e-14 (s + 2)/((s + 1)(s + 3))
        const StateSpace<2, 1, 1> sys{
            .A = {{0.0, 1.0}, {-3.0, -4.0}},
            .B = {{0.0}, {1.0}},
            .C = {{2e-14, 1e-14}},
        };
        const auto z = to_zpk(sys);
        REQUIRE(z.success);
        REQUIRE(z.n_zeros == 1);
        CHECK(z.zeros[0].real() == doctest::Approx(-2.0));
        CHECK(z.gain == doctest::Approx(1e-14));
    }

    TEST_CASE("ZPK dcgain is defined for a slow but nonzero pole") {
        // H = 1/(s + 1e-13): DC gain 1e13
        const ZPK<0, 1> sys{.poles = {damp::complex<double>{-1e-13, 0.0}}, .gain = 1.0};
        const auto      dc = sys.dcgain();
        REQUIRE(dc.has_value());
        CHECK(*dc == doctest::Approx(1e13));
    }

    TEST_CASE("ZPK cancel_matching is relative to the root magnitude") {
        // A right-half-plane zero must not cancel a left-half-plane pole just
        // because both are near the origin.
        const ZPK<1, 1> near_origin{
            .zeros = {damp::complex<double>{1e-13, 0.0}},
            .poles = {damp::complex<double>{-1e-13, 0.0}},
        };
        CHECK(near_origin.cancel_matching().n_cancelled == 0);
        // Identical to 12 digits at 1 MHz-scale: a numerical common factor.
        const ZPK<1, 1> fast{
            .zeros = {damp::complex<double>{-1e6 * (1.0 + 1e-12), 0.0}},
            .poles = {damp::complex<double>{-1e6, 0.0}},
        };
        CHECK(fast.cancel_matching().n_cancelled == 1);
        // The roundoff floor scales with the fastest root, but stays far below a 20% gap between
        // two slow roots in float.
        const ZPK<1, 2, float> slow_pair{
            .zeros = {damp::complex<float>{-1.0F, 0.0F}},
            .poles = {damp::complex<float>{-1.2F, 0.0F}, damp::complex<float>{-1e4F, 0.0F}},
        };
        CHECK(slow_pair.cancel_matching().n_cancelled == 0);
    }

    TEST_CASE("cancel_matching cancels a computed common factor at the origin") {
        //! 1/(s + 1) with an unobservable integrator, in coordinates rotated by 30°: the zero of
        //! the common factor s comes out at roundoff level, not at 0.
        const double       c = damp::cos(damp::numbers::pi_v<double> / 6.0);
        const double       s = damp::sin(damp::numbers::pi_v<double> / 6.0);
        const Matrix<2, 2> R{{c, -s}, {s, c}};
        const Matrix<2, 2> J{{0.0, 0.0}, {0.0, -1.0}};
        const StateSpace   sys{.A = R * J * R.transpose(), .B = R * Matrix<2, 1>{{1.0}, {1.0}}, .C = Matrix<1, 2>{{0.0, 1.0}} * R.transpose()};
        const auto         factors = to_zpk(sys).as_zpk<1>();
        REQUIRE(factors.has_value());
        const auto cancelled = factors->cancel_matching();
        CHECK(cancelled.n_cancelled == 1);
        CHECK(damp::abs(cancelled.zpk.poles[0] - damp::complex<double>{-1.0, 0.0}) <= 1e-12);
    }

    TEST_CASE("cancel_matching cancels a double common factor with a tolerance above its spread") {
        //! (s+1)²/((s+1)²(s+2)): the double pole comes out split by about √ε.
        const auto art = to_zpk(TransferFunction<3, 4>{.num = {1.0, 2.0, 1.0}, .den = {1.0, 4.0, 5.0, 2.0}});
        REQUIRE(art.success);
        CHECK(art.zpk.cancel_matching(damp::cbrt(std::numeric_limits<double>::epsilon())).n_cancelled == 2);
    }

    TEST_CASE("to_zpk refuses a constant denominator of zero") {
        CHECK_FALSE(to_zpk(TransferFunction<1, 1>{.num = {1.0}, .den = {0.0}}).success);
        CHECK(to_zpk(TransferFunction<1, 1>{.num = {0.0}, .den = {1.0}}).success); // G = 0 is a model
    }

    TEST_CASE("SS → TF keeps its digits when the poles spread over five decades, in double and float") {
        //! A = diag(−0.1, −1, …, −1e4), B = C = 1: den = ∏(s + pᵢ), num = Σᵢ ∏_{j≠i}(s + pⱼ), sums of
        //! positive terms that the oracle forms without cancellation.
        constexpr damp::array<double, 6> p{0.1, 1.0, 10.0, 100.0, 1e3, 1e4};
        damp::array<double, 7>           den{1.0};
        damp::array<double, 7>           num{};
        for (size_t i = 0; i < 6; ++i) {
            for (size_t k = i + 1; k > 0; --k) {
                den[k] += p[i] * den[k - 1];
            }
            damp::array<double, 6> q{1.0};
            size_t                 degree = 0;
            for (size_t j = 0; j < 6; ++j) {
                if (j != i) {
                    for (size_t k = degree + 1; k > 0; --k) {
                        q[k] += p[j] * q[k - 1];
                    }
                    ++degree;
                }
            }
            for (size_t k = 0; k < 6; ++k) {
                num[k + 1] += q[k];
            }
        }
        const auto check = [&]<typename T>(T den_tol, T num_tol) {
            StateSpace<6, 1, 1, T> sys{};
            for (size_t i = 0; i < 6; ++i) {
                sys.A(i, i) = static_cast<T>(-p[i]);
                sys.B(i, 0) = T{1};
                sys.C(0, i) = T{1};
            }
            const auto tf = to_transfer_function(sys);
            CHECK(tf.num[0] == T{0});
            for (size_t k = 0; k < 7; ++k) {
                CHECK(static_cast<double>(tf.den[k]) == doctest::Approx(den[k]).epsilon(static_cast<double>(den_tol)));
                if (k > 0) {
                    CHECK(static_cast<double>(tf.num[k]) == doctest::Approx(num[k]).epsilon(static_cast<double>(num_tol)));
                }
            }
        };
        check(1e-13, 1e-11);
        check(1e-6F, 1e-3F);
    }

    TEST_CASE("SS → TF keeps a Tustin PI + lead's integral gain in float") {
        //! C(s) = (s + 300)(s + 3000)/(s(s + 2e4)) at 20 kHz. Its Tustin realization has B ~ 1e-9 and
        //! C ~ 1e6, where the float transmission zeros move by percents; num(1) = Σ num carries the
        //! integral gain, about +1.5e-3, and a numerator from those zeros flips its sign.
        const TransferFunction<3, 3, float> pi_lead{.num = {1.0F, 3.3e3F, 9.0e5F}, .den = {1.0F, 2.0e4F, 0.0F}};
        const auto                          realization = pi_lead.to_state_space();
        REQUIRE(realization.has_value());
        const auto sys = discretize(*realization, 5e-5F, DiscretizationMethod::Tustin);
        REQUIRE(sys);
        const auto tf = to_transfer_function(*sys);
        const auto ref = to_transfer_function(sys->as<double>()); // the same float model, in double
        double     num_at_1 = 0.0;
        double     ref_at_1 = 0.0;
        for (size_t k = 0; k < 3; ++k) {
            CHECK(static_cast<double>(tf.num[k]) == doctest::Approx(ref.num[k]).epsilon(1e-5));
            CHECK(static_cast<double>(tf.den[k]) == doctest::Approx(ref.den[k]).epsilon(1e-5));
            num_at_1 += static_cast<double>(tf.num[k]);
            ref_at_1 += ref.num[k];
        }
        CHECK(ref_at_1 == doctest::Approx(1.5e-3).epsilon(1e-2));
        CHECK(num_at_1 / ref_at_1 == doctest::Approx(1.0).epsilon(1e-3));
    }
}
