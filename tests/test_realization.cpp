// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <limits>

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/realization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
// A well-scaled 3-state model, then the same plant with states in wildly different units
// (x_bad = S·x, S = diag(1e6, 1e-3, 1)).
constexpr StateSpace<3, 1, 1> make_nice() {
    return StateSpace<3, 1, 1>{
        .A = Matrix<3, 3>{{-1.0, 2.0, 0.0}, {-3.0, -4.0, 1.0}, {0.5, 0.0, -2.0}},
        .B = Matrix<3, 1>{{1.0}, {0.5}, {0.0}},
        .C = Matrix<1, 3>{{1.0, 0.0, 2.0}},
        .D = Matrix<1, 1>{{0.25}},
    };
}

constexpr StateSpace<3, 1, 1> make_badly_scaled() {
    const auto          nice = make_nice();
    const ColVec<3>     s{1e6, 1e-3, 1.0};
    StateSpace<3, 1, 1> bad = nice;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            bad.A(i, j) = s[i] * nice.A(i, j) / s[j];
        }
        bad.B(i, 0) = s[i] * nice.B(i, 0);
        bad.C(0, i) = nice.C(0, i) / s[i];
    }
    return bad;
}

// Spread of the nonzero off-diagonal / B / C magnitudes: max / min.
template<size_t NX>
double spread(const StateSpace<NX, 1, 1>& sys) {
    double     lo = std::numeric_limits<double>::max();
    double     hi = 0.0;
    const auto take = [&](double v) {
        const double a = damp::abs(v);
        if (a > 0.0) {
            lo = damp::min(lo, a);
            hi = damp::max(hi, a);
        }
    };
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NX; ++j) {
            if (i != j) {
                take(sys.A(i, j));
            }
        }
        take(sys.B(i, 0));
        take(sys.C(0, i));
    }
    return hi / lo;
}

// Markov parameter C·Aᵏ·B.
template<size_t NX>
double markov(const StateSpace<NX, 1, 1>& sys, size_t k) {
    Matrix<NX, 1> v = sys.B;
    for (size_t i = 0; i < k; ++i) {
        v = Matrix<NX, 1>(sys.A * v);
    }
    return (sys.C * v)(0, 0);
}
} // namespace

TEST_SUITE("State balancing (prescale / ssbal)") {
    TEST_CASE("balance_states undoes a bad choice of units, exactly and without changing the I/O map") {
        constexpr auto bad = make_badly_scaled();
        constexpr auto bal = balance_states(bad);
        CHECK(spread(bad) > 1e8);
        CHECK(spread(bal.sys) < 1e2);

        //! Powers of two: undoing the scaling gives back every entry bit for bit.
        const auto s = bal.scale();
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(bal.sys.A(i, j) * s[i] / s[j] == bad.A(i, j)); // A = D·Ã·D⁻¹
            }
            CHECK(bal.sys.B(i, 0) * s[i] == bad.B(i, 0));
            CHECK(bal.sys.C(0, i) / s[i] == bad.C(0, i));
        }
        CHECK(bal.sys.D == bad.D);

        //! Same transfer function: the Markov parameters agree.
        for (size_t k = 0; k < 6; ++k) {
            CHECK(markov(bal.sys, k) == doctest::Approx(markov(bad, k)).epsilon(1e-12));
        }
    }

    TEST_CASE("balance_states is constexpr, float-ready, and leaves a balanced or non-finite model alone") {
        constexpr auto bal = balance_states(make_badly_scaled());
        static_assert(bal.exponent[0] > 10 && bal.exponent[1] < -5); // x_bad = D·x̃ recovers S = diag(1e6 ≈ 2^20, 1e-3 ≈ 2^-10, 1)
        constexpr auto bal_f = balance_states(make_badly_scaled().as<float>());
        static_assert(bal_f.exponent[0] == bal.exponent[0]);
        CHECK(bal.as<float>().sys.A(0, 1) == doctest::Approx(bal_f.sys.A(0, 1)));

        constexpr StateSpace<2, 1, 1> even{
            .A = Matrix<2, 2>{{-1.0, 1.0}, {-1.0, -1.0}},
            .B = Matrix<2, 1>{{1.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 1.0}}
        };
        constexpr auto same = balance_states(even);
        static_assert(same.exponent[0] == 0 && same.exponent[1] == 0 && same.sys.A == even.A);

        auto broken = make_badly_scaled();
        broken.A(0, 1) = std::numeric_limits<double>::quiet_NaN();
        const auto untouched = balance_states(broken);
        CHECK(untouched.exponent[0] == 0);
        CHECK(untouched.exponent[2] == 0);
        CHECK(untouched.sys.B == broken.B);
    }

    TEST_CASE("balance_states counts each state's own dynamics: a dominant diagonal is left alone") {
        //! Couplings of 100 and 0.01 beside poles at −1000: each state's row and column already
        //! agree once the diagonal counts, as in xGEBAL since LAPACK 3.5.0. Without it the
        //! couplings alone set the scale: exponents +2 and −2.
        constexpr StateSpace<2, 1, 1> stiff{
            .A = Matrix<2, 2>{{-1000.0, 100.0}, {0.01, -1000.0}},
            .B = Matrix<2, 1>{{1.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 1.0}}
        };
        constexpr auto same = balance_states(stiff);
        static_assert(same.exponent[0] == 0 && same.exponent[1] == 0 && same.sys.A == stiff.A);
    }

    TEST_CASE("matlab::prescale and matlab::ssbal forward to balance_states") {
        constexpr auto bad = make_badly_scaled();
        static_assert(matlab::prescale(bad).exponent == balance_states(bad).exponent);
        static_assert(matlab::ssbal(bad).sys.A == balance_states(bad).sys.A);
    }
}

namespace {
// Non-normal 4-state model: one complex pair and two real modes, badly scaled on purpose.
constexpr StateSpace<4, 1, 2, double, 1, 1> make_modal_plant() {
    return StateSpace<4, 1, 2, double, 1, 1>{
        .A = Matrix<4, 4>{{-1.0, 2.0, 0.5, 0.0}, {-3.0, -1.0, 1.0, 0.2}, {0.3, 0.0, -4.0, 10.0}, {0.0, 0.1, 0.0, -0.5}},
        .B = Matrix<4, 1>{{1.0}, {0.0}, {2.0}, {-1.0}},
        .C = Matrix<2, 4>{{1.0, 0.0, 0.0, 1.0}, {0.0, 1.0, 0.5, 0.0}},
        .D = Matrix<2, 1>{{0.0}, {0.1}},
        .G = Matrix<4, 1>{{0.0}, {1.0}, {0.0}, {0.0}},
        .H = Matrix<2, 1>{{1.0}, {1.0}},
    };
}
} // namespace

TEST_SUITE("Modal canonical form (canon)") {
    TEST_CASE("modal_form: block-diagonal A, the same modes, the same I/O map") {
        constexpr auto sys = make_modal_plant();
        constexpr auto res = modal_form(sys);
        static_assert(res.success && res.status == design::DesignStatus::ok);
        const auto& Am = res.sys.A;

        //! Exactly block diagonal, pairs as [σ ω; −ω σ] with ω > 0.
        size_t pairs = 0;
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                const bool in_pair = (i + 1 == j && Am(i, j) != 0.0) || (j + 1 == i && Am(j, i) != 0.0);
                if (i != j && !in_pair) {
                    CHECK(Am(i, j) == 0.0);
                }
            }
            if (i + 1 < 4 && Am(i, i + 1) != 0.0) {
                ++pairs;
                CHECK(Am(i, i + 1) > 0.0);
                CHECK(Am(i + 1, i) == -Am(i, i + 1));
                CHECK(Am(i, i) == Am(i + 1, i + 1));
            }
        }
        CHECK(pairs == 1);

        //! Same eigenvalues: the modal blocks' trace and determinant match A's.
        const auto ev = mat::compute_eigenvalues(sys.A);
        double     tr = 0.0;
        for (size_t i = 0; i < 4; ++i) {
            tr += ev.values[i].real();
        }
        CHECK(Am(0, 0) + Am(1, 1) + Am(2, 2) + Am(3, 3) == doctest::Approx(tr).epsilon(1e-12));
        CHECK(mat::det(Am) == doctest::Approx(mat::det(sys.A)).epsilon(1e-10));

        //! z = T·x: T·A = A_modal·T, and the Markov parameters C Aᵏ B and C Aᵏ G agree.
        const Matrix<4, 4> lhs = res.transform * sys.A;
        const Matrix<4, 4> rhs = Am * res.transform;
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                CHECK(lhs(i, j) == doctest::Approx(rhs(i, j)).epsilon(1e-9).scale(1.0));
            }
        }
        Matrix<4, 1> v = sys.B;
        Matrix<4, 1> vm = res.sys.B;
        Matrix<4, 1> g = sys.G;
        Matrix<4, 1> gm = res.sys.G;
        for (size_t k = 0; k < 5; ++k) {
            const Matrix<2, 1> y = sys.C * v;
            const Matrix<2, 1> ym = res.sys.C * vm;
            const Matrix<2, 1> yg = sys.C * g;
            const Matrix<2, 1> ygm = res.sys.C * gm;
            for (size_t r = 0; r < 2; ++r) {
                CHECK(ym(r, 0) == doctest::Approx(y(r, 0)).epsilon(1e-9).scale(1.0));
                CHECK(ygm(r, 0) == doctest::Approx(yg(r, 0)).epsilon(1e-9).scale(1.0));
            }
            v = sys.A * v;
            vm = Am * vm;
            g = sys.A * g;
            gm = Am * gm;
        }
        CHECK(res.sys.D == sys.D);
        CHECK(res.sys.Ts == sys.Ts);
    }

    TEST_CASE("modal_form: float, repeated modes, and a defective A") {
        constexpr auto res_f = modal_form(make_modal_plant().as<float>());
        static_assert(res_f.success);
        constexpr auto res_d = modal_form(make_modal_plant());
        CHECK(mat::det(res_f.sys.A) == doctest::Approx(mat::det(res_d.sys.A)).epsilon(1e-4));

        //! A repeated eigenvalue with a full set of eigenvectors is already modal.
        constexpr StateSpace<2, 1, 1> repeated{.A = Matrix<2, 2>{{-2.0, 0.0}, {0.0, -2.0}}, .B = Matrix<2, 1>{{1.0}, {1.0}}, .C = Matrix<1, 2>{{1.0, 0.0}}};
        constexpr auto                rep = modal_form(repeated);
        static_assert(rep.success && rep.sys.A(0, 1) == 0.0 && rep.sys.A(1, 0) == 0.0);

        //! A double integrator is a Jordan block: no modal form.
        constexpr StateSpace<2, 1, 1> jordan{.A = Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}}, .B = Matrix<2, 1>{{0.0}, {1.0}}, .C = Matrix<1, 2>{{1.0, 0.0}}};
        constexpr auto                defective = modal_form(jordan);
        static_assert(!defective.success && defective.status == design::DesignStatus::singular);

        //! One state: the model itself.
        constexpr StateSpace<1, 1, 1> one{.A = Matrix<1, 1>{{-3.0}}, .B = Matrix<1, 1>{{2.0}}, .C = Matrix<1, 1>{{1.0}}};
        static_assert(modal_form(one).success && modal_form(one).sys.A(0, 0) == -3.0);
    }

    TEST_CASE("matlab::canon forwards to modal_form") {
        constexpr auto sys = make_modal_plant();
        static_assert(matlab::canon(sys).sys.A == modal_form(sys).sys.A);
        static_assert(matlab::canon(sys).as<float>().success);
    }
}

namespace {
// (s + 0.5)(2s + 1) / ((s + 1)(s + 2)(s + 3)): three real poles, a full numerator.
constexpr TransferFunction<3, 4, double> kTf{.num = {2.0, 2.0, 0.5}, .den = {1.0, 6.0, 11.0, 6.0}};

// An invertible, non-orthogonal change of basis that scrambles any canonical structure.
constexpr Matrix<3, 3> kScramble{{2.0, 1.0, 0.0}, {0.0, 1.0, -1.0}, {1.0, 0.0, 3.0}}; // det = 5

// The TF's controllable canonical realization in scrambled coordinates.
constexpr StateSpace<3, 1, 1> make_scrambled() {
    return similarity_transform(*kTf.to_controllable_canonical_form(), kScramble).sys;
}

// x₃ is decoupled from the input: (A, b) has a 2-state controllable part.
constexpr StateSpace<3, 1, 1> make_uncontrollable() {
    const StateSpace<3, 1, 1> plain{
        .A = Matrix<3, 3>{{-1.0, 0.0, 0.0}, {1.0, -2.0, 0.0}, {0.0, 0.0, -3.0}},
        .B = Matrix<3, 1>{{1.0}, {0.0}, {0.0}},
        .C = Matrix<1, 3>{{1.0, 1.0, 1.0}},
        .D = Matrix<1, 1>{{0.5}},
    };
    return similarity_transform(plain, kScramble).sys;
}

// Its dual: x₃ invisible to the output, a 2-state observable part.
constexpr StateSpace<3, 1, 1> make_unobservable() {
    const auto c = make_uncontrollable();
    return StateSpace<3, 1, 1>{.A = c.A.transpose(), .B = c.C.transpose(), .C = c.B.transpose(), .D = c.D};
}

template<size_t R, size_t C>
void check_close(const Matrix<R, C>& a, const Matrix<R, C>& b, double tol) {
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            CHECK(a(i, j) == doctest::Approx(b(i, j)).epsilon(tol).scale(1.0));
        }
    }
}
} // namespace

TEST_SUITE("Similarity, canonical and staircase forms (ss2ss, canon companion, ctrbf, obsvf)") {
    TEST_CASE("similarity_transform: T·A·T⁻¹, T·B, C·T⁻¹ by solves; the I/O map is unchanged") {
        constexpr auto sys = make_nice();
        constexpr auto z = similarity_transform(sys, kScramble);
        static_assert(z.success && z.status == design::DesignStatus::ok);
        check_close(z.sys.A * kScramble, kScramble * sys.A, 1e-12); // Ā·T = T·A
        check_close(z.sys.B, kScramble * sys.B, 1e-15);
        check_close(z.sys.C * kScramble, sys.C, 1e-12); // C̄·T = C
        CHECK(z.sys.D == sys.D);
        for (size_t k = 0; k < 6; ++k) {
            CHECK(markov(z.sys, k) == doctest::Approx(markov(sys, k)).epsilon(1e-12));
        }
        const auto back = similarity_transform(z.sys, *mat::lu_solve(kScramble, Matrix<3, 3>::identity()));
        REQUIRE(back.success);
        check_close(back.sys.A, sys.A, 1e-12);

        constexpr Matrix<3, 3> rank_two{{1.0, 2.0, 3.0}, {2.0, 4.0, 6.0}, {0.0, 0.0, 1.0}};
        static_assert(similarity_transform(sys, rank_two).status == design::DesignStatus::singular);
        static_assert(matlab::ss2ss(sys, kScramble).sys.A == z.sys.A);
    }

    TEST_CASE("controllable / observable canonical forms match the TransferFunction realizations") {
        constexpr auto scrambled = make_scrambled();
        constexpr auto ccf = controllable_canonical_form(scrambled);
        static_assert(ccf.success);
        constexpr auto ccf_tf = *kTf.to_controllable_canonical_form();
        check_close(ccf.sys.A, ccf_tf.A, 1e-10);
        CHECK(ccf.sys.B == ccf_tf.B); // eₙ, exact by construction
        check_close(ccf.sys.C, ccf_tf.C, 1e-10);
        CHECK(ccf.sys.D == ccf_tf.D);
        check_close(ccf.transform * scrambled.A, ccf.sys.A * ccf.transform, 1e-10); // T·A = A_c·T

        constexpr auto ocf = observable_canonical_form(scrambled);
        static_assert(ocf.success);
        constexpr auto ocf_tf = *kTf.to_observable_canonical_form();
        check_close(ocf.sys.A, ocf_tf.A, 1e-10);
        check_close(ocf.sys.B, ocf_tf.B, 1e-10);
        CHECK(ocf.sys.C == ocf_tf.C); // eₙᵀ, exact by construction
        for (size_t k = 0; k < 6; ++k) {
            CHECK(markov(ocf.sys, k) == doctest::Approx(markov(scrambled, k)).epsilon(1e-10));
        }
    }

    TEST_CASE("companion_form is MATLAB's canon(sys, 'companion'): B = e₁, −α in A's last column") {
        constexpr auto scrambled = make_scrambled();
        constexpr auto comp = companion_form(scrambled);
        static_assert(comp.success);
        CHECK(comp.sys.B == Matrix<3, 1>{{1.0}, {0.0}, {0.0}});
        CHECK(comp.sys.A(1, 0) == 1.0);
        CHECK(comp.sys.A(2, 1) == 1.0);
        CHECK(comp.sys.A(0, 0) == 0.0);
        CHECK(comp.sys.A(2, 0) == 0.0);
        for (size_t i = 0; i < 3; ++i) {
            CHECK(comp.sys.A(i, 2) == doctest::Approx(-kTf.den[3 - i] / kTf.den[0]).epsilon(1e-10)); // sⁿ + Σ αᵢ sⁱ, αᵢ = den[3 − i]
        }
        for (size_t k = 0; k < 6; ++k) {
            CHECK(markov(comp.sys, k) == doctest::Approx(markov(scrambled, k)).epsilon(1e-10));
        }
        static_assert(matlab::canon(scrambled, "companion").sys.A == comp.sys.A);
        static_assert(matlab::canon(scrambled, "modal").sys.A == modal_form(scrambled).sys.A);
        static_assert(matlab::canon(scrambled, "jordan").status == design::DesignStatus::bad_weights);
        static_assert(matlab::canon(make_modal_plant(), "companion").status == design::DesignStatus::bad_weights); // MIMO
    }

    TEST_CASE("canonical and companion forms refuse an uncontrollable or unobservable SISO model") {
        static_assert(controllable_canonical_form(make_uncontrollable()).status == design::DesignStatus::singular);
        static_assert(companion_form(make_uncontrollable()).status == design::DesignStatus::singular);
        static_assert(observable_canonical_form(make_unobservable()).status == design::DesignStatus::singular);
        static_assert(controllable_canonical_form(make_scrambled().as<float>()).success);
    }

    TEST_CASE("controllability_staircase splits off the uncontrollable mode, in MATLAB's ctrbf layout") {
        constexpr auto sys = make_uncontrollable();
        constexpr auto st = controllability_staircase(sys);
        static_assert(st.success && st.rank == 2);
        //! [A_nc 0; A21 A_c], B̄ = [0; B_c]: the decoupling is exact.
        CHECK(st.sys.A(0, 1) == 0.0);
        CHECK(st.sys.A(0, 2) == 0.0);
        CHECK(st.sys.B(0, 0) == 0.0);
        CHECK(st.sys.A(0, 0) == doctest::Approx(-3.0).epsilon(1e-12)); // the uncontrollable mode
        check_close(st.transform * st.transform.transpose(), Matrix<3, 3>::identity(), 1e-14);
        for (size_t k = 0; k < 6; ++k) {
            CHECK(markov(st.sys, k) == doctest::Approx(markov(sys, k)).epsilon(1e-12));
        }
        CHECK(st.sys.D == sys.D);

        constexpr auto full = controllability_staircase(make_modal_plant()); // MIMO outputs, controllable
        static_assert(full.success && full.rank == 4);
        static_assert(matlab::ctrbf(sys).rank == 2);
    }

    TEST_CASE("observability_staircase splits off the unobservable mode, in MATLAB's obsvf layout") {
        constexpr auto sys = make_unobservable();
        constexpr auto st = observability_staircase(sys);
        static_assert(st.success && st.rank == 2);
        //! [A_no A12; 0 A_o], C̄ = [0 C_o].
        CHECK(st.sys.A(1, 0) == 0.0);
        CHECK(st.sys.A(2, 0) == 0.0);
        CHECK(st.sys.C(0, 0) == 0.0);
        CHECK(st.sys.A(0, 0) == doctest::Approx(-3.0).epsilon(1e-12)); // the unobservable mode
        check_close(st.transform * st.transform.transpose(), Matrix<3, 3>::identity(), 1e-14);
        for (size_t k = 0; k < 6; ++k) {
            CHECK(markov(st.sys, k) == doctest::Approx(markov(sys, k)).epsilon(1e-12));
        }
        static_assert(matlab::obsvf(sys).rank == 2);
        static_assert(observability_staircase(make_scrambled()).rank == 3);
    }

    TEST_CASE("canonical and companion forms of a stiff model are refused or keep its DC gain") {
        //! Modes from 1 to 1e6 (and to 1e8), coupled, in rotated coordinates: the forward Markov
        //! parameters see only the fastest mode, the inverse moments the slowest.
        const auto stiff = [](const damp::array<double, 5>& p) {
            Matrix<5, 5> J{};
            for (size_t i = 0; i < 5; ++i) {
                J(i, i) = -p[i];
                if (i + 1 < 5) {
                    J(i, i + 1) = 0.5 * p[i + 1];
                }
            }
            const ColVec<5>    v{1.0, 2.0, 3.0, 4.0, 5.0};
            const Matrix<5, 5> Q = Matrix<5, 5>::identity() - ((v * v.transpose()) * (2.0 / 55.0)); // a reflection
            return StateSpace<5, 1, 1>{.A = Q * J * Q, .B = {{1.0}, {2.0}, {3.0}, {1.0}, {2.0}}, .C = {{1.0, 0.75, 0.5, 0.25, 0.0}}};
        };
        const auto dc_gain = [](const StateSpace<5, 1, 1>& s) { return -(s.C * *mat::lu_solve(s.A, s.B))(0, 0); };
        const auto check = [&](const StateSpace<5, 1, 1>& sys, const StateTransformResult<5, 1, 1>& form) {
            if (form.success) {
                CHECK(dc_gain(form.sys) == doctest::Approx(dc_gain(sys)).epsilon(1e-7));
            } else {
                CHECK((form.status == design::DesignStatus::not_converged || form.status == design::DesignStatus::singular)); // refused: Markov certificate or W_c
            }
        };
        for (const auto& p : {damp::array<double, 5>{1.0, 30.0, 1e3, 3e4, 1e6}, damp::array<double, 5>{1.0, 1e2, 1e4, 1e6, 1e8}}) {
            const auto sys = stiff(p);
            check(sys, controllable_canonical_form(sys));
            check(sys, observable_canonical_form(sys));
            check(sys, companion_form(sys));
        }
        const auto mild = stiff({1.0, 2.0, 3.0, 4.0, 5.0}); // still certified
        CHECK(controllable_canonical_form(mild).success);
        CHECK(observable_canonical_form(mild).success);
        CHECK(companion_form(mild).success);
    }

    TEST_CASE("the staircase rank does not depend on the units of the input") {
        //! det[B AB] = 1e11 with B = 1e9: controllable, whatever the input's scale.
        for (const double b : {1e9, 1.0, 1e-9}) {
            CAPTURE(b);
            const StateSpace<2, 1, 1> sys{.A = {{-1.0, 0.0}, {1e-7, -2.0}}, .B = {{b}, {0.0}}, .C = {{0.0, 1.0}}};
            const auto                st = controllability_staircase(sys);
            CHECK(st.success);
            CHECK(st.rank == 2);
        }
    }
}
