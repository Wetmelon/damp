// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/minreal.hpp"
#include "damp/design/stability.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::design;

namespace {

/// Continuous DC gain G(0) = D − C A^{-1} B via solve.
template<size_t NX, size_t NU, size_t NY>
Matrix<NY, NU> cont_dc(const StateSpace<NX, NU, NY>& sys) {
    const auto X = mat::solve(-sys.A, sys.B);
    REQUIRE(X.has_value());
    return sys.C * (*X) + sys.D;
}

/// Max-abs entry-wise difference.
template<size_t R, size_t C>
double max_abs_diff(const Matrix<R, C>& A, const Matrix<R, C>& B) {
    double m = 0.0;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            m = damp::max(m, std::abs(A(i, j) - B(i, j)));
        }
    }
    return m;
}

/// Fully controllable + observable second-order plant.
StateSpace<2, 1, 1> coupled_minimal() {
    return StateSpace<2, 1, 1>{
        .A = Matrix<2, 2>{{-2.0, 1.0}, {0.0, -3.0}},
        .B = Matrix<2, 1>{{1.0}, {1.0}},
        .C = Matrix<1, 2>{{1.0, 0.0}},
        .D = Matrix<1, 1>{{0.0}},
        .Ts = 0.0
    };
}

} // namespace

TEST_SUITE("minreal") {

    TEST_CASE("minimal plant is unchanged (order preserved, I/O map)") {
        const auto sys = coupled_minimal();
        REQUIRE(stability::is_controllable(sys.A, sys.B));
        REQUIRE(stability::is_observable(sys.A, sys.C));

        const auto mr = minreal(sys);
        REQUIRE(mr.success);
        CHECK(mr.order == 2);
        CHECK(mr.eliminated == 0);

        // Already minimal ⇒ identity return (bit-exact matrices).
        CHECK(max_abs_diff(mr.system.A, sys.A) == doctest::Approx(0.0));
        CHECK(max_abs_diff(mr.system.B, sys.B) == doctest::Approx(0.0));
        CHECK(max_abs_diff(mr.system.C, sys.C) == doctest::Approx(0.0));

        const auto red = mr.extract<2>();
        REQUIRE(red.has_value());
        CHECK(max_abs_diff(red->A, sys.A) == doctest::Approx(0.0));

        // Descriptive alias agrees.
        const auto mr2 = minimal_realization(sys);
        CHECK(mr2.success);
        CHECK(mr2.order == mr.order);
    }

    TEST_CASE("uncontrollable mode is removed") {
        // A = diag(-1,-2), B drives only state 0; C sees both.
        // I/O map is 1/(s+1); state 1 never enters from u.
        const StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{-1.0, 0.0}, {0.0, -2.0}},
            .B = Matrix<2, 1>{{1.0}, {0.0}},
            .C = Matrix<1, 2>{{1.0, 1.0}},
            .D = Matrix<1, 1>{{0.0}},
        };
        REQUIRE_FALSE(stability::is_controllable(sys.A, sys.B));
        REQUIRE(stability::is_observable(sys.A, sys.C));

        const auto mr = minreal(sys);
        REQUIRE(mr.success);
        CHECK(mr.order == 1);
        CHECK(mr.eliminated == 1);

        const auto red = mr.extract<1>();
        REQUIRE(red.has_value());
        REQUIRE(stability::is_controllable(red->A, red->B));
        REQUIRE(stability::is_observable(red->A, red->C));

        // DC gain of original: C (-A)^{-1} B = [1 1] * diag(1, 1/2) * [1;0] = 1
        const auto dc_full = cont_dc(sys);
        const auto dc_red = cont_dc(*red);
        CHECK(dc_red(0, 0) == doctest::Approx(dc_full(0, 0)).epsilon(1e-10));
        CHECK(dc_red(0, 0) == doctest::Approx(1.0).epsilon(1e-10));

        // Pole of reduced system ≈ -1
        CHECK(red->A(0, 0) == doctest::Approx(-1.0).epsilon(1e-8));
    }

    TEST_CASE("unobservable mode is removed") {
        // A = diag(-1,-2), both states driven; C sees only state 0.
        const StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{-1.0, 0.0}, {0.0, -2.0}},
            .B = Matrix<2, 1>{{1.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>{{0.0}},
        };
        REQUIRE(stability::is_controllable(sys.A, sys.B));
        REQUIRE_FALSE(stability::is_observable(sys.A, sys.C));

        const auto mr = minreal(sys);
        REQUIRE(mr.success);
        CHECK(mr.order == 1);
        CHECK(mr.eliminated == 1);

        const auto red = mr.extract<1>();
        REQUIRE(red.has_value());
        REQUIRE(stability::is_controllable(red->A, red->B));
        REQUIRE(stability::is_observable(red->A, red->C));

        const auto dc_full = cont_dc(sys);
        const auto dc_red = cont_dc(*red);
        CHECK(dc_red(0, 0) == doctest::Approx(dc_full(0, 0)).epsilon(1e-10));
        CHECK(red->A(0, 0) == doctest::Approx(-1.0).epsilon(1e-8));
    }

    TEST_CASE("series interconnection cancellation (pole-zero cancel)") {
        // G1 = (s+2)/(s+1), G2 = 1/(s+2)  ⇒  G1*G2 = 1/(s+1)
        // Companion SS series has order 2 with a non-minimal cancelled mode.
        const TransferFunction<2, 2> tf1{{2.0, 1.0}, {1.0, 1.0}}; // (2 + s)/(1 + s)
        const TransferFunction<1, 2> tf2{{1.0}, {2.0, 1.0}};      // 1/(2 + s)
        const auto                   ss1 = tf1.to_state_space().value();
        const auto                   ss2 = tf2.to_state_space().value();
        const auto                   series_sys = *(ss1 * ss2); // order 2

        CHECK(series_sys.A.rows() == 2);

        const auto mr = minreal(series_sys);
        REQUIRE(mr.success);
        CHECK(mr.order == 1);
        CHECK(mr.eliminated == 1);

        const auto red = mr.extract<1>();
        REQUIRE(red.has_value());

        // Expected minimal TF: 1/(s+1) → DC = 1, pole = -1
        const auto dc_full = cont_dc(series_sys);
        const auto dc_red = cont_dc(*red);
        CHECK(dc_full(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(dc_red(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(red->A(0, 0) == doctest::Approx(-1.0).epsilon(1e-6));

        // Residue of 1/(s+1) ⇒ product of input/output maps has magnitude 1
        CHECK(std::abs(red->B(0, 0) * red->C(0, 0)) == doctest::Approx(1.0).epsilon(1e-6));
    }

    TEST_CASE("both uncontrollable and unobservable modes") {
        // State 0: CO (pole -1). State 1: unctrl. State 2: unobs.
        const StateSpace<3, 1, 1> sys{
            .A = Matrix<3, 3>{{-1.0, 0.0, 0.0}, {0.0, -2.0, 0.0}, {0.0, 0.0, -3.0}},
            .B = Matrix<3, 1>{{1.0}, {0.0}, {1.0}},
            .C = Matrix<1, 3>{{1.0, 1.0, 0.0}},
            .D = Matrix<1, 1>{{0.0}},
        };

        const auto mr = minreal(sys);
        REQUIRE(mr.success);
        CHECK(mr.order == 1);
        CHECK(mr.eliminated == 2);

        const auto red = mr.extract<1>();
        REQUIRE(red.has_value());
        CHECK(cont_dc(*red)(0, 0) == doctest::Approx(cont_dc(sys)(0, 0)).epsilon(1e-9));
        CHECK(red->A(0, 0) == doctest::Approx(-1.0).epsilon(1e-6));
    }

    TEST_CASE("MIMO: remove unobservable state") {
        // 2-state, 2-in, 1-out: second state driven but not seen
        const StateSpace<2, 2, 1> sys{
            .A = Matrix<2, 2>{{-1.0, 0.0}, {0.0, -4.0}},
            .B = Matrix<2, 2>{{1.0, 0.0}, {0.0, 1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 2>{{0.0, 0.0}},
        };
        REQUIRE(stability::is_controllable(sys.A, sys.B));
        REQUIRE_FALSE(stability::is_observable(sys.A, sys.C));

        const auto mr = minreal(sys);
        REQUIRE(mr.success);
        CHECK(mr.order == 1);

        const auto red = mr.extract<1>();
        REQUIRE(red.has_value());
        CHECK(max_abs_diff(cont_dc(sys), cont_dc(*red)) < 1e-9);
    }

    TEST_CASE("discrete-time non-minimal plant") {
        const StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{0.5, 0.0}, {0.0, 0.2}},
            .B = Matrix<2, 1>{{1.0}, {0.0}},
            .C = Matrix<1, 2>{{1.0, 1.0}},
            .D = Matrix<1, 1>{{0.0}},
            .Ts = 0.01,
        };
        const auto mr = minreal(sys);
        REQUIRE(mr.success);
        CHECK(mr.order == 1);
        CHECK(mr.system.Ts == doctest::Approx(0.01));
        const auto red = mr.extract<1>();
        REQUIRE(red.has_value());
        CHECK(red->Ts == doctest::Approx(0.01));
        CHECK(red->A(0, 0) == doctest::Approx(0.5).epsilon(1e-8));
    }

    TEST_CASE("extract rejects wrong order") {
        const auto sys = coupled_minimal();
        const auto mr = minreal(sys);
        REQUIRE(mr.success);
        CHECK_FALSE(mr.extract<1>().has_value());
        CHECK(mr.extract<2>().has_value());
    }

    TEST_CASE("as<float> converts result") {
        const auto sys = coupled_minimal();
        const auto mr = minreal(sys);
        REQUIRE(mr.success);
        const auto mrf = mr.as<float>();
        CHECK(mrf.success);
        CHECK(mrf.order == 2);
        CHECK(static_cast<double>(mrf.system.A(0, 0)) == doctest::Approx(sys.A(0, 0)));
    }

    TEST_CASE("constexpr minreal on uncontrollable pair") {
        constexpr StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{-1.0, 0.0}, {0.0, -2.0}},
            .B = Matrix<2, 1>{{1.0}, {0.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
        };
        constexpr auto mr = minreal(sys);
        static_assert(mr.success);
        static_assert(mr.order == 1);
        static_assert(mr.eliminated == 1);
    }
}
