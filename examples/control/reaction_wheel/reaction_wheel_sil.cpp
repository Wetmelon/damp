// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file reaction_wheel_sil.cpp
 * @brief Reaction-wheel gallery — print gains + one tick each (host)
 *
 * Not closed-loop plant SIL. Flashable cascade: reaction_wheel_sketch.cpp.
 */

#include "damp/backend.hpp"
#include "damp/controllers/adrc.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/ofmpc.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/controllers/smc.hpp"
#include "damp/controllers/stsmc.hpp"
#include "damp/matrix/colvec.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "reaction_wheel_controller.hpp"

using namespace damp;
using namespace damp::examples_reaction_wheel;

static constinit PIController          rate_pi = rate_pi_design.discretize(Ts).as<float>();
static constinit PIController          angle_pi = angle_pi_design.discretize(Ts).as<float>();
static constinit Law<adrc_res, float>  adrc;
static constinit SMCController         smc{smc_res.as<float>()};
static constinit STSMCController       stsmc{stsmc_res.as<float>()};
static constinit Law<place_res, float> place_ctrl;
// Each design runs as the law named by its result, on its matrices as constants: LQR on the measured
// state (design::place, design::lqr), LQI with integrators, LQG with the Kalman estimator, LQGI with both,
// and OFMPC with its filter on constants.
static constinit Law<lqr_res, float>  lqr;
static constinit Law<lqi_res, float>  lqi;
static constinit Law<lqg_res, float>  lqg;
static constinit Law<lqgi_res, float> lqgi;
static constinit Law<mpc_art, float>  of_mpc;

int main() {
    fmt::print("===== Reaction-wheel roll — controller synthesis =====\n\n");
    fmt::print("Plant:  θ̈ = ({:.3g})² θ + ({:.3g}) u ,  Ts = {:.4g} s\n", w0, b_u, Ts);
    fmt::print("Cascade:  w_angle={:.3g} rad/s → w_rate={:.3g} rad/s , |u|≤{}\n\n", w_angle, w_rate, u_max);

    fmt::print(
        "cascade PI:  angle Kp={:.4g} Ki={:.4g}  |  rate Kp={:.4g} Ki={:.4g}\n",
        angle_pi_design.Kp,
        angle_pi_design.Ki,
        rate_pi_design.Kp,
        rate_pi_design.Ki
    );
    fmt::print("ADRC (NX=2):  wc={:.3g}  wo={:.3g}  b0={:.3g}\n", wc_adrc, wo_adrc, b_u);
    fmt::print("SMC:  λ={:.3g}  k={:.3g}  b0={:.3g}\n", smc_res.lambda, smc_res.k, smc_res.b0);
    fmt::print("STSMC:  k1={:.4g}  k2={:.4g}  λ={:.3g}\n", stsmc_res.k1, stsmc_res.k2, stsmc_res.lambda);

    fmt::print(
        "place_sampled (s=[{:.3g},{:.3g}], Ts={:.3g}) K = [{:.4g}, {:.4g}]\n",
        place_poles_s[0].real(),
        place_poles_s[1].real(),
        Ts,
        place_ctrl.K(0, 0),
        place_ctrl.K(0, 1)
    );
    fmt::print("LQR   K = [{:.4g}, {:.4g}]\n", lqr.K(0, 0), lqr.K(0, 1));
    fmt::print("LQI   K = [{:.4g}, {:.4g} | Ki={:.4g}]\n", lqi.K(0, 0), lqi.K(0, 1), lqi.K(0, 2));
    fmt::print("LQG   K = [{:.4g}, {:.4g}]\n", lqg.K(0, 0), lqg.K(0, 1));
    fmt::print("LQGI  Kx = [{:.4g}, {:.4g}]  Ki={:.4g}\n", lqgi.K(0, 0), lqgi.K(0, 1), lqgi.K(0, 2));
    fmt::print("OF-MPC:  success={}  NP=10  NC=4  |u|≤{}\n", mpc_art.success, u_max);

    const ColVec<NX, float> x0{0.05f, 0.0f};
    const float             theta0 = x0(0);
    const float             omega0 = x0(1);
    const float             r0 = 0.0f;

    // Every law's command goes through the actuator limit; the laws that track the applied input
    // (ADRC's observer, the state-feedback integrators and estimators) are told what went out.
    const auto clamp_u = [](float u) { return damp::clamp(u, -u_max, u_max); };
    const auto applied = [&](auto& law, const ColVec<NU, float>& u_req) {
        const float u = clamp_u(u_req(0));
        law.back_calculate(u);
        return u;
    };

    const float u_cascade = cascade_period(angle_pi, rate_pi, r0, theta0, omega0); // clamps internally
    const float u_adrc = clamp_u(adrc.control(r0, theta0));
    adrc.back_calculate(u_adrc);
    const float u_smc = clamp_u(smc.control(r0, theta0));
    const float u_st = clamp_u(stsmc.control(r0, theta0));
    const float u_place = applied(place_ctrl, place_ctrl.step(x0));
    const float u_lqr = applied(lqr, lqr.step(x0));
    const float u_lqi = applied(lqi, lqi.step(r0, x0));
    const float u_lqg = applied(lqg, lqg.step(theta0));
    const float u_lqgi = applied(lqgi, lqgi.step(r0, theta0));
    const float u_mpc = applied(of_mpc, of_mpc.step(r0, theta0));

    fmt::print("\nOne tick from x=[{:.3g}, {:.3g}], r=0, |u| ≤ {} (I/O differs — not a bake-off):\n", theta0, omega0, u_max);
    fmt::print("  cascade={:.4g}  adrc={:.4g}  smc={:.4g}  stsmc={:.4g}\n", u_cascade, u_adrc, u_smc, u_st);
    fmt::print(
        "  place={:.4g}  lqr={:.4g}  lqi={:.4g}  lqg={:.4g}  lqgi={:.4g}  of_mpc={:.4g}\n",
        u_place,
        u_lqr,
        u_lqi,
        u_lqg,
        u_lqgi,
        u_mpc
    );

    return 0;
}
