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

#include "damp/controllers/adrc.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/offset_free_mpc.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/controllers/smc.hpp"
#include "damp/controllers/stsmc.hpp"
#include "damp/matrix/colvec.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "reaction_wheel_controller.hpp"

using namespace damp;
using namespace damp::examples_reaction_wheel;

static constinit PIController<float>                     rate_pi = rate_pi_design.discretize(Ts).as<float>();
static constinit PIController<float>                     angle_pi = angle_pi_design.discretize(Ts).as<float>();
static constinit ADRCController<2, float>                adrc{design::adrc<2>(wc_adrc, wo_adrc, b_u).as<float>(), static_cast<float>(Ts)};
static constinit SMCController<float>                    smc{smc_res.as<float>(), static_cast<float>(Ts), 0.05f};
static constinit STSMCController<float>                  stsmc{stsmc_res.as<float>(), static_cast<float>(Ts)};
static constinit StateFeedback<NX, NU, float>            place_ctrl{K_place->as<float>()};
static constinit LQR<NX, NU, float>                      lqr{lqr_res.as<float>()};
static constinit LQI<NX, NU, NY, float>                  lqi{lqi_res.as<float>()};
static constinit LQG<NX, NU, NY, float, NX, NY>          lqg{lqg_res.as<float>()};
static constinit LQGI<NX, NU, NY, float, NX, NY>         lqgi{lqgi_res.as<float>()};
static constinit OffsetFreeMPC<NX, NU, NY, 10, 4, float> of_mpc{mpc_art.as<float>()};

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
        place_poles_s[0],
        place_poles_s[1],
        Ts,
        place_ctrl.K(0, 0),
        place_ctrl.K(0, 1)
    );
    fmt::print("LQR   K = [{:.4g}, {:.4g}]\n", lqr.K(0, 0), lqr.K(0, 1));
    fmt::print("LQI   K = [{:.4g}, {:.4g} | Ki={:.4g}]\n", lqi.K(0, 0), lqi.K(0, 1), lqi.K(0, 2));
    fmt::print("LQG   K = [{:.4g}, {:.4g}]\n", lqg.lqr.K(0, 0), lqg.lqr.K(0, 1));
    fmt::print("LQGI  Kx = [{:.4g}, {:.4g}]  Ki={:.4g}\n", lqgi.lqi.K(0, 0), lqgi.lqi.K(0, 1), lqgi.lqi.K(0, 2));
    fmt::print("OF-MPC:  success={}  NP=10  NC=4  |u|≤{}\n", mpc_art.success, u_max);

    const ColVec<NX, float> x0{0.05f, 0.0f};
    const float             theta0 = x0(0);
    const float             omega0 = x0(1);
    const float             r0 = 0.0f;

    const float u_cascade = cascade_period(angle_pi, rate_pi, r0, theta0, omega0);
    const float u_adrc = adrc.control(r0, theta0);
    const float u_smc = smc.control(r0, theta0);
    const float u_st = stsmc.control(r0, theta0);
    const float u_place = place_ctrl.control(x0)(0);
    const float u_lqr = lqr.control(x0)(0);
    const float u_lqi = lqi.control(ColVec<NY, float>{r0}, ColVec<NY, float>{theta0}, x0)(0);
    const float u_lqg = lqg.step(ColVec<NY, float>{theta0})(0);
    const float u_lqgi = lqgi.step(ColVec<NY, float>{r0}, ColVec<NY, float>{theta0})(0);
    const float u_mpc = of_mpc.control(ColVec<NY, float>{r0}, ColVec<NY, float>{theta0})(0);

    fmt::print("\nOne tick from x=[{:.3g}, {:.3g}], r=0 (I/O differs — not a bake-off):\n", theta0, omega0);
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
