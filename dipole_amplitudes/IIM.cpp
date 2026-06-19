#include <iostream>
#include <fstream>
#include <chrono>

#include "../other/utils.hpp"
#include "../other/ctes.hpp"
#include "../calculations/wavefunctions.hpp"
#include "../other/integration.hpp"
#include "../other/correcs.hpp"
#include "IIM.hpp"


param_IIM IIM_S("IIM_S", 26.25 * mb_to_gev2, 0.6194, 0.2131e-4, 0.2545, 9.9, 0.7);
param_IIM IIM_RS("IIM_RS", 21.85 * mb_to_gev2, 0.762, 6.226e-5, 0.2319, 9.9, 0.7);


namespace IIM{

double N_p(double r,  double x, const param_IIM& param)
{

    double Qs  = std::pow(param.x_0/x, param.lambda/2.0);
    double rQs = r * Qs;

    rQs = std::max(rQs, 1e-12);      // proteção numérica

    if (rQs < 2.0) {

        double Lx = std::log(1.0/x);

        double exp =
            2.0 * (param.gamma_s +
            std::log(2.0/rQs) /
            (param.kappa * param.lambda * Lx));

        return param.N_0 * std::pow(rQs/2.0, exp)*low_x_factor(x, 6.0);

    } else {

        double a = -param.N_0*param.N_0*param.gamma_s*param.gamma_s /
                   ((1-param.N_0)*(1-param.N_0)*std::log(1.0-param.N_0));

        double b = 0.5 * std::pow(1.0-param.N_0,
                   -(1.0-param.N_0)/(param.N_0*param.gamma_s));

        double ln = std::log(b * rQs);

        return 1.0 - std::exp(-a * ln * ln)*low_x_factor(x, 6.0);
    }
}

double sigma_qq(double r, double x, const param_IIM& param){return param.sigma0 * N_p(r, x, param);}

double amplitude_p(double x, double Q2, const Meson& M, param_IIM param, bool fc)
{
    auto amp_r = [x, Q2, &M, param, fc](double r) {
        double Ov = sigma_qq(r, x, param) * overlap_r(r, Q2, M, fc);// * sqrt_fc; // r de d²r = 2π r dr
        return r * Ov;
    };
    double amp = 0.5 * integrate_simpson(amp_r, 1e-4, 10.0, 300);
    return amp;
}














}