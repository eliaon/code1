#include "../other/utils.h"
#include "../other/ctes.h"
#include "../calculations/wavefunctions.h"
#include "../calculations/nuclear.hpp"
#include "../other/integration.hpp"
#include "../other/correcs.h"
#include "../other/plot.h"
#include "../calculations/nuclear.hpp"

#include <sstream>
#include <fstream>
#include <cmath>
#include <iostream>
#include <format>
#include <iomanip>
#include <string>
#include <chrono>
#include <boost/math/special_functions/bessel.hpp>


parametros_GBW gbw("GBW_old", 23*mb_to_gev2, 3e-4, 0.29);


// ----------- parametros GBW NEW 2018 Saturation model of DIS: an update https://doi.org/10.1007/JHEP03(2018)102
// sequência: sigma_0 , x_0, lambda
parametros_GBW gbw_5("GBW(new)_5", 28.18*mb_to_gev2, 0.31e-4, 0.237);

parametros_GBW gbw_10("GBW(new)_10", 27.43*mb_to_gev2, 0.40e-4, 0.248);

parametros_GBW gbw_20("GBW(new)_20", 26.60*mb_to_gev2, 0.53e-4, 0.259);

parametros_GBW gbw_50("GBW(new)_50", 25.21*mb_to_gev2, 0.80e-4, 0.281);





namespace GBW {
double QS2( double x, parametros_GBW params)
{
    return pow(params.x0 / x, params.lambda);
}   

double N_p(double r,  double x, parametros_GBW params)
{
    double Qs2 = QS2(x, params);
    double arg = (r * r) * Qs2 / 4.0;
    return (1.0 - exp(-arg))*std::pow(1.0-x, 5.26);
}

double sigma_qq_p(double r, double x, parametros_GBW params)
{
    return params.sigma0 * N_p(r, x, params);
}



double amplitude_p(double x, double Q2, const Meson& M,
                const parametros_GBW& gbw)
{
    auto amp_r = [x, Q2, M, gbw](double r) {
        double Ov = overlap_r(r, Q2, M);
        double sigma_dip = sigma_qq_p(r, x, gbw);
        //double sqrt_fc = std::sqrt(f_c(r, -0.979599, 0.403569));
        return 2.0 * M_PI * r * Ov * sigma_dip; //* sqrt_fc; // r de d²r = 2π r dr
    };
    double amp = integrate_simpson(amp_r, rmin, rmax, Nr);
    //cout << "x: " << x << ", Q2: " << Q2 << ", amp: " << amp << "\n";
    return amp / (4.0 * M_PI); // normaliza pela área da seção de choque total
}
// --------------- modelos para parametrizar B e omega ------------

double amplitude_model(double x,  double B, double omega,
                        const Meson& M, double Q2)
{
    auto amp_r = [x, Q2, M, B, omega](double r) {
        double Ov = overlap_r(r, Q2, M);
        double sigma_qq = GBW::sigma_qq_p(r, x, gbw);
        double sqrt_fc = std::sqrt(f_c(r, B, omega));
        return  2*M_PI *r * Ov * sigma_qq * sqrt_fc; // r de d²r = 2π r dr
    };
    double amp = integrate_simpson(amp_r, rmin, rmax, Nr);
    return amp / (4.0 * M_PI); // normaliza pela área da seção de choque total
}

double sigma_model(double x, double B, double omega, const Meson& M, double Q2)
    {
     parametros_GBW params = gbw; // ou escolha outro conjunto de parâmetros se desejar
     double amp = GBW::amplitude_model(x, B, omega, M, Q2);
     double B_slope_val = B_slope(x, Q2, M);
     double lambda_e = calculate_lambda(x, 0.0, Q2, M, "GBW");
     double RG_val = RG(x, Q2, lambda_e, M);
     double beta_val = beta(x, Q2, lambda_e, M);
    
     double correction_factor = RG_val * RG_val * (1.0 + beta_val * beta_val);
     double sigma_gev = correction_factor * (amp * amp) / (16.0 * M_PI * B_slope_val);
     double sigma_nb = sigma_gev * GeV2_to_nb;
     return sigma_nb;
    }



//------------------ caso nuclear -----------------------

double N_A(double r, double x, double b, parametros_GBW params, const TA_Table& table)
{
    double sigmaqq_p = params.sigma0 * N_p(r, x, params);
    double TA = interpolate_TA(b, table); 
    double arg = 0.5 * TA * sigmaqq_p;
    return (1.0 - exp(-arg));   
}

double sigma_qq_A(double r, double x, double b, parametros_GBW params, const TA_Table& table)
{
    auto b_integrand = [r, x, params, &table](double b) {
        return 2.0 * M_PI * b * N_A(r, x, b, params, table);
    };
    double bmax = table.b_vals.back(); // limite de integração em b
    int nb = table.b_vals.size(); // número de pontos para integração em b
    return integrate_simpson(b_integrand, 0.0, bmax, nb);
}



}