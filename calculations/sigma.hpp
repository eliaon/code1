#ifndef SIGMA_HPP
#define SIGMA_HPP

#include "../other/utils.h"
#include "nuclear.hpp"

namespace gamma_p {

    double sigma(double W, double Q2, const Meson& M, std::string modelo);

    double dsigma_dt(double t, double W, double Q2, const Meson& M, std::string modelo);

    double sigma_integrado(double W, double Q2, const Meson& M, std::string modelo);
}

namespace gamma_A{

    double sigma_qq_A(double r, double x, double b, double Delta ,std::string modelo, const TA_Table& table);

    double sigma(double W, double Q2, const Meson&M, std::string modelo, const TA_Table& table);

    double dN_domega(double omega, double sqrt_s, nucleous&Nucleo);

    double d_sigma_dy_AA(double y, double sqrt_s, double Q2, const Meson& M, std::string model,
                    TA_Table& table, nucleous& Nucleo);

    double d_sigma_dy_AB(double y, double sqrt_s, double Q2, const Meson& M, std::string model,
                    TA_Table& tableA, TA_Table& tableB,
                    nucleous& NucleoA, nucleous& NucleoB);
}








#endif // SIGMA_HPP