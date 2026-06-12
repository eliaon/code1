#ifndef SIGMA_HPP
#define SIGMA_HPP

#include "../other/utils.hpp"
#include "nuclear.hpp"

namespace gamma_p {

    double sigma(double W, double Q2, const Meson& M, std::string modelo);

    double dsigma_dt(double t, double W, double Q2, const Meson& M, std::string modelo);

    double sigma_integrado(double W, double Q2, const Meson& M, std::string modelo);
}

namespace gamma_A{

    double N_qq_A(double r, double x, double b, double Delta ,std::string modelo, Nucleus& Nucleo, const TA_Table& table);

    double sigma(double W, double Q2, const Meson&M, std::string modelo, const TA_Table& table, Nucleus& Nucleo);

    double dN_domega(double omega, double sqrt_s, Nucleus&Nucleo, double bmin);

    double d_sigma_dy_AA(double y, double sqrt_s, double Q2, const Meson& M, std::string model,
                    TA_Table& table, Nucleus& Nucleo);

    double d_sigma_dy_AB(double y, double sqrt_s, double Q2, const Meson& M, std::string model,
                    TA_Table& tableA, TA_Table& tableB,
                    Nucleus& NucleoA, Nucleus& NucleoB);
}








#endif // SIGMA_HPP