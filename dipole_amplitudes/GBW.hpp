#include <cmath>
#include "../other/utils.hpp"
#include "../calculations/nuclear.hpp"

#ifndef GBW_H
#define GBW_H

// ---------------- Classe de parâmetros GBW --------------
class parametros_GBW {
    public:
        std::string name; // nome do conjunto de parâmetros
        double sigma0; // GeV^-2
        double x0;
        double lambda;

        parametros_GBW(std::string name_, double sigma0_, double x0_, double lambda_)
            : name(name_), sigma0(sigma0_), x0(x0_), lambda(lambda_) {}
};

extern parametros_GBW gbw;

extern parametros_GBW gbw_10; 


namespace GBW {

double QS2( double x, parametros_GBW params);

double N_p(double r,  double x, parametros_GBW params);

double sigma_qq_p(double r, double x, parametros_GBW params);


double amplitude_p(double x, double Q2, const Meson& M,
                const parametros_GBW& gbw, bool fc = false);


double amplitude_model(double x,  double B, double omega,
                        const Meson& M, double Q2 = 0);

double sigma_model(double x, double B, double omega, const Meson& M, double Q2 = 0);

}
#endif// GBW_H