#ifndef BCGC_H
#define BCGC_H

#include "../other/utils.hpp"

namespace bCGC {

class parametros_bCGC {
public:
    std::string name;
    double N0;
    double gamma_s;
    double lambda;
    double x0;
    double Bcgc;
    double kappa;

    parametros_bCGC(std::string name_, double N0_, double gamma_s_,
                    double lambda_, double x0_, double Bcgc_, double kappa_)
        : name(name_), N0(N0_), gamma_s(gamma_s_), lambda(lambda_),
          x0(x0_), Bcgc(Bcgc_), kappa(kappa_) {}
};

extern const parametros_bCGC bg_fitted;

double Qs(double x, double b, const parametros_bCGC& params = bg_fitted);

double N_p(double r, double x, double b,
           const parametros_bCGC& params = bg_fitted);

double sigma_qq_p(double r, double x, double Delta,
                  const parametros_bCGC& params = bg_fitted);

double amplitude_p(double x, double Delta, double Q2, const Meson& M,
                   const parametros_bCGC& params = bg_fitted);

} // namespace bCGC

#endif // BCGC_H
