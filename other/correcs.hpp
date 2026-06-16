#ifndef CORRECS_H
#define CORRECS_H

#include "../dipole_amplitudes/GBW.hpp"
#include "../libraries/mantysaari/dipoleamplitude.hpp"
#include "utils.hpp"

#include <string>

// Correções de skewedness (Shuvaev et al.) e parte real (relações de dispersão).
// Referências: Gonçalves, Machado, Meneses EPJC 68 (2010) 133;
//              Kowalski & Teaney PRD 68 (2003); Mäntysaari & Lappi PRD 60 (1999);
//              Mäntysaari & Zurita PRD 98 (2018) Sec. V, Ref. [47].

struct SkewCorrection {
    double lambda_e = 0.0;
    double Rg       = 1.0;
    double beta     = 0.0;
    double factor   = 1.0; // Rg^2 * (1 + beta^2)
};

double calculate_lambda(double x, double Delta, double Q2, const Meson& M,
                        std::string dipolemodel);

double RG(double x, double Q2, double lambda_e, const Meson& M);

double beta(double x, double Q2, double lambda_e, const Meson& M);

SkewCorrection compute_skew_correction(double x, double Delta, double Q2,
                                       const Meson& M, const std::string& dipolemodel);

void debug_correc(std::string model);

double f_c(double r, double B, double omega, double R = 6.8);

double xg_p(double x, double Q);
double xg_Pb(double x, double Q);
double S_Pb(double x, double Q);

#endif // CORRECS_H
