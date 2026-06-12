#ifndef WAVEFUNCTIONS_H
#define WAVEFUNCTIONS_H
#include "../other/utils.hpp"
#include "../other/ctes.hpp"



double laplacian_phi_L(double r, double z, const Meson& M);

double psi_Vpsi_T(double z, double r, double Q2, const Meson& M);

double psi_Vpsi_L(double z, double r, double Q2, const Meson& M);

double phi_T(double r, double z, const Meson& M);

double phi_L(double r, double z, const Meson& M);

double overlap_r(double r, double Q2, const Meson& M);

void overlap_csv(void);

void overlap_csv_fc(void);

#endif// WAVEFUNCTIONS_H