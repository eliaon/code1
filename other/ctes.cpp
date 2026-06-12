#include "ctes.hpp"
#include <cmath>

// --------- PARÂMETROS DE INTEGRAÇÃO ---------------

const double rmin = 1e-4;
const double rmax = 20;
const double Nr = 200;

const double tmin = 0.0;
const double tmax = 2.5;
const double Nt = 100;

const double zmin = 1e-6;
const double zmax = 1.0 - 1e-6;
const double Nz = 200;

const double alfem = 1.0 / 137.0;
const double Nc    = 3.0;
const double Nf = 4.0;
const double lambda = 0.29;
const double gamma_s = 0.46;
const double CFAC = 5.07;
const double sigma0 = 23.0 / 0.3894;
const double gev2_to_mb = 0.3894;
const double mb_to_gev2 = 1/gev2_to_mb;

const double m_p = 0.938; // massa do próton em GeV
// =========================================================
// Cargas efetivas dos quarks
// =========================================================

// charm
const double qJ = 2.0/3.0;

// strange
const double qS = 1.0/3.0;

// bottom
const double qB = 1.0/3.0;

// light quarks (rho)
// e_f = 1/sqrt(2)
const double qU2 = 1.0/sqrt(2.0);

// =========================================================
// Massas dos mésons [GeV]
// =========================================================

const double massa_psi     = 3.097;
const double massa_phi     = 1.019;
const double massa_rho     = 0.776;
const double massa_upsilon = 9.460;

// =========================================================
// Massas dos quarks [GeV]
// =========================================================

const double muq_GBW = 0.14;
const double ms_GBW  = 0.14;

const double mc_GBW  = 1.3528;

const double mb_GBW  = 4.5;


// Mäntysaari & Zurita 2018 (Phys. Rev. D 98, 036002), Table I (IPsat, HERA I)
const double mc_ipsat = 1.3528;
const double ms_ipsat = 0.03;

// ------------ PARAMETROS DGLAP ------------------

const double sigma0_dglap = 22.60*mb_to_gev2;

const double R2_psi = 1.5070 * 1.5070;
const double R2_phi = 3.3922 * 3.3922;

const double GeV2_to_nb = 3.89379e5;