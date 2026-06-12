#include <chrono>
#include <omp.h>
#include <fstream>
#include <sstream>

#include "../other/integration.hpp"
#include "../libraries/mantysaari/dipoleamplitude.hpp"
#include "../libraries/mantysaari/dglap_cpp/AlphaStrong.h"
#include "../libraries/mantysaari/dglap_cpp/EvolutionLO_nocoupling.h"
#include "../other/utils.hpp"
#include "../other/ctes.hpp"
#include "../calculations/wavefunctions.hpp"
#include "../other/correcs.hpp"
#include <boost/math/special_functions/bessel.hpp>

#include "ipsat.hpp"

using namespace MZ_ipsat;

namespace IPSAT{

// ----------------- função de calcular a seção de choque de dipolo -----------------
double sigma_qq_p(double r, double x, double Delta, DipoleAmplitude& dipole)
{
    auto N_b = [&](double b) {
        return 2.0 * M_PI * b * dipole.N(r, x, b) * boost::math::cyl_bessel_j(0, Delta * b);
    };
    double bmax = 10.0;
    double Nval = integrate_simpson(N_b, 0.0, bmax, 200);
    return 2 * Nval;
}





//----------- amplitude -----------
double amplitude_p(double x, double Delta, double Q2, const Meson& M,
                 int Nr, int Nz,
                 double rmin, double rmax)
{
    DipoleAmplitude dipole(MZ_IPSAT);
    dipole.EnableLookupTable();
    auto fr = [&](double r) {
        double Ov = overlap_r(r, Q2, M);
        double sigma_qq = sigma_qq_p(r, x, Delta, dipole);
        return 0.5 * r * Ov * sigma_qq; // r de d²r = 2π r dr
    };
    double Ir = integrate_simpson(fr, rmin, rmax, Nr);
    return Ir;
}





}


