#ifndef IPSAT_H
#define IPSAT_H

#include "../other/utils.hpp"
#include "../other/ctes.hpp"
#include "../calculations/wavefunctions.hpp"
#include "../other/correcs.hpp"
#include "../other/integration.hpp"
#include "../libraries/mantysaari/dglap_cpp/EvolutionLO_nocoupling.h"

#include <cmath>

namespace MZ_ipsat {
    class DipoleAmplitude;
}

namespace IPSAT{

double sigma_qq_p(double r, double x, double Delta, MZ_ipsat::DipoleAmplitude& dipole);

double amplitude_p(double x, double Delta, double Q2, const Meson& M, bool fc = false);
};

#endif //IPSAT_H