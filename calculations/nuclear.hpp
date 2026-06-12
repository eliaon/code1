#ifndef NUCLEAR_HPP
#define NUCLEAR_HPP

#include <vector>
#include <iostream>


class Nucleus {
    public:
        std::string name; // nome do núcleo (ex: Pb)
        double Z; // número atômico
        double R; // raio nuclear GeV^-1
        double a; // diffusividade GeV^-1
        double rho0; // densidade central Gev^-3

        Nucleus(std::string name_, double Z_, double R_, double a_, double rho0_)
            : name(name_), Z(Z_), R(R_), a(a_), rho0(rho0_) {}
};

extern Nucleus Pb208; // parâmetros típicos para Pb-208: R=6.62 fm, a=0.546 fm, rho0=0.1603 fm^-3

extern Nucleus proton; // parâmetros típicos para próton: R=0.84 fm, a=0.0 (sem diffusividade), rho0=0.17 fm^-3 (ajustado para A=1)

struct TA_Table {
    std::vector<double> b_vals; // valores de b
    std::vector<double> TA_vals; // valores de T_A(b)
};

double rho_WS(double r, const Nucleus& N);

double interpolate_TA(double b, const TA_Table& table);

TA_Table precompute_TA(int Nb, double bmax, const Nucleus& N);

double integral_rho(const Nucleus& N);

double compute_rho0(int A, const Nucleus& N);

TA_Table get_TA_table(const Nucleus& nuc,
                      int Nb = 200,
                      double bmax = 50.0);








#endif // NUCLEAR_HPP