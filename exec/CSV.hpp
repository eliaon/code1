#ifndef CSV_HPP
#define CSV_HPP

#include <iostream>
#include "../calculations/nuclear.hpp"

namespace CSV {

    void N_p(std::string model);

    void sigma_gamma_p(std::string model, bool fc = false);

    void sigma_gamma_A(std::string model, Nucleus&Nucleo, bool fc = false);

    void rapidez_AA(double sqrt_s, std::string model, Nucleus &Nucleo1, Nucleus& Nucleo2, bool fc = false);

}
void Shadowing_factor_csv(double Q);

namespace Compare_all {
    const std::vector<std::string> modelos;

    void N(double x);

    void sigma_gamma_p(bool fc = false);

    void rapidez_AA(double sqrt_s, Nucleus& Nucleo, bool fc = false);

    void overlap(bool fc);

}


#endif // CSV_HPP