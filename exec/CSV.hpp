#ifndef CSV_HPP
#define CSV_HPP

#include <iostream>
#include "../calculations/nuclear.hpp"

namespace CSV {

    void N_p(std::string model);

    void sigma_gamma_p(std::string model);

    void sigma_gamma_A(std::string model, nucleous&Nucleo);

    void rapidez_AA(double sqrt_s, std::string model, nucleous &Nucleo1, nucleous& Nucleo2);

}



#endif // CSV_HPP