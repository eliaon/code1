#include <iostream>
#include "exec/CSV.hpp"
#include "calculations/nuclear.hpp"

int main() {
    //Modelos disponíveis(?): "GBW(new)", "GBW(old)", "IPSAT", "LHAnPDF"
    //Átomos disponíveis: Pb208, proton
    //sqrt(s):5.36e3 -> CMS
            //13.6e3 -> run 3 LHCb

    CSV::N_p("LHAnPDF");

    //CSV::sigma_gamma_p("GBW(new)"); 

    //CSV::sigma_gamma_A("GBW(new)", Pb208);

    //CSV::rapidez_AA(5.36e3, "GBW(new)", Pb208, Pb208);
    
    return 0;
}