#include <iostream>
#include "exec/CSV.hpp"
#include "calculations/nuclear.hpp"
#include "other/plots/plot_other.hpp"
//#include "other/plots/plot_gnuplot.hpp"

int main() {
    //Modelos disponíveis(?): "GBW(new)", "GBW(old)", "bCGC", "IPSAT", "LHAnPDF", "IIM_S", "IIM_RS"
    //Átomos disponíveis: Pb208, proton
    //sqrt(s):5.36e3 -> CMS
            //13.6e3 -> run 3 LHCb

    /*
    std::vector<std::string> modelos = {"GBW(new)", "GBW(old)", "bCGC", "IPSAT", "LHAnPDF", "IIM_S", "IIM_RS"};
    
    for (const auto& modelo : modelos) {
        CSV::sigma_gamma_p(modelo);
    }*/
    

    //CSV::N_p("IPSAT");

    //CSV::sigma_gamma_p("GBW(new)"); 

    // para colisões A: GBW_Shadowing está disponível!

    //CSV::sigma_gamma_A("GBW(new)", proton, true);

    //Compare_all::N(1e-5);
    Compare_all::sigma_gamma_p(false); 

    Compare_all::sigma_gamma_p(true); 

    //Compare_all::rapidez_AA(3e3, Pb208, false);
    //Compare_all::rapidez_AA(3e3, Pb208, true);

    //Compare_all::overlap(true);
    //Compare_all::overlap(false);

    //CSV::rapidez_AA(5.36e3, "GBW(new)", Pb208, Pb208);
    

    //Shadowing_factor_csv(3.16227);


    return 0;
}
