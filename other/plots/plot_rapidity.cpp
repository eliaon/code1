#include "plot_io.hpp"
#include "../../libraries/matplotlib-cpp/matplotlibcpp.h"
#include "../utils.hpp"
#include <chrono>
#include <sstream>
#include <iomanip>

namespace plt = matplotlibcpp;

// -------------- PLOTA CURVAS DE RAPIDEZ PARA UM SQRT S ESCOLHIDO



void plot_rapidity()
{
    std::string meson = get_meson();
    std::string filename =
        "out/plots/Rapidez/" + meson +"/" + meson + "_rapidez_5.26e+03GeV.csv";

    std::vector<double> Y, rap_GLC, rap_BG;

    read_csv(filename, Y, rap_GLC, rap_BG);

    plt::figure_size(800,600);

    plt::plot(Y, rap_GLC,
        {{"label","GLC"},
         {"color","red"},
         {"linestyle","-"},
         {"linewidth","1.8"}});

    plt::plot(Y, rap_BG,
        {{"label","BG"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.2"}});

    plt::xlim(-8,8);

    plt::xlabel("Y");
    plt::ylabel("dσ/dY [nb]");

    double s = 5.26e3;

    std::string meson_label = meson;
    if(meson == "Jpsi")
        meson_label = "psi";

    std::stringstream title;
    title << "Distribuição de rapidez do $" << meson_label
          << "$ em pp a "
          << std::fixed << std::setprecision(2)
          << s/1e3 << " TeV";

    plt::title(title.str());

    plt::grid(true);
    plt::legend();

    std::string out =
        "out/plots/Rapidez/"+meson+"/" + meson +
        "_rapidez_5p26TeV_" +
        timestamp() + ".pdf";

    plt::save(out);
    plt::show();
}

void plot_rapidez_PbPb_Jpsi(std::string filename, std::string model, double sqrt_s)
{
    std::vector<double> Y, rap_GLC, rap_BG;

    read_csv(filename, Y, rap_GLC, rap_BG);

    plt::figure_size(800,600);

    plt::plot(Y, rap_GLC,
        {{"label","GLC"},
         {"color","red"},
         {"linestyle","-"},
         {"linewidth","1.8"}});

    plt::plot(Y, rap_BG,
        {{"label","BG"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.2"}});

    plt::xlim(-8,8);

    plt::xlabel("Y");
    plt::ylabel("dσ/dY [nb]");

    double sqrt_s_TeV = sqrt_s / 1e3;

    std::stringstream title;
    title << "Distribuição de rapidez do $J/\\psi$ em Pb-Pb a " << sqrt_s_TeV << " TeV";

    plt::title(title.str());

    plt::grid(true);
    plt::legend();


    std::string out =
        "out/plots/Rapidez/psi/PbPb-Jpsi_rapidez("+ model +")_+" + doubleParaString(sqrt_s_TeV) + "TeV_" +
        timestamp() + ".pdf";

    plt::save(out);
    plt::show();
}

void plot_rapidez_PbPb_phi(std::string csv, std::string model, double sqrt_s)
{
    std::vector<double> Y, rap_GLC, rap_BG;


    read_csv(csv, Y, rap_GLC, rap_BG);

    plt::figure_size(800,600);

    // --- dados experimentais ---
    std::vector<double> y, dsigma_dy, error;

    if(sqrt_s == 5.36e3) // só tem dados para 5.36 TeV
    {
    read_rapidity_hepdata(
        "out/csv/Rapidez/phi/data/rapidity_PbPb-phi_5360TeV.csv", // Observation of coherent ϕ(1020) meson photoproduction
        y,                                                   //    in ultraperipheral PbPb collisions at √sNN = 5.36 TeV
        dsigma_dy,                                           //http://dx.doi.org/10.1103/2ssw-wwyy
        error
    );


    plt::errorbar(y, dsigma_dy, error,
        {{"fmt","o"},
         {"color","blue"},
         {"label","CMS (2025)"}});

        plt::plot(Y, rap_GLC,
        {{"label","GLC "+ model},
         {"color","orange"},
         {"linestyle","-"},
         {"linewidth","1.8"}});

    plt::plot(Y, rap_BG,
        {{"label","BG "+ model},
         {"color","orange"},
         {"linestyle","--"},
         {"linewidth","1.2"}});
    }

    

    plt::ylim(10,130);
    plt::xlim(0.0,1.5);

    PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gcf()\n"   // garante que usa a figura atual
    "plt.yscale('log')\n"
    "plt.grid(True)\n"
);

    plt::xlabel("Y");
    plt::ylabel("dσ/dY [mb]");

    double sqrt_s_TeV = sqrt_s / 1e3;
    std::stringstream title;
    title << "Distribuição de rapidez do $\\phi$ em Pb-Pb a " << sqrt_s_TeV << " TeV";

    plt::title(title.str());

    plt::grid(true);
    plt::legend();

    std::string out =
        "out/plots/Rapidez/PbPb-phi_rapidez_+" + doubleParaString(sqrt_s_TeV) + "TeV_" +
        timestamp() + ".pdf";

    plt::save(out);
    plt::show();
}

void plot_rapidez_generico(std::string csv, std::string model, std::string meson, std::string collision, double sqrt_s)
{
    std::vector<double> Y, rap_GLC, rap_BG;

    read_csv(csv, Y, rap_GLC, rap_BG);

    plt::figure_size(800,600);

    plt::plot(Y, rap_GLC,
        {{"label","GLC"},
         {"color","red"},
         {"linestyle","-"},
         {"linewidth","1.8"}});

    plt::plot(Y, rap_BG,
        {{"label","BG"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.2"}});

    plt::xlim(-8,8);

    plt::xlabel("Y");
    plt::ylabel("dσ/dY [nb]");

    double sqrt_s_TeV = sqrt_s / 1e3;

    std::stringstream title;
    title << "Distribuição de rapidez do $"+ meson +"$ em " + collision +" a " << sqrt_s_TeV << " TeV no modelo" + model+".";

    plt::title(title.str());

    plt::grid(true);
    plt::legend();


    std::string out =
        "out/plots/Rapidez/"+ meson +"/"+ collision +"-"+ meson +"_rapidez("+ model +")_+" + doubleParaString(sqrt_s_TeV) + "TeV_" +
        timestamp() + ".pdf";

    plt::save(out);
    plt::show();
}

void plot_rapidez(std::string csv, std::string model, std::string meson, std::string collision,double sqrt_s)
{
    if(meson == "psi" and collision == "PbPb"){
        plot_rapidez_PbPb_Jpsi(csv, model, sqrt_s);
    } else if (meson == "phi" and collision == "PbPb"){
        plot_rapidez_PbPb_phi(csv, model, sqrt_s);
    } else{
        plot_rapidez_generico(csv, model, meson, collision, sqrt_s);
    }
}






















