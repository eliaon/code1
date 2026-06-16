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
    std::cout << "Entrou na função plot_rapidez_PbPb_Jpsi\n";
    std::cout << "csv = " << filename << '\n';

    std::cout << "Antes extrair_nome_base\n";
    std::string plotname = extrair_nome_base(filename);
    std::cout << "Depois extrair_nome_base\n";

    // =========================================================
    // dados teóricos
    // =========================================================

    std::vector<double> Y_th, rap_GLC, rap_BG;
    std::cout << "Antes read_csv\n";
    try {
        read_csv(
            filename,
            Y_th,
            rap_GLC,
            rap_BG
        );
    }
    catch (const std::exception& e) {
        std::cerr << "Error reading theoretical data: " << e.what() << std::endl;
        return;
    }
    std::cout << "Depois read_csv\n";

    std::cout << "Y_th = " << Y_th.size() << '\n';
    std::cout << "rap_GLC = " << rap_GLC.size() << '\n';
    std::cout << "rap_BG = " << rap_BG.size() << '\n';

    // =========================================================
    // figura
    // =========================================================
    std::cout << "Antes figure_size\n";
    plt::figure_size(800,600);
    std::cout << "Depois figure_size\n";

    // =========================================================
    // plotar dados teóricos
    // =========================================================
    std::cout << "Antes plot\n";
    fflush(stdout);

    // Convert vectors to Python lists using PyRun_SimpleString to avoid matplotlib-cpp bugs
    std::string Y_th_py = vec_to_pylist(Y_th);
    std::string rap_GLC_py = vec_to_pylist(rap_GLC);
    std::string rap_BG_py = vec_to_pylist(rap_BG);

    PyRun_SimpleString(
        ("import matplotlib.pyplot as plt\n"
         "Y_th = " + Y_th_py + "\n"
         "rap_GLC = " + rap_GLC_py + "\n"
         "rap_BG = " + rap_BG_py + "\n"
         "plt.plot(Y_th, rap_GLC, label='GLC', color='red', linestyle='-', linewidth=1.8)\n"
         "plt.plot(Y_th, rap_BG, label='BG', color='red', linestyle='--', linewidth=1.2)\n"
        ).c_str()
    );

    std::cout << "Depois plot\n";

    plt::ylabel("$d\\sigma/dY$ [nb]");

    std::stringstream title;
    title << "$J/\\psi$ fotoprodução coerente em PbPb "
          << "($\\sqrt{s_{NN}}="
          << sqrt_s/1e3
          << "\\,\\mathrm{TeV}$) - "
          << model;

    plt::title(title.str());
    plt::xlabel("$Y$");
    plt::legend();

    // =========================================================
    // ajustes finais
    // =========================================================

    PyRun_SimpleString(
        "import matplotlib.pyplot as plt\n"
        "ax = plt.gca()\n"
        "ax.set_yscale('log')\n"
        "ax.grid(True, which='both', linestyle='--', alpha=0.6)\n"
        "ax.set_xlim(-2.5, 2.5)\n"
    );

    // =========================================================
    // salvar
    // =========================================================

    std::string out =
        "out/plots/Rapidez/psi/"
        + plotname + "_"
        + timestamp()
        + ".pdf";

    std::string label =
        "Arquivo: " + plotname;

    PyRun_SimpleString(
        ("import matplotlib.pyplot as plt\n"
         "plt.figtext(0.01,0.01,'"
         + label +
         "', fontsize=8, alpha=0.7)\n").c_str()
    );

    std::cout << "Salvando figura em: "
              << out << std::endl;

    plt::save(out);
    plt::show();
}

void plot_rapidez_PbPb_phi(std::string csv,
                           std::string model,
                           double sqrt_s)
{

std::cout << "csv = " << csv << '\n';


    std::string plotname = extrair_nome_base(csv);

    // =========================================================
    // dados teóricos
    // =========================================================

    std::vector<double> Y_th, rap_GLC, rap_BG;

    try {
        read_csv(
            csv,
            Y_th,
            rap_GLC,
            rap_BG
        );
    }
    catch (const std::exception& e) {
        std::cerr << "Error reading theoretical data: " << e.what() << std::endl;
        return;
    }

    // =========================================================
    // dados experimentais
    // =========================================================
    std::vector<double>
        Y_exp,
        dsigma_exp,
        error_exp;


    if(sqrt_s == 5.36e3)
    {
        fflush(stdout);
        read_rapidity_hepdata(
            "out/csv/Rapidez/phi/data/rapidity_PbPb-phi_5360TeV.csv",
            Y_exp,
            dsigma_exp,
            error_exp
        );
        fflush(stdout);
    }
    else
    {
        std::cout << "Warning: sqrt_s = " << sqrt_s << " does not match experimental data (5.36e3)" << std::endl;
    }

    fflush(stdout);
    for(size_t i=0;i<Y_exp.size();++i)
{

 
}
    fflush(stdout);

    // =========================================================
    // figura
    // =========================================================

    plt::figure_size(800,600);

    // =========================================================
    // painel superior
    // =========================================================

    plt::subplot2grid(2,1,0);

    // curvas teóricas

fflush(stdout);

    // Convert vectors to Python lists using PyRun_SimpleString to avoid matplotlib-cpp bugs
    std::string Y_th_py = vec_to_pylist(Y_th);
    std::string rap_GLC_py = vec_to_pylist(rap_GLC);
    std::string rap_BG_py = vec_to_pylist(rap_BG);

    PyRun_SimpleString(
        ("import matplotlib.pyplot as plt\n"
         "Y_th = " + Y_th_py + "\n"
         "rap_GLC = " + rap_GLC_py + "\n"
         "rap_BG = " + rap_BG_py + "\n"
         "plt.plot(Y_th, rap_GLC, label='GLC', color='orange', linestyle='-', linewidth=1.8)\n"
         "plt.plot(Y_th, rap_BG, label='BG', color='orange', linestyle='--', linewidth=1.2)\n"
        ).c_str()
    );



    // dados experimentais

    if(!Y_exp.empty())
    {
        std::string Y_exp_py = vec_to_pylist(Y_exp);
        std::string dsigma_py = vec_to_pylist(dsigma_exp);
        std::string error_py = vec_to_pylist(error_exp);

        PyRun_SimpleString(
            ("import matplotlib.pyplot as plt\n"
             "Y_exp = " + Y_exp_py + "\n"
             "dsigma = " + dsigma_py + "\n"
             "error = " + error_py + "\n"
             "plt.errorbar(Y_exp, dsigma, error, fmt='o', color='blue', label='CMS (2025)')\n"
            ).c_str()
        );
    }

    plt::ylabel("$d\\sigma/dY$ [mb]");

    std::stringstream title;

    title << "$\\phi$ fotoprodução coerente em PbPb "
          << "($\\sqrt{s_{NN}}="
          << sqrt_s/1e3
          << "\\,\\mathrm{TeV}$) - "
          << model;

    plt::title(title.str());

    plt::legend();

    // =========================================================
    // painel inferior : desvio relativo
    // =========================================================

    plt::subplot2grid(2,1,1);

    std::vector<double> dev_GLC;
    std::vector<double> dev_BG;

    if(!Y_exp.empty())
    {
        for(size_t i=0; i<Y_exp.size(); ++i)
        {
            double Yd = Y_exp[i];

            double theo_GLC = 0.0;
            double theo_BG  = 0.0;

            double best = 1e99;

            for(size_t j=0; j<Y_th.size(); ++j)
            {
                double diff = std::abs(Y_th[j] - Yd);

                if(diff < best)
                {
                    best = diff;

                    theo_GLC = rap_GLC[j];
                    theo_BG  = rap_BG[j];
                }
            }

            if(dsigma_exp[i] > 0.0)
            {
                dev_GLC.push_back(
                    (theo_GLC - dsigma_exp[i])
                    / dsigma_exp[i]
                );

                dev_BG.push_back(
                    (theo_BG - dsigma_exp[i])
                    / dsigma_exp[i]
                );
            }
        }
    }

    // GLC

    if(!Y_exp.empty() && !dev_GLC.empty())
    {
        std::string Y_exp_py = vec_to_pylist(Y_exp);
        std::string dev_GLC_py = vec_to_pylist(dev_GLC);
        std::string dev_BG_py = vec_to_pylist(dev_BG);

        PyRun_SimpleString(
            ("import matplotlib.pyplot as plt\n"
             "Y_exp = " + Y_exp_py + "\n"
             "dev_GLC = " + dev_GLC_py + "\n"
             "dev_BG = " + dev_BG_py + "\n"
             "plt.plot(Y_exp, dev_GLC, marker='o', linestyle='none', color='orange', label='GLC')\n"
             "plt.plot(Y_exp, dev_BG, marker='x', linestyle='none', color='darkred', label='BG')\n"
            ).c_str()
        );
    }

    // linha horizontal

    PyRun_SimpleString(
        "import matplotlib.pyplot as plt\n"
        "plt.plot([0.0, 1.5], [0.0, 0.0], color='black', linestyle='--')\n"
    );

    plt::xlabel("$Y$");
    plt::ylabel("$\\Delta$");

    // =========================================================
    // ajustes finais
    // =========================================================

    if(!Y_exp.empty())
    {
        PyRun_SimpleString(
            "import matplotlib.pyplot as plt\n"

            "axes = plt.gcf().axes\n"

            // eixo superior
            "axes[0].set_position([0.12,0.32,0.83,0.63])\n"

            // eixo inferior
            "axes[1].set_position([0.12,0.10,0.83,0.16])\n"

            // escalas
            "axes[0].set_yscale('log')\n"

            // grids
            "axes[0].grid(True, which='both', linestyle='--', alpha=0.6)\n"
            "axes[1].grid(True, linestyle='--', alpha=0.5)\n"

            // limites
            "axes[0].set_xlim(0.0,1.5)\n"
            "axes[0].set_ylim(10,130)\n"

            "axes[1].set_xlim(0.0,1.5)\n"

            // remover labels superiores
            "axes[0].tick_params(labelbottom=False)\n"
        );
    }

    // =========================================================
    // salvar
    // =========================================================

    std::string out =
        "out/plots/Rapidez/phi/"
        + plotname + "_"
        + timestamp()
        + ".pdf";

    std::string label =
        "Arquivo: " + plotname;

    PyRun_SimpleString(
        ("import matplotlib.pyplot as plt\n"
         "plt.figtext(0.01,0.01,'"
         + label +
         "', fontsize=8, alpha=0.7)\n").c_str()
    );

    std::cout << "Salvando figura em: "
              << out << std::endl;

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






















