#include "plot_io.hpp"
#include "../utils.hpp"
#include "../../calculations/nuclear.hpp"

#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "../../libraries/matplotlib-cpp/matplotlibcpp.h"

namespace plt = matplotlibcpp;

namespace {
bool is_proton_target(const Nucleus& nucleus)
{
    return nucleus.name == "p" || nucleus.name == "proton";
}
} // namespace

void plot_sigma_placeholder(std::string csv, std::string dipolemodel, std::string meson, Nucleus& nucleus);

// ------------- PLOTA CURVA PARA dσ/dt PARA UMA DADA ENERGIA EM GeV

void dsigma_dt_curve(const std::string& W)
{
    std::string filename = "out/csv/dsigmadt/psi/Jpsi_dsigma_dt_W=" + W + "GeV.csv";

    std::vector<double> t, glc, bg;

    read_csv(filename, t, glc, bg);

    if(W == "100")
    {
        for(auto& v : glc) v *= 5.0;
        for(auto& v : bg)  v *= 5.0;

        plt::plot(t, glc, {{"label","GLC"}});
        plt::plot(t, bg, {{"label","BG"}});
    }
    else
    {
        plt::plot(t, glc,
            {{"label","GLC (W=" + W + " GeV)"}});
        plt::plot(t, bg,
            {{"label","BG (W=" + W + " GeV)"},{"linestyle","--"}});
    }
}


// ------------- PLOTA DADOS EXPERIMENTAIS PARA dσ/dt DE FOTOPRODUÇÃO DE φ

void plot_dsigma_exp_data_phi()
{
    std::ifstream file("out/csv/dsigmadt/phi/data/HEPData-dsigdt_gammap_phi.csv");

    std::vector<double> t;
    std::vector<double> dsdt;
    std::vector<double> err;

    std::string line;

    while(std::getline(file, line))
    {
        // ignorar comentários HEPData
        if(line.empty() || line[0] == '#')
            continue;

        // ignorar header textual
        if(!std::isdigit(line[0]))
            continue;

        std::stringstream ss(line);
        std::string value;

        std::vector<double> row;

        while(std::getline(ss, value, ','))
            row.push_back(std::stod(value));

        if(row.size() < 5)
            continue;

        double t_val = row[0];
        double sigma = row[3] * 1000.0;   // µb → nb
        double error = std::abs(row[4]) * 1000.0;

        t.push_back(t_val);
        dsdt.push_back(sigma);
        err.push_back(error);
    }

    plt::errorbar(t, dsdt, err,
    {{"fmt","o"},
     {"color","black"},
     {"label","HERA (W≈70 GeV)"}});
}

// ------------- PLOTAR AS CURVAS dσ/dt PARA TODOS OS W ESCOLHIDOS

void plot_dsigma_dt(std::string meson)
{
    std::vector<std::string> w_values = {"70"};

    plt::figure_size(700,500);

    for(const auto& W : w_values)
        dsigma_dt_curve(W);

    // plota dados experimentais apenas para o méson phi
    if(meson == "phi")
        plot_dsigma_exp_data_phi();

    plt::xlim(0.0, 2.5);
    plt::ylim(1e-3, 1e4);

    PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.xscale('linear')\n"
    "plt.yscale('log')\n"
    );

    plt::xlabel("$|t| (GeV^2)$");
    plt::ylabel("$d\\sigma/dt$ (nb)");

    
    plt::grid(true);
    plt::legend();

    plt::save("out/plots/dsigmadt/" + meson + "/"+ meson + "_dsigma_dt_" + timestamp() + ".pdf");
    plt::show();
}














// ------------------ PLOT SEÇÕES DE CHOQUE INTEGRAIS



void plot_sigma_Jpsi(std::string csv_file, std::string dipolemodel)
{
    int Q2 = 0;

   
    // =========================================================
    // carregar curvas teóricas
    // =========================================================

    std::vector<double> W_th, sigma_GLC, sigma_BG;

    try
    {
        read_csv(
            csv_file,
            W_th,
            sigma_GLC,
            sigma_BG);
    }
    catch(...)
    {
        std::cerr << "Falha ao carregar curva teórica\n";
        return;
    }

    // =========================================================
    // FIGURA
    // =========================================================

    plt::figure_size(800,600);

    // =========================================================
    // estilos
    // =========================================================
 std::vector<int> dataset;
    std::vector<double> W_exp, sigma_exp, err_exp;

    read_sigma_exp(
        "out/csv/sigma/psi/data/sigma_gammap_jpsi.csv",
        dataset,
        W_exp,
        sigma_exp,
        err_exp);

    std::map<int,std::string> exp_map = {
        {0,"H1"},
        {1,"H1"},
        {2,"ALICE"},
        {3,"LHCb"}
    };

    std::map<std::string,std::string> colors = {
        {"H1","blue"},
        {"ALICE","black"},
        {"LHCb","purple"}
    };

    std::map<std::string,std::string> markers = {
        {"H1","o"},
        {"ALICE","s"},
        {"LHCb","^"}
    };

    // =========================================================
    // PAINEL SUPERIOR
    // =========================================================
    plt::subplot2grid(2,1,0);
    std::map<int,std::vector<int>> groups;

    for(size_t i=0;i<dataset.size();++i)
        groups[dataset[i]].push_back(i);

    for(auto& g : groups)
    {
        int d = g.first;
        std::string exp = exp_map[d];

        std::vector<double> W,s,e;

        for(int idx : g.second)
        {
            W.push_back(W_exp[idx]);
            s.push_back(sigma_exp[idx]);
            e.push_back(err_exp[idx]);
        }

        plt::errorbar(
            W,
            s,
            e,
            {{"label",exp},
             {"marker",markers[exp]},
             {"color",colors[exp]},
             {"linestyle","none"}}
        );
    }

    // curvas teóricas

    plt::plot(
        W_th,
        sigma_GLC,
        {{"label","GLC"},
         {"color","red"},
         {"linewidth","1.8"}}
    );

    plt::plot(
        W_th,
        sigma_BG,
        {{"label","BG"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.2"}}
    );

    plt::ylabel("$\\sigma$ [nb]");
    plt::title("$J/\\psi$ produção exclusiva ($Q^2=0$) - " + dipolemodel );

    plt::legend();

    // =========================================================
    // PAINEL INFERIOR : DESVIO
    // =========================================================

    plt::subplot2grid(2,1,1);

    std::vector<double> dev_GLC;
    std::vector<double> dev_BG;

    for(size_t i=0; i<W_exp.size(); ++i)
    {
        double Wd = W_exp[i];

        double theo_GLC = 0.0;
        double theo_BG  = 0.0;

        double best = 1e99;

        for(size_t j=0; j<W_th.size(); ++j)
        {
            double diff = std::abs(W_th[j] - Wd);

            if(diff < best)
            {
                best = diff;

                theo_GLC = sigma_GLC[j];
                theo_BG  = sigma_BG[j];
            }
        }

        dev_GLC.push_back(
            (theo_GLC - sigma_exp[i]) / sigma_exp[i]
        );

        dev_BG.push_back(
            (theo_BG - sigma_exp[i]) / sigma_exp[i]
        );
    }

    // GLC

    plt::plot(
        W_exp,
        dev_GLC,
        {{"marker","o"},
         {"linestyle","none"},
         {"color","red"},
         {"label","GLC"}}
    );

    // BG

    plt::plot(
        W_exp,
        dev_BG,
        {{"marker","x"},
         {"linestyle","none"},
         {"color","darkred"},
         {"label","BG"}}
    );

    // linha horizontal

    plt::plot(
        std::vector<double>{20,10000},
        std::vector<double>{0,0},
        {{"color","black"},
         {"linestyle","--"}}
    );

    plt::xlabel("$W$ [GeV]");
    plt::ylabel("$\\Delta$");

    // =========================================================
    // AJUSTES FINAIS
    // =========================================================

    
    PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"

    "axes = plt.gcf().axes\n"

    // eixo superior
    "axes[0].set_position([0.12,0.32,0.83,0.63])\n"

    // eixo inferior
    "axes[1].set_position([0.12,0.10,0.83,0.16])\n"
    "axes[1].set_xlim(30, 3000)\n"
// escalas
        "axes[0].set_xscale('log')\n"
        "axes[0].set_yscale('log')\n"
);

    // =========================================================
    // salvar
    // =========================================================

    std::string plotname = extrair_nome_base(csv_file);

    std::string out =
        "out/plots/sigma/psi/" +
        plotname + "_" +
        timestamp() + ".pdf";

    std::string label = "Arquivo: " + plotname;

    PyRun_SimpleString(
        ("import matplotlib.pyplot as plt\n"
         "plt.figtext(0.01, 0.01, '" +
         label +
         "', fontsize=8, alpha=0.7)\n").c_str()
    );

    plt::save(out);

    plt::show();
}

void plot_sigma_phi(std::string csv, std::string dipolemodel)
{
    int Q2 = 0;

    std::string plotname = extrair_nome_base(csv);

    // =========================================================
    // dados teóricos
    // =========================================================

    std::vector<double> W_th, sigma_GLC, sigma_BG;

    read_csv(
        csv,
        W_th,
        sigma_GLC,
        sigma_BG
    );

    // =========================================================
    // dados experimentais
    // =========================================================

    std::vector<int> dataset_fixedpoint;

    std::vector<double>
        W_fixedpoint,
        sigma_fixedpoint,
        error_fixedpoint;

    read_sigma_exp(
        "out/csv/sigma/phi/data/phi_fixedpoint_data(nb).csv",
        dataset_fixedpoint,
        W_fixedpoint,
        sigma_fixedpoint,
        error_fixedpoint
    );

    std::vector<int> dataset_ZEUS;

    std::vector<double>
        W_ZEUS,
        sigma_ZEUS,
        error_ZEUS;

    read_sigma_exp(
        "out/csv/sigma/phi/data/phi_sigma_expdata_ZEUS(1994).csv",
        dataset_ZEUS,
        W_ZEUS,
        sigma_ZEUS,
        error_ZEUS
    );

    // =========================================================
    // figura
    // =========================================================

    plt::figure_size(800,600);

    // =========================================================
    // painel superior
    // =========================================================

    plt::subplot2grid(2,1,0);

    // curvas teóricas

    plt::plot(
        W_th,
        sigma_GLC,
        {{"label","GLC"},
         {"color","red"},
         {"linestyle","-"},
         {"linewidth","1.8"}}
    );

    plt::plot(
        W_th,
        sigma_BG,
        {{"label","BG"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.2"}}
    );

    // dados experimentais

    plt::errorbar(
        W_fixedpoint,
        sigma_fixedpoint,
        error_fixedpoint,
        {{"fmt","o"},
         {"color","blue"},
         {"label","Fixed Point"}}
    );

    plt::errorbar(
        W_ZEUS,
        sigma_ZEUS,
        error_ZEUS,
        {{"fmt","s"},
         {"color","green"},
         {"label","ZEUS (1994)"}}
    );

    plt::ylabel("$\\sigma$ [nb]");

    std::stringstream title;

    title << "$\\phi$ fotoprodução exclusiva "
          << "($\\gamma p \\to \\phi p$) - " << dipolemodel;

    plt::title(title.str());

    plt::legend();

    // =========================================================
    // painel inferior : desvio relativo
    // =========================================================

    plt::subplot2grid(2,1,1);

    std::vector<double> W_exp;
    std::vector<double> sigma_exp;

    // juntar datasets experimentais

    for(size_t i=0;i<W_fixedpoint.size();++i)
    {
        W_exp.push_back(W_fixedpoint[i]);
        sigma_exp.push_back(sigma_fixedpoint[i]);
    }

    for(size_t i=0;i<W_ZEUS.size();++i)
    {
        W_exp.push_back(W_ZEUS[i]);
        sigma_exp.push_back(sigma_ZEUS[i]);
    }

    // desvios

    std::vector<double> dev_GLC;
    std::vector<double> dev_BG;

    for(size_t i=0; i<W_exp.size(); ++i)
    {
        double Wd = W_exp[i];

        double theo_GLC = 0.0;
        double theo_BG  = 0.0;

        double best = 1e99;

        for(size_t j=0; j<W_th.size(); ++j)
        {
            double diff = std::abs(W_th[j] - Wd);

            if(diff < best)
            {
                best = diff;

                theo_GLC = sigma_GLC[j];
                theo_BG  = sigma_BG[j];
            }
        }

        dev_GLC.push_back(
            (theo_GLC - sigma_exp[i]) / sigma_exp[i]
        );

        dev_BG.push_back(
            (theo_BG - sigma_exp[i]) / sigma_exp[i]
        );
    }

    // GLC

    plt::plot(
        W_exp,
        dev_GLC,
        {{"marker","o"},
         {"linestyle","none"},
         {"color","red"},
         {"label","GLC"}}
    );

    // BG

    plt::plot(
        W_exp,
        dev_BG,
        {{"marker","x"},
         {"linestyle","none"},
         {"color","darkred"},
         {"label","BG"}}
    );

    // linha horizontal

    plt::plot(
        std::vector<double>{8,100},
        std::vector<double>{0,0},
        {{"color","black"},
         {"linestyle","--"}}
    );

    plt::xlabel("$W$ [GeV]");
    plt::ylabel("$\\Delta$");

    // =========================================================
    // ajustes finais
    // =========================================================

    PyRun_SimpleString(
        "import matplotlib.pyplot as plt\n"

        "axes = plt.gcf().axes\n"

     "axes[0].set_position([0.12,0.32,0.83,0.63])\n"

    // eixo inferior
    "axes[1].set_position([0.12,0.10,0.83,0.16])\n"
        // escalas
        "axes[0].set_xscale('log')\n"
        "axes[0].set_yscale('log')\n"

        "axes[1].set_xscale('log')\n"

        // grids
        "axes[0].grid(True, which='both', linestyle='--', alpha=0.6)\n"
        "axes[1].grid(True, linestyle='--', alpha=0.5)\n"

        // limites
        "axes[0].set_xlim(8,100)\n"
        "axes[0].set_ylim(100,10000)\n"

        "axes[1].set_xlim(8,100)\n"

        // remover labels superiores
        "axes[0].tick_params(labelbottom=False)\n"

    );

    // =========================================================
    // salvar
    // =========================================================

    std::string out =
        "out/plots/sigma/phi/" +
        plotname +"_" +
        timestamp() + ".pdf";

    std::string label = "Arquivo: " + plotname;

    PyRun_SimpleString(
        ("import matplotlib.pyplot as plt\n"
         "plt.figtext(0.01, 0.01, '" +
         label +
         "', fontsize=8, alpha=0.7)\n").c_str()
    );

    std::cout << "Salvando figura em: " << out << std::endl;
    plt::save(out);

    plt::show();
}


void plot_sigma_rho(std::string csv, std::string dipolemodel)
{
    int Q2 = 0;

    std::string plotname = extrair_nome_base(csv);

    // =========================================================
    // dados teóricos
    // =========================================================

    std::vector<double> W_th, sigma_GLC, sigma_BG;

    read_csv(
        csv,
        W_th,
        sigma_GLC,
        sigma_BG
    );

    // =========================================================
    // dados experimentais
    // =========================================================

    std::vector<int> dataset_ALICE;

    std::vector<double>
        W_ALICE,
        sigma_ALICE,
        error_ALICE;

        read_xyscan_csv("out/csv/sigma/rho/data/rho_photoproduction_ZEUS(abramovic2000).csv",
                        W_ALICE,
                        sigma_ALICE,
                        error_ALICE);


    // =========================================================
    // figura
    // =========================================================

    plt::figure_size(800,600);

    // =========================================================
    // painel superior
    // =========================================================

    plt::subplot2grid(2,1,0);

    // curvas teóricas

    plt::plot(
        W_th,
        sigma_GLC,
        {{"label","GLC"},
         {"color","red"},
         {"linestyle","-"},
         {"linewidth","1.8"}}
    );

    plt::plot(
        W_th,
        sigma_BG,
        {{"label","BG"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.2"}}
    );

    // dados experimentais

    plt::errorbar(
        W_ALICE,
        sigma_ALICE, 
        error_ALICE ,
        {{"fmt","o"},
         {"color","blue"},
         {"label","ALICE (2000)"}}
    );


    plt::ylabel("$\\sigma$ [nb]");

    std::stringstream title;

    title << "$\\rho$ fotoprodução exclusiva "
          << "($\\gamma p \\to \\rho p$) - " << dipolemodel ;

    plt::title(title.str());

    plt::legend();

    // =========================================================
    // painel inferior : desvio relativo
    // =========================================================

    plt::subplot2grid(2,1,1);

    std::vector<double> W_exp;
    std::vector<double> sigma_exp;

    // juntar datasets experimentais

    for(size_t i=0;i<W_ALICE.size();++i)
    {
        W_exp.push_back(W_ALICE[i]);
        sigma_exp.push_back(sigma_ALICE[i]);
    }

    // desvios

    std::vector<double> dev_GLC;
    std::vector<double> dev_BG;

    for(size_t i=0; i<W_exp.size(); ++i)
    {
        double Wd = W_exp[i];

        double theo_GLC = 0.0;
        double theo_BG  = 0.0;

        double best = 1e99;

        for(size_t j=0; j<W_th.size(); ++j)
        {
            double diff = std::abs(W_th[j] - Wd);

            if(diff < best)
            {
                best = diff;

                theo_GLC = sigma_GLC[j];
                theo_BG  = sigma_BG[j];
            }
        }

        dev_GLC.push_back(
            (theo_GLC - sigma_exp[i]) / sigma_exp[i]
        );

        dev_BG.push_back(
            (theo_BG - sigma_exp[i]) / sigma_exp[i]
        );
    }

    // GLC

    plt::plot(
        W_exp,
        dev_GLC,
        {{"marker","o"},
         {"linestyle","none"},
         {"color","red"},
         {"label","GLC"}}
    );

    // BG

    plt::plot(
        W_exp,
        dev_BG,
        {{"marker","x"},
         {"linestyle","none"},
         {"color","darkred"},
         {"label","BG"}}
    );

    // linha horizontal

    plt::plot(
        std::vector<double>{8,1000},
        std::vector<double>{0,0},
        {{"color","black"},
         {"linestyle","--"}}
    );

    plt::xlabel("$W$ [GeV]");
    plt::ylabel("$\\Delta$");

    // =========================================================
    // ajustes finais
    // =========================================================

    PyRun_SimpleString(
        "import matplotlib.pyplot as plt\n"

        "axes = plt.gcf().axes\n"

          // eixo superior
"axes[0].set_position([0.12,0.32,0.83,0.63])\n"

    // eixo inferior
    "axes[1].set_position([0.12,0.10,0.83,0.16])\n"

        // escalas
        "axes[0].set_xscale('log')\n"
        "axes[0].set_yscale('log')\n"

        "axes[1].set_xscale('log')\n"

        // grids
        "axes[0].grid(True, which='both', linestyle='--', alpha=0.6)\n"
        "axes[1].grid(True, linestyle='--', alpha=0.5)\n"

        // limites
        "axes[0].set_xlim(8,1000)\n"

        "axes[1].set_xlim(8,1000)\n"

        // remover labels superiores
        "axes[0].tick_params(labelbottom=False)\n"
    );

    // =========================================================
    // salvar
    // =========================================================

    std::string out =
        "out/plots/sigma/rho/" +
        plotname + "_" +
        timestamp() + ".pdf";

    std::string label = "Arquivo: " + plotname;

    PyRun_SimpleString(
        ("import matplotlib.pyplot as plt\n"
         "plt.figtext(0.01, 0.01, '" +
         label +
         "', fontsize=8, alpha=0.7)\n").c_str()
    );

    std::cout << "Salvando figura em: " << out << std::endl;
    plt::save(out);

    plt::show();
}

void plot_sigma_placeholder(std::string csv, std::string dipolemodel, std::string meson)
{
    plot_sigma_placeholder(csv, dipolemodel, meson, proton);
}

void plot_sigma_placeholder(std::string csv, std::string dipolemodel, std::string meson, Nucleus& nucleus)
{
        plt::clf();
        plt::figure_size(800,600);

        std::vector<double> W, sigma_GLC, sigma_BG;
        read_csv(
            csv,
            W,
            sigma_GLC,
            sigma_BG
        );
        

        if(meson == "psi" and nucleus.name == "Pb")
        {
            std::vector<double> W_CMS, sigma_CMS, error_CMS;
            read_xyscan_csv("out/csv/sigma/psi/data/sigma_psi_gammaPb_CMS.csv",
                            W_CMS, sigma_CMS, error_CMS);

            plt::errorbar(W_CMS, sigma_CMS, error_CMS,
                {{"fmt","o"}, {"color","blue"}, {"label","CMS"}});// CMS Collaboration. Probing Small Bjorken-x Nuclear Gluonic Structure via
                                                                //Coherent J/ψ Photoproduction in Ultraperipheral Pb-Pb Collisions at √sNN = 5.02 TeV.
                                                                //Phys. Rev. Lett., American Physical Society, 2023.

            std::vector<double> W_ALICE_2023, sigma_ALICE_2023, error_ALICE_2023;
            read_xyscan_csv("out/csv/sigma/psi/data/sigma_psi_gammaPb_ALICE_2023.csv",  
                            W_ALICE_2023, sigma_ALICE_2023, error_ALICE_2023);

            plt::errorbar(W_ALICE_2023, sigma_ALICE_2023, error_ALICE_2023,      //ALICE Collaboration. Energy dependence of coherent photonuclear production of
                                                                        //J/ψ mesons in ultra-peripheral Pb-Pb collisions at √sNN = 5.02 TeV. 2023. 
                {{"fmt","s"}, {"color","green"}, {"label","ALICE(2023)"}});

            std::vector<double> W_ALICE_2024, sigma_ALICE_2024, error_ALICE_2024;
            read_xyscan_csv("out/csv/sigma/psi/data/sigma_psi_gammaPb_ALICE_2024.csv",
                            W_ALICE_2024, sigma_ALICE_2024, error_ALICE_2024);  // SHATAT, A. Charmonium photoproduction in Pb-Pb collisions with nuclear overlap
                                                                        //measured with ALICE at the LHC. Tese (Theses) — Université Paris-Saclay, 2024. https://theses.hal.science/tel-04797642
    
            plt::errorbar(W_ALICE_2024, sigma_ALICE_2024, error_ALICE_2024,
                {{"fmt","^"}, {"color","orange"}, {"label","ALICE(2024)"}});
        }

        if (meson == "psi" and is_proton_target(nucleus))
        {
            std::vector<int> dataset;
            std::vector<double> W_exp, sigma_exp, err_exp;

            read_sigma_exp(
                "out/csv/sigma/psi/data/sigma_gammap_jpsi.csv",
                dataset,
                W_exp,
                sigma_exp,
                err_exp);

            std::map<int,std::string> exp_map = {
                {0,"H1"},
                {1,"H1"},
                {2,"ALICE"},
                {3,"LHCb"}
            };

            std::map<std::string,std::string> colors = {
                {"H1","blue"},
                {"ALICE","black"},
                {"LHCb","purple"}
            };

            std::map<std::string,std::string> markers = {
                {"H1","o"},
                {"ALICE","s"},
                {"LHCb","^"}
            };

            std::map<int,std::vector<int>> groups;

            for(size_t i=0;i<dataset.size();++i)
                groups[dataset[i]].push_back(i);

            for(auto& g : groups)
            {
                int d = g.first;
                std::string exp = exp_map[d];

                std::vector<double> W,s,e;

                for(int idx : g.second)
                {
                    W.push_back(W_exp[idx]);
                    s.push_back(sigma_exp[idx]);
                    e.push_back(err_exp[idx]);
                }

                plt::errorbar(
                    W,
                    s,
                    e,
                    {{"label",exp},
                     {"marker",markers[exp]},
                     {"color",colors[exp]},
                     {"linestyle","none"}}
                );
            }
        }

        if (meson == "phi" and is_proton_target(nucleus))
        {
            std::vector<int> dataset_fixedpoint;

            std::vector<double>
                W_fixedpoint,
                sigma_fixedpoint,
                error_fixedpoint;

            read_sigma_exp(
                "out/csv/sigma/phi/data/phi_fixedpoint_data(nb).csv",
                dataset_fixedpoint,
                W_fixedpoint,
                sigma_fixedpoint,
                error_fixedpoint
            );

            plt::errorbar(
                W_fixedpoint,
                sigma_fixedpoint,
                error_fixedpoint,
                {{"fmt","o"},
                 {"color","blue"},
                 {"label","Fixed Point"}}
            );
        }

        if (meson == "rho" and is_proton_target(nucleus))
        {
            std::vector<int> dataset_ALICE;

            std::vector<double>
                W_ALICE,
                sigma_ALICE,
                error_ALICE;

                read_xyscan_csv("out/csv/sigma/rho/data/rho_photoproduction_ZEUS(abramovic2000).csv",
                                W_ALICE,
                                sigma_ALICE,
                                error_ALICE);

            plt::errorbar(
                W_ALICE,
                sigma_ALICE, 
                error_ALICE ,
                {{"fmt","o"},
                 {"color","blue"},
                 {"label","ALICE (2000)"}}
            );
        }

        plt::plot(
            W,
            sigma_GLC,
            {{"label","GLC"},
             {"color","red"},
             {"linestyle","-"},
             {"linewidth","1.8"}}
        );
        plt::plot(
            W,
            sigma_BG,
            {{"label","BG"},
             {"color","red"},
             {"linestyle","--"},
             {"linewidth","1.2"}}
        );

        plt::xlabel("$W$ [GeV]");
        plt::ylabel("$\\sigma$ [nb]");
        plt::title("Produção exclusiva de $\\" + meson + "$ em γ-" + nucleus.name + " (" + dipolemodel + ")");

        std::string out =
            "out/plots/sigma/" + meson + "/" +
            meson + "_" +
            dipolemodel + "_" + "gamma" + nucleus.name + "_" +
            timestamp() + ".pdf";

        std::filesystem::create_directories("out/plots/sigma/" + meson);
        std::cout << "Salvando figura em: " << out << std::endl;
        plt::save(out);
        plt::show();
}




void plot_sigma(std::string meson, std::string csv, std::string dipolemodel)
{
    if(meson == "psi")
        plot_sigma_Jpsi(csv, dipolemodel);
    else if(meson == "phi")
        plot_sigma_phi(csv, dipolemodel);
    else if(meson == "rho")
        plot_sigma_rho(csv, dipolemodel);
    else
        plot_sigma_placeholder(csv, dipolemodel, meson);
}

void plot_sigma_Jpsi_gammaPb(std::string csv_file, std::string dipolemodel)
{
    std::vector<double> W, sigma_GLC, sigma_BG;
    read_csv(csv_file, W, sigma_GLC, sigma_BG);

    plt::clf();

    // =========================
    // Dados experimentais
    // =========================
    std::vector<double> W_CMS, sigma_CMS, error_CMS;
    read_xyscan_csv("out/csv/sigma/psi/data/sigma_psi_gammaPb_CMS.csv",
                    W_CMS, sigma_CMS, error_CMS);

    plt::errorbar(W_CMS, sigma_CMS, error_CMS,
        {{"fmt","o"}, {"color","blue"}, {"label","CMS"}});// CMS Collaboration. Probing Small Bjorken-x Nuclear Gluonic Structure via
                                                        //Coherent J/ψ Photoproduction in Ultraperipheral Pb-Pb Collisions at √sNN = 5.02 TeV.
                                                        //Phys. Rev. Lett., American Physical Society, 2023.

    std::vector<double> W_ALICE_2023, sigma_ALICE_2023, error_ALICE_2023;
    read_xyscan_csv("out/csv/sigma/psi/data/sigma_psi_gammaPb_ALICE_2023.csv",  
                    W_ALICE_2023, sigma_ALICE_2023, error_ALICE_2023);

    plt::errorbar(W_ALICE_2023, sigma_ALICE_2023, error_ALICE_2023,      //ALICE Collaboration. Energy dependence of coherent photonuclear production of
                                                                        //J/ψ mesons in ultra-peripheral Pb-Pb collisions at √sNN = 5.02 TeV. 2023. 
        {{"fmt","s"}, {"color","green"}, {"label","ALICE(2023)"}});

    std::vector<double> W_ALICE_2024, sigma_ALICE_2024, error_ALICE_2024;
    read_xyscan_csv("out/csv/sigma/psi/data/sigma_psi_gammaPb_ALICE_2024.csv",
                    W_ALICE_2024, sigma_ALICE_2024, error_ALICE_2024);  // SHATAT, A. Charmonium photoproduction in Pb-Pb collisions with nuclear overlap
                                                                        //measured with ALICE at the LHC. Tese (Theses) — Université Paris-Saclay, 2024. https://theses.hal.science/tel-04797642
    
    plt::errorbar(W_ALICE_2024, sigma_ALICE_2024, error_ALICE_2024,
        {{"fmt","^"}, {"color","orange"}, {"label","ALICE(2024)"}});

    // =========================
    // Banda sombreada
    // =========================

    std::vector<double> sigma_min(W.size()), sigma_max(W.size());

for (size_t i = 0; i < W.size(); ++i) {
    sigma_min[i] = std::min(sigma_GLC[i], sigma_BG[i]);
    sigma_max[i] = std::max(sigma_GLC[i], sigma_BG[i]);
}

// converter para string Python
std::string W_py   = vec_to_pylist(W);
std::string min_py = vec_to_pylist(sigma_min);
std::string max_py = vec_to_pylist(sigma_max);

// plot banda
std::string cmd =
    "import matplotlib.pyplot as plt\n"
    "W = " + W_py + "\n"
    "y1 = " + min_py + "\n"
    "y2 = " + max_py + "\n"
    "plt.fill_between(W, y1, y2, alpha=0.3, color='orange', label='" + dipolemodel + " band')\n";

PyRun_SimpleString(cmd.c_str());


    // =========================
    // Estética
    // =========================

    plt::xlabel("W (GeV)");
    plt::ylabel("σ(γ Pb → J/ψ Pb) [nb]");
    plt::title("Produção exclusiva de J/ψ em γ-Pb (" + dipolemodel + ")");

    plt::xlim(25,1000);
    plt::ylim(1e3,1e5);

    plt::legend();

    PyRun_SimpleString(
        "ax = plt.gca()\n"
        "ax.set_xscale('log')\n"
        "ax.set_yscale('log')\n"
        "ax.grid(True)\n"
    );

    std::filesystem::create_directories("out/plots/sigma/psi");
    plt::save("out/plots/sigma/psi/Jpsi_gammaPb_" + dipolemodel + "_" + timestamp() + ".pdf");

    plt::show();
}

void plot_sigma_phi_gammaPb(std::string csv_file, std::string dipolemodel)
{
    std::vector<double> W, sigma_GLC, sigma_BG;

    read_csv(csv_file, W, sigma_GLC, sigma_BG);

    plt::clf();
    plt::figure_size(800,600);

    plt::named_plot(dipolemodel + " (GLC)", W, sigma_GLC, "-");
    plt::named_plot(dipolemodel + " (Boosted Gaussian)", W, sigma_BG, "--");

    plt::xlabel("W (GeV)");
    plt::ylabel("σ(γ Pb → φ Pb) [nb]");
    plt::title("Produção exclusiva de φ em γ-Pb (" + dipolemodel + ")");

    plt::legend();

    PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gcf()\n"   // garante que usa a figura atual
    "plt.xscale('log')\n"
    "plt.yscale('log')\n"
    "plt.grid(True)\n"
);

    std::filesystem::create_directories("out/plots/sigma/phi");
    plt::save("out/plots/sigma/phi/phi_gammaPb_" + dipolemodel + "_" + timestamp() + ".pdf");

    plt::show();
}


void plot_sigma_gammaA(std::string meson, std::string csv, 
                        std::string modelo, Nucleus& nucleus)
{
    if(is_proton_target(nucleus)){plot_sigma(meson, csv, modelo);}
    else{
    if (meson == "psi" and nucleus.name == "Pb"){
        plot_sigma_Jpsi_gammaPb(csv, modelo);
    } else if (meson == "phi" and nucleus.name == "Pb"){
        plot_sigma_phi_gammaPb(csv, modelo);
    } else {
        plot_sigma_placeholder(csv, modelo, meson, nucleus);
    }
    }
    
}
void plot_sigma_models(const std::vector<std::pair<std::string,std::string>>& files, const Meson& M,
                   bool fc)
{
    std::string pyfile =
        "out/plots/sigma/tmp_plot_sigma.py";

    std::ofstream py(pyfile);

    

    if(!py.is_open())
    {
        std::cerr << "Erro ao criar script Python.\n";
        return;
    }

    py << "import matplotlib.pyplot as plt\n";
    py << "import csv\n\n";
    letra_slide(py);

    py << "plt.figure(figsize=(9,6))\n\n";
if(M.meson == "phi"){
    std::vector<int> dataset_fixedpoint;

    std::vector<double>
        W_fixedpoint,
        sigma_fixedpoint,
        error_fixedpoint;

    read_sigma_exp(
        "out/csv/sigma/phi/data/phi_fixedpoint_data(nb).csv",
        dataset_fixedpoint,
        W_fixedpoint,
        sigma_fixedpoint,
        error_fixedpoint
    );

    std::vector<int> dataset_ZEUS;

    std::vector<double>
        W_ZEUS,
        sigma_ZEUS,
        error_ZEUS;

    read_sigma_exp(
        "out/csv/sigma/phi/data/phi_sigma_expdata_ZEUS(1994).csv",
        dataset_ZEUS,
        W_ZEUS,
        sigma_ZEUS,
        error_ZEUS
    );
     write_python_vector(py,"W_fixedpoint",W_fixedpoint);
    write_python_vector(py,"sigma_fixedpoint",sigma_fixedpoint);
    write_python_vector(py,"err_fixedpoint",error_fixedpoint);

    py << "plt.errorbar("
          "W_fixedpoint,"
          "sigma_fixedpoint,"
          "yerr=err_fixedpoint,"
          "fmt='o',"
          "capsize=3,"
          "color='black',"
          "label='Fixed target')\n\n";

    write_python_vector(py,"W_ZEUS",W_ZEUS);
    write_python_vector(py,"sigma_ZEUS",sigma_ZEUS);
    write_python_vector(py,"err_ZEUS",error_ZEUS);

    py << "plt.errorbar("
          "W_ZEUS,"
          "sigma_ZEUS,"
          "yerr=err_ZEUS,"
          "fmt='s',"
          "capsize=3,"
          "color='blue',"
          "label='ZEUS (1994)')\n\n";
    } else if(M.meson == "psi"){ {
std::vector<int> dataset;
    std::vector<double> W_exp, sigma_exp, err_exp;

    read_sigma_exp(
        "out/csv/sigma/psi/data/sigma_gammap_jpsi.csv",
        dataset,
        W_exp,
        sigma_exp,
        err_exp);

    std::map<int,std::string> exp_map = {
        {0,"H1"},
        {1,"H1"},
        {2,"ALICE"},
        {3,"LHCb"}
    };

    std::map<std::string,std::string> colors = {
        {"H1","blue"},
        {"ALICE","black"},
        {"LHCb","purple"}
    };

    std::map<std::string,std::string> markers = {
        {"H1","o"},
        {"ALICE","s"},
        {"LHCb","^"}
    };
    for(const auto& exp : {"H1","ALICE","LHCb"})
    {
        py << "W=[]\n";
        py << "S=[]\n";
        py << "E=[]\n";

        for(size_t i=0;i<W_exp.size();++i)
        {
            if(exp_map[dataset[i]] == exp)
            {
                py << "W.append(" << W_exp[i] << ")\n";
                py << "S.append(" << sigma_exp[i] << ")\n";
                py << "E.append(" << err_exp[i] << ")\n";
            }
        }

        py << "plt.errorbar("
              "W,S,"
              "yerr=E,"
              "fmt='" << markers[exp] << "',"
              "color='" << colors[exp] << "',"
              "capsize=3,"
              "label='" << exp << "')\n\n";
    };}
    } else if(M.meson == "rho"){
        std::vector<int> dataset_ALICE;

        std::vector<double>
            W_ALICE,
            sigma_ALICE,
            error_ALICE;

            read_xyscan_csv("out/csv/sigma/rho/data/rho_photoproduction_ZEUS(abramovic2000).csv",
                            W_ALICE,
                            sigma_ALICE,
                            error_ALICE);

    write_python_vector(py,"W_ALICE",W_ALICE);
    write_python_vector(py,"sigma_ALICE",sigma_ALICE);
    write_python_vector(py,"err_ALICE",error_ALICE);

    py << "plt.errorbar("
          "W_ALICE,"
          "sigma_ALICE,"
          "yerr=err_ALICE,"
          "fmt='^',"
          "capsize=3,"
          "color='black',"
          "label='ZEUS')\n\n";                
        } else {
        std::cerr << "Meson sem dados experimentais: " << M.meson << std::endl;
        return;
    }
    for(const auto& [fname,label] : files)
{
    py << "W=[]\n";
    py << "sigma=[]\n";

    py << "with open(r'" << fname << "') as f:\n";
    py << "    reader = csv.reader(f)\n";
    py << "    next(reader)\n";
    py << "    for row in reader:\n";
    py << "        W.append(float(row[0]))\n";
    py << "        sigma.append(float(row[1]))\n";
    py << "print('" << label << "', len(W))\n";
    py << "plt.plot("
          "W,"
          "sigma,"
          "lw=2,"
          "label=r'" << label << "')\n\n";
}

    py << "plt.xlabel(r'$W$ (GeV)')\n";
    py << "plt.ylabel(r'$\\sigma$ (nb)')\n";

    py << "plt.title(r'Seção de choque de $\\gamma p \\to \\" << M.meson << " p$')\n";

    py << "plt.xscale('log')\n";
    py << "plt.yscale('log')\n";

    py << "plt.xlim(8.0,3000.0)\n";
    //py << "plt.ylim(,1.2)\n";

    py << "plt.grid(True)\n";
    py << "plt.legend()\n";

    std::string filename =
        "out/plots/sigma/"
        + M.meson+"/" + M.nome +
        + "_sigma_models_"
        + timestamp()
        + ".pdf";

    py << "plt.savefig(r'" << filename << "')\n";
    py << "plt.show()\n";

    py.close();

    std::cout << "Executando script Python...\n";

    std::string cmd = "python3 \"" + pyfile + "\"";

    int ret = std::system(cmd.c_str());

    if(ret != 0)
        std::cerr << "Erro ao executar Python.\n";
    else
        std::cout << "Figura salva em: "
                  << filename
                  << std::endl;
}
