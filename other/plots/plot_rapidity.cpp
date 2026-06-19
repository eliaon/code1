#include "plot_io.hpp"
#include "../../libraries/matplotlib-cpp/matplotlibcpp.h"
#include "../utils.hpp"
#include <chrono>
#include <sstream>
#include <fstream>
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
    std::string pyfile =
        "out/plots/Rapidez/tmp_plot_rapidez.py";

    std::ofstream py(pyfile);

    if(!py.is_open())
    {
        std::cerr << "Erro ao criar script Python.\n";
        return;
    }

    //=========================================================
    // leitura dos dados teóricos
    //=========================================================

    std::vector<double> Y_th, rap_GLC, rap_BG;

    try
    {
        read_csv(csv,Y_th,rap_GLC,rap_BG);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        return;
    }

    //=========================================================
    // dados experimentais
    //=========================================================

    std::vector<double>
        Y_exp,
        dsigma_exp,
        error_exp;

    if(sqrt_s == 5.36e3)
    {
        read_rapidity_hepdata(
            "out/csv/Rapidez/phi/data/rapidity_PbPb-phi_5360TeV.csv",
            Y_exp,
            dsigma_exp,
            error_exp
        );
    }

    //=========================================================
    // desvios relativos
    //=========================================================

    std::vector<double>
        dev_GLC,
        dev_BG,
        Y_dev;

    for(size_t i=0;i<Y_exp.size();++i)
    {
        double Yd = Y_exp[i];

        double theo_GLC = 0.0;
        double theo_BG  = 0.0;

        double best = 1e99;

        for(size_t j=0;j<Y_th.size();++j)
        {
            double diff = std::abs(Y_th[j]-Yd);

            if(diff < best)
            {
                best = diff;

                theo_GLC = rap_GLC[j];
                theo_BG  = rap_BG[j];
            }
        }

        if(dsigma_exp[i] > 0.0)
        {
            Y_dev.push_back(Yd);

            dev_GLC.push_back(
                (theo_GLC-dsigma_exp[i])
                / dsigma_exp[i]
            );

            dev_BG.push_back(
                (theo_BG-dsigma_exp[i])
                / dsigma_exp[i]
            );
        }
    }

    //=========================================================
    // script Python
    //=========================================================

    py << "import matplotlib.pyplot as plt\n";
    py << "import numpy as np\n\n";

    write_python_vector(py,"Y_th",Y_th);
    write_python_vector(py,"rap_GLC",rap_GLC);
    write_python_vector(py,"rap_BG",rap_BG);

    write_python_vector(py,"Y_exp",Y_exp);
    write_python_vector(py,"dsigma_exp",dsigma_exp);
    write_python_vector(py,"error_exp",error_exp);

    write_python_vector(py,"Y_dev",Y_dev);
    write_python_vector(py,"dev_GLC",dev_GLC);
    write_python_vector(py,"dev_BG",dev_BG);

    py << "\n";

    py << "fig = plt.figure(figsize=(8,6))\n";

    py << "ax1 = plt.subplot2grid((2,1),(0,0))\n";

    py << "ax1.plot(Y_th,rap_GLC,"
          "color='orange',"
          "lw=1.8,"
          "label='GLC')\n";

    py << "ax1.plot(Y_th,rap_BG,"
          "color='orange',"
          "ls='--',"
          "lw=1.2,"
          "label='BG')\n";

    py << "if len(Y_exp)>0:\n";
    py << "    ax1.errorbar("
          "Y_exp,"
          "dsigma_exp,"
          "yerr=error_exp,"
          "fmt='o',"
          "color='blue',"
          "capsize=3,"
          "label='CMS (2025)')\n";

    py << "ax1.set_ylabel(r'$d\\sigma/dY$ [mb]')\n";

    std::stringstream title;

title
<< "Fotoprodução coerente de "
<< "$\\phi$ em PbPb "
<< "($\\sqrt{s_{NN}}="
<< sqrt_s/1e3
<< "\\,\\mathrm{TeV}$)"
<< " - "
<< model;

py << "ax1.set_title(r'" << title.str() << "')\n";

    py << "ax1.legend()\n";

    py << "ax1.set_yscale('log')\n";

    py << "ax1.grid(True,which='both',"
          "linestyle='--',alpha=0.6)\n";

    py << "ax1.set_xlim(0.0,1.5)\n";
    py << "ax1.set_ylim(10,200)\n";

    py << "ax1.tick_params(labelbottom=False)\n";

    py << "\n";

    py << "ax2 = plt.subplot2grid((2,1),(1,0))\n";

    py << "if len(Y_dev)>0:\n";
    py << "    ax2.plot("
          "Y_dev,"
          "dev_GLC,"
          "marker='o',"
          "linestyle='none',"
          "color='orange',"
          "label='GLC')\n";

    py << "    ax2.plot("
          "Y_dev,"
          "dev_BG,"
          "marker='x',"
          "linestyle='none',"
          "color='darkred',"
          "label='BG')\n";

    py << "ax2.axhline(0,"
          "color='black',"
          "linestyle='--')\n";

    py << "ax2.grid(True,"
          "linestyle='--',"
          "alpha=0.5)\n";

    py << "ax2.set_xlim(0.0,1.5)\n";

    py << "ax2.set_xlabel(r'$Y$')\n";

    py << "ax2.set_ylabel(r'$\\Delta$')\n";

    
    py << "fig.tight_layout()\n";
    py << "fig.subplots_adjust(hspace=0.03)\n";

    std::string out =
        "out/plots/Rapidez/phi/"
        + extrair_nome_base(csv)
        + "_"
        + timestamp()
        + ".pdf";

    py << "fig.text("
          "0.01,"
          "0.01,"
          "r'Arquivo: "
       << extrair_nome_base(csv)
       << "',"
          "fontsize=8,"
          "alpha=0.7)\n";

    py << "plt.savefig(r'" << out << "')\n";
    py << "plt.show()\n";

    py.close();

    std::string cmd =
        "python3 \"" + pyfile + "\"";

    int ret = std::system(cmd.c_str());

    if(ret != 0)
    {
        std::cerr << "Erro ao executar Python.\n";
    }
    else
    {
        std::cout
            << "Figura salva em: "
            << out
            << std::endl;
    }
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

void plot_rapidezAA_models(const std::vector<std::pair<std::string,std::string>>& files, double sqrt_s, const Meson& M,
                   bool fc)
{
    std::string pyfile =
        "out/plots/Rapidez/tmp_plot_rapidezAA.py";

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

    // =========================================================
    // dados experimentais
    // =========================================================
    std::vector<double>
        Y_exp,
        dsigma_exp,
        error_exp;


        fflush(stdout);
        read_rapidity_hepdata(
            "out/csv/Rapidez/phi/data/rapidity_PbPb-phi_5360TeV.csv",
            Y_exp,
            dsigma_exp,
            error_exp
        );
        fflush(stdout);
    
     write_python_vector(py,"Y_exp",Y_exp);
    write_python_vector(py,"dsigma_exp",dsigma_exp);
    write_python_vector(py,"error_exp",error_exp);

    py << "plt.errorbar("
          "Y_exp,"
          "dsigma_exp,"
          "yerr=error_exp,"
          "fmt='o',"
          "capsize=3,"
          "color='black',"
          "label='Fixed target')\n\n";
} else {
        std::cerr << "Meson sem dados experimentais: " << M.meson << std::endl;
        return;
    }
    for(const auto& [fname,label] : files)
{
    py << "Y=[]\n";
    py << "dsigmadY=[]\n";

    py << "with open(r'" << fname << "') as f:\n";
    py << "    reader = csv.reader(f)\n";
    py << "    next(reader)\n";
    py << "    for row in reader:\n";
    py << "        Y.append(float(row[0]))\n";
    py << "        dsigmadY.append(float(row[1]))\n";
    py << "print('" << label << "', len(Y))\n";
    py << "plt.plot("
          "Y,"
          "dsigmadY,"
          "lw=2,"
          "label=r'" << label << "')\n\n";
}

    py << "plt.xlabel(r'$Y$')\n";
    py << "plt.ylabel(r'$\\frac{d\\sigma}{dY}$ (mb)')\n";

    py << "plt.title(r'Distribuição de rapidez do $PbPb \\to \\" << M.meson << " PbPb$')\n";

    //py << "plt.xscale('log')\n";
    //py << "plt.yscale('log')\n";

    py << "plt.xlim(0.0,1.5)\n";
    py << "plt.ylim(10,130)\n";

    py << "plt.grid(True)\n";
    py << "plt.legend()\n";

    std::string filename =
        "out/plots/Rapidez/"
        + M.meson +"/" + M.nome +
        + "_Rapidez_models_"
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




















