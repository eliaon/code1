#include "plot_io.hpp"
#include "../../libraries/matplotlib-cpp/matplotlibcpp.h"
#include "../utils.hpp"
#include "../../calculations/nuclear.hpp"
#include "../../dipole_amplitudes/LHAPDF.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <utility>
#include <vector>


namespace plt = matplotlibcpp;

namespace {
std::string py_string(const std::string& value)
{
    std::string escaped = "'";
    for (char c : value) {
        if (c == '\\' || c == '\'') escaped += '\\';
        escaped += c;
    }
    escaped += "'";
    return escaped;
}
}

void plot_N_multi(
    const std::vector<std::string>& filenames,
    const std::vector<double>& x_vals,
    const std::string& modelo
)
{
    if (!x_vals.empty() && x_vals.size() != filenames.size()) {
        throw std::runtime_error("x_vals e filenames devem ter o mesmo tamanho");
    }

    plt::figure();

    std::ostringstream py;
    py << "import matplotlib.pyplot as plt\n"
       << "import numpy as np\n"
       << "files=[";

    for (size_t i = 0; i < filenames.size(); ++i) {
        if (i > 0) py << ",";
        py << py_string(filenames[i]);
    }

    py << "]\n"
       << "labels=[";

    for (size_t i = 0; i < filenames.size(); ++i) {
        if (i > 0) py << ",";
        if (x_vals.empty())
            py << py_string("");
        else
            py << py_string("x=" + doubleParaString(x_vals[i]));
    }

    py << "]\n"
       << "for f, label in zip(files, labels):\n"
       << "    data = np.genfromtxt(f, delimiter=',', names=True)\n"
       << "    plt.plot(data[data.dtype.names[0]], data[data.dtype.names[1]], label=label)\n"
       << "plt.xscale('log')\n"
       << "plt.yscale('log')\n";

    for (size_t i = 0; i < filenames.size(); ++i)
    {
        std::vector<double> r, N;

        read_two_columns(filenames[i], r, N);
    }

    PyRun_SimpleString(py.str().c_str());

    plt::xlim(1e-3, 10.0);
    plt::ylim(1e-4, 1.2);

    plt::xlabel("r [fm]");
    plt::ylabel("N(r)");

    if (!x_vals.empty())
        plt::legend();

    plt::grid(true);
    plt::title("Função de dipolo N(r) - " + modelo);
    

    plt::save("out/plots/N/N_"+modelo+"_" + timestamp() + ".pdf");
    plt::show();
    plt::close();
}

void plot_N_dglap(std::string csv_file1, std::string csv_file2)
{
    std::vector<double> r2_x1, N_dglap_x1;
    std::vector<double> r2_x2, N_dglap_x2;

    read_two_columns(csv_file1, r2_x1, N_dglap_x1);
    read_two_columns(csv_file2, r2_x2, N_dglap_x2);


    plt::figure_size(800,600);
    plt::plot(r2_x1, N_dglap_x1, {{"label","N(r^2) DGLAP x=1e-4"}});
    plt::plot(r2_x2, N_dglap_x2, {{"label","N(r^2) DGLAP x=1e-2"}});

    PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gcf()\n"   // garante que usa a figura atual
    "plt.xscale('log')\n"
    "plt.yscale('log')\n"
    );

    plt::xlim(1e-3,100.0);
    plt::ylim(1e-4,1.2);

    plt::xlabel("$r^2$");
    plt::ylabel("$N_p$");

    plt::title("Função de dipolo N(r^2) do modelo DGLAP");

    plt::grid(true);
    plt::legend();

    std::string out =
        "out/plots/N/N_dglap_" +
        timestamp() + ".pdf";

    plt::save(out);
    plt::show();
}

// -------------- PLOTA AS CURVAS DE N

void plot_N_models(const std::vector<std::pair<std::string,std::string>>& files, double x)
{
    plt::figure_size(700,500);

    for (const auto& [fname, label] : files)
    {
        std::vector<double> r2, N;

        std::ifstream fin(fname);
        if (!fin.is_open()) {
            std::cerr << "Erro ao abrir: " << fname << std::endl;
            continue;
        }

        std::string line;
        std::getline(fin, line); // pula header

        while (std::getline(fin, line))
        {
            std::stringstream ss(line);
            std::string val1, val2;

            std::getline(ss, val1, ',');
            std::getline(ss, val2, ',');

            r2.push_back(std::stod(val1));
            N.push_back(std::stod(val2));
        }

        plt::plot(r2, N, {{"label", label}});
    }

    plt::xlabel("$r$");
    plt::ylabel("$N_p$");

    plt::title("Amplitude de dipolo $N_p(r)$ para diferentes modelos ");

    plt::xlim(0.001,10.0);
    plt::ylim(0.0001,1.2);

    PyRun_SimpleString( "import matplotlib.pyplot as plt\n"
                        "plt.yscale('log')\n"
                        "plt.xscale('log')\n"

    );

    plt::grid(true);
    plt::legend();

    std::string filename =
        "out/plots/N/N_models_" + timestamp() + "_x=" + doubleParaString(x) + ".pdf";

    plt::save(filename);
    plt::show();
}

void compare_N_models(double x)
{
    std::vector<std::pair<std::string,std::string>> files = {
        {N_file(x, "IPSAT"), "IPsat"},
        {N_file(x, "bCGC"), "bCGC"},
        {N_file(x, "GBW(new)"), "GBW(new)"}
    };

    plot_N_models(files, x);
}

void plot_xf_multiQ2()
{
    LHAPDF::setPaths("/home/elian/local/share/LHAPDF");


    std::string pdf_name = "NNPDF30_nlo_as_0118";

    LHAnPDF lhapdf("NNPDF30_nlo_as_0118");
    const std::vector<double> Q2_values = {4.0};

    namespace fs = std::filesystem;
    fs::create_directories("out/plots/xf");

    std::vector<int> flavors = {21}; //  // 21 para xg, 1: xd, 2: xu, 3: xs, 4: xc, 5: xb, 6: xt
    
std::string pyfile = "out/plots/xf/plot_xf.py";
std::ofstream py(pyfile);

py << "import matplotlib.pyplot as plt\n";
py << "import numpy as np\n\n";

for (int flavor : flavors)
{
    std::string sflavor = flavorName(flavor);
    std::string xf = flavorToString(flavor);

    py << "plt.figure()\n";

    for (double Q2 : Q2_values)
    {
        std::string fname = lhapdf.xf_vs_x(Q2, flavor);

        py << "data = np.loadtxt('" << fname << "', delimiter=',', skiprows=1)\n";
        py << "x = data[:,0]\n";
        py << "xf = data[:,1]\n";
        py << "plt.plot(x, xf, label='Q²=" << Q2 << " GeV²')\n\n";
    }

    py << "plt.xscale('log')\n";
    py << "plt.yscale('log')\n";
    py << "plt.xlabel('x')\n";
    py << "plt.ylabel('" << xf << "(x,Q²)')\n";
    py << "plt.xlim(1e-5, 1.0)\n";
    py << "plt.title('" << sflavor << " PDF - " << pdf_name << "')\n";
    py << "plt.legend()\n";
    py << "plt.grid(True, which='both', ls='--')\n";
    py << "plt.show()\n";
    py << "plt.savefig('out/plots/xf/" << xf << "_multiQ2_" << pdf_name << "_" << timestamp() << ".pdf', dpi=300)\n";
    py << "plt.close()\n\n";  // fecha figura corretamente
}

py.close(); 

int ret = std::system(("python3 " + pyfile).c_str());
    
}

void plot_TA_b(){

    TA_Table table = precompute_TA(200, 50.0, Pb208); // pré-calcula T_A(b) para Pb-208

    std::vector<double> b = table.b_vals;
    std::vector<double> TA = table.TA_vals;

    plt::figure_size(800,600);
    plt::plot(b, TA, {{"label","T_A(b)"}});
    plt::xlabel("b (GeV^-1)");
    plt::ylabel("T_A(b) (GeV^2)");
    plt::title("Função de perfil nuclear T_A(b) para Pb-208");
    plt::grid(true);
    plt::legend();
    plt::save("out/plots/other/TA_b(Pb-208)_" + timestamp() + ".pdf");
    plt::show();
    
}




// ------------------ PĹOT OVERLAPS

void plot_Jpsi_overlap()
{

    std::string csv_file_psi = "out/csv/overlap/Jpsi_overlap_r.csv";
    std::vector<double> r_psi, overlap_GLC_psi, overlap_BG_psi;

    read_csv(csv_file_psi, r_psi, overlap_GLC_psi, overlap_BG_psi);

    plt::figure_size(1000,600);

    plt::plot(r_psi, overlap_GLC_psi,
        {{"label","GLC"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.8"}});

    plt::plot(r_psi, overlap_BG_psi,
        {{"label","BG"},
         {"color","black"},
         {"linestyle","-"},
         {"linewidth","1.2"}});

    std::string csv_file_Y = "out/csv/overlap/Upsilon_overlap_r.csv";
    std::vector<double> r_Y, overlap_GLC_Y, overlap_BG_Y;

    read_csv(csv_file_Y, r_Y, overlap_GLC_Y, overlap_BG_Y);

    plt::plot(r_Y, overlap_GLC_Y,
        {{"label","GLC (Upsilon)"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.8"}});
    
    plt::plot(r_Y, overlap_BG_Y,
        {{"label","BG (Upsilon)"},
         {"color","black"},
         {"linestyle","-"},
         {"linewidth","1.2"}});

    plt::xlim(0.0001,1.0);
    plt::ylim(0.0,0.025);

    plt::xlabel("$r$ [fm]");
    plt::ylabel("Overlap $r \\Psi_V \\Psi_{\\gamma}$");

    // escala log em r
    PyRun_SimpleString(
        "import matplotlib.pyplot as plt\n"
        "plt.xscale('log')\n"
    );

    plt::title("Sobreposição de funções de onda de $\\psi$ e $\\Upsilon$");

    plt::grid(true);
    plt::legend();

    std::string out =
        "out/plots/overlap/Jpsi_Upsilon_overlap_" +
        timestamp() + ".pdf";

    plt::save(out);
    plt::show();
}

void plot_phi_overlap()
{

    std::string csv_file_phi = "out/csv/overlap/phi_overlap_r.csv";
    std::vector<double> r_phi, overlap_GLC_phi, overlap_BG_phi;

    read_csv(csv_file_phi, r_phi, overlap_GLC_phi, overlap_BG_phi);

    plt::figure_size(1000,600);

    plt::plot(r_phi, overlap_GLC_phi,
        {{"label","GLC (phi)"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.8"}});

    plt::plot(r_phi, overlap_BG_phi,
        {{"label","BG (phi)"},
         {"color","black"},
         {"linestyle","-"},
         {"linewidth","1.2"}});

    std::string csv_file_rho = "out/csv/overlap/rho_overlap_r.csv";
    std::vector<double> r_rho, overlap_GLC_rho, overlap_BG_rho;

    read_csv(csv_file_rho, r_rho, overlap_GLC_rho, overlap_BG_rho);
    plt::plot(r_rho, overlap_GLC_rho,
        {{"label","GLC (rho)"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.8"}});
    plt::plot(r_rho, overlap_BG_rho,
        {{"label","BG (rho)"},
         {"color","black"},
         {"linestyle","-"},
         {"linewidth","1.2"}});
    
    plt::xlim(0.01,3.0);
    plt::ylim(0.0,0.01);

    plt::xlabel("$r$ [fm]");
    plt::ylabel("Overlap $r \\phi_V \\phi_{\\gamma}$");

    // escala log em r
    PyRun_SimpleString(
        "import matplotlib.pyplot as plt\n"
        "plt.xscale('log')\n"
    );

    plt::title("Sobreposição de funções de onda de $\\phi$ e $\\rho$");

    plt::grid(true);
    plt::legend();

    std::string out =
        "out/plots/overlap/phi-rho_overlap_" +
        timestamp() + ".pdf";

    plt::save(out);
    plt::show();
}

void plot_overlap( std::string meson)
{
    if(meson == "Jpsi" || meson == "Upsilon")
        plot_Jpsi_overlap();

    else if(meson == "phi" || meson == "rho")
        plot_phi_overlap();

    else
        throw std::runtime_error("Meson desconhecido: " + meson);
}


void plot_overlap_fc(std::string csv_file, std::string csv_file_fc, std::string meson)
{
    std::vector<double> r1, r2;
    std::vector<double> overlap_BG, overlap_GLC;
    std::vector<double> overlap_BG_fc, overlap_GLC_fc;

    read_csv(csv_file, r1, overlap_BG, overlap_GLC);
    read_csv(csv_file_fc, r2, overlap_BG_fc, overlap_GLC_fc);

    if(r1.empty() || r2.empty()) {
    throw std::runtime_error("Vetores vazios — falha no read_csv");
}
if(meson == "Jpsi") {
    // -------- SEM FC --------

plt::plot(r1, overlap_BG, {{"label","BG"}});
plt::plot(r1, overlap_GLC, {{"label","GLC"}});

PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gca().set_xscale('log')\n"
);

plt::title("Sem fator de correção");
plt::xlabel("$r$ [fm]");
plt::ylabel("Overlap $r \\Psi_V \\Psi_{\\gamma}$");


plt::xlim(0.0001, 1.0);
plt::ylim(0.0, 0.025);

plt::legend();
plt::grid(true);



plt::save("out/plots/overlap/sem_fc_" + meson + "_" + timestamp() + ".pdf");

// -------- COM FC --------
plt::figure();

plt::plot(r2, overlap_BG_fc, {{"label","BG (fc)"}});
plt::plot(r2, overlap_GLC_fc, {{"label","GLC (fc)"}});

PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gca().set_xscale('log')\n"
);

plt::title("Com fator de correção");
plt::xlabel("$r$ [fm]");
plt::ylabel("Overlap $r \\Psi_V \\Psi_{\\gamma}$");
plt::legend();
plt::grid(true);

plt::xlim(0.0001, 1.0);
plt::ylim(0.0, 0.025);


plt::save("out/plots/overlap/com_fc_" + meson + "_" + timestamp() + ".pdf");

plt::show();
}
else if(meson == "phi") {
    // -------- SEM FC --------
plt::figure();

plt::plot(r1, overlap_BG, {{"label","BG"}});
plt::plot(r1, overlap_GLC, {{"label","GLC"}});
PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gca().set_xscale('log')\n"
);

plt::title("Sem fator de correção");
plt::xlabel("$r$ [fm]");
plt::xlim(0.01, 3.0);
plt::ylim(0.0, 0.01);
plt::ylabel("Overlap $r \\Psi_V \\Psi_{\\gamma}$");


plt::legend();
plt::grid(true);

plt::save("out/plots/overlap/sem_fc_" + meson + "_" + timestamp() + ".pdf");

// -------- COM FC --------
plt::figure();

plt::plot(r2, overlap_BG_fc, {{"label","BG (fc)"}});
plt::plot(r2, overlap_GLC_fc, {{"label","GLC (fc)"}});

PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gca().set_xscale('log')\n"
);
plt::title("Com fator de correção");
plt::xlabel("$r$ [fm]");
plt::ylabel("Overlap $r \\Psi_V \\Psi_{\\gamma}$");
plt::xlim(0.01, 3.0);
plt::ylim(0.0, 0.01);


plt::legend();
plt::grid(true);

plt::save("out/plots/overlap/com_fc_" + meson + "_" + timestamp() + ".pdf");

plt::show();
}
}



void plot_XY(std::vector<double> x, std::vector<double> y, std::string xlabel, std::string ylabel, std::string title)
{
    std::cout << "plot_XY: x.size() = " << x.size() << ", y.size() = " << y.size() << std::endl;
    plt::figure_size(800,600);
    plt::plot(x, y, {{"label","Dados"}});

    std::cout << "Iniciando pyrun...\n";
    PyRun_SimpleString(
        "import matplotlib.pyplot as plt\n"
        "plt.gca().set_xscale('log')\n"
    );

    plt::xlabel(xlabel);
    plt::ylabel(ylabel);
    plt::title(title);
    plt::legend();
    plt::grid(true);

    std::string out =
        "out/plots/other/" + title + "_" + timestamp() + ".pdf";

    plt::save(out);
    plt::show();
}

