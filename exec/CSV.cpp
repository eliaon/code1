#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <filesystem>

#include "../calculations/sigma.hpp"
#include "../other/utils.h"
#include "../other/plot.h"
#include "../calculations/nuclear.hpp"

namespace CSV {

    void N_p(std::string model){
        std::vector<double> x = {1e-2, 1e-3, 1e-4, 1e-5};
        std::vector<std::string> filenames;

        

        for (double x_val : x) {
            std::string filename = N_file(x_val, model);
            
            filenames.push_back(filename);
        }
    plot_N_multi(filenames, x, model);
    }

    void sigma_gamma_p(std::string model)
{
    double Q2 = 0.0;

    const Meson& M_GLC = input_meson("GBW");
    const Meson& M_BG  = meson_modelsGBW.find(M_GLC.meson)->second.M_BG;

    std::cout << "Calculando sigma_p para " << M_GLC.meson
              << " com modelo " << model << "\n";

    std::filesystem::create_directories("out/csv/sigma/" + M_GLC.meson);

    std::string basepath = "out/csv/sigma/" + M_GLC.meson + "/" +
                           M_GLC.meson + "_sigma_" + model;

    std::string filename_csv = basepath + ".csv";
    std::string filename_dat = basepath + ".dat";

    std::ofstream fout_csv(filename_csv);
    std::ofstream fout_dat(filename_dat);

    if (!fout_csv.is_open() || !fout_dat.is_open()) {
        std::cerr << "Erro ao abrir arquivos de saída.\n";
        return;
    }

    // cabeçalho csv
    fout_csv << "W,sigma_GLC,sigma_BG\n";

    // comentário opcional no .dat (xmgrace ignora linhas com #)
    fout_dat << "# W sigma_GLC sigma_BG\n";

    const int Nw = 10;
    double Wmin = x_to_W(1e-2, M_GLC);
    double Wmax = 3e3;

    using clock = std::chrono::steady_clock;
    auto start = clock::now();

    std::vector<double> Wv(Nw), sigma_GLC_v(Nw), sigma_BG_v(Nw);

    std::cout << "W (GeV), sigma_GLC (nb), sigma_BG (nb)\n";

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < Nw; ++i) {
        double frac = static_cast<double>(i) / (Nw - 1);
        double W = Wmin * std::pow(Wmax / Wmin, frac);

        Wv[i] = W;
        sigma_GLC_v[i] = gamma_p::sigma(W, Q2, M_GLC, model) * GeV2_to_nb;
        sigma_BG_v[i]  = gamma_p::sigma(W, Q2, M_BG, model) * GeV2_to_nb;

        std::cout << Wv[i] << " "
                  << sigma_GLC_v[i] << " "
                  << sigma_BG_v[i] << "\n";
    }

    for (int i = 0; i < Nw; ++i) {
        // CSV
        fout_csv << Wv[i] << ","
                 << sigma_GLC_v[i] << ","
                 << sigma_BG_v[i] << "\n";

        // DAT para xmgrace
        fout_dat << Wv[i] << " "
                 << sigma_GLC_v[i] << " "
                 << sigma_BG_v[i] << "\n";
    }

    fout_csv.close();
    fout_dat.close();

    auto end = clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Tempo de execução: " << duration << " ms\n";
    std::cout << "Arquivo CSV: " << filename_csv << "\n";
    std::cout << "Arquivo DAT: " << filename_dat << "\n";

    plot_sigma(M_GLC.meson, filename_csv, model);
}

    void sigma_gamma_A(std::string model, nucleous&Nucleo)
    {
    double Q2=0.0;

    const Meson& M_GLC = input_meson("GBW"); //importa o glc do meson escolhido
    const Meson& M_BG  = meson_modelsGBW.find(M_GLC.meson)->second.M_BG; //pega o bg correspondente ao meson escolhido

    std::string Q2_str = std::format("{:.3g}", Q2);
    std::string path = "out/csv/sigma/" + M_GLC.meson + "/"+M_GLC.meson + "_sigma-gamma"+
                            Nucleo.name + "_" + model ;
    std::string filename_csv = path + ".csv";
    std::string filename_dat = path + ".dat";
    std::ofstream fout_csv(filename_csv);
    std::ofstream fout_dat(filename_dat);
    // cabeçalho csv
    fout_csv << "W,sigma_GLC,sigma_BG\n";

    // comentário opcional no .dat (xmgrace ignora linhas com #)
    fout_dat << "# W sigma_GLC sigma_BG\n";

    const int Nw = 100;
    double Wmin = x_to_W(1e-2, M_GLC);
    double Wmax = 3e3;

    using clock = std::chrono::steady_clock;
    auto start = clock::now();
    
    TA_Table table = precompute_TA(200, 50.0, Nucleo); // pré-calcula T_A(b) para Pb-208

    std::vector<double> Wv(Nw), sigma_GLC_v(Nw), sigma_BG_v(Nw);


    std::cout << "W (Gev), sigma_GLC (nb), sigma_BG (nb)\n";
    #pragma omp parallel for schedule(dynamic)
for (int i = 0; i < Nw; ++i) {

    double frac = static_cast<double>(i) / (Nw - 1);
    double W = Wmin * pow(Wmax / Wmin, frac);
    double x = (M_GLC.MV * M_GLC.MV) / (W * W);

    Wv[i] = W;

    sigma_GLC_v[i] = gamma_A::sigma(W, Q2, M_GLC, model, table) * GeV2_to_nb;
    sigma_BG_v[i]  = gamma_A::sigma(W, Q2, M_BG, model, table) * GeV2_to_nb;
    std::cout << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << "\n";
}
//for (int i = 0; i < Nw; ++i) {
//    std::cout << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << "\n";
//}
    for (int i = 0; i < Nw; ++i) {
    fout_csv << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << "\n";
    fout_dat << Wv[i] << " " << sigma_GLC_v[i] << " " << sigma_BG_v[i] << "\n";
}

    fout_csv.close();
    fout_dat.close();

    auto end = clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Tempo de execução: " << duration << " ms" << std::endl;
    std::cout << "Arquivo '" << path << "'(csv e dat) gerado em " << duration << " ms" << std::endl;
    plot_sigma_gammaA(M_GLC.meson, filename_csv, model);
}


    void rapidez_AA(double sqrt_s, std::string model, nucleous &Nucleo1, nucleous& Nucleo2)
{
    //double sqrt_s = 2.76e3; // GeV

    const Meson& M_GLC = input_meson("GBW");
    const Meson& M_BG  = meson_modelsGBW.find(M_GLC.meson)->second.M_BG;

    std::string sqrt_s_str = std::format("{:.3g}", sqrt_s);
    std::string path = "out/csv/Rapidez/" + M_GLC.meson + "/"+M_GLC.meson + "_Rapidez_"+
                            Nucleo1.name + "-" + Nucleo2.name +"_" + model ;

    std::string filename_csv = path + ".csv";
    std::string filename_dat = path + ".dat";

    double Q2 = 0.0;


    TA_Table tableA = get_TA_table(Nucleo1); // pré-calcula T_A(b) para o núcleo A

    TA_Table tableB =
    (Nucleo1.name == Nucleo2.name)? tableA: get_TA_table(Nucleo2);
    

    const int Ny = 100;
    double ymax = std::log(sqrt_s / M_GLC.MV);

    int Ntot = 2*Ny + 1;

    std::vector<double> Y(Ntot), GLC(Ntot), BG(Ntot);

    #pragma omp parallel for schedule(dynamic)
for (int i = 0; i <= 2*Ny; ++i)
{
    int k = i - Ny;
    double y = (double(k) / Ny) * ymax;

    double dsdy_GLC, dsdy_BG;

    if (Nucleo1.name == Nucleo2.name)
    {
        dsdy_GLC = gamma_A::d_sigma_dy_AA(y, sqrt_s, Q2, M_GLC, model, tableA, Nucleo1) * gev2_to_mb;
        dsdy_BG  = gamma_A::d_sigma_dy_AA(y, sqrt_s, Q2, M_BG,  model, tableA, Nucleo1) * gev2_to_mb;
    }
    else
    {
        dsdy_GLC = gamma_A::d_sigma_dy_AB(y, sqrt_s, Q2, M_GLC, model, tableA, tableB, Nucleo1, Nucleo2) * gev2_to_mb;
        dsdy_BG  = gamma_A::d_sigma_dy_AB(y, sqrt_s, Q2, M_BG,  model, tableA, tableB, Nucleo1, Nucleo2) * gev2_to_mb;
    }

    int idx = Ny + k;

    Y[idx] = y;
    GLC[idx] = dsdy_GLC;
    BG[idx] = dsdy_BG;
}

    // ✔ escrita segura
    std::ofstream fout_csv(filename_csv);
    fout_csv << "y,d_sigma_dy_GLC,d_sigma_dy_BG\n";
    std::ofstream fout_dat(filename_dat);
    fout_dat << "#y d_sigma_dy_GLC d_sigma_dy_BG\n";


    for (int i = 0; i < Ntot; ++i)
    {
        fout_csv << Y[i] << "," << GLC[i] << "," << BG[i] << "\n";

        fout_dat << Y[i] << " " << GLC[i] << " " << BG[i] << "\n";
    }

    fout_csv.close();
    fout_dat.close();

    std::cout << "Arquivo '" << path << "'(csv e dat) gerados.\n";

    plot_rapidez(filename_csv, model, M_GLC.meson, Nucleo1.name+Nucleo2.name, sqrt_s);

    std::cout << "Plot gerado para " << M_GLC.meson << " para a colisão " <<Nucleo1.name << Nucleo2.name 
                <<" a "<< sqrt_s_str <<" TeV.\n";
}
}