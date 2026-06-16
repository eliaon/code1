#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <filesystem>

#include "../calculations/sigma.hpp"
#include "../other/utils.hpp"
#include "../calculations/nuclear.hpp"
#include "../other/correcs.hpp"

#include "../other/plots/plot_sigma.hpp"
#include "../other/plots/plot_rapidity.hpp"
#include "../other/plots/plot_other.hpp"



namespace CSV {
    namespace {
        bool is_proton_target(const Nucleus& nucleo)
        {
            return nucleo.name == "p" || nucleo.name == "proton";
        }
    }

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

    const Meson& M_GLC = input_meson(model);
    const Meson& M_BG  = get_meson_bg(M_GLC.meson, model);
    warn_ipsat_light_meson(M_GLC, model);

    std::cout << "Calculando sigma_p para " << M_GLC.meson
              << " com modelo " << model << "\n";

    std::filesystem::create_directories("out/csv/sigma/" + M_GLC.meson);

    std::string basepath = "out/csv/sigma/" + M_GLC.meson + "/" +
                           M_GLC.meson + "_sigma-gammap_" + model;

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

    using clock = std::chrono::steady_clock;
    auto start = clock::now();

    std::vector<double> sigma_GLC_v(Nw), sigma_BG_v(Nw);
    std::vector<double> Wv = W_space(Nw, M_GLC);

    std::cout << "W (GeV), sigma_GLC (nb), sigma_BG (nb)\n";

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < Nw; ++i) {
        sigma_GLC_v[i] = gamma_p::sigma(Wv[i], Q2, M_GLC, model) * GeV2_to_nb;
        sigma_BG_v[i]  = gamma_p::sigma(Wv[i], Q2, M_BG, model) * GeV2_to_nb;
    }

    for (int i = 0; i < Nw; ++i) {
        std::cout << Wv[i] << " "
                  << sigma_GLC_v[i] << " "
                  << sigma_BG_v[i] << "\n";

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
void dsigma_dt_csv(std::string model)
{
    double W = 70.0; // GeV

    const Meson& M_GLC = input_meson(model);
    const Meson& M_BG  = get_meson_bg(M_GLC.meson, model);
    warn_ipsat_light_meson(M_GLC, model);

    double Q2 = 0.0;
    std::string W_str = doubleParaString(W);
    std::string filename ="csv/" + M_GLC.meson + "_dsigma_dt_W=" + W_str + "GeV.csv";
    std::ofstream fout(filename);
    fout << "t,dsigma_dt_GLC,dsigma_dt_BG\n";

    const int Npoints = 150;
    double tmin = 0.0, tmax = 4.0; // GeV^2
    std::vector<double> t_v(Npoints), dsdt_glc_v(Npoints), dsdt_bg_v(Npoints);

    //esse for escolhe os valores de t e calcula a seção de choque para cada modelo, salvando no arquivo csv
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < Npoints; ++i) {
        double frac = static_cast<double>(i) / (Npoints - 1);
        double t = tmin + frac * (tmax - tmin);

        t_v[i] = t;
        dsdt_glc_v[i] = gamma_p::dsigma_dt(t, W, Q2, M_GLC, model) * GeV2_to_nb; // converte para nb/GeV^2
        dsdt_bg_v[i]  = gamma_p::dsigma_dt(t, W, Q2, M_BG, model) * GeV2_to_nb;  // converte para nb/GeV^2
    }

    for (int i = 0; i < Npoints; ++i) {
        std::cout << t_v[i] << "," << dsdt_glc_v[i] << "," << dsdt_bg_v[i] << "\n";
        fout << t_v[i] << "," << dsdt_glc_v[i] << "," << dsdt_bg_v[i] << "\n";
    }
    fout.close();
    std::cout << "Arquivo '" << filename << "' gerado." << std::endl;

    plot_dsigma_dt(M_GLC.meson);
}


    void sigma_gamma_A(std::string model, Nucleus&Nucleo)
    {
    double Q2=0.0;

    const Meson& M_GLC = input_meson("GBW"); //importa o glc do meson escolhido
    const Meson& M_BG  = meson_modelsGBW.find(M_GLC.meson)->second.M_BG; //pega o bg correspondente ao meson escolhido

    std::string Q2_str = std::format("{:.3g}", Q2);
    std::string path = "out/csv/sigma/" + M_GLC.meson + "/"+M_GLC.meson + "_sigma-gamma"+
                            Nucleo.name + "_" + model + "_fc";
    std::filesystem::create_directories("out/csv/sigma/" + M_GLC.meson);
    std::string filename_csv = path + ".csv";
    std::string filename_dat = path + ".dat";
    std::ofstream fout_csv(filename_csv);
    std::ofstream fout_dat(filename_dat);
    // cabeçalho csv
    fout_csv << "W,sigma_GLC,sigma_BG\n";

    // comentário opcional no .dat (xmgrace ignora linhas com #)
    fout_dat << "# W sigma_GLC sigma_BG\n";

    const int Nw = 100;


    std::vector<double> Wv = W_space(Nw, M_GLC);

    using clock = std::chrono::steady_clock;
    auto start = clock::now();
    
    TA_Table table = precompute_TA(200, 50.0, Nucleo); // pré-calcula T_A(b) para núcleos

    std::vector<double> sigma_GLC_v(Nw), sigma_BG_v(Nw);


    std::cout << "W (Gev), sigma_GLC (nb), sigma_BG (nb)\n";
    #pragma omp parallel for schedule(dynamic)
for (int i = 0; i < Nw; ++i) {
        sigma_GLC_v[i] = gamma_A::sigma(Wv[i], Q2, M_GLC, model, table, Nucleo) * GeV2_to_nb;
        sigma_BG_v[i]  = gamma_A::sigma(Wv[i], Q2, M_BG, model, table, Nucleo) * GeV2_to_nb;
    }

    for (int i = 0; i < Nw; ++i) {
    std::cout << Wv[i] << " " << sigma_GLC_v[i] << " " << sigma_BG_v[i] << "\n";
    fout_csv << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << "\n";
    fout_dat << Wv[i] << " " << sigma_GLC_v[i] << " " << sigma_BG_v[i] << "\n";
}

    fout_csv.close();
    fout_dat.close();

    auto end = clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Tempo de execução: " << duration << " ms" << std::endl;
    std::cout << "Arquivo '" << path << "'(csv e dat) gerado em " << duration << " ms" << std::endl;
    plot_sigma_gammaA(M_GLC.meson, filename_csv, model, Nucleo);
}


    void rapidez_AA(double sqrt_s, std::string model, Nucleus &Nucleo1, Nucleus& Nucleo2)
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
    

    const int Ny = 25;
    double ymax = 4.0;//std::log(sqrt_s / M_GLC.MV);

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
    std::cout << "Calculado dσ/dy ="<< dsdy_GLC << " ," << dsdy_BG << " para y = " << y << "...\n";

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

void create_TA_b_csv(std::string meson)
{
    std::vector<double> W, TA_GLC, TA_BG;

    read_csv("out/csv/other/" + meson + "_TA_b.csv", W, TA_GLC, TA_BG);

    std::ofstream out("out/csv/other/" + meson + "_TA_b_processed.csv");
    out << "W,TA_GLC,TA_BG\n";
    for(size_t i=0;i<W.size();++i)
        out << W[i] << "," << TA_GLC[i] << "," << TA_BG[i] << "\n";
    out.close();
}

void Shadowing_factor_csv()
{
    std::vector<double> x = linspace(1e-5, 0.5, 100);
    std::vector<double> xg_p_vals(x.size()), xg_Pb_vals(x.size()), S(x.size());

    for(size_t i = 0; i < x.size(); ++i){
        xg_p_vals[i] = xg_p(x[i], 1.0);
        xg_Pb_vals[i] = xg_Pb(x[i], 1.0);
        S[i] = S_Pb(x[i], 1.0);
        cout<< "x = " << x[i] << ", xg_p = " << xg_p_vals[i] << ", xg_Pb = " 
        << xg_Pb_vals[i] << ", S_Pb = " << S[i] << endl;
    }

    cout << "Plotando xg do proton\n";
    plot_XY(x, xg_p_vals, "x", "xg_p(x)", "Gluon PDF do próton");
    cout << "Plotando xg do chumbo\n";

    plot_XY(x, xg_Pb_vals, "x", "xg_Pb(x)", "Gluon PDF do Pb");
    cout << "Plotando fator de shadowing\n";
    plot_XY(x, S, "x", "S_Pb(x)", "Fator de shadowing S_Pb(x)");
}