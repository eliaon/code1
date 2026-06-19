#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <filesystem>

#include "../calculations/sigma.hpp"
#include "../calculations/wavefunctions.hpp"
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

    namespace string_files{
    std::string overlap_file(const Meson&M, bool fc){
        std::string model_str = std::string("GBW(new)") + (fc ? "_fc" : "");
        std::string filename = "out/csv/overlap_" + M.nome + "_" + model_str + ".csv";
        std::ofstream fout(filename);
        fout << "r,overlap\n";
        std::vector<double> r = logspace(1e-6, 15.0, 500);
        std::vector<double> overlap(r.size());
        const int Nr = 200;
        double rmin = 1e-4, rmax = 15.0;
        for (size_t i = 0; i < r.size(); ++i)
{
    overlap[i] = overlap_r(r[i], 0.0, M, fc);

    fout << r[i]/CFAC
         << ","
         << 0.5*r[i]*overlap[i]
         << "\n";
}
        std::cout << "Arquivo CSV gerado: " << filename << std::endl;
        return filename;
    }
string sigma_file(const Meson& M, string modelo, bool fc){
    string modelo_str = modelo + (fc ? "_fc" : "");

    std::string filename = "out/csv/curves/sigma_" + modelo_str + "_" + M.nome + ".csv";
    std::ofstream fout(filename);
    fout << "W, sigma\n";
    std::vector<double> W_values = W_space(100, M);
    std::vector<double> sigma_values(W_values.size());
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < 100; ++i){
        double sigma = gamma_p::sigma(W_values[i], 0.0, M, modelo, fc); 
        //std::cout << "W: " << W_values[i] << " GeV, sigma: " << sigma * GeV2_to_nb << " nb\n";
        sigma_values[i] = sigma * GeV2_to_nb;
    }
    for(size_t i=0; i<W_values.size(); ++i){
        fout << W_values[i] << "," << sigma_values[i] << "\n";
        //std::cout << "Escrevendo no CSV: W = " << W_values[i] << " GeV, sigma = " << sigma_values[i] << " nb\n";
    }
    fout.close();
    return filename;
}

string rapidezAA_file(double sqrt_s, const Meson& M, string modelo, Nucleus& Nucleo, bool fc){
    std::string modelo_str = modelo + (fc ? "_fc" : "");
    std::string filename = "out/csv/curves/rapidez_"+ M.nome + "_" + modelo_str + "_" + Nucleo.name + Nucleo.name + ".csv";
    std::ofstream fout(filename);
    fout << "y,rapidez\n";

    TA_Table tableA = get_TA_table(Nucleo); // pré-calcula T_A(b) para o núcleo A


    std::vector<double> y_values = linspace(0.0, 2.0, 20);
    std::vector<double> rap_values(y_values.size());
    #pragma omp parallel for schedule(dynamic)
    for(size_t i=0; i<y_values.size(); ++i){
        double rap = gamma_A::d_sigma_dy_AA(y_values[i], sqrt_s, 0.0, M, modelo, tableA, Nucleo, fc) * gev2_to_mb;
        rap_values[i] = rap;
    }
    for(size_t i=0; i<y_values.size(); ++i){
        fout << y_values[i] << "," << rap_values[i] << "\n";
    }
    fout.close();
    return filename;

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


    void sigma_gamma_p(std::string model, bool fc)
{
    double Q2 = 0.0;

    const Meson& M_GLC = input_meson(model);
    const Meson& M_BG  = get_meson_bg(M_GLC.meson, model);
    warn_ipsat_light_meson(M_GLC, model);
    std::string model_str = model + (fc ? "_fc" : "");

    std::cout << "Calculando sigma_p para " << M_GLC.meson
              << " com modelo " << model_str << "\n";

    std::filesystem::create_directories("out/csv/sigma/" + M_GLC.meson);

    std::string basepath = "out/csv/sigma/" + M_GLC.meson + "/" +
                           M_GLC.meson + "_sigma-gammap_" + model_str;

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
        sigma_GLC_v[i] = gamma_p::sigma(Wv[i], Q2, M_GLC, model, fc) * GeV2_to_nb;
        sigma_BG_v[i]  = gamma_p::sigma(Wv[i], Q2, M_BG, model, fc) * GeV2_to_nb;
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

    plot_sigma(M_GLC.meson, filename_csv, model_str);
}

void dsigma_dt_csv(std::string model, bool fc)
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
        dsdt_glc_v[i] = gamma_p::dsigma_dt(t, W, Q2, M_GLC, model, fc) * GeV2_to_nb; // converte para nb/GeV^2
        dsdt_bg_v[i]  = gamma_p::dsigma_dt(t, W, Q2, M_BG, model, fc) * GeV2_to_nb;  // converte para nb/GeV^2
    }

    for (int i = 0; i < Npoints; ++i) {
        std::cout << t_v[i] << "," << dsdt_glc_v[i] << "," << dsdt_bg_v[i] << "\n";
        fout << t_v[i] << "," << dsdt_glc_v[i] << "," << dsdt_bg_v[i] << "\n";
    }
    fout.close();
    std::cout << "Arquivo '" << filename << "' gerado." << std::endl;

    plot_dsigma_dt(M_GLC.meson);
}


    void sigma_gamma_A(std::string model, Nucleus&Nucleo, bool fc)
    {
    double Q2=0.0;

    const Meson& M_GLC = input_meson("GBW"); //importa o glc do meson escolhido
    const Meson& M_BG  = meson_modelsGBW.find(M_GLC.meson)->second.M_BG; //pega o bg correspondente ao meson escolhido

    std::string model_str = model + (fc ? "_fc" : "");

    std::string Q2_str = std::format("{:.3g}", Q2);
    std::string path = "out/csv/sigma/" + M_GLC.meson + "/"+M_GLC.meson + "_sigma-gamma"+
                            Nucleo.name + "_" + model_str ;
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
        sigma_GLC_v[i] = gamma_A::sigma(Wv[i], Q2, M_GLC, model, table, Nucleo, fc) * GeV2_to_nb;
        sigma_BG_v[i]  = gamma_A::sigma(Wv[i], Q2, M_BG, model, table, Nucleo, fc) * GeV2_to_nb;
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
    plot_sigma_gammaA(M_GLC.meson, filename_csv, model_str, Nucleo);
}


    void rapidez_AA(double sqrt_s, std::string model, Nucleus &Nucleo1, Nucleus& Nucleo2, bool fc)
{
    //double sqrt_s = 2.76e3; // GeV

    const Meson& M_GLC = input_meson("GBW");
    const Meson& M_BG  = meson_modelsGBW.find(M_GLC.meson)->second.M_BG;
    const string model_str = model + (fc ? "_fc" : "");

    std::string sqrt_s_str = std::format("{:.3g}", sqrt_s);
    std::string path = "out/csv/Rapidez/" + M_GLC.meson + "/"+M_GLC.meson + "_Rapidez_"+
                            Nucleo1.name + "-" + Nucleo2.name +"_" + model_str ;

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
        dsdy_GLC = gamma_A::d_sigma_dy_AA(y, sqrt_s, Q2, M_GLC, model, tableA, Nucleo1, fc) * gev2_to_mb;
        dsdy_BG  = gamma_A::d_sigma_dy_AA(y, sqrt_s, Q2, M_BG,  model, tableA, Nucleo1, fc) * gev2_to_mb;
    }
    else
    {
        dsdy_GLC = gamma_A::d_sigma_dy_AB(y, sqrt_s, Q2, M_GLC, model, tableA, tableB, Nucleo1, Nucleo2, fc) * gev2_to_mb;
        dsdy_BG  = gamma_A::d_sigma_dy_AB(y, sqrt_s, Q2, M_BG,  model, tableA, tableB, Nucleo1, Nucleo2, fc) * gev2_to_mb;
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

    plot_rapidez(filename_csv, model_str, M_GLC.meson, Nucleo1.name+Nucleo2.name, sqrt_s);

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

void Shadowing_factor_csv(double Q)
{
    std::vector<double> x = logspace(1e-5, 0.5, 100);
    std::vector<double> xg_p_vals(x.size()), xg_Pb_vals(x.size()), S(x.size());

    for(size_t i = 0; i < x.size(); ++i){
        xg_p_vals[i] = xg_p(x[i], Q);
        xg_Pb_vals[i] = xg_Pb(x[i], Q);
        S[i] = S_Pb(x[i], sqrt(Q));
        cout<< "x = " << x[i] << ", xg_p = " << xg_p_vals[i] << ", xg_Pb = " 
        << xg_Pb_vals[i] << ", S_Pb = " << S[i] << endl;
    }

    cout << "Plotando fator de shadowing\n";
    plot_XY(x, S, "x", "S_Pb", "Fator de shadowing S_Pb");
    //cout << "Plotando xg do proton\n";
    plot_XY(x, xg_p_vals, "x", "xg_p", "Gluon PDF do próton");
    //cout << "Plotando xg do chumbo\n";
    plot_XY(x, xg_Pb_vals, "x", "xg_Pb", "Gluon PDF do Pb");
}

namespace Compare_all {

    const std::vector<std::string> modelos = 
    {"GBW(new)", 
        //"GBW(old)", 
        "bCGC", 
        //"IIM_S", 
        "IPSAT",
         "IIM_RS"};
    
    void N(double x){

        std::vector<std::pair<std::string, std::string>> filenames;

        for (const auto& modelo : modelos) {
            std::cout << "Gerando CSV para N(r) com modelo " << modelo << "...\n";
            std::string filename = N_file(x, modelo);
            filenames.push_back({filename, modelo});
            std::cout << "CSV gerado: " << filename << "\n";
        }

        plot_N_models(filenames, x);
    }

    void sigma_gamma_p(bool fc)
    {
        std::vector<std::pair<std::string, std::string>> filenames;

        const Meson& M_GLC = input_meson("GBW");
        const Meson& M_BG  = get_meson_bg(M_GLC.meson, "GBW");
        for(const auto& modelo: modelos){
            std::cout << "Gerando CSV para sigma_gamma_p com modelo " << modelo << "...\n";
            std::string filename = CSV::string_files::sigma_file(M_BG, modelo, fc);
            filenames.push_back({filename, modelo});
            std::cout << "CSV gerado: " << filename << "\n";
        }

        plot_sigma_models(filenames, M_BG, fc);
    }

    void rapidez_AA(double sqrt_s, Nucleus& Nucleo, bool fc)
    {
        std::vector<std::pair<std::string, std::string>> filenames;
        const Meson& M_GLC = input_meson("GBW");
        const Meson& M_BG  = get_meson_bg(M_GLC.meson, "GBW");
        for(const auto& modelo: modelos){
            std::cout << "Gerando CSV para rapidez_AA com modelo " << modelo << "...\n";
            std::string filename = CSV::string_files::rapidezAA_file(sqrt_s, M_GLC, modelo, Nucleo, fc);
            filenames.push_back({filename, modelo});
            std::cout << "CSV gerado: " << filename << "\n";
        }

        plot_rapidezAA_models(filenames, sqrt_s, M_BG, fc);
    }

    void overlap(bool fc)
    {
        std::vector<std::pair<std::string, std::string>> filenames;
        const Meson& M_GLC = input_meson("GBW");
        const Meson& M_BG  = get_meson_bg(M_GLC.meson, "GBW");
        
        filenames.push_back({CSV::string_files::overlap_file(M_GLC, fc), "GLC"});
        filenames.push_back({CSV::string_files::overlap_file(M_BG, fc), "BG"});

        plot_overlap_models(filenames, M_GLC, fc);
    }
}












