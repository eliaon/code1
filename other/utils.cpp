#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <functional>
#include <map>
#include <algorithm>
#include <filesystem>

#include "ctes.h"
#include "utils.h"
#include "integration.hpp"
#include "plot.h"

#include "../dipole_amplitudes/GBW.h"
#include "../dipole_amplitudes/LHAPDF.hpp"
#include "../dipole_amplitudes/ipsat.h"
#include "../libraries/mantysaari/dipoleamplitude.hpp"
#include "../calculations/nuclear.hpp"




std::string extrair_nome_base(const std::string& caminho)
{
    return std::filesystem::path(caminho).stem().string();
}

Meson::Meson(std::string m, std::string n, double MV_, double mf_, double ef_,
             double NT_, double R2T_, double NL_, double R2L_)
    : meson(m), nome(n), MV(MV_), mf(mf_), ef(ef_),
      NT(NT_), NL(NL_), R2T(R2T_), R2L(R2L_), R2(0.0), isGLC(true) {}

Meson::Meson(std::string m, std::string n, double MV_, double mf_, double ef_,
             double NT_, double NL_, double R2_)
    : meson(m), nome(n), MV(MV_), mf(mf_), ef(ef_),
      NT(NT_), NL(NL_), R2T(0.0), R2L(0.0), R2(R2_), isGLC(false) {}


// =========================================================
// Gaus-LC
// =========================================================

// J/psi
Meson Jpsi_GLC_GBW(
    "psi",
    "Jpsi_GLC",
    massa_psi,
    mc_GBW,
    qJ,
    1.23,
    6.5,
    0.83,
    3.0
);

// phi
Meson phi_GLC_GBW(
    "phi",
    "phi_GLC",
    massa_phi,
    ms_GBW,
    qS,
    4.75,
    16.0,
    1.41,
    9.7
);

// rho
Meson rho_GLC_GBW(
    "rho",
    "rho_GLC",
    massa_rho,
    muq_GBW,
    qU2,
    4.47,
    21.9,
    1.79,
    10.4
);

// Upsilon (m_b = 4.5 GeV)
Meson Upsilon_GLC_GBW_45(
    "upsilon45",
    "Upsilon_GLC_45",
    massa_upsilon,
    mb_GBW,
    qB,
    0.67,
    2.16,
    0.47,
    1.01
);

// Upsilon (m_b = 4.2 GeV)
Meson Upsilon_GLC_GBW_42(
    "upsilon42",
    "Upsilon_GLC_42",
    massa_upsilon,
    4.2,
    qB,
    0.76,
    1.91,
    0.47,
    1.01
);


// =========================================================
// Boosted Gaussian
// =========================================================

// J/psi
Meson Jpsi_BG_GBW(
    "psi",
    "Jpsi_BG",
    massa_psi,
    mc_GBW,
    qJ,
    0.578,
    0.575,
    2.3
);

// phi
Meson phi_BG_GBW(
    "phi",
    "phi_BG",
    massa_phi,
    ms_GBW,
    qS,
    0.919,
    0.825,
    11.2
);

// rho
Meson rho_BG_GBW(
    "rho",
    "rho_BG",
    massa_rho,
    muq_GBW,
    qU2,
    0.911,
    0.853,
    12.9
);

// Upsilon (m_b = 4.5 GeV)
Meson Upsilon_BG_GBW_45(
    "upsilon45",
    "Upsilon_BG_45",
    massa_upsilon,
    mb_GBW,
    qB,
    0.469,
    0.469,
    0.55
);

// Upsilon (m_b = 4.2 GeV)
Meson Upsilon_BG_GBW_42(
    "upsilon42",
    "Upsilon_BG_42",
    massa_upsilon,
    4.2,
    qB,
    0.481,
    0.480,
    0.57
);
Meson Jpsi_GLC_ipsat("psi", "Jpsi_GLC", massa_psi, mc_ipsat, qJ, 1.23, 6.5, 0.83, 3.0);
Meson phi_GLC_ipsat ("phi",  "phi_GLC",  massa_phi, ms_ipsat, qS, 4.75, 16.0, 1.41, 9.7);

Meson Jpsi_BG_ipsat ("psi", "Jpsi_BG",  massa_psi, mc_ipsat, qJ, 0.578, 0.575, 2.3);//R2_psi);
Meson phi_BG_ipsat  ("phi",  "phi_BG",   massa_phi, ms_ipsat, qS, 0.919, 0.825, 11.2);//R2_phi);




std::string doubleParaString(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(1) << x;
    return oss.str();
}

// Derivada por extrapolação de Richardson
double dfridr(
    const std::function<double(double)>& func,
    double x,
    double h,
    double& err
) {
    constexpr int NTAB = 10;
    constexpr double CON  = 1.4;
    constexpr double CON2 = CON * CON;
    constexpr double BIG  = 1.0e30;
    constexpr double SAFE = 2.0;

    if (h == 0.0) {
        throw std::runtime_error("h must be nonzero in dfridr");
    }

    double a[NTAB][NTAB];

    double hh = h;
    a[0][0] = (func(x + hh) - func(x - hh)) / (2.0 * hh);

    err = BIG;
    double dfridr_val = a[0][0];

    for (int i = 1; i < NTAB; ++i) {
        hh /= CON;
        a[0][i] = (func(x + hh) - func(x - hh)) / (2.0 * hh);

        double fac = CON2;

        for (int j = 1; j <= i; ++j) {
            a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);
            fac *= CON2;

            double errt = std::max(
                std::abs(a[j][i] - a[j - 1][i]),
                std::abs(a[j][i] - a[j - 1][i - 1])
            );

            if (errt <= err) {
                err = errt;
                dfridr_val = a[j][i];
            }
        }

        if (std::abs(a[i][i] - a[i - 1][i - 1]) >= SAFE * err) {
            return dfridr_val;
        }
    }

    return dfridr_val;
}
double derivative_richardson(const std::function<double(double)>& f,
                             double x,
                             double h)
{
    double D1 = (f(x + h) - f(x - h)) / (2.0*h);
    double D2 = (f(x + h/2) - f(x - h/2)) / h;

    return (4.0*D2 - D1) / 3.0;
}

double derivative_poly5(
    const std::function<double(double)>& f,
    double x)
{
    double h = 1e-10;

    double f1=f(x-2*h);
    double f2=f(x-h);
    double f3=f(x+h);
    double f4=f(x+2*h);

    return (-f4 + 8*f3 - 8*f2 + f1)/(12*h);
}

// ----------------- conversão W <-> x ----------------

double x_to_W(double x, const Meson& M) {
    return M.MV / std::sqrt(x);
}

double W_to_x(double W, const Meson& M) {
    return (M.MV * M.MV) / (W * W);
}



std::map<std::string, MesonModelsGBW> meson_modelsGBW = {
    {"psi", {Jpsi_GLC_GBW, Jpsi_BG_GBW}},
    {"phi", {phi_GLC_GBW, phi_BG_GBW}},
    {"rho", {rho_GLC_GBW, rho_BG_GBW}},
    {"upsilon45", {Upsilon_GLC_GBW_45, Upsilon_BG_GBW_45}},
    {"upsilon42", {Upsilon_GLC_GBW_42, Upsilon_BG_GBW_42}}
};

std::map<std::string, MesonModelsipsat> meson_modelsipsat = {
    {"psi", {Jpsi_GLC_ipsat, Jpsi_BG_ipsat}},
    {"phi", {phi_GLC_ipsat, phi_BG_ipsat}},
};

// ----------------- função escolhe meson ----------
Meson input_meson(const std::string model) // devolve o meson glc
{
    std::string meson_input;
    std::cout << "Insira o meson (psi, phi, rho, upsilon42, upsilon45): ";
    std::cin >> meson_input;

    // -------- NORMALIZAÇÃO --------
    // transforma tudo em minúsculo
    std::transform(meson_input.begin(), meson_input.end(),
                   meson_input.begin(), ::tolower);

    // padroniza nomes
    if (meson_input == "psi") meson_input = "psi";
    else if (meson_input == "phi") meson_input = "phi";
    else if (meson_input == "rho") meson_input = "rho";
    else if (meson_input == "upsilon45") meson_input = "upsilon45";
    else if (meson_input == "upsilon42") meson_input = "upsilon42";
    else {
        std::cerr << "Meson invalido. Usando psi por padrao.\n";
        meson_input = "psi";
    }

    // -------- SELEÇÃO DO MODELO --------
    if (model == "GBW") {
        auto it = meson_modelsGBW.find(meson_input);
        if (it == meson_modelsGBW.end()) {
            std::cerr << "Meson nao encontrado no modelo GBW. Usando psi.\n";
            return meson_modelsGBW.at("psi").M_GLC;
        }
        return it->second.M_GLC;
    }
    else if (model == "ipsat") {
        auto it = meson_modelsipsat.find(meson_input);
        if (it == meson_modelsipsat.end()) {
            std::cerr << "Meson nao encontrado no modelo ipsat. Usando psi.\n";
            return meson_modelsipsat.at("psi").M_GLC;
        }
        return it->second.M_GLC;
    }
    else {
        std::cerr << "Modelo desconhecido: " << model
                  << ". Usando GBW por padrao.\n";

        auto it = meson_modelsGBW.find(meson_input);
        if (it == meson_modelsGBW.end()) {
            return meson_modelsGBW.at("psi").M_GLC;
        }
        return it->second.M_GLC;
    }
}
void perfil(const Meson& meson){
    if (meson.isGLC){
        std::cout << "Perfil do méson Gaus-LC:\n";
        std::cout << "Méson: " << meson.meson << "\n";
        std::cout << "Massa do méson (GeV): " << meson.MV << "\n";
        std::cout << "Massa do quark (GeV): " << meson.mf << "\n";
        std::cout << "Carga efetiva: " << meson.ef << "\n";
        std::cout << "Normalização transversa NT: " << meson.NT << "\n";
        std::cout << "Parâmetro transverso R2T (GeV^-2): " << meson.R2T << "\n";
        std::cout << "Normalização longitudinal NL: " << meson.NL << "\n";
        std::cout << "Parâmetro longitudinal R2L (GeV^-2): " << meson.R2L << "\n";
    } else {
        std::cout << "Perfil do méson Boosted Gaussian:\n";
        std::cout << "Méson: " << meson.meson << "\n";
        std::cout << "Massa do méson (GeV): " << meson.MV << "\n";
        std::cout << "Massa do quark (GeV): " << meson.mf << "\n";
        std::cout << "Carga efetiva: " << meson.ef << "\n";
        std::cout << "Normalização transversa NT: " << meson.NT << "\n";
        std::cout << "Normalização longitudinal NL: " << meson.NL << "\n";
        std::cout << "Parâmetro R2 (GeV^-2): " << meson.R2 << "\n";
    }
}

// ---------------- slope B(Q2) ----------------
double B_slope(double x, double Q2, const Meson& M) {
    double W = x_to_W(x, M);
    if (M.meson == "psi" || M.meson == "upsilon42" || M.meson == "upsilon45"){
        double B1 = 4.80 + 4.0* 0.133 *log(W/90.0); //valores do lhcb dados pelo haimon xdxd
        return B1;
    }
        else if (M.meson == "phi" || M.meson == "rho"){
            double B2 = 0.55 * (14.0 / pow((Q2 + M.MV*M.MV), 0.2) + 1.0);
            return B2;}
            else {
                std::cerr << "Méson desconhecido para cálculo de B: " << M.meson << std::endl;
                return 0.0;
            }
}







std::string timestamp(void)
{
    std::time_t now = std::time(nullptr);
    std::tm* local = std::localtime(&now);

    char buffer[32];
    std::strftime(buffer, sizeof(buffer), "%Y%m%d_%H%M%S", local);
    return std::string(buffer);
}


// ------------ LÊ CSV COM 2 COLUNAS, PARA N 
void read_two_columns(
    const std::string& filename,
    std::vector<double>& x,
    std::vector<double>& y)
{
    std::ifstream file(filename);
    std::cout << "Reading file: " << filename << std::endl;
    if(!file)
        throw std::runtime_error("Cannot open file: " + filename);

    std::string line;

    std::getline(file,line); // header

    while(std::getline(file,line))
    {
        if(line.empty()) continue;

        std::stringstream ss(line);
        std::string a,b;

        std::getline(ss,a,',');
        std::getline(ss,b,',');

        x.push_back(std::stod(a));
        y.push_back(std::stod(b));
    }
}

void read_sigma_exp(
    const std::string& filename,
    std::vector<int>& dataset,
    std::vector<double>& W,
    std::vector<double>& sigma,
    std::vector<double>& err)
{
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line))
    {
        // remove espaços iniciais
        line.erase(0, line.find_first_not_of(" \t"));

        if (line.empty() || line[0] == '#')
            continue;

        std::vector<double> cols;

        // --- tenta parsing como CSV ---
        {
            std::stringstream ss(line);
            std::string token;

            while (std::getline(ss, token, ','))
            {
                try {
                    cols.push_back(std::stod(token));
                }
                catch (...) {
                    // ignora tokens não numéricos (ex: "-")
                }
            }
        }

        // --- fallback: separação por espaço/tab ---
        if (cols.size() < 3)
        {
            cols.clear();
            std::stringstream ss(line);
            std::string token;

            while (ss >> token)
            {
                try {
                    cols.push_back(std::stod(token));
                }
                catch (...) {
                    // ignora lixo
                }
            }
        }

        if (cols.size() < 3)
            continue;

        // --- caso 1: formato simples ---
        if (cols.size() == 3)
        {
            int d = 0;
            dataset.push_back(d);
            W.push_back(cols[0]);
            sigma.push_back(cols[1]);
            err.push_back(cols[2]);
        }

        // --- caso 2: erro assimétrico (4 colunas) ---
        else if (cols.size() >= 4 && cols.size() < 8)
        {
            int d = 0;

            double W_val = cols[0];
            double s_val = cols[1];
            double dy_minus = std::abs(cols[2]);
            double dy_plus  = std::abs(cols[3]);

            double e_val = 0.5 * (dy_minus + dy_plus);

            dataset.push_back(d);
            W.push_back(W_val);
            sigma.push_back(s_val);
            err.push_back(e_val);
        }

        // --- caso 3: HEPData completo ---
        else if (cols.size() >= 8)
        {
            int d = static_cast<int>(cols[0]);

            double W_val     = cols[1];
            double sigma_val = cols[3];

            double stat_p = cols[4];
            double stat_m = std::abs(cols[5]);
            double sys_p  = cols[6];
            double sys_m  = std::abs(cols[7]);

            double stat = 0.5 * (stat_p + stat_m);
            double sys  = 0.5 * (sys_p + sys_m);

            double err_val = std::sqrt(stat*stat + sys*sys);

            dataset.push_back(d);
            W.push_back(W_val);
            sigma.push_back(sigma_val * 1000.0); // μb → nb
            err.push_back(err_val * 1000.0);
        }
    }
}


//-------------- LEITOR DE DADOS EXPERIMENTAIS DO XYSCAN ----------------------
void read_xyscan_csv(const std::string& filename,
                     std::vector<double>& W,
                     std::vector<double>& sigma,
                     std::vector<double>& error)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erro ao abrir arquivo: " << filename << "\n";
        return;
    }

    std::string line;

    while (std::getline(file, line)) {

        // ignora comentários e linhas vazias
        if (line.empty() || line[0] == '#')
            continue;

        std::stringstream ss(line);
        std::string col;
        std::vector<double> values;

        while (std::getline(ss, col, ',')) {
            if (col.empty()) continue; // ignora vírgula final
            values.push_back(std::stod(col));
        }

        if (values.size() < 3)
    continue;

double w   = values[0];
double sig = values[1];
double err;

// caso 1: CSV com 3 colunas (erro já simétrico)
if (values.size() == 3) {
    err = values[2];
}
// caso 2: CSV com 4 colunas (erro assimétrico)
else {
    double dy_minus = values[2];
    double dy_plus  = values[3];
    err = 0.5 * (dy_minus + dy_plus);
}

W.push_back(w);
sigma.push_back(sig);
error.push_back(err);
    }

    file.close();
}

// ------------ LÊ CSV'S COM 1 INDEPENDENTE E DUAS DEPENDENTES, IDEAL PARA COMPARAÇÕES
// ------------- BOOSTED GAUSSIAN - GAUSSIAN LIGHT CONE

void read_csv(
    const std::string& filename,
    std::vector<double>& x,
    std::vector<double>& glc,
    std::vector<double>& bg)
{
    std::ifstream file(filename);

    if(!file)
        throw std::runtime_error("Cannot open file: " + filename);

    x.clear();
    glc.clear();
    bg.clear();

    std::string line;

    std::getline(file, line); // header
    
    while (std::getline(file, line))
    {
        if(line.empty()) continue;

        std::stringstream ss(line);
        std::string a,b,c;

        if (!std::getline(ss, a, ',')) continue;
        if (!std::getline(ss, b, ',')) continue;
        if (!std::getline(ss, c, ',')) continue;

        try {
            x.push_back(std::stod(a));
            glc.push_back(std::stod(b));
            bg.push_back(std::stod(c));
        }
        catch (...) {
            // ignora linha mal formatada
            continue;
        }
    }

    // sanity check forte
    if (x.size() != glc.size() || x.size() != bg.size()) {
        throw std::runtime_error("CSV columns size mismatch in: " + filename);
    }
}

void read_rapidity_hepdata(
    const std::string& filename,
    std::vector<double>& y,
    std::vector<double>& dsdy,
    std::vector<double>& err)
{
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line))
    {
        if (line.empty()) continue;

        //  remove espaços iniciais
        line.erase(0, line.find_first_not_of(" \t"));

        if (line.empty()) continue;

        // ignora headers
        if (line[0] == '#' || line[0] == '|')
            continue;
        std::stringstream ss(line);
        std::string token;
        std::vector<double> cols;

        while (std::getline(ss, token, ','))
        {
            try {
                cols.push_back(std::stod(token));
            } catch (...) {}
        }

        if (cols.size() < 8) continue;

        double y_val = cols[0];
        double sigma = cols[3];

        double stat_p = cols[4];
        double stat_m = std::abs(cols[5]);
        double sys_p  = cols[6];
        double sys_m  = std::abs(cols[7]);

        double stat = 0.5 * (stat_p + stat_m);
        double sys  = 0.5 * (sys_p + sys_m);

        double error = std::sqrt(stat*stat + sys*sys);

        y.push_back(y_val);
        dsdy.push_back(sigma);
        err.push_back(error);
    }
}

// -------------- CARREGA SET DE VALORES DE CURVAS DE N

bool load_set(
    const std::string& xstr,
    std::vector<double>& r2_ipsat,
    std::vector<double>& N_ipsat,
    std::vector<double>& r2_gbw,
    std::vector<double>& N_gbw,
    std::vector<double>& r2_iim,
    std::vector<double>& N_iim)
{
    std::cout << "Loading x=" << xstr << std::endl;
    read_two_columns("csv/N_ipsat_x=" + xstr + ".csv", r2_ipsat, N_ipsat);
    read_two_columns("csv/N_GBW_x=" + xstr + ".csv",  r2_gbw,  N_gbw);
    read_two_columns("csv/N_IIM_x=" + xstr + ".csv",  r2_iim,  N_iim);
    if(r2_ipsat.empty() || r2_gbw.empty() || r2_iim.empty())
{
    std::cerr << "Empty dataset for x=" << xstr << std::endl;
    return false;
}
return true;
}



std::string get_meson()
{
    std::cout << "Escreva o meson (psi, phi): ";
    std::string meson;
    std::cin >> meson;

    if (meson != "psi" && meson != "phi")
    {
        std::cerr << "Meson inválido. Use 'psi' ou 'phi'." << std::endl;
        exit(1);
    }

    return meson;
}


std::string vec_to_pylist(const std::vector<double>& v)
{
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        oss << v[i];
        if (i != v.size() - 1) oss << ",";
    }
    oss << "]";
    return oss.str();
}

std::string flavorToString(int flavor) {
    switch (flavor) {
        case 1: return "xd";
        case 2: return "xu";
        case 3: return "xs";
        case 4: return "xc";
        case 5: return "xb";
        case 6: return "xt";
        case 21: return "xg";
        default: return "unknown";
    }
}

std::string flavorName(int flavor) {
    switch (flavor) {
        case 1: return "down";
        case 2: return "up";
        case 3: return "strange";
        case 4: return "charm";
        case 5: return "bottom";
        case 6: return "top";
        case 21: return "gluon";
        default: return "unknown";
    }
}





double get_amplitude_p(double x, double Delta,double Q2, const Meson& M, std::string modelo)
    {
        if (modelo == "GBW(old)") {
            double amp = GBW::amplitude_p(x, Q2, M, gbw);
          //  cout <<"x = " << x << " Amplitude GBW: " << amp << std::endl;
            return amp;
        } else if (modelo == "GBW(new)") {
            double amp = GBW::amplitude_p(x, Q2, M, gbw_10);
          //  cout <<"x = " << x << " Amplitude GBW: " << amp << std::endl;
            return amp;
        }
        else if (modelo == "LHAnPDF") {
            // create an instance with a PDF set name and call it
            // Provide a default PDF set name since LHAPDF has no default constructor
            LHAnPDF lhapdf("CT14lo");
            return lhapdf.amplitude_p(x, Q2, M);
        }
        else if(modelo == "IPSAT") {
            return IPSAT::amplitude_p(x, Delta, Q2, M);
        } else {
            std::cerr << "Modelo desconhecido: " << modelo << std::endl;
            return 0.0;
        }
    }

double get_Np(double r, double x, std::string modelo, double b)
{
    if (modelo == "GBW(old)") {
        return GBW::N_p(r, x, gbw);
    }
    else if (modelo == "GBW(new)") {
        return GBW::N_p(r, x, gbw_10);
    }
    else if (modelo == "LHAnPDF") {
        LHAnPDF lhapdf("CT14lo");
        return lhapdf.N_p(r, x);
    }
    else if (modelo == "IPSAT") {
        DipoleAmplitude dipole(MZ_IPSAT);
        return dipole.N(r, x, b);
    }
    else {
        std::cerr << "Modelo desconhecido: " << modelo << std::endl;
        return 0.0;
    }
}

double get_dipolo_p(double r, double x, double Delta, std::string modelo)
{
    if (modelo == "GBW(new)"){
        return GBW::sigma_qq_p(r, x, gbw_10);
    } else if (modelo =="GBW(old)"){
        return GBW::sigma_qq_p(r, x, gbw);
    } else if (modelo == "IPSAT"){
        DipoleAmplitude dipole(MZ_IPSAT);
        return IPSAT::sigma_qq_p(r, x, Delta, dipole);
    } else if(modelo == "LHAnPDF"){
        LHAnPDF lhapdf("CT14lo");
        return lhapdf.sigma_qq_p(x, r);
    } else {
        std::cerr << "Modelo desconhecido ou não implementado:" << modelo <<std::endl;
        return 0.0;
    }
}

string N_file(double x, std::string modelo)
{
    std::string filename = "out/csv/curves/N_" + modelo + "_x=" + doubleParaString(x) + ".csv";
    std::ofstream fout(filename);
    fout << "r,N\n";


    const int Npoints = 5000;
    double rmin = 1e-4, rmax = 10.0;

    for (int i = 0; i < Npoints; ++i)
    {
        double frac = (double)i / (Npoints - 1);
        double r = rmin * std::pow(rmax / rmin, frac)*CFAC; // escala log
        double N_val = get_Np(r, x, modelo); // exemplo para x=1e-4 e x0=1e-2
        fout << r/CFAC << "," << N_val << "\n"; // converte r para fm
    }
    return filename;
}

