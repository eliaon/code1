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



#include "ctes.hpp"
#include "correcs.hpp"
#include "utils.hpp"
#include "integration.hpp"

#include "../dipole_amplitudes/GBW.hpp"
#include "../dipole_amplitudes/LHAPDF.hpp"
#include "../dipole_amplitudes/bCGC.hpp"
#include "../dipole_amplitudes/ipsat.hpp"
#include "../dipole_amplitudes/IIM.hpp"
#include "../libraries/mantysaari/dipoleamplitude.hpp"
#include "../calculations/nuclear.hpp"


namespace {
LHAnPDF& ct14lo_pdf()
{
    thread_local LHAnPDF pdf("CT14lo");
    return pdf;
}

MZ_ipsat::DipoleAmplitude& ipsat_dipole()
{
    thread_local MZ_ipsat::DipoleAmplitude dipole(MZ_ipsat::MZ_IPSAT);
    thread_local bool lookup_enabled = [] {
        dipole.EnableLookupTable();
        return true;
    }();
    (void)lookup_enabled;
    return dipole;
}
} // namespace

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
// Mäntysaari & Zurita 2018, Table II (IPsat, boosted Gaussian)
namespace {
constexpr double R2_Jpsi_ipsat = 1.5070 * 1.5070;
constexpr double R2_phi_ipsat  = 3.3922 * 3.3922;
constexpr double R2_rho_ipsat  = 3.6376 * 3.6376;

bool uses_ipsat_mesons(const std::string& dipole_model)
{
    std::string lower = dipole_model;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return lower == "ipsat";
}
} // namespace

// GLC: sem parâmetros na Table II; mantidos para comparação (massas IPsat)
Meson Jpsi_GLC_ipsat("psi", "Jpsi_GLC", massa_psi, mc_ipsat, qJ, 1.23, 6.5, 0.83, 3.0);
Meson phi_GLC_ipsat ("phi",  "phi_GLC",  massa_phi, ms_ipsat, qS, 4.75, 16.0, 1.41, 9.7);
Meson rho_GLC_ipsat ("rho",  "rho_GLC",  massa_rho, ms_ipsat, qU2, 4.47, 21.9, 1.79, 10.4);

Meson Jpsi_BG_ipsat("psi", "Jpsi_BG", massa_psi, mc_ipsat, qJ, 0.5890, 0.5860, R2_Jpsi_ipsat);
Meson phi_BG_ipsat ("phi", "phi_BG",  massa_phi, ms_ipsat, qS, 0.9950, 0.8400, R2_phi_ipsat);
Meson rho_BG_ipsat ("rho", "rho_BG",  massa_rho, ms_ipsat, qU2, 0.9942, 0.8926, R2_rho_ipsat);




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
    {"rho", {rho_GLC_ipsat, rho_BG_ipsat}},
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
    if (uses_ipsat_mesons(model)) {
        auto it = meson_modelsipsat.find(meson_input);
        if (it == meson_modelsipsat.end()) {
            std::cerr << "Meson nao encontrado no modelo IPSAT. Usando psi.\n";
            return meson_modelsipsat.at("psi").M_GLC;
        }
        return it->second.M_GLC;
    }

    if (model == "GBW" || model == "bCGC" || model == "BCGC" || model == "bcgc" ||
        model == "GBW(new)" || model == "GBW(old)" || model == "IIM_S" ||
        model == "IIM_RS" || model == "LHAnPDF") {
        auto it = meson_modelsGBW.find(meson_input);
        if (it == meson_modelsGBW.end()) {
            std::cerr << "Meson nao encontrado no modelo " << model << ". Usando psi.\n";
            return meson_modelsGBW.at("psi").M_GLC;
        }
        return it->second.M_GLC;
    }

    std::cerr << "Modelo desconhecido: " << model
              << ". Usando GBW por padrao.\n";

    auto it = meson_modelsGBW.find(meson_input);
    if (it == meson_modelsGBW.end()) {
        return meson_modelsGBW.at("psi").M_GLC;
    }
    return it->second.M_GLC;
}

const Meson& get_meson_bg(const std::string& meson_key, const std::string& dipole_model)
{
    if (uses_ipsat_mesons(dipole_model)) {
        return meson_modelsipsat.at(meson_key).M_BG;
    }
    return meson_modelsGBW.at(meson_key).M_BG;
}

void warn_ipsat_light_meson(const Meson& M, const std::string& dipole_model)
{
    if (!uses_ipsat_mesons(dipole_model)) {
        return;
    }
    if (M.meson == "phi" || M.meson == "rho") {
        std::cerr << "AVISO [IPSAT]: producao exclusiva de " << M.meson
                  << " nao e confiavel com m_l=0.03 GeV (Mantysaari & Zurita 2018, Sec. V). "
                  << "Dipolos grandes dominam; interprete os resultados com cautela.\n";
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
double B_regge(double W, const Meson& M) // cálculo de B usando a fórmula de Regge, com parâmetros ajustados para rho e phi
{                                        // Light vector meson photoproduction in ultraperipheral heavy ion
    if (M.meson == "rho"){               // collisions at the LHC within the Reggeometric Pomeron approach (https://arxiv.org/abs/2301.05136v2)
        double a = 0.60, b = 0.9, alpha = 0.21, W0 = 1.0;
        double B0 = 4.0 * (a/std::pow(M.MV, 2) + b / (2.0*0.938));
        double B = B0 + 4.0 * alpha * std::log(W / W0);
        cout << "B de Regge para rho: " << B << endl;
        return B;
    } else if (M.meson == "phi"){
        double a = 0.0, b = 1.34, alpha = 0.17, W0 = 1.0;
        double B0 = 4.0 * (a/std::pow(M.MV, 2) + b / (2.0*0.938));
        double B = B0 + 4.0 * alpha * std::log(W / W0);
        cout << "B de Regge para phi: " << B << endl;
        return B;
    } else {
        std::cerr << "Méson desconhecido para cálculo de B: " << M.meson << std::endl;
        return 0.0;
    }
}


double B_slope(double x, double Q2, const Meson& M) {
    double W = x_to_W(x, M);
    if (M.meson == "psi" || M.meson == "upsilon42" || M.meson == "upsilon45"){
        double B1 = 4.80 + 4.0* 0.133 *log(W/90.0); //valores do lhcb dados pelo haimon xdxd
        return B1;
    }
        else if (M.meson == "phi" || M.meson == "rho"){
            double B2 = 0.55 * (14.0 / pow((Q2 + M.MV*M.MV), 0.2) + 1.0);
            return B2;}
            //return B_regge(W, M);}
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




// -------------- CARREGA SET DE VALORES DE CURVAS DE N


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
            return ct14lo_pdf().amplitude_p(x, Q2, M);
        }
        else if(modelo == "IPSAT") {
            return IPSAT::amplitude_p(x, Delta, Q2, M);
        }
        else if(modelo == "bCGC" || modelo == "BCGC" || modelo == "bcgc") {
            return bCGC::amplitude_p(x, Delta, Q2, M);
        } else if(modelo == "IIM_S" || modelo == "iim_s"){
        return IIM::amplitude_p(x, Q2, M, IIM_S);

        }else if(modelo == "IIM_RS" || modelo == "iim_rs"){
        return IIM::amplitude_p(x, Q2, M, IIM_RS);
        }else {
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
        return ct14lo_pdf().N_p(x, r);
    }
    else if (modelo == "IPSAT") {;
        return ipsat_dipole().N(r, x, b);
    }
    else if (modelo == "bCGC" || modelo == "BCGC" || modelo == "bcgc") {
        return bCGC::N_p(r, x, b);
    }
    else if(modelo == "IIM_S" || modelo == "iim_s"){
        return IIM::N_p(r, x, IIM_S);

    }else if(modelo == "IIM_RS" || modelo == "iim_rs"){
        return IIM::N_p(r, x, IIM_RS);
    }else {
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
        return IPSAT::sigma_qq_p(r, x, Delta, ipsat_dipole());
    } else if (modelo == "bCGC" || modelo == "BCGC" || modelo == "bcgc"){
        return bCGC::sigma_qq_p(r, x, Delta);
    } else if(modelo == "LHAnPDF"){
        return ct14lo_pdf().sigma_qq_p(x, r);
    } else if (modelo == "IIM_S" || modelo == "iim_s"){
        return IIM::sigma_qq(r, x, IIM_S);

    }else if(modelo == "IIM_RS" || modelo == "iim_rs"){
        return IIM::sigma_qq(r, x, IIM_RS);
    }else if(modelo == "GBW_Shadowing"){
        return GBW::sigma_qq_p(r, x, gbw_10) * S_Pb(x, 1.0);}
    else {
        std::cerr << "Modelo desconhecido ou não implementado (get_dipolo_p):" << modelo <<std::endl;
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

std::vector<double> W_space(double N, const Meson& M)
{
    double Wmin = x_to_W(1e-2, M);
    double Wmax = 3e3;
    int n_points = std::max(1, static_cast<int>(std::round(N)));
    std::vector<double> values(n_points);

    if (n_points == 1) {
        values.front() = Wmin;
        return values;
    }

    const double log_min = std::log(Wmin);
    const double log_max = std::log(Wmax);

    for (int i = 0; i < n_points; ++i) {
        const double frac = static_cast<double>(i) / (n_points - 1);
        values[i] = std::exp(log_min + frac * (log_max - log_min));
    }

    return values;
}

double low_x_factor(double x, double exp)
{
    const double xc = std::clamp(x, 0.0, 1.0);
    return std::pow(1.0 - xc, exp);
}





std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result;
    if (num == 0) return result;
    if (num == 1) {
        result.push_back(start);
        return result;
    }

    double delta = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result.push_back(start + delta * i);
    }
    return result;
}