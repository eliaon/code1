#ifndef UTILS_H
#define UTILS_H

#include <functional>
#include <map>
#include <string>

#include "ctes.hpp"

std::string extrair_nome_base(const std::string& caminho);

class Meson {
public:
    std::string meson;
    std::string nome;
    double MV;
    double mf;
    double ef;
    double NT;
    double NL;
    double R2T;
    double R2L;
    double R2;
    bool isGLC;

    Meson(std::string m, std::string n, double MV_, double mf_, double ef_,
          double NT_, double R2T_, double NL_, double R2L_);

    Meson(std::string m, std::string n, double MV_, double mf_, double ef_,
          double NT_, double NL_, double R2_);
};

// 🔹 Só declaração (extern)
extern Meson Jpsi_GLC_GBW;
extern Meson phi_GLC_GBW;
extern Meson rho_GLC_GBW;
extern Meson Upsilon_GLC_GBW_45;
extern Meson Upsilon_GLC_GBW_42;
extern Meson Jpsi_BG_GBW;
extern Meson phi_BG_GBW;
extern Meson rho_BG_GBW;
extern Meson Upsilon_BG_GBW_45;
extern Meson Upsilon_BG_GBW_42;

extern Meson Jpsi_GLC_ipsat;
extern Meson phi_GLC_ipsat;
extern Meson rho_GLC_ipsat;
extern Meson Jpsi_BG_ipsat;
extern Meson phi_BG_ipsat;
extern Meson rho_BG_ipsat;

struct MesonModelsGBW{
    const Meson& M_GLC;
    const Meson& M_BG;
};

struct MesonModelsipsat{
    const Meson& M_GLC;
    const Meson& M_BG;
};

extern std::map<std::string, MesonModelsGBW> meson_modelsGBW;
extern std::map<std::string, MesonModelsipsat> meson_modelsipsat;


// funções
std::string doubleParaString(double x);

double dfridr(const std::function<double(double)>& func,
              double x,
              double h,
              double& err);

double derivative_richardson(const std::function<double(double)>& f,
                             double x,
                             double h = 1e-4);

double derivative_poly5(const std::function<double(double)>& f,
                        double x);

double x_to_W(double x, const Meson& M);
double W_to_x(double W, const Meson& M);

Meson input_meson(const std::string model);

const Meson& get_meson_bg(const std::string& meson_key, const std::string& dipole_model);

void warn_ipsat_light_meson(const Meson& M, const std::string& dipole_model);

void perfil(const Meson& M);

double B_slope(double x, double Q2, const Meson& M);


void read_sigma_exp(const std::string& filename,
                    std::vector<int>& dataset,
                    std::vector<double>& W,
                    std::vector<double>& sigma,
                    std::vector<double>& error);

std::string timestamp();

void read_two_columns(const std::string& filename,
                      std::vector<double>& col1,
                      std::vector<double>& col2);

void read_csv(const std::string& filename,
              std::vector<double>& col1,
              std::vector<double>& col2,
              std::vector<double>& col3);


void read_xyscan_csv(const std::string& filename,
                    std::vector<double>& x,
                    std::vector<double>& y,
                    std::vector<double>& error);



std::string get_meson();

std::string vec_to_pylist(const std::vector<double>& vec);
std::string flavorToString(int flavor);
std::string flavorName(int flavor);

void read_rapidity_hepdata(
    const std::string& filename,
    std::vector<double>& y,
    std::vector<double>& dsdy,
    std::vector<double>& err);


double get_amplitude_p(double x, double Delta, double Q2, const Meson& M, std::string modelo);

double get_Np(double r, double x, std::string modelo, double b= 2.0);

double get_dipolo_p(double r, double x, double Delta, std::string modelo);

std::string N_file(double x, std::string modelo);

std::vector<double> W_space(double N, const Meson& M);

double low_x_factor(double x, double exp = 5.0);

#endif