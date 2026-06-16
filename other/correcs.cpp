
#include "correcs.hpp"
#include "integration.hpp"
#include "utils.hpp"
#include "ctes.hpp"
#include "../calculations/wavefunctions.hpp"
#include "../calculations/sigma.hpp"
#include "../dipole_amplitudes/GBW.hpp"
#include "../libraries/mantysaari/dipoleamplitude.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <LHAPDF/LHAPDF.h>


using namespace MZ_ipsat;

double lnA(double y, double Delta, double Q2, const Meson& M,
           std::string dipolemodel);

namespace {

constexpr double lambda_fortran_bcgc   = 0.28;
constexpr double lambda_derivative_step = 1e-3;
constexpr double lambda_boundary_x      = 1e-2;
constexpr double lambda_min             = 0.05;
constexpr double lambda_max             = 0.45;

std::string normalize_model_name(std::string model)
{
    std::transform(model.begin(), model.end(), model.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return model;
}

bool uses_fortran_fixed_lambda(const std::string& dipolemodel)
{
    return normalize_model_name(dipolemodel) == "bcgc";
}

bool has_gbw_lambda(const std::string& dipolemodel)
{
    const std::string m = normalize_model_name(dipolemodel);
    return m == "gbw(old)" || m == "gbw(new)" || m == "gbw";
}

double gbw_lambda_for_model(const std::string& dipolemodel)
{
    const std::string m = normalize_model_name(dipolemodel);
    if (m == "gbw(old)" || m == "gbw") {
        return gbw.lambda;
    }
    if (m == "gbw(new)") {
        return gbw_10.lambda;
    }
    return gbw_10.lambda;
}

double clamp_lambda(double lambda)
{
    return std::clamp(lambda, lambda_min, lambda_max);
}

bool is_valid_lambda(double lambda)
{
    return std::isfinite(lambda) && lambda > 0.0;
}

// λ_e = ∂ ln|A| / ∂ ln(1/x), com y = ln(1/x) = -ln(x).
double lambda_from_amplitude(double x, double Delta, double Q2, const Meson& M,
                             const std::string& dipolemodel)
{
    double err = 0.0;
    const double x_eval = (x > lambda_boundary_x) ? lambda_boundary_x : x;
    const double y = -std::log(x_eval);

    auto f_lnA = [&](double y_var) {
        return lnA(y_var, Delta, Q2, M, dipolemodel);
    };

    return dfridr(f_lnA, y, lambda_derivative_step, err);
}

} // namespace

double lnA(double y, double Delta, double Q2, const Meson& M,
           std::string dipolemodel)
{
    const double x = std::exp(-y);
    const double amp = get_amplitude_p(x, Delta, Q2, M, dipolemodel);

    if (amp == 0.0) {
        std::cerr << "Amp = 0 em x=" << x << std::endl;
        return -1e30;
    }

    return std::log(std::abs(amp));
}

double calculate_lambda(double x, double Delta, double Q2, const Meson& M,
                        std::string dipolemodel)
{
    if (uses_fortran_fixed_lambda(dipolemodel)) {
        return lambda_fortran_bcgc;
    }

    const double lambda_amp = lambda_from_amplitude(x, Delta, Q2, M, dipolemodel);
    if (is_valid_lambda(lambda_amp)) {
        return clamp_lambda(lambda_amp);
    }

    if (has_gbw_lambda(dipolemodel)) {
        return clamp_lambda(gbw_lambda_for_model(dipolemodel));
    }

    std::cerr << "AVISO: lambda_e invalido para modelo " << dipolemodel
              << "; usando fallback 0.28.\n";
    return lambda_fortran_bcgc;
}

// R_g(λ) = 2^(2λ+3) / sqrt(π) * Γ(λ+5/2) / Γ(λ+4)  [Shuvaev; Gonçalves et al.]
double RG(double /*x*/, double /*Q2*/, double lambda_e, const Meson& /*M*/)
{
    const double lam = clamp_lambda(lambda_e);
    return std::pow(2.0, 2.0 * lam + 3.0) * tgamma(lam + 2.5)
         / (std::sqrt(M_PI) * tgamma(lam + 4.0));
}

// Re A / Im A ≈ tan(π λ_e / 2) via relações de dispersão [Machado EPJC 2009].
double beta(double /*x*/, double /*Q2*/, double lambda_e, const Meson& /*M*/)
{
    const double lam = clamp_lambda(lambda_e);
    return std::tan(M_PI * lam / 2.0);
}

SkewCorrection compute_skew_correction(double x, double Delta, double Q2,
                                       const Meson& M, const std::string& dipolemodel)
{
    SkewCorrection out;
    out.lambda_e = calculate_lambda(x, Delta, Q2, M, dipolemodel);
    out.Rg       = RG(x, Q2, out.lambda_e, M);
    out.beta     = beta(x, Q2, out.lambda_e, M);
    out.factor   = out.Rg * out.Rg * (1.0 + out.beta * out.beta);
    return out;
}

void debug_correc(std::string model)
{
    double x = 1e-4;
    double Q2 = 0.0;
    const Meson& M = Jpsi_GLC_GBW;

    for (int i = 0; i < 120; ++i) {
        const double xi = x + i * 1e-4;
        const auto corr = compute_skew_correction(xi, 0.0, Q2, M, model);
        std::cout << "x: " << xi
                  << "  lambda_e: " << corr.lambda_e
                  << "  Rg: " << corr.Rg
                  << "  beta: " << corr.beta
                  << "  factor: " << corr.factor << std::endl;
    }
}

double f_c(double r, double B, double omega, double R)
{
    const double omega2 = omega * omega;
    const double fc_num = 1.0 + B * std::exp(-omega2 * (r - R) * (r - R));
    const double fc_den = 1.0 + B * std::exp(-omega2 * R * R);
    return fc_num / fc_den;
}


namespace {
    LHAPDF::PDF* & nnpdf_proton_pdf() {
        thread_local LHAPDF::PDF* pdf = LHAPDF::mkPDF("NNPDF30_nlo_as_0118", 0);
        return pdf;
    }

    LHAPDF::PDF* & nnpdf_pb_pdf() {
        thread_local LHAPDF::PDF* pdf = LHAPDF::mkPDF("nNNPDF30_nlo_as_0118_A208_Z82", 0);
        return pdf;
    }
}

double xg_p(double x, double Q) {
    return nnpdf_proton_pdf()->xfxQ(21, x, Q); // Retorna xg(x, Q^2) = x * g(x, Q^
}

double xg_Pb(double x, double Q){
    return nnpdf_pb_pdf()->xfxQ(21, x, Q); 
}

// ------- shadowing ----------------

double S_Pb(double x, double Q){
    double xg_p_val = xg_p(x, Q);
    double xg_Pb_val = xg_Pb(x, Q);

    return xg_Pb_val / (xg_p_val);
}













