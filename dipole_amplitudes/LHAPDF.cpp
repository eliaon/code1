#include "LHAPDF.hpp"
#include <LHAPDF/LHAPDF.h>
#include "../other/utils.hpp"
#include "../other/ctes.hpp"
#include "../other/plots/plot_other.hpp"
#include "../other/integration.hpp"
#include "../calculations/wavefunctions.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

LHAnPDF::LHAnPDF(const std::string& setname, int member)
{
    pdf = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(setname, member));

    // BGBK/BGK-like defaults. sigma0 must be in GeV^-2 because r is in GeV^-1.
    sigma0 = 23.0 * mb_to_gev2;
    C = 0.29;
    mu0_2 = 1.85; // GeV^2
}

LHAnPDF::~LHAnPDF() = default;

// --- PDF ---
double LHAnPDF::gluon(double x, double Q2) const {
    double Q = std::sqrt(Q2);
    return pdf->xfxQ(21, x, Q);
}

double LHAnPDF::alpha_s(double Q2) const {
    return pdf->alphasQ(std::sqrt(Q2));
}

// --- escala ---
double LHAnPDF::mu2(double r) const {
    return C/(r*r) + mu0_2;
}

// --- dipolo ---
double LHAnPDF::sigma_qq_p(double x, double r) const {
    double mu2_val = mu2(r);

    double xg = gluon(x, mu2_val);
    double as = alpha_s(mu2_val);

    double arg = (M_PI * M_PI / (3.0 * sigma0)) * r * r * as * xg;

    return sigma0 * (1.0 - std::exp(-arg));
}

double LHAnPDF::N_p(double x, double r) const {
    return sigma_qq_p(x, r) / sigma0;
}

void LHAnPDF::N_csv()
{

    const int Npoints = 5000;
    double rmin = 1e-4, rmax = 10.0;

    std::vector<double> x_values = {1e-5, 1e-4, 1e-3, 1e-2};

    std::vector<std::string> filenames;

    for (double x : x_values)
    {
        std::string filename = "csv/N_LHAPDF_x_" + doubleParaString(x) + ".csv";
        std::ofstream fout(filename);

        for (int i = 0; i < Npoints; ++i)
        {
            double frac = (double)i / (Npoints - 1);
            double r = rmin * std::pow(rmax / rmin, frac); // log spacing

            double N_val = N_p(x, r);

            fout << r * CFAC << "," << N_val << "\n"; // r em fm
        }
        filenames.push_back(filename);
        std::cout << "Gerado csv para N(r) em x = " << x << ": " << filename << std::endl;
        fout.close();
    }
    plot_N_multi(filenames, x_values, "LHAPDF - CT14lo");

}

std::string LHAnPDF::xf_vs_x(double Q2, int flavor)
{
    int nflavor = flavor; // 21 para xg, 1: xd, 2: xu, 3: xs, 4: xc, 5: xb, 6: xt
    std::string fstring = flavorName(nflavor);
    std::string xf_string = flavorToString(flavor);

    std::string filename = "csv/"+fstring+"_LHAPDF_Q2=" + doubleParaString(Q2) + ".csv";
    std::ofstream fout(filename);
    fout << "x," << xf_string << "\n";
    for(double x = 1e-5; x < 1.0; x *= 1.2) {
    double Q = sqrt(Q2);
    double xf = pdf->xfxQ(nflavor, x, Q); 
    fout << x << "," << xf << "\n";
    std::cout << x << " " << xf << std::endl;
}
    return filename;
}

double LHAnPDF::amplitude_p(double x, double Q2, const Meson& M, bool fc)
{
    // Implementação simplificada da amplitude de difração
    // A amplitude real envolveria a convolução da função de onda do méson com a função de dipolo N(r)

    auto amp_r = [&](double r) {
        double Ov = overlap_r(r, Q2, M, fc);
        double sigma_qq = sigma_qq_p(x, r);
        return 0.5 * r * Ov * sigma_qq; // r de d²r = 2π r dr/4π
    };

    const int Nr = 1000;
    double rmin = 1e-4, rmax = 10.0;
    double amp = integrate_simpson(amp_r, rmin, rmax, Nr);

    return amp;
}

// --- setters ---
void LHAnPDF::set_sigma0(double val) { sigma0 = val; }
void LHAnPDF::set_C(double val) { C = val; }
void LHAnPDF::set_mu0_2(double val) { mu0_2 = val; }
