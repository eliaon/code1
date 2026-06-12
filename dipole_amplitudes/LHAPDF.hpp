
#ifndef LHAPDF_HPP
#define LHAPDF_HPP

#include <memory>
#include <string>
#include <LHAPDF/LHAPDF.h> // Make sure LHAPDF is installed and your includePath is set correctly

struct Meson;

class LHAnPDF {
private:
    std::unique_ptr<LHAPDF::PDF> pdf;

    double sigma0;
    double C;
    double mu0_2;

public:
    LHAnPDF(const std::string& setname, int member = 0);
    ~LHAnPDF();

    double gluon(double x, double Q2) const;
    double alpha_s(double Q2) const;
    double mu2(double r) const;
    double sigma_qq_p(double x, double r) const;
    double N_p(double x, double r) const;

    void N_csv();
    std::string xf_vs_x(double Q2, int flavor);
    double amplitude_p(double x, double Q2, const Meson& M);

    void set_sigma0(double val);
    void set_C(double val);
    void set_mu0_2(double val);
};

#endif