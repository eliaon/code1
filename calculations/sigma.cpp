#include <cmath>
#include <iostream>

#include "../dipole_amplitudes/GBW.hpp"
#include "../dipole_amplitudes/LHAPDF.hpp"
#include "../other/correcs.hpp"
#include "../other/utils.hpp"
#include "../other/integration.hpp"

#include <boost/math/special_functions/bessel.hpp>

#include "nuclear.hpp"
#include "wavefunctions.hpp"


namespace gamma_p {

    double amplitude_p(double x, double Delta, double Q2, const Meson& M, std::string modelo)
    {
        auto amp_r = [x, Delta, Q2, &M, &modelo](double r){
            double Ov = overlap_r(r, Q2, M);
            double sigma_qq = get_dipolo_p(r, x, Delta, modelo);
            return 2* M_PI * r * Ov * sigma_qq;
        };
    double amp = integrate_simpson(amp_r, rmin, rmax, Nr);
    return amp/(4.0 * M_PI);
    }


    double sigma(double W, double Q2, const Meson& M, std::string modelo)
    {
        double x = W_to_x(W, M);
        double amp = get_amplitude_p(x, 0.0, Q2, M, modelo);

        double B_val = B_slope(x, Q2, M);
        const SkewCorrection skew = compute_skew_correction(x, 0.0, Q2, M, modelo);

        return skew.factor * (amp * amp) / (16.0 * M_PI * B_val);
    }

    double dsigma_dt(double t, double W, double Q2, const Meson& M, std::string modelo)
    {
        double Delta = std::sqrt(t);
        double x = W_to_x(W, M);
        double amp = get_amplitude_p(x, Delta, Q2, M, modelo);
        const SkewCorrection skew = compute_skew_correction(x, Delta, Q2, M, modelo);

        return (amp * amp) / (16.0 * M_PI) * skew.factor;
    }

    double sigma_integrado(double W, double Q2, const Meson&M, std::string modelo)
    {
        auto ds_dt = [&](double t){
            return dsigma_dt(t, W, Q2, M, modelo);
        };
        return integrate_simpson(ds_dt, tmin, tmax, Nt);
    }
}

// ------------- CASO NUCLEAR ---------------------------
// -------------- MODELO BGBK ---------------------------


namespace gamma_A{
    namespace {
        bool is_proton_target(const Nucleus& nucleo)
        {
            return nucleo.name == "p" || nucleo.name == "proton";
        }
    }

    double N_qq_A(double r, double x, double b, double Delta ,std::string modelo, 
                    Nucleus& Nucleo, const TA_Table& table)
    {
        if (is_proton_target(Nucleo)){return get_Np(r, x, modelo, b);}
        double sigmaqq_p = get_dipolo_p(r, x, Delta, modelo);
        double TA = interpolate_TA(b, table); 
        double arg = 0.5 * TA * sigmaqq_p;
        //cout << "sigma_qq_p = " << sigmaqq_p <<", TA = " <<TA <<", arg = " <<arg <<endl;
        return (1.0 - exp(-arg));   
    }

    double amplitude_A(double x, double b, double Delta, double Q2, const Meson&M,
                     std::string modelo, Nucleus& Nucleo, const TA_Table& table)
    {
        auto amp_r = [x, Q2, b, Delta, modelo, M, &Nucleo, &table](double r) {
        double sqrt_fc = std::sqrt(f_c(r, -0.979599, 0.403569)); // valores de B e omega fitados no minuit para o proton
        double Ov = overlap_r(r, Q2, M);
        double N_qq = N_qq_A(r, x, b, Delta, modelo, Nucleo, table); 
        //cout << "Overlap = " << Ov << ", sigma_qq = " << sigma_qq << endl;
        return 2 * M_PI * r * Ov * N_qq * sqrt_fc; // r de d²r = 2π r dr => 2π r * Ov * N_qq para integrar em dr depois
    };
    return integrate_simpson(amp_r, rmin, rmax, Nr)/(4 * M_PI);
    }

    double sigma(double W, double Q2, const Meson&M, std::string modelo, const TA_Table& table, Nucleus& Nucleo)
    {
        double x = W_to_x(W, M);
        auto sig_b = [&](double b){
            double amp = amplitude_A(x, b, 0.0, Q2, M, modelo, Nucleo, table);
            //cout << "Amplitude A = " << amp << endl;
            double amp2 = amp * amp;
            return 2 * M_PI * b * amp2;
        };
    double bmax = table.b_vals.back(); // limite de integração em b
    int nb = table.b_vals.size(); // número de pontos para integração em b
    //cout << "bmax = "<< bmax <<", nb =" << nb << endl;
    //const SkewCorrection skew = compute_skew_correction(x, 0.0, Q2, M, modelo);
    return integrate_simpson(sig_b, 0.0, bmax, nb); //* skew.factor;
    }

    double dN_domega(double omega, double sqrt_s, Nucleus&Nucleo, double bmin)
{
    double Z = Nucleo.Z; // número atômico do átomo
    double Z2 = Z * Z;
    double alpha = 1.0 / 137.0; // constante de estrutura fina
    double gamma = sqrt_s / (2.0 * m_p); // fator de Lorentz para o núcleo (massa em GeV)
    double ksi = (omega * bmin) / gamma;
    double K0 = boost::math::cyl_bessel_k(0, ksi);
    double K1 = boost::math::cyl_bessel_k(1, ksi);
    double coeff = (2 * Z2 * alpha) / (M_PI * omega);
    double flux = coeff * (ksi * K0 * K1 - (ksi * ksi / 2.0) * (K1 * K1 - K0 * K0));
    //std::cout << "gamma = " << gamma << std::endl;
    //std::cout << "ksi = " << ksi << std::endl;
    return flux;
    
}

double d_sigma_dy_AA(double y, double sqrt_s, double Q2, const Meson& M, std::string model,
                    TA_Table& table, Nucleus& Nucleo)
{

    double omega_plus = (M.MV / 2.0) * exp(y);
    double omega_minus = (M.MV / 2.0) * exp(-y);

    double W_plus = std::sqrt(2 * omega_plus * sqrt_s);
    double W_minus = std::sqrt(2 * omega_minus * sqrt_s);

    double bmin = 2 * Nucleo.R; // b mínimo para evitar sobreposição dos núcleos

    double n_plus  = omega_plus  * dN_domega(omega_plus, sqrt_s, Nucleo, bmin);
    double n_minus = omega_minus * dN_domega(omega_minus, sqrt_s, Nucleo, bmin);

    //std::cout<<"omega = " << omega_plus << "n(omega)=" << n_plus << std::endl;

    double sigma_plus  = sigma(W_plus, Q2, M, model, table, Nucleo);
    double sigma_minus = sigma(W_minus, Q2, M, model, table, Nucleo);
    

    return n_plus * sigma_plus + n_minus * sigma_minus;
}

double d_sigma_dy_AB(double y, double sqrt_s, double Q2, const Meson& M, std::string model,
                    TA_Table& tableA, TA_Table& tableB,
                    Nucleus& NucleoA, Nucleus& NucleoB)
{

    double omega_plus = (M.MV / 2.0) * exp(y);
    double omega_minus = (M.MV / 2.0) * exp(-y);

    double W_plus = std::sqrt(2 * omega_plus * sqrt_s);
    double W_minus = std::sqrt(2 * omega_minus * sqrt_s);

    double bmin = NucleoA.R + NucleoB.R; // soma dos raios dos núcleos como b mínimo para evitar sobreposição



    double n_plus  = omega_plus  * dN_domega(omega_plus, sqrt_s, NucleoB, bmin);
    double n_minus = omega_minus * dN_domega(omega_minus, sqrt_s, NucleoA, bmin);

    //std::cout<<"omega = " << omega_plus << "n(omega)=" << n_plus << std::endl;

    double sigma_plus  = sigma(W_plus, Q2, M, model, tableB, NucleoB);
    double sigma_minus = sigma(W_minus, Q2, M, model, tableA, NucleoA);
    

    return n_plus * sigma_plus + n_minus * sigma_minus;
}

}
