#include "../other/utils.hpp"
#include "../other/ctes.hpp"
#include "../other/integration.hpp"
#include "../other/correcs.hpp"

#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>


// ---------------- phi_T e derivada ----------------
double phi_T(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        double zz = z * (1.0 - z);
        return M.NT * (zz * zz) * exp(-(r * r) / (2.0 * M.R2T));
    } else {
        double zz = z * (1.0 - z);
        double part1 = M.NT * zz;
        double arg1  = -(M.mf*M.mf * M.R2) / (8.0 * zz);
        double arg2  = -(2.0 * zz * r * r) / M.R2;
        double arg3  = (M.mf*M.mf * M.R2) / 2.0;
        return part1 * exp(arg1 + arg2 + arg3);
    }
}
double phi_L(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        double zz = z * (1.0 - z);
        return M.NL * zz * exp(-(r * r) / (2.0 * M.R2L));
    } else {
        double zz = z * (1.0 - z);
        double part1 = M.NL * zz;
        double arg1  = -(M.mf*M.mf * M.R2) / (8.0 * zz);
        double arg2  = -(2.0 * zz * r * r) / M.R2;
        double arg3  = (M.mf*M.mf * M.R2) / 2.0;
        return part1 * exp(arg1 + arg2 + arg3);
    }
}

double dphiT_dr(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        return -(r / M.R2T) * phi_T(r, z, M);
    } else {
        double zz = z * (1.0 - z);
        return -((4.0 * r * zz) / M.R2) * phi_T(r, z, M);
    }
}
double dphiL_dr(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        return -(r / M.R2L) * phi_L(r, z, M);
    } else {
        double zz = z * (1.0 - z);
        return -((4.0 * r * zz) / M.R2) * phi_L(r, z, M);
    }
}

double laplacian_phi_L(double r, double z, const Meson& M)
{
    const double dr = 1e-4;

    if (r < 1e-6)
    {
        double phi_p = phi_L(dr, z, M);
        double phi_0 = phi_L(0.0, z, M);

        double d2phi = (phi_p - phi_0) / (dr * dr);
        return 2.0 * d2phi;
    }

    double phi_p = phi_L(r + dr, z, M);
    double phi_m = phi_L(r - dr, z, M);
    double phi_0 = phi_L(r,       z, M);

    double dphi_dr   = (phi_p - phi_m) / (2.0 * dr);
    double d2phi_dr2 = (phi_p - 2.0 * phi_0 + phi_m) / (dr * dr);

    return d2phi_dr2 + (1.0 / r) * dphi_dr;
}



// ----------------psi_V psi_T ----------------
double psi_Vpsi_T(double z, double r,double Q2,  const Meson& M)
{
    double Mf2 = M.mf * M.mf;
    double EPS2 = z * (1.0 - z) * Q2 + Mf2;
    double EPS  = sqrt(EPS2);

    double K0 = boost::math::cyl_bessel_k(0, EPS * r);
    double K1 = boost::math::cyl_bessel_k(1, EPS * r);

    double PHIT  = phi_T(r, z, M);
    double DPHIT = dphiT_dr(r, z, M);

    double ZZZ = z*z + (1.0 - z)*(1.0 - z);
    double ANORM = M.ef * sqrt(4.0 * M_PI * alfem) * Nc / (M_PI * z * (1.0 - z));

    return ANORM * ((Mf2 * K0 * PHIT) - (ZZZ * EPS * K1 * DPHIT));
}

// ----------------psi_V*psi_L ----------------
double psi_Vpsi_L(double z, double r, double Q2, const Meson& M)
{
    const double mf2 = M.mf * M.mf;
    const double eps2 = z * (1.0 - z) * Q2 + mf2;
    const double eps  = std::sqrt(eps2);

    const double K0 = boost::math::cyl_bessel_k(0, eps * r);

    const double phiL = phi_L(r, z, M);
    const double lap_phiL = laplacian_phi_L(r, z, M);

    const double pref =
        M.ef * std::sqrt(4.0 * M_PI * alfem) * Nc / M_PI;

    const double bracket =
        M.MV * phiL
        +
        (mf2 * phiL - lap_phiL) / (M.MV * z * (1.0 - z));

    return pref
           * 2.0 * std::sqrt(Q2)
           * z * (1.0 - z)
           * K0
           * bracket;
}


// ---------------- overlap integrado em z ----------------
double overlap_r(double r, double Q2, const Meson& M) {
    auto fz = [r, Q2, &M](double z) {
        return psi_Vpsi_T(z, r, Q2, M); //+ psi_Vpsi_L(z, r, Q2, M);
    };
    return integrate_simpson( fz, zmin, zmax, Nz);
}

