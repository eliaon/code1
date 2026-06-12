#include "bCGC.hpp"

#include "../calculations/wavefunctions.hpp"
#include "../other/ctes.hpp"
#include "../other/integration.hpp"

#include <boost/math/special_functions/bessel.hpp>
#include <algorithm>
#include <cmath>

namespace bCGC {

const parametros_bCGC bg_fitted(
    "bCGC_BG",
    0.3159,
    0.6924,
    0.2039,
    0.001849,
    5.746,
    9.9
);

namespace {
constexpr double bmax = 10.0;
constexpr int Nb = 200;

} // namespace

double Qs(double x, double b, const parametros_bCGC& params)
{
    const double x_safe = std::max(x, 1e-12);
    const double qsfac1 = std::pow(params.x0 / x_safe, params.lambda / 2.0);
    const double qsfac2 = std::pow(
        std::exp(-(b * b) / (2.0 * params.Bcgc)),
        1.0 / (2.0 * params.gamma_s)
    );
    return qsfac1 * qsfac2;
}

double N_p(double r, double x, double b, const parametros_bCGC& params)
{
    const double x_safe = std::max(x, 1e-12);
    const double y = std::log(1.0 / x_safe);

    if (y <= 0.0) {
        return 0.0;
    }

    const double qs = Qs(x_safe, b, params);
    const double rqs = std::max(r * qs, 1e-300);

    const double A = -std::pow(params.N0 * params.gamma_s, 2.0)
        / (std::pow(1.0 - params.N0, 2.0) * std::log(1.0 - params.N0));
    const double B = 0.5 * std::pow(
        1.0 - params.N0,
        -(1.0 - params.N0) / (params.N0 * params.gamma_s)
    );

    double N = 0.0;
    if (rqs > 2.0) {
        N = 1.0 - std::exp(-A * std::pow(std::log(B * rqs), 2.0));
    } else {
        const double exponent = 2.0 * (
            params.gamma_s
            + std::log(2.0 / rqs) / (params.kappa * params.lambda * y)
        );
        N = params.N0 * std::pow(rqs / 2.0, exponent);
    }

    return N * low_x_factor(x_safe, 6.0);
}

double sigma_qq_p(double r, double x, double Delta, const parametros_bCGC& params)
{
    auto b_integrand = [r, x, Delta, &params](double b) {
        return 2.0 * M_PI * b
            * N_p(r, x, b, params)
            * boost::math::cyl_bessel_j(0, Delta * b);
    };

    return 2.0 * integrate_simpson(b_integrand, 0.0, bmax, Nb);
}

double amplitude_p(double x, double Delta, double Q2, const Meson& M,
                   const parametros_bCGC& params)
{
    auto r_integrand = [x, Delta, Q2, &M, &params](double r) {
        const double Ov = overlap_r(r, Q2, M);
        const double sigma_dip = sigma_qq_p(r, x, Delta, params);
        return 0.5 * r * Ov * sigma_dip;
    };

    return integrate_simpson(r_integrand, rmin, rmax, Nr);
}

} // namespace bCGC
