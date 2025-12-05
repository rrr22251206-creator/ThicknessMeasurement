#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include "FilmThickness.hpp"

//
// =====================================================
//               SpectrumModel
// =====================================================
double SpectrumModel::gaussian(double lambda, double mu, double sigma) {
    return exp(-(lambda - mu) * (lambda - mu) / (2 * sigma * sigma));
}

double SpectrumModel::interference(double lambda,
    double A, double B, double C,
    double mu, double eta,
    double phi0, double z)
{
    double k = 2.0 * M_PI / lambda;      // k = 2π/λ
    double env = C * gaussian(lambda, mu, eta);
    return (A + B * cos(k * z + phi0)) * env;   // Eq.(5)
}

//
// =====================================================
//               EMD (High-pass substitute)
// =====================================================
std::vector<double> EMD::highpass(const std::vector<double>& x)
{
    double alpha = 0.95;
    std::vector<double> y(x.size());
    y[0] = x[0];
    for (size_t i = 1; i < x.size(); i++)
        y[i] = alpha * (y[i - 1] + x[i] - x[i - 1]);
    return y;
}

std::vector<double> EMD::extractIMF1(const std::vector<double>& x, EMDMode)
{
    return highpass(x);  // Approximate IMF1
}

//
// =====================================================
//           Lomb–Scargle (Eq.(2) & Eq.(3))
// =====================================================
static double computePower(
    double z,
    const std::vector<double>& sigma,
    const std::vector<double>& I)
{
    double sum_c = 0, sum_s = 0, sum_cc = 0, sum_ss = 0;

    for (size_t i = 0; i < sigma.size(); i++) {
        double w = 2 * M_PI * sigma[i] * z;     // 2πσz
        double c = cos(w);
        double s = sin(w);
        sum_c += I[i] * c;
        sum_s += I[i] * s;
        sum_cc += c * c;
        sum_ss += s * s;
    }

    // Power P(z) — Eq.(2)
    return (sum_c * sum_c) / sum_cc + (sum_s * sum_s) / sum_ss;
}

double LombScargle::findPeak(
    const std::vector<double>& sigma,
    const std::vector<double>& I,
    double z_min, double z_max, double step)
{
    double best_z = 0;
    double best_power = -1;

    for (double z = z_min; z <= z_max; z += step) {
        double Pz = computePower(z, sigma, I);

        if (Pz > best_power) {
            best_power = Pz;
            best_z = z;
        }
    }
    return best_z;
}

//
// =====================================================
//               Peak Refinement (NEW)
//     Using 3-point quadratic interpolation
// =====================================================
static double refinePeak(
    double z0,
    const std::vector<double>& sigma,
    const std::vector<double>& I)
{
    const double dz = 1e-8;  // same as main scan

    double P1 = computePower(z0 - dz, sigma, I);
    double P2 = computePower(z0, sigma, I);
    double P3 = computePower(z0 + dz, sigma, I);

    double denom = (P1 - 2 * P2 + P3);
    if (fabs(denom) < 1e-20)
        return z0;

    // Parabola vertex formula
    double refined = z0 + (dz * (P1 - P3)) / (2 * denom);

    return refined;
}

//
// =====================================================
//               FilmThicknessSolver
// =====================================================
FilmThicknessSolver::FilmThicknessSolver(double refrIndex)
    : n1(refrIndex) {
}

double FilmThicknessSolver::computeThickness(
    const std::vector<double>& lambda,
    const std::vector<double>& intensity,
    EMDMode mode)
{
    // Step 1: Extract IMF1
    auto imf1 = EMD::extractIMF1(intensity, mode);

    // Step 2: Convert λ → σ = 1/λ
    std::vector<double> sigma(lambda.size());
    for (size_t i = 0; i < lambda.size(); i++)
        sigma[i] = 1.0 / lambda[i];

    // Step 3: Lomb–Scargle coarse scan
    double z_peak = LombScargle::findPeak(
        sigma, imf1,
        0e-6, 20e-6, 1e-8
    );

    // Step 4: Peak refinement (NEW)
    z_peak = refinePeak(z_peak, sigma, imf1);

    // Step 5: Convert z → thickness d = z / n1 (Eq.(4))
    return z_peak / n1;
}
