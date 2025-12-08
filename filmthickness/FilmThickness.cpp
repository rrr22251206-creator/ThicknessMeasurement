#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <iostream>
#include "FilmThickness.hpp"

// ======================================================================
// SpectrumModel
// ======================================================================
double SpectrumModel::gaussian(double lambda, double mu, double sigma)
{
    return exp(-(lambda - mu) * (lambda - mu) / (2 * sigma * sigma));
}

double SpectrumModel::interference(
    double lambda,
    double A, double B, double C,
    double mu, double eta,
    double phi0,
    double n1,
    double d
)
{
    double phase = 4.0 * M_PI * n1 * d / lambda;
    double env = C * gaussian(lambda, mu, eta);
    return (A + B * cos(phase + phi0)) * env;
}

// ======================================================================
// EMD IMF1 (工程版：零相位平滑 + 去基线)
// ======================================================================
std::vector<double> EMD::highpass(const std::vector<double>& x)
{
    int N = x.size();
    const int w = 15; // 31 点平滑窗（零相位）
    std::vector<double> y(N);

    for (int i = 0; i < N; i++) {
        double sum = 0;
        int count = 0;

        for (int k = -w; k <= w; k++) {
            int idx = i + k;
            if (idx >= 0 && idx < N) {
                sum += x[idx];
                count++;
            }
        }

        double smooth = sum / count;
        y[i] = x[i] - smooth;  // 高频成分 = 原信号 − 平滑
    }

    return y;
}

std::vector<double> EMD::extractIMF1(const std::vector<double>& x, EMDMode)
{
    return highpass(x);
}

// ======================================================================
// Lomb–Scargle
// ======================================================================
static double computePower(
    double z,
    const std::vector<double>& sigma,
    const std::vector<double>& I)
{
    double sum_c = 0, sum_s = 0, sum_cc = 0, sum_ss = 0;

    for (size_t i = 0; i < sigma.size(); i++) {
        double w = 2 * M_PI * sigma[i] * z;
        double c = cos(w);
        double s = sin(w);

        sum_c += I[i] * c;
        sum_s += I[i] * s;

        sum_cc += c * c;
        sum_ss += s * s;
    }

    return (sum_c * sum_c) / sum_cc + (sum_s * sum_s) / sum_ss;
}

double LombScargle::findPeak(
    const std::vector<double>& sigma,
    const std::vector<double>& I,
    double z_min, double z_max, double step)
{
    double best_z = 0;
    double best_P = -1;

    for (double z = z_min; z <= z_max; z += step) {
        double Pz = computePower(z, sigma, I);
        if (Pz > best_P) {
            best_P = Pz;
            best_z = z;
        }
    }
    return best_z;
}

// ======================================================================
// Peak refinement
// ======================================================================
static double refinePeak(
    double z0,
    const std::vector<double>& sigma,
    const std::vector<double>& I)
{
    const double dz = 1e-8;

    double P1 = computePower(z0 - dz, sigma, I);
    double P2 = computePower(z0, sigma, I);
    double P3 = computePower(z0 + dz, sigma, I);

    double denom = (P1 - 2 * P2 + P3);
    if (fabs(denom) < 1e-20)
        return z0;

    return z0 + (dz * (P1 - P3)) / (2 * denom);
}

// ======================================================================
// FilmThicknessSolver
// ======================================================================
FilmThicknessSolver::FilmThicknessSolver(double refrIndex)
    : n1(refrIndex) {
}

double FilmThicknessSolver::computeThickness(
    const std::vector<double>& lambda,
    const std::vector<double>& intensity,
    EMDMode mode)
{
    auto imf1 = EMD::extractIMF1(intensity, mode);

    std::vector<double> sigma(lambda.size());
    for (size_t i = 0; i < lambda.size(); i++)
        sigma[i] = 1.0 / lambda[i];

    double z_peak = LombScargle::findPeak(
        sigma, imf1,
        0.0, 60e-6, 1e-8);

    std::cout << "[DEBUG] z_peak coarse = " << z_peak << "\n";

    z_peak = refinePeak(z_peak, sigma, imf1);

    std::cout << "[DEBUG] z_peak refined = " << z_peak << "\n";

    return z_peak / (2.0 * n1);
}
