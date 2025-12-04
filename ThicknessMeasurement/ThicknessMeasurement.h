#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Global parameter configuration (consistent with original simulation/experiment parameters, can be adjusted according to actual needs)
const double MU = 660.0;         // Gaussian light source center wavelength (nm)
const double ETA = 80.0;         // Light source wavelength standard deviation (nm)
const double C = 18000.0;        // Maximum spectral intensity amplitude
const double PHI0 = M_PI_2;      // Initial phase (rad) (using M_PI_2 instead of M_PI/2)
const double DEFAULT_REFRACTIVE = 1.5; // Default film refractive index (e.g., PVC coating)
const double LAMBDA_MIN = 500.0; // Spectrometer minimum wavelength (nm)
const double LAMBDA_MAX = 1000.0;// Spectrometer maximum wavelength (nm)
const int PIXEL_NUM = 512;       // Spectrometer pixel count (wavelength sampling points)
const int EMD_MAX_ITER = 100;    // Maximum iterations for EMD sifting
const double Z_SEARCH_MIN = 0.0; // LSP optical path difference search minimum (μm)
const double Z_SEARCH_MAX = 100.0;// LSP optical path difference search maximum (μm)
const double Z_SEARCH_STEP = 0.01;// LSP optical path difference search step (μm)

/**
 * @brief Simulate interference intensity signal (Equation 5)
 * @param A Reflection coefficient parameter A (0<A≤1)
 * @param B Reflection coefficient parameter B (0<B≤1)
 * @param z Optical path difference (μm)
 * @param lambda Wavelength array (nm)
 * @return Intensity array I_r
 */
VectorXd simulateInterferenceIntensity(double A, double B, double z, const VectorXd& lambda);

/**
 * @brief EMD filtering: Extract first IMF (high-frequency interference components, formula-related processing)
 * @param signal Original intensity signal
 * @return IMF1: First intrinsic mode function
 */
VectorXd emdFilter(const VectorXd& signal);

/**
 * @brief Wavelength to wavenumber transformation (λ→σ=1/λ, solving non-uniform sampling problem)
 * @param lambda Wavelength array (nm)
 * @param imf1 Wavelength domain IMF1 signal
 * @param sigma Output wavenumber array (1/μm)
 * @param imf1_sigma Output wavenumber domain IMF1 signal
 */
void lambdaToSigma(const VectorXd& lambda, const VectorXd& imf1, VectorXd& sigma, VectorXd& imf1_sigma);

/**
 * @brief Lomb-Scargle periodogram analysis (Equations 2, 3), extract optical path difference z_peak corresponding to power spectrum peak
 * @param sigma Wavenumber array (1/μm)
 * @param imf1_sigma Wavenumber domain IMF1 signal
 * @return z_peak: Optical path difference corresponding to power spectrum peak (μm)
 */
double lombScargle(const VectorXd& sigma, const VectorXd& imf1_sigma);

/**
 * @brief Calculate film thickness (Equation 4: d = z_peak / n)
 * @param z_peak Optical path difference corresponding to LSP peak (μm)
 * @param refractive Film refractive index (default 1.5, can be modified according to material)
 * @return d: Film thickness (μm)
 */
double calculateThickness(double z_peak, double refractive = DEFAULT_REFRACTIVE);

// Function implementations

VectorXd simulateInterferenceIntensity(double A, double B, double z, const VectorXd& lambda) {
    Eigen::Index n = lambda.size();
    VectorXd I(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        double lambda_μm = lambda(i) * 1e-3; // Convert wavelength to μm
        double k = 2 * M_PI / lambda_μm;     // Wave number (1/μm)
        double cos_term = cos(k * z + PHI0); // Cosine interference term

        // Gaussian light source intensity distribution G(λ)
        double gauss = exp(-pow(lambda(i) - MU, 2) / (2 * pow(ETA, 2)))
            / (ETA * sqrt(2 * M_PI));

        // Calculate reflected interference intensity (Equation 5)
        I(i) = (A + B * cos_term) * C * gauss;
    }
    return I;
}

VectorXd emdFilter(const VectorXd& signal) {
    Eigen::Index n = signal.size();
    VectorXd x = signal;
    VectorXd imf1(n);

    for (int iter = 0; iter < EMD_MAX_ITER; ++iter) {
        // 1. Extract local maximum and minimum indices
        vector<Eigen::Index> max_idx, min_idx;
        for (Eigen::Index i = 1; i < n - 1; ++i) {
            if (x(i) > x(i - 1) && x(i) > x(i + 1)) max_idx.push_back(i);
            if (x(i) < x(i - 1) && x(i) < x(i + 1)) min_idx.push_back(i);
        }
        if (max_idx.size() < 2 || min_idx.size() < 2) break; // Insufficient extrema, terminate sifting

        // 2. Linear interpolation to fit upper and lower envelopes (simplified implementation, can be replaced with cubic spline)
        VectorXd upper_env(n), lower_env(n);
        // Upper envelope interpolation
        for (size_t i = 0; i < max_idx.size() - 1; ++i) {
            Eigen::Index left = max_idx[i];
            Eigen::Index right = max_idx[i + 1];
            double slope = (x(right) - x(left)) / (right - left);
            for (Eigen::Index j = left; j <= right; ++j) {
                upper_env(j) = x(left) + slope * (j - left);
            }
        }
        // Lower envelope interpolation
        for (size_t i = 0; i < min_idx.size() - 1; ++i) {
            Eigen::Index left = min_idx[i];
            Eigen::Index right = min_idx[i + 1];
            double slope = (x(right) - x(left)) / (right - left);
            for (Eigen::Index j = left; j <= right; ++j) {
                lower_env(j) = x(left) + slope * (j - left);
            }
        }
        // Boundary filling (using values at first and last extrema)
        for (Eigen::Index j = 0; j < max_idx[0]; ++j) upper_env(j) = x(max_idx[0]);
        for (Eigen::Index j = max_idx.back() + 1; j < n; ++j) upper_env(j) = x(max_idx.back());
        for (Eigen::Index j = 0; j < min_idx[0]; ++j) lower_env(j) = x(min_idx[0]);
        for (Eigen::Index j = min_idx.back() + 1; j < n; ++j) lower_env(j) = x(min_idx.back());

        // 3. Calculate envelope mean and sift
        VectorXd m = (upper_env + lower_env) / 2;
        VectorXd h = x - m;
        imf1 = h;
        x = m; // Remaining signal as input for next sifting iteration
    }
    return imf1;
}

void lambdaToSigma(const VectorXd& lambda, const VectorXd& imf1, VectorXd& sigma, VectorXd& imf1_sigma) {
    Eigen::Index n = lambda.size();
    sigma.resize(n);
    imf1_sigma.resize(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        double lambda_μm = lambda(i) * 1e-3; // Convert wavelength to μm
        sigma(i) = 1.0 / lambda_μm;          // Wavenumber σ=1/λ (1/μm)
        imf1_sigma(i) = imf1(i);             // Intensity value mapped to wavenumber domain
    }
}

double lombScargle(const VectorXd& sigma, const VectorXd& imf1_sigma) {
    Eigen::Index z_num = static_cast<Eigen::Index>((Z_SEARCH_MAX - Z_SEARCH_MIN) / Z_SEARCH_STEP + 1);
    VectorXd z_list(z_num);
    VectorXd P(z_num); // LSP power spectrum

    // Iterate over all candidate optical path differences z, calculate power spectrum
    for (Eigen::Index z_idx = 0; z_idx < z_num; ++z_idx) {
        double z = Z_SEARCH_MIN + z_idx * Z_SEARCH_STEP;
        z_list(z_idx) = z;

        // Calculate time-shift invariant parameter δ (Equation 3)
        double sum_sin = 0.0, sum_cos = 0.0;
        for (Eigen::Index i = 0; i < sigma.size(); ++i) {
            double arg = 4 * M_PI * z * sigma(i);
            sum_sin += sin(arg);
            sum_cos += cos(arg);
        }
        
        double delta = 0.0;
        if (z != 0.0) {
            delta = atan2(sum_sin, sum_cos) / (4 * M_PI * z);
        } else {
            delta = 0.0; // Avoid division by zero when z=0
        }

        // Calculate LSP power spectrum P(z) (Equation 2)
        double sum_cos_term = 0.0, sum_sin_term = 0.0;
        double sum_cos2 = 0.0, sum_sin2 = 0.0;
        for (Eigen::Index i = 0; i < sigma.size(); ++i) {
            double arg = 2 * M_PI * z * (sigma(i) - delta);
            double cos_arg = cos(arg);
            double sin_arg = sin(arg);

            sum_cos_term += imf1_sigma(i) * cos_arg;
            sum_sin_term += imf1_sigma(i) * sin_arg;
            sum_cos2 += cos_arg * cos_arg;
            sum_sin2 += sin_arg * sin_arg;
        }

        // Avoid division by zero (theoretically impossible, adding protection)
        sum_cos2 = max(sum_cos2, 1e-10);
        sum_sin2 = max(sum_sin2, 1e-10);
        double P_cos = pow(sum_cos_term, 2) / sum_cos2;
        double P_sin = pow(sum_sin_term, 2) / sum_sin2;
        P(z_idx) = 0.5 * (P_cos + P_sin);
        
        // Ensure no NaN values
        if (isnan(P(z_idx)) || isinf(P(z_idx))) {
            P(z_idx) = -1.0; // Set invalid values to a very low value
        }
    }

    // Extract z_peak corresponding to power spectrum peak
    Eigen::Index peak_idx;
    P.maxCoeff(&peak_idx);
    return z_list(peak_idx);
}

double calculateThickness(double z_peak, double refractive) {
    return z_peak / (2 * refractive);
}
