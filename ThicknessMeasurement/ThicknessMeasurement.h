#pragma once
#include <Eigen/Dense>

using Eigen::VectorXd;

// ================= 常量声明 =================
extern const double MU;                  // Gaussian center wavelength (nm)
extern const double ETA;                 // Gaussian std deviation (nm)
extern const double C;                   // Intensity amplitude
extern const double PHI0;                // Initial phase (rad)
extern const double DEFAULT_REFRACTIVE;  // Default refractive index

extern const double LAMBDA_MIN;           // nm
extern const double LAMBDA_MAX;           // nm
extern const int    PIXEL_NUM;

extern const int    EMD_MAX_ITER;
extern const double Z_SEARCH_MIN;         // μm
extern const double Z_SEARCH_MAX;         // μm
extern const double Z_SEARCH_STEP;        // μm

// ================= 算法接口 =================

VectorXd simulateInterferenceIntensity(
    double A, double B, double z,
    const VectorXd& lambda);

VectorXd emdFilter(
    const VectorXd& signal);

void lambdaToSigma(
    const VectorXd& lambda,
    const VectorXd& imf1,
    VectorXd& sigma,
    VectorXd& imf1_sigma);

double lombScargle(
    const VectorXd& sigma,
    const VectorXd& imf1_sigma);

double calculateThickness(
    double z_peak,
    double refractive = DEFAULT_REFRACTIVE);
