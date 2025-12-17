#define _USE_MATH_DEFINES
#include "ThicknessMeasurement.h"
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Eigen;

// ================= 常量定义（只在此文件） =================
const double MU = 660.0;
const double ETA = 80.0;
const double C = 18000.0;
const double PHI0 = M_PI_2;
const double DEFAULT_REFRACTIVE = 1.5;

const double LAMBDA_MIN = 500.0;
const double LAMBDA_MAX = 1000.0;
const int    PIXEL_NUM = 512;

const int    EMD_MAX_ITER = 100;
const double Z_SEARCH_MIN = 0.0;
const double Z_SEARCH_MAX = 100.0;
const double Z_SEARCH_STEP = 0.01;

// ================= 光谱仿真（公式 5） =================
VectorXd simulateInterferenceIntensity(
    double A, double B, double z,
    const VectorXd& lambda)
{
    VectorXd I(lambda.size());

    for (Eigen::Index i = 0; i < lambda.size(); ++i) {
        double lambda_um = lambda(i) * 1e-3;
        double k = 2.0 * M_PI / lambda_um;

        double gauss =
            std::exp(-std::pow(lambda(i) - MU, 2)
                / (2.0 * ETA * ETA))
            / (ETA * std::sqrt(2.0 * M_PI));

        I(i) = (A + B * std::cos(k * z + PHI0)) * C * gauss;
    }
    return I;
}

// ================= EMD（简化 IMF1） =================
VectorXd emdFilter(const VectorXd& signal)
{
    VectorXd x = signal;
    VectorXd imf1(signal.size());

    for (int iter = 0; iter < EMD_MAX_ITER; ++iter) {
        std::vector<int> max_idx, min_idx;

        for (int i = 1; i < x.size() - 1; ++i) {
            if (x(i) > x(i - 1) && x(i) > x(i + 1)) max_idx.push_back(i);
            if (x(i) < x(i - 1) && x(i) < x(i + 1)) min_idx.push_back(i);
        }
        if (max_idx.size() < 2 || min_idx.size() < 2) break;

        VectorXd upper = x, lower = x;
        for (size_t i = 0; i + 1 < max_idx.size(); ++i)
            upper.segment(max_idx[i], max_idx[i + 1] - max_idx[i])
            .setLinSpaced(
                max_idx[i + 1] - max_idx[i],
                x(max_idx[i]),
                x(max_idx[i + 1]));

        for (size_t i = 0; i + 1 < min_idx.size(); ++i)
            lower.segment(min_idx[i], min_idx[i + 1] - min_idx[i])
            .setLinSpaced(
                min_idx[i + 1] - min_idx[i],
                x(min_idx[i]),
                x(min_idx[i + 1]));

        VectorXd m = 0.5 * (upper + lower);
        imf1 = x - m;
        x = m;
    }
    return imf1;
}

// ================= λ → σ =================
void lambdaToSigma(
    const VectorXd& lambda,
    const VectorXd& imf1,
    VectorXd& sigma,
    VectorXd& imf1_sigma)
{
    sigma.resize(lambda.size());
    imf1_sigma.resize(lambda.size());

    for (Eigen::Index i = 0; i < lambda.size(); ++i) {
        double lambda_um = lambda(i) * 1e-3;
        sigma(i) = 1.0 / lambda_um;
        imf1_sigma(i) = imf1(i);
    }
}

// ================= Lomb–Scargle =================
double lombScargle(
    const VectorXd& sigma,
    const VectorXd& imf1_sigma)
{
    int z_num = int((Z_SEARCH_MAX - Z_SEARCH_MIN) / Z_SEARCH_STEP) + 1;
    double best_z = 0.0;
    double best_P = -1.0;

    for (int i = 0; i < z_num; ++i) {
        double z = Z_SEARCH_MIN + i * Z_SEARCH_STEP;
        double sum_c = 0, sum_s = 0, sum_cc = 0, sum_ss = 0;

        for (int k = 0; k < sigma.size(); ++k) {
            double w = 2.0 * M_PI * z * sigma(k);
            double c = std::cos(w);
            double s = std::sin(w);

            sum_c += imf1_sigma(k) * c;
            sum_s += imf1_sigma(k) * s;
            sum_cc += c * c;
            sum_ss += s * s;
        }

        double P = (sum_c * sum_c) / sum_cc
            + (sum_s * sum_s) / sum_ss;

        if (P > best_P) {
            best_P = P;
            best_z = z;
        }
    }
    return best_z;
}

// ================= 厚度计算（公式 4） =================
double calculateThickness(double z_peak, double refractive)
{
    // z = 2 n d  →  d = z / (2 n)
    return z_peak / (2.0 * refractive);
}
