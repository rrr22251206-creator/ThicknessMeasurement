#define _USE_MATH_DEFINES
#include "ThicknessMeasurement.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace Eigen;

// ================= 常量定义 =================
const double MU = 1315.0;     // 中心波长 (nm)
const double ETA = 25.0;      // 带宽 (nm)
const double C = 18000.0;     // 强度系数
const double PHI0 = M_PI_2;   // 相位偏移

// 【关键修改】将默认折射率设为 1.0，使算法输出光程差(OPD)
const double DEFAULT_REFRACTIVE = 1.0;

// 【关键修改】搜索范围设置
const double Z_SEARCH_MIN = 100.0;
const double Z_SEARCH_MAX = 1000.0;
const double Z_SEARCH_STEP = 0.1; // 建议设为 0.1 以提高精度

// ================= 光谱仿真 =================
VectorXd simulateInterferenceIntensity(
    double A, double B, double z_opd_um,
    const VectorXd& lambda_nm)
{
    VectorXd I(lambda_nm.size());
    for (Eigen::Index i = 0; i < lambda_nm.size(); ++i) {
        double lambda_um = lambda_nm(i) * 1e-3;
        double k = 2.0 * M_PI / lambda_um;
        double delta_lambda = lambda_nm(i) - MU;
        double gauss = std::exp(-(delta_lambda * delta_lambda) / (2.0 * ETA * ETA)) / (ETA * std::sqrt(2.0 * M_PI));
        I(i) = (A + B * std::cos(k * z_opd_um + PHI0)) * C * gauss;
    }
    return I;
}

// ================= EMD滤波 =================
VectorXd emdFilter(const VectorXd& signal)
{
    VectorXd x = signal;
    VectorXd imf1(signal.size());
    const int EMD_MAX_ITER = 20;

    for (int iter = 0; iter < EMD_MAX_ITER; ++iter) {
        std::vector<int> max_idx, min_idx;
        for (int i = 1; i < x.size() - 1; ++i) {
            if (x(i) >= x(i - 1) && x(i) > x(i + 1)) max_idx.push_back(i);
            else if (x(i) <= x(i - 1) && x(i) < x(i + 1)) min_idx.push_back(i);
        }

        if (max_idx.size() < 2 || min_idx.size() < 2) {
            imf1 = x;
            break;
        }

        VectorXd upper = x;
        VectorXd lower = x;

        for (size_t i = 0; i + 1 < max_idx.size(); ++i) {
            int len = max_idx[i + 1] - max_idx[i];
            upper.segment(max_idx[i], len + 1).setLinSpaced(len + 1, x(max_idx[i]), x(max_idx[i + 1]));
        }
        for (size_t i = 0; i + 1 < min_idx.size(); ++i) {
            int len = min_idx[i + 1] - min_idx[i];
            lower.segment(min_idx[i], len + 1).setLinSpaced(len + 1, x(min_idx[i]), x(min_idx[i + 1]));
        }

        VectorXd m = 0.5 * (upper + lower);
        imf1 = x - m;
        x = imf1;
    }
    return imf1;
}

// ================= 坐标变换 =================
void lambdaToSigma(const VectorXd& lambda_nm, const VectorXd& imf1, VectorXd& sigma, VectorXd& imf1_sigma)
{
    sigma.resize(lambda_nm.size());
    imf1_sigma.resize(lambda_nm.size());
    for (Eigen::Index i = 0; i < lambda_nm.size(); ++i) {
        double lambda_um = lambda_nm(i) * 1e-3;
        if (lambda_um > 1e-6) sigma(i) = 1.0 / lambda_um;
        else sigma(i) = 0.0;
        imf1_sigma(i) = imf1(i);
    }
}

// ================= 计算功率谱 =================
SpectrumData calculatePowerSpectrum(const VectorXd& sigma, const VectorXd& imf1_sigma)
{
    int z_num = int((Z_SEARCH_MAX - Z_SEARCH_MIN) / Z_SEARCH_STEP) + 1;

    SpectrumData result;
    result.z_axis.resize(z_num);
    result.power.resize(z_num);

    for (int i = 0; i < z_num; ++i) {
        double z = Z_SEARCH_MIN + i * Z_SEARCH_STEP;
        double sum_c = 0, sum_s = 0, sum_cc = 0, sum_ss = 0;

        for (int k = 0; k < sigma.size(); ++k) {
            double w = 2.0 * M_PI * z * sigma(k);
            double c = std::cos(w);
            double s = std::sin(w);
            double val = imf1_sigma(k);

            sum_c += val * c;
            sum_s += val * s;
            sum_cc += c * c;
            sum_ss += s * s;
        }

        double P = 0.0;
        if (sum_cc > 1e-9 && sum_ss > 1e-9) {
            P = (sum_c * sum_c) / sum_cc + (sum_s * sum_s) / sum_ss;
        }

        result.z_axis(i) = z;
        result.power(i) = P;
    }
    return result;
}

// ================= 【核心修改】寻找多峰 (带距离抑制) =================
std::vector<PeakInfo> findPeaks(const SpectrumData& spectrum, double thresholdRatio)
{
    std::vector<PeakInfo> raw_peaks;

    // 1. 找到全局最大值
    double max_p = spectrum.power.maxCoeff();
    double threshold = max_p * thresholdRatio;

    // 2. 找出所有局部最大值
    for (int i = 1; i < spectrum.power.size() - 1; ++i) {
        double current = spectrum.power(i);
        double prev = spectrum.power(i - 1);
        double next = spectrum.power(i + 1);

        if (current > prev && current > next && current > threshold) {
            raw_peaks.push_back({ spectrum.z_axis(i), current });
        }
    }

    // 3. 按强度降序排序
    std::sort(raw_peaks.begin(), raw_peaks.end(), [](const PeakInfo& a, const PeakInfo& b) {
        return a.power > b.power;
        });

    // 4. 距离抑制 (过滤旁瓣)
    std::vector<PeakInfo> final_peaks;
    const double MIN_DIST = 50.0; // 设置为50um以去除旁瓣

    for (const auto& candidate : raw_peaks) {
        bool keep = true;
        for (const auto& existing : final_peaks) {
            if (std::abs(candidate.z_opd - existing.z_opd) < MIN_DIST) {
                keep = false;
                break;
            }
        }
        if (keep) {
            final_peaks.push_back(candidate);
        }
    }

    // 5. 按 OPD 从小到大排序返回
    std::sort(final_peaks.begin(), final_peaks.end(), [](const PeakInfo& a, const PeakInfo& b) {
        return a.z_opd < b.z_opd;
        });

    return final_peaks;
}

// ================= 厚度换算 =================
double opdToThickness(double opd, double refractive)
{
    if (refractive <= 0.0) return 0.0;
    return opd / (2.0 * refractive);
}