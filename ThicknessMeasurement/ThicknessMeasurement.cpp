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
const double DEFAULT_REFRACTIVE = 3.7; // 折射率

// 搜索范围设置
// 注意：Z代表光程差(OPD)，OPD ≈ 2 * n * d
// 如果测厚范围是 10um ~ 200um，n=1.5，则 OPD 范围应为 30um ~ 600um
const double Z_SEARCH_MIN = 100.0;
const double Z_SEARCH_MAX = 6000.0;    // 扩大到 6000 um (也就是对应最大厚度 2000 um)
const double Z_SEARCH_STEP = 0.5;      // 搜索范围变大后，步长不能太小，否则计算太慢。

// ================= 光谱仿真 =================
// 输入：lambda (nm), z (um, 注意这里 z 代表光程差 OPD)
VectorXd simulateInterferenceIntensity(
    double A, double B, double z_opd_um,
    const VectorXd& lambda_nm)
{
    VectorXd I(lambda_nm.size());

    for (Eigen::Index i = 0; i < lambda_nm.size(); ++i) {
        // 关键：将 nm 转换为 um，以便与 z_opd_um 单位一致
        double lambda_um = lambda_nm(i) * 1e-3;

        // k = 2 * PI / lambda
        double k = 2.0 * M_PI / lambda_um;

        // 高斯包络 (注意这里必须用 nm 计算，因为 MU 和 ETA 是 nm)
        double delta_lambda = lambda_nm(i) - MU;
        double gauss = std::exp(-(delta_lambda * delta_lambda)
            / (2.0 * ETA * ETA))
            / (ETA * std::sqrt(2.0 * M_PI));

        // 信号生成: cos(k * OPD)
        I(i) = (A + B * std::cos(k * z_opd_um + PHI0)) * C * gauss;
    }
    return I;
}

// ================= EMD滤波 (提取高频干涉项) =================
VectorXd emdFilter(const VectorXd& signal)
{
    VectorXd x = signal;
    VectorXd imf1(signal.size());
    const int EMD_MAX_ITER = 20; // 通常10-20次迭代足矣，100次太慢

    for (int iter = 0; iter < EMD_MAX_ITER; ++iter) {
        std::vector<int> max_idx, min_idx;

        // 寻找极值点
        for (int i = 1; i < x.size() - 1; ++i) {
            if (x(i) >= x(i - 1) && x(i) > x(i + 1)) max_idx.push_back(i);
            else if (x(i) <= x(i - 1) && x(i) < x(i + 1)) min_idx.push_back(i);
        }

        // 极值点过少，无法继续分解
        if (max_idx.size() < 2 || min_idx.size() < 2) {
            imf1 = x; // 返回当前残差
            break;
        }

        VectorXd upper = x;
        VectorXd lower = x;

        // 简单的线性插值包络
        for (size_t i = 0; i + 1 < max_idx.size(); ++i) {
            int len = max_idx[i + 1] - max_idx[i];
            upper.segment(max_idx[i], len + 1)
                .setLinSpaced(len + 1, x(max_idx[i]), x(max_idx[i + 1]));
        }
        for (size_t i = 0; i + 1 < min_idx.size(); ++i) {
            int len = min_idx[i + 1] - min_idx[i];
            lower.segment(min_idx[i], len + 1)
                .setLinSpaced(len + 1, x(min_idx[i]), x(min_idx[i + 1]));
        }

        VectorXd m = 0.5 * (upper + lower);
        imf1 = x - m;
        x = imf1; // 下一次迭代针对残差
    }
    return imf1;
}

// ================= 坐标变换: λ(nm) -> σ(um^-1) =================
void lambdaToSigma(
    const VectorXd& lambda_nm,
    const VectorXd& imf1,
    VectorXd& sigma,
    VectorXd& imf1_sigma)
{
    sigma.resize(lambda_nm.size());
    imf1_sigma.resize(lambda_nm.size());

    for (Eigen::Index i = 0; i < lambda_nm.size(); ++i) {
        // 关键：将 nm 转为 um
        double lambda_um = lambda_nm(i) * 1e-3;

        // 波数 σ = 1 / λ (单位: um^-1)
        if (lambda_um > 1e-6) {
            sigma(i) = 1.0 / lambda_um;
        }
        else {
            sigma(i) = 0.0;
        }
        imf1_sigma(i) = imf1(i);
    }
}

// ================= Lomb–Scargle Periodogram =================
// 返回峰值对应的 z (即 OPD，单位 um)
double lombScargle(
    const VectorXd& sigma,
    const VectorXd& imf1_sigma)
{
    // 计算搜索步数
    int z_num = int((Z_SEARCH_MAX - Z_SEARCH_MIN) / Z_SEARCH_STEP);

    double best_z = 0.0;
    double best_P = -1.0;

    // 预计算数据的方差（可选，为了标准归一化，此处简化略去）

    // 遍历搜索 Z (OPD)
    for (int i = 0; i <= z_num; ++i) {
        double z = Z_SEARCH_MIN + i * Z_SEARCH_STEP;

        double sum_c = 0, sum_s = 0;
        double sum_cc = 0, sum_ss = 0, sum_cs = 0;

        // 针对每一个波数点求和
        for (int k = 0; k < sigma.size(); ++k) {
            // 相位 = 2 * PI * OPD * sigma
            // sigma 是 1/lambda，单位 um^-1
            // z 是 OPD，单位 um
            // 结果是弧度，量级正确
            double w = 2.0 * M_PI * z * sigma(k);

            double c = std::cos(w);
            double s = std::sin(w);
            double val = imf1_sigma(k);

            sum_c += val * c;
            sum_s += val * s;
            sum_cc += c * c;
            sum_ss += s * s;
        }

        // 简化的 Lomb-Scargle 功率谱公式 (假设均值为0)
        // 避免除以零
        if (sum_cc < 1e-9 || sum_ss < 1e-9) continue;

        double P = (sum_c * sum_c) / sum_cc + (sum_s * sum_s) / sum_ss;

        if (P > best_P) {
            best_P = P;
            best_z = z;
        }
    }

    return best_z; // 返回 OPD (um)
}

// ================= 厚度计算 =================
double calculateThickness(double z_peak_opd, double refractive)
{
    // 物理公式: OPD = 2 * n * d
    // 所以: d = OPD / (2 * n)
    // z_peak_opd 单位是 um，结果 d 单位也是 um
    if (refractive <= 0.0) return 0.0;
    return z_peak_opd / (2.0 * refractive);
}