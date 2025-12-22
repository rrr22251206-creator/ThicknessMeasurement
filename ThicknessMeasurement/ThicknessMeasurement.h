#pragma once
#include <Eigen/Dense>
#include <vector>

using Eigen::VectorXd;

// ================= 结构体定义 =================
// 存储功率谱数据 (用于绘图或分析)
struct SpectrumData {
    VectorXd z_axis; // 横坐标：光程差 (um)
    VectorXd power;  // 纵坐标：信号强度
};

// 存储峰值信息
struct PeakInfo {
    double z_opd;    // 峰值对应的光程差 (um)
    double power;    // 峰值强度
};

// ================= 常量声明 =================
extern const double MU;
extern const double ETA;
extern const double C;
extern const double PHI0;
extern const double DEFAULT_REFRACTIVE;

extern const double LAMBDA_MIN;
extern const double LAMBDA_MAX;
extern const int    PIXEL_NUM;

extern const int    EMD_MAX_ITER;
extern const double Z_SEARCH_MIN;
extern const double Z_SEARCH_MAX;
extern const double Z_SEARCH_STEP;

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

// 返回完整的功率谱数据
SpectrumData calculatePowerSpectrum(
    const VectorXd& sigma,
    const VectorXd& imf1_sigma);

// 从谱图中寻找多个峰值
std::vector<PeakInfo> findPeaks(
    const SpectrumData& spectrum,
    double thresholdRatio = 0.3); // 阈值系数，低于最大峰 30% 的峰会被忽略

// 简单的物理厚度换算辅助函数
double opdToThickness(
    double opd,
    double refractive);