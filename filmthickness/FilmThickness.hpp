#pragma once
#include <vector>

/**
 * @brief IMF1 提取模式（保留接口用于扩展）
 */
enum class EMDMode {
    FastHighPass    ///< 使用零相位平滑去基线代替 EMD IMF1
};

/**
 * @class SpectrumModel
 * @brief 光谱仿真模块，对应论文公式 (5)
 *
 * 公式 (5):
 *  I(λ) = (A + B cos(4π n1 d / λ + φ0)) * C * exp(-(λ - μ)^2 / (2η²))
 */
class SpectrumModel {
public:
    static double gaussian(double lambda, double mu, double sigma);

    static double interference(
        double lambda,
        double A, double B, double C,
        double mu, double eta,
        double phi0,
        double n1,
        double d
    );
};

/**
 * @class EMD
 * @brief IMF1 提取，对应论文 Fig.3
 *
 * 此工程版采用：
 *   IMF1 ≈ 原信号 − 零相位平滑基线
 * 完全零相位，不改变主频位置，效果等价于 EMD 的 IMF1。
 */
class EMD {
public:
    static std::vector<double> extractIMF1(
        const std::vector<double>& x,
        EMDMode mode = EMDMode::FastHighPass);

private:
    static std::vector<double> highpass(const std::vector<double>& x);
};

/**
 * @class LombScargle
 * @brief Lomb–Scargle 功率谱，对应论文公式 (2)(3)
 */
class LombScargle {
public:
    static double findPeak(
        const std::vector<double>& sigma,
        const std::vector<double>& I,
        double z_min, double z_max, double step);
};

/**
 * @class FilmThicknessSolver
 * @brief PPS 主流程，对应论文公式 (4)
 *
 * 厚度公式：
 *      d = z_peak / (2 n1)
 */
class FilmThicknessSolver {
public:
    FilmThicknessSolver(double refrIndex);

    double computeThickness(
        const std::vector<double>& lambda,
        const std::vector<double>& intensity,
        EMDMode mode = EMDMode::FastHighPass);

private:
    double n1;
};
