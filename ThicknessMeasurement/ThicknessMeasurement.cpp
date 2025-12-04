#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

// 全局参数配置（与原文仿真/实验参数一致，可根据实际需求调整）
const double MU = 660.0;         // 高斯光源中心波长(nm)
const double ETA = 80.0;         // 光源波长标准差(nm)
const double C = 18000.0;        // 光谱最大强度振幅
const double PHI0 = M_PI_2;      // 初始相位(rad)（使用M_PI_2替代M_PI/2）
const double DEFAULT_REFRACTIVE = 1.5; // 默认薄膜折射率（如PVC涂层）
const double LAMBDA_MIN = 500.0; // 光谱仪最小波长(nm)
const double LAMBDA_MAX = 1000.0;// 光谱仪最大波长(nm)
const int PIXEL_NUM = 512;       // 光谱仪像素数（波长采样点数量）
const int EMD_MAX_ITER = 100;    // EMD筛分最大迭代次数
const double Z_SEARCH_MIN = 0.0; // LSP光程差搜索最小值(μm)
const double Z_SEARCH_MAX = 100.0;// LSP光程差搜索最大值(μm)
const double Z_SEARCH_STEP = 0.01;// LSP光程差搜索步长(μm)

/**
 * @brief 仿真反射干涉光强信号（公式5）
 * @param A 反射系数参数A (0<A≤1)
 * @param B 反射系数参数B (0<B≤1)
 * @param z 光程差(μm)
 * @param lambda 波长数组(nm)
 * @return 光强数组I_r
 */
VectorXd simulateInterferenceIntensity(double A, double B, double z, const VectorXd& lambda) {
    Eigen::Index n = lambda.size();
    VectorXd I(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        double lambda_μm = lambda(i) * 1e-3; // 波长转换为μm
        double k = 2 * M_PI / lambda_μm;     // 波数(1/μm)
        double cos_term = cos(k * z + PHI0); // 余弦干涉项

        // 高斯光源强度分布G(λ)
        double gauss = exp(-pow(lambda(i) - MU, 2) / (2 * pow(ETA, 2)))
            / (ETA * sqrt(2 * M_PI));

        // 反射干涉光强计算（公式5）
        I(i) = (A + B * cos_term) * C * gauss;
    }
    return I;
}

/**
 * @brief EMD滤波：提取第一阶IMF（高频干涉分量，公式相关处理）
 * @param signal 原始光强信号
 * @return IMF1：第一阶固有模态函数
 */
VectorXd emdFilter(const VectorXd& signal) {
    Eigen::Index n = signal.size();
    VectorXd x = signal;
    VectorXd imf1(n);

    for (int iter = 0; iter < EMD_MAX_ITER; ++iter) {
        // 1. 提取局部极大值和极小值索引
        vector<Eigen::Index> max_idx, min_idx;
        for (Eigen::Index i = 1; i < n - 1; ++i) {
            if (x(i) > x(i - 1) && x(i) > x(i + 1)) max_idx.push_back(i);
            if (x(i) < x(i - 1) && x(i) < x(i + 1)) min_idx.push_back(i);
        }
        if (max_idx.size() < 2 || min_idx.size() < 2) break; // 极值点不足，终止筛分

        // 2. 线性插值拟合上、下包络线（简化实现，实际可替换为三次样条）
        VectorXd upper_env(n), lower_env(n);
        // 上包络插值
        for (size_t i = 0; i < max_idx.size() - 1; ++i) {
            Eigen::Index left = max_idx[i];
            Eigen::Index right = max_idx[i + 1];
            double slope = (x(right) - x(left)) / (right - left);
            for (Eigen::Index j = left; j <= right; ++j) {
                upper_env(j) = x(left) + slope * (j - left);
            }
        }
        // 下包络插值
        for (size_t i = 0; i < min_idx.size() - 1; ++i) {
            Eigen::Index left = min_idx[i];
            Eigen::Index right = min_idx[i + 1];
            double slope = (x(right) - x(left)) / (right - left);
            for (Eigen::Index j = left; j <= right; ++j) {
                lower_env(j) = x(left) + slope * (j - left);
            }
        }
        // 边界填充（使用首尾极值点值）
        for (Eigen::Index j = 0; j < max_idx[0]; ++j) upper_env(j) = x(max_idx[0]);
        for (Eigen::Index j = max_idx.back() + 1; j < n; ++j) upper_env(j) = x(max_idx.back());
        for (Eigen::Index j = 0; j < min_idx[0]; ++j) lower_env(j) = x(min_idx[0]);
        for (Eigen::Index j = min_idx.back() + 1; j < n; ++j) lower_env(j) = x(min_idx.back());

        // 3. 计算包络均值并筛分
        VectorXd m = (upper_env + lower_env) / 2;
        VectorXd h = x - m;
        imf1 = h;
        x = m; // 剩余信号作为下一次筛分输入
    }
    return imf1;
}

/**
 * @brief 波长域到波数域的变换（λ→σ=1/λ，解决非等间隔采样问题）
 * @param lambda 波长数组(nm)
 * @param imf1 波长域IMF1信号
 * @param sigma 输出波数数组(1/μm)
 * @param imf1_sigma 输出波数域IMF1信号
 */
void lambdaToSigma(const VectorXd& lambda, const VectorXd& imf1,
    VectorXd& sigma, VectorXd& imf1_sigma) {
    Eigen::Index n = lambda.size();
    sigma.resize(n);
    imf1_sigma.resize(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        double lambda_μm = lambda(i) * 1e-3; // 波长转换为μm
        sigma(i) = 1.0 / lambda_μm;          // 波数σ=1/λ (1/μm)
        imf1_sigma(i) = imf1(i);             // 光强值对应到波数域
    }
}

/**
 * @brief Lomb-Scargle周期图分析（公式2、3），提取功率谱峰值对应的光程差z_peak
 * @param sigma 波数数组(1/μm)
 * @param imf1_sigma 波数域IMF1信号
 * @return z_peak：功率谱峰值对应的光程差(μm)
 */
double lombScargle(const VectorXd& sigma, const VectorXd& imf1_sigma) {
    Eigen::Index z_num = static_cast<Eigen::Index>((Z_SEARCH_MAX - Z_SEARCH_MIN) / Z_SEARCH_STEP + 1);
    VectorXd z_list(z_num);
    VectorXd P(z_num); // LSP功率谱

    // 遍历所有候选光程差z，计算功率谱
    for (Eigen::Index z_idx = 0; z_idx < z_num; ++z_idx) {
        double z = Z_SEARCH_MIN + z_idx * Z_SEARCH_STEP;
        z_list(z_idx) = z;

        // 计算时移不变参数δ（公式3）
        double sum_sin = 0.0, sum_cos = 0.0;
        for (Eigen::Index i = 0; i < sigma.size(); ++i) {
            double arg = 4 * M_PI * z * sigma(i);
            sum_sin += sin(arg);
            sum_cos += cos(arg);
        }
        double delta = atan2(sum_sin, sum_cos) / (4 * M_PI * z);

        // 计算LSP功率谱P(z)（公式2）
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

        // 避免分母为0（理论上不会发生，添加保护）
        sum_cos2 = max(sum_cos2, 1e-10);
        sum_sin2 = max(sum_sin2, 1e-10);
        double P_cos = pow(sum_cos_term, 2) / sum_cos2;
        double P_sin = pow(sum_sin_term, 2) / sum_sin2;
        P(z_idx) = 0.5 * (P_cos + P_sin);
    }

    // 提取功率谱峰值对应的z_peak
    Eigen::Index peak_idx;
    P.maxCoeff(&peak_idx);
    return z_list(peak_idx);
}

/**
 * @brief 计算薄膜厚度（公式4：d = z_peak / n）
 * @param z_peak LSP峰值对应的光程差(μm)
 * @param refractive 薄膜折射率（默认1.5，可根据材料修改）
 * @return d：薄膜厚度(μm)
 */
double calculateThickness(double z_peak, double refractive = DEFAULT_REFRACTIVE) {
    return z_peak / refractive;
}

/**
 * @brief 主函数：串联完整流程，测试薄膜厚度测量
 */
int main() {
    // 1. 生成波长数组（光谱仪像素对应的波长，500-1000nm均匀采样）
    VectorXd lambda(PIXEL_NUM);
    for (Eigen::Index i = 0; i < PIXEL_NUM; ++i) {
        lambda(i) = LAMBDA_MIN + (LAMBDA_MAX - LAMBDA_MIN) * i / (PIXEL_NUM - 1);
    }

    // 2. 仿真输入：假设真实厚度d_true=9μm（原文仿真参数）
    double d_true = 9.0;
    double refractive = 1.5;
    double z_true = 2 * d_true * refractive; // 往返光程差
    VectorXd I_raw = simulateInterferenceIntensity(0.5, 0.3, z_true, lambda); // A=0.5, B=0.3

    // 3. 波数域变换（λ→σ） - 跳过EMD滤波，直接对原始信号进行分析
    VectorXd sigma, imf1_sigma;
    lambdaToSigma(lambda, I_raw, sigma, imf1_sigma);

    // 4. LSP分析：提取光程差z_peak
    double z_peak = lombScargle(sigma, imf1_sigma);

    // 5. 计算薄膜厚度
    double d_calc = calculateThickness(z_peak, refractive);

    // 7. 输出测量结果
    cout << "========================================" << endl;
    cout << "          薄膜厚度测量结果              " << endl;
    cout << "========================================" << endl;
    cout << "真实厚度:    " << fixed << setprecision(3) << d_true << " μm" << endl;
    cout << "计算厚度:    " << fixed << setprecision(3) << d_calc << " μm" << endl;
    cout << "测量误差:    " << fixed << setprecision(3) << abs(d_calc - d_true) << " μm" << endl;
    cout << "是否满足精度: " << (abs(d_calc - d_true) <= 0.1 ? "是（≤0.1μm）" : "否") << endl;
    cout << "========================================" << endl;

    return 0;
}