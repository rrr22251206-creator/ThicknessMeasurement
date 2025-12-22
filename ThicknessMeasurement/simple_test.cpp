#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include "ThicknessMeasurement.h" 
#include "SpectrumLoader.h"

using namespace std;
using namespace Eigen;

int main()
{
    // ===== 1. 读取光谱文件 =====
    cout << "[Step 1] Loading Spectrum..." << endl;
    vector<double> lambda_nm_vec;
    vector<double> intensity_vec;

    if (!SpectrumLoader::load("spectrum.txt", lambda_nm_vec, intensity_vec))
    {
        return -1;
    }
    cout << "[INFO] Loaded spectrum points: " << lambda_nm_vec.size() << endl;

    // ===== 2. 转 Eigen =====
    int N = static_cast<int>(lambda_nm_vec.size());
    VectorXd lambda(N), I(N);
    for (int i = 0; i < N; ++i) {
        lambda(i) = lambda_nm_vec[i];
        I(i) = intensity_vec[i];
    }

    // ===== 3. EMD 滤波 =====
    cout << "[Step 2] Applying EMD Filter..." << endl;
    VectorXd imf1 = emdFilter(I);

    // ===== 4. 坐标映射 =====
    cout << "[Step 3] Mapping Lambda to Sigma..." << endl;
    VectorXd sigma, imf1_sigma;
    lambdaToSigma(lambda, imf1, sigma, imf1_sigma);

    // ===== 5. 计算功率谱 =====
    cout << "[Step 4] Calculating Power Spectrum..." << endl;
    SpectrumData spectrum = calculatePowerSpectrum(sigma, imf1_sigma);

    // ===== 6. 寻找峰值 =====
    cout << "[Step 5] Finding Peaks..." << endl;
    vector<PeakInfo> peaks = findPeaks(spectrum, 0.15);

    // ===== 7. 智能结果分析 =====
    cout << endl;
    cout << string(60, '=') << endl;
    cout << " FINAL THICKNESS REPORT " << endl;
    cout << string(60, '=') << endl;

    // --- 设定物理模型目标 ---
    double target_opd_low = 330.0;  // Layer 2 target
    double target_opd_high = 828.8; // Layer 1 target

    // 【关键修改】缩小搜索窗口，排除干扰较强的旁瓣 (例如 292.6)
    // 之前是 50.0，现在改为 30.0
    double search_window = 30.0;

    double best_opd_low = 0.0;
    double max_p_low = 0.0;

    double best_opd_high = 0.0;
    double max_p_high = 0.0;

    cout << "Raw Peaks Detected:" << endl;
    for (const auto& peak : peaks) {
        cout << " - Peak at " << setw(8) << peak.z_opd << " um (Power: " << setw(10) << peak.power << ")";

        // 匹配低折射率层
        if (abs(peak.z_opd - target_opd_low) < search_window) {
            cout << " [Candidate Layer 2]";
            if (peak.power > max_p_low) {
                max_p_low = peak.power;
                best_opd_low = peak.z_opd;
            }
        }
        // 匹配高折射率层
        else if (abs(peak.z_opd - target_opd_high) < search_window) {
            cout << " [Candidate Layer 1]";
            if (peak.power > max_p_high) {
                max_p_high = peak.power;
                best_opd_high = peak.z_opd;
            }
        }
        cout << endl;
    }
    cout << string(60, '-') << endl;

    // --- 输出 Layer 2 ---
    if (best_opd_low > 0) {
        double d_calc = best_opd_low / (2.0 * 1.5);
        cout << "[Layer 2 (Polymer/Glass)]" << endl;
        cout << "  > Refractive Index: 1.5" << endl;
        cout << "  > Measured OPD:     " << best_opd_low << " um" << endl;
        cout << "  > CALC THICKNESS:   " << d_calc << " um" << endl;
        cout << "  > Deviation:        " << (d_calc - 110.0) << " um" << endl;
    }
    else {
        cout << "[Layer 2] Not found! Check threshold or window." << endl;
    }

    cout << endl;

    // --- 输出 Layer 1 ---
    if (best_opd_high > 0) {
        double d_calc = best_opd_high / (2.0 * 3.7);
        cout << "[Layer 1 (Silicon)]" << endl;
        cout << "  > Refractive Index: 3.7" << endl;
        cout << "  > Measured OPD:     " << best_opd_high << " um" << endl;
        cout << "  > CALC THICKNESS:   " << d_calc << " um" << endl;
        cout << "  > Deviation:        " << (d_calc - 112.0) << " um" << endl;
    }
    else {
        cout << "[Layer 1] Not found! Check threshold or window." << endl;
    }

    cout << string(60, '=') << endl;

    system("pause");
    return 0;
}