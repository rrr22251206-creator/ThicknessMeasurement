#include <iostream>
#include <cmath>    // For math operations
#include <iomanip>  // For output formatting

using namespace std;

/*
 * ======================================================================================
 * 核心公式推导 (Derivation of the Constraint Formula)
 * ======================================================================================
 *
 * 1. 干涉原理 (Interference Principle):
 *    在红外LSP（低相干光谱干涉）测厚中，薄膜上下表面的反射光发生干涉。
 *    干涉极值条件满足光程差公式：
 *    OPD = 2 * n * d = m * λ
 *    其中：
 *      n: 折射率 (Refractive index)
 *      d: 物理厚度 (Physical thickness)
 *      m: 干涉级次 (Interference order, 整数)
 *      λ: 波长 (Wavelength)
 *
 * 2. 光谱周期性 (Spectral Periodicity):
 *    在波长域中，相邻两个干涉条纹（级次 m 和 m+1）之间的波长间隔 Δλ (自由光谱范围) 近似为：
 *    Δλ_period ≈ λ^2 / (2 * n * d)
 *    这意味着：厚度 d 越大，光谱上的条纹越密集（即 Δλ 越小）。
 *
 * 3. 奈奎斯特采样定理 (Nyquist Sampling Theorem):
 *    为了能够分辨出正弦形式的干涉条纹，光谱仪的采样频率必须至少是信号频率的2倍。
 *    在光谱仪中，采样间隔就是光学分辨率 δλ (Resolution)。
 *    根据奈奎斯特判据，分辨率必须小于等于条纹周期的一半：
 *    δλ ≤ Δλ_period / 2
 *
 * 4. 最终公式推导 (Final Calculation):
 *    将第2步的 Δλ_period 代入第3步：
 *    δλ ≤ (λ^2 / (2 * n * d)) / 2
 *    δλ ≤ λ^2 / (4 * n * d)
 *
 *    反转公式以求最大厚度 d_max：
 *    4 * n * d * δλ ≤ λ^2
 *
 *    ==> d_max = λ^2 / (4 * n * δλ)
 *
 * ======================================================================================
 */

 /**
  * @brief Calculate the theoretical maximum measurable thickness based on spectrometer parameters.
  * @param center_wavelength_nm Center Wavelength lambda (nm)
  * @param resolution_nm Spectrometer Resolution delta_lambda (nm)
  * @param refractive_index Refractive Index n (default 1.0)
  * @return double Theoretical Max Thickness (um)
  */
double calculate_max_thickness(double center_wavelength_nm, double resolution_nm, double refractive_index = 1.0) {

    // --- Step A: Unit Conversion (nm -> m) ---
    // Use standard SI units (meters) to avoid order of magnitude errors.
    double lambda_m = center_wavelength_nm * 1e-9;
    double res_m = resolution_nm * 1e-9;

    // --- Step B: Apply Formula ---
    // d_max = lambda^2 / (4 * n * delta_lambda)
    double numerator = lambda_m * lambda_m;
    double denominator = 4.0 * refractive_index * res_m;

    double d_max_m = numerator / denominator;

    // --- Step C: Result Conversion (m -> um) ---
    double d_max_um = d_max_m * 1e6;

    return d_max_um;
}

int main() {
    // English Output Mode
    cout << "=================================================" << endl;
    cout << "   IR LSP Thickness: Resolution vs Limit Verify" << endl;
    cout << "=================================================" << endl << endl;

    // 1. Define Input Parameters (Based on PDF Page 13)
    double lambda_0 = 1150.0;   // Center Wavelength (nm)
    double res_req = 0.066;     // Required Resolution derived in PDF (nm)
    double n_air = 1.0;         // Refractive Index of Air

    // 2. Perform Calculation
    double limit_thickness = calculate_max_thickness(lambda_0, res_req, n_air);

    // 3. Print Detailed Info
    cout << fixed << setprecision(4);
    cout << "[Parameters Setup]" << endl;
    cout << "  Center Wavelength (lambda)      : " << lambda_0 << " nm" << endl;
    cout << "  Spectral Resolution (delta_lambda): " << res_req << " nm" << endl;
    cout << "  Refractive Index (n)            : " << n_air << endl;
    cout << "-------------------------------------------------" << endl;
    cout << "[Calculation: d_max = lambda^2 / (4 * n * delta_lambda)]" << endl;
    cout << "  Theoretical Max Thickness       : " << limit_thickness << " um" << endl;
    cout << "-------------------------------------------------" << endl;

    // 4. Verification Logic
    double target_thickness = 5000.0;
    cout << "[Verification Conclusion]" << endl;
    if (limit_thickness >= target_thickness) {
        cout << "  [PASS] Current resolution (" << res_req << " nm) meets the "
            << (int)target_thickness << " um thickness requirement." << endl;
        cout << "  (Theoretical Margin: " << (limit_thickness - target_thickness) << " um)" << endl;
    }
    else {
        cout << "  [FAIL] Resolution insufficient. Cannot measure " << (int)target_thickness << " um." << endl;
    }

    cout << endl;
    system("pause"); // Keep console open
    return 0;
}