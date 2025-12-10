#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include "FilmThickness.hpp"

/**
 * @brief 主程序：使用论文公式 (5) 仿真光谱并使用 PPS 求厚度
 */
int main()
{
    const int N = 2048;
    std::vector<double> lambda(N);
    std::vector<double> I(N);

    // 光源参数（论文公式 5）
    double A = 0.5;
    double B = 0.2;
    double C = 18000;
    double mu = 660e-9;
    double eta = 80e-9;
    double phi0 = M_PI / 2;

    // 真实薄膜
    double true_d = 9e-6;  // 9 µm
    double n1 = 1.5;

    // 生成反射干涉光谱
    for (int i = 0; i < N; i++) {
        lambda[i] = 500e-9 + (800e-9 - 500e-9) * i / (N - 1);
        I[i] = SpectrumModel::interference(
            lambda[i],
            A, B, C,
            mu, eta,
            phi0,
            n1,
            true_d
        );
    }

    // PPS 求解
    FilmThicknessSolver solver(n1);
    double d_est = solver.computeThickness(lambda, I);

    // 输出
    std::cout << "True thickness:      " << true_d << " m\n";
    std::cout << "Estimated thickness: " << d_est << " m\n";
    std::cout << "Error:               "
        << (d_est - true_d) * 1e9 << " nm\n";

    return 0;
}
