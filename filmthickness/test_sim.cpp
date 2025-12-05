#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "FilmThickness.hpp"

int main() {
    const int N = 2048;
    std::vector<double> lambda(N), I(N);

    double A = 0.5, B = 0.2, C = 18000;
    double mu = 660e-9, eta = 80e-9;
    double phi0 = M_PI / 2;

    double true_d = 9e-6;   // 9 ¦Ìm
    double n1 = 1.5;
    double z = true_d * n1;

    for (int i = 0; i < N; i++) {
        lambda[i] = 500e-9 + (800e-9 - 500e-9) * i / (N - 1);
        I[i] = SpectrumModel::interference(lambda[i],
            A, B, C,
            mu, eta,
            phi0, z);
    }

    FilmThicknessSolver solver(n1);
    double d_est = solver.computeThickness(lambda, I);

    std::cout << "True thickness: " << true_d << "\n";
    std::cout << "Estimated thickness: " << d_est << "\n";
    std::cout << "Error: " << (d_est - true_d) * 1e9 << " nm\n";

    return 0;
}
