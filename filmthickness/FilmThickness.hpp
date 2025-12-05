#pragma once
#include <vector>
#include <cmath>

//
// ======= SpectrumModel（含 PDF Eq.(5)）=======
//
class SpectrumModel {
public:
    static double gaussian(double lambda, double mu, double sigma);

    // PDF Eq.(5)：I(λ) = (A + B cos(k z + φ0)) * C * G(λ)
    static double interference(double lambda,
        double A, double B, double C,
        double mu, double eta,
        double phi0, double z);
};

//
// ======= EMD（使用高通滤波替代 IMF1）=======
//
enum class EMDMode {
    FastHighPass
};

class EMD {
public:
    static std::vector<double> extractIMF1(
        const std::vector<double>& x,
        EMDMode mode = EMDMode::FastHighPass
    );

private:
    static std::vector<double> highpass(const std::vector<double>& x);
};

//
// ======= Lomb–Scargle（含 PDF Eq.(2), Eq.(3)）=======
//
class LombScargle {
public:
    static double findPeak(
        const std::vector<double>& sigma,
        const std::vector<double>& I,
        double z_min, double z_max, double step);
};

//
// ======= PPS 主算法（含 PDF Eq.(4)）======
//
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
