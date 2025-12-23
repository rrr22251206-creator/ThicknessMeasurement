#pragma once
#include <Eigen/Dense>
#include <vector>

using Eigen::VectorXd;

// ================= Data Structures =================

struct SpectrumData {
    VectorXd z_axis;
    VectorXd power;
};

struct PeakInfo {
    double z_opd;
    double power;
};

// Configuration struct for runtime tuning
struct AnalysisConfig {
    double z_search_min = 50.0;    // Min OPD (um)
    double z_search_max = 1600.0;  // Max OPD (um) - Extended to cover total sum
    double z_search_step = 0.2;    // Step size (um)
    double peak_threshold = 0.05;  // Relative threshold (0.0~1.0). Low value to catch weak layers.
    double min_peak_dist = 40.0;   // Min distance between peaks (um) to suppress side lobes
    double additivity_tol = 30.0;  // Tolerance for A+B=C check (um)
};

// Result structure
struct LayerStructure {
    bool is_double_layer;
    double opd_layer_1; // Signal A
    double opd_layer_2; // Signal B
    double opd_total;   // Signal C (Sum)
    double error;       // |(A+B) - C|
    double score;       // Confidence score
};

// ================= Constants =================
extern const double MU;
extern const double ETA;
extern const double C;
extern const double PHI0;

// ================= Algorithms =================

// Preprocessing: Remove DC bias
VectorXd emdFilter(const VectorXd& signal);

// Domain conversion: Wavelength -> Wavenumber
void lambdaToSigma(
    const VectorXd& lambda,
    const VectorXd& imf1,
    VectorXd& sigma,
    VectorXd& imf1_sigma);

// Power Spectrum (Lomb-Scargle / FFT)
SpectrumData calculatePowerSpectrum(
    const VectorXd& sigma,
    const VectorXd& imf1_sigma,
    const AnalysisConfig& config);

// Peak Finding with suppression
std::vector<PeakInfo> findPeaks(
    const SpectrumData& spectrum,
    const AnalysisConfig& config);

// Blind Topology Solver (Auto-detect single/double layers)
LayerStructure solveLayerTopology(
    const std::vector<PeakInfo>& peaks,
    const AnalysisConfig& config);