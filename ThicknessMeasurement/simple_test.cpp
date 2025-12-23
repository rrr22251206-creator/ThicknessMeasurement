#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm> 
#include "ThicknessMeasurement.h"
#include "SpectrumLoader.h"

using namespace std;
using namespace Eigen;

struct Material {
    string name;
    double refractive_index;
};

int main() {
    // ==========================================
    // 0. CONFIGURATION
    // ==========================================
    AnalysisConfig config;
    config.z_search_min = 50.0;
    config.z_search_max = 1600.0;   // Extended range to see the Total Sum peak
    config.z_search_step = 0.2;
    config.min_peak_dist = 40.0;

    // Crucial: Low threshold to catch the weak Layer 2 signal (Power ~43k vs Max ~300k)
    config.peak_threshold = 0.05;

    // Crucial: Tolerance for A+B=C check
    config.additivity_tol = 30.0;

    vector<Material> material_db = {
        { "Polymer", 1.5 },
        { "Silicon", 3.7 },
        { "Glass", 1.45 },
        { "Air", 1.0 }
    };

    // ==========================================
    // 1. Load Data
    // ==========================================
    cout << "[INFO] Loading spectrum..." << endl;
    vector<double> lambda_vec, intensity_vec;
    if (!SpectrumLoader::load("spectrum.txt", lambda_vec, intensity_vec)) return -1;

    int N = (int)lambda_vec.size();
    VectorXd lambda(N), I(N);
    for (int i = 0; i < N; ++i) { lambda(i) = lambda_vec[i]; I(i) = intensity_vec[i]; }

    // ==========================================
    // 2. Signal Processing
    // ==========================================
    VectorXd imf1 = emdFilter(I);
    VectorXd sigma, imf1_sigma;
    lambdaToSigma(lambda, imf1, sigma, imf1_sigma);

    cout << "[INFO] Running Spectral Analysis..." << endl;
    SpectrumData spectrum = calculatePowerSpectrum(sigma, imf1_sigma, config);
    vector<PeakInfo> peaks = findPeaks(spectrum, config);

    cout << "--- Detected Raw Peaks (OPD) ---" << endl;
    for (const auto& p : peaks) {
        cout << "OPD: " << setw(8) << p.z_opd << " um | Power: " << setw(10) << p.power << endl;
    }
    cout << "--------------------------------" << endl;

    // ==========================================
    // 3. Blind Topology Solve
    // ==========================================
    LayerStructure result = solveLayerTopology(peaks, config);

    cout << string(60, '=') << endl;
    cout << " ANALYSIS REPORT " << endl;
    cout << string(60, '=') << endl;

    if (result.is_double_layer) {
        cout << "[RESULT] Double-Layer Structure Detected!" << endl;

        double opd1 = result.opd_layer_1;
        double opd2 = result.opd_layer_2;
        if (opd1 > opd2) std::swap(opd1, opd2);

        cout << "  Layer Signal A (OPD): " << opd1 << " um" << endl;
        cout << "  Layer Signal B (OPD): " << opd2 << " um" << endl;
        cout << "  Total Sum Peak (OPD): " << result.opd_total << " um" << endl;
        cout << "  Verification Error:   " << result.error << " um" << endl;
        cout << "  Confidence Score:     " << fixed << setprecision(0) << result.score << endl;

        cout << "\n[ANALYSIS] Material Matching Matrix:" << endl;
        cout << setw(12) << "Signal" << " | ";
        for (const auto& m : material_db) cout << setw(10) << m.name << " | ";
        cout << endl << string(60, '-') << endl;

        auto print_row = [&](double opd, string label) {
            cout << setw(12) << label << " | ";
            for (const auto& m : material_db) {
                double d = opd / (2.0 * m.refractive_index);
                cout << setw(10) << fixed << setprecision(1) << d << " | ";
            }
            cout << endl;
            };

        print_row(opd1, "Signal A");
        print_row(opd2, "Signal B");

        cout << endl << string(60, '-') << endl;
        cout << "Auto-Interpretation (Assuming High-n & Low-n mix):" << endl;

        // Scenario 1: opd1 is low index, opd2 is high index
        double d_poly = opd1 / (2.0 * 1.5);
        double d_si = opd2 / (2.0 * 3.7);
        cout << "  > Scenario 1: Polymer (" << d_poly << " um) + Silicon (" << d_si << " um)" << endl;

        // Scenario 2: opd1 is high index, opd2 is low index
        double d_si_2 = opd1 / (2.0 * 3.7);
        double d_poly_2 = opd2 / (2.0 * 1.5);
        cout << "  > Scenario 2: Silicon (" << d_si_2 << " um) + Polymer (" << d_poly_2 << " um)" << endl;

    }
    else if (!peaks.empty()) {
        cout << "[RESULT] Single-Layer Structure Detected (Most likely)." << endl;
        double best_opd = result.opd_layer_1;

        cout << "  Dominant Signal (OPD): " << best_opd << " um" << endl;
        cout << "  Physical Thickness Candidates:" << endl;
        for (const auto& m : material_db) {
            cout << "    - If " << setw(8) << m.name << " (n=" << m.refractive_index << "): "
                << (best_opd / (2.0 * m.refractive_index)) << " um" << endl;
        }
    }
    else {
        cout << "[RESULT] No clear structure found." << endl;
    }

    cout << string(60, '=') << endl;
    system("pause");
    return 0;
}