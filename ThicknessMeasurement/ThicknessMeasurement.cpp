#define _USE_MATH_DEFINES
#include "ThicknessMeasurement.h"
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace Eigen;

// Physics Constants
const double MU = 1315.0;
const double ETA = 25.0;
const double C = 18000.0;
const double PHI0 = M_PI_2;

// ================= EMD Filter =================
VectorXd emdFilter(const VectorXd& signal) {
    VectorXd x = signal;
    VectorXd imf1(signal.size());
    const int max_iter = 20;

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<int> max_idx, min_idx;
        for (int i = 1; i < x.size() - 1; ++i) {
            if (x(i) >= x(i - 1) && x(i) > x(i + 1)) max_idx.push_back(i);
            else if (x(i) <= x(i - 1) && x(i) < x(i + 1)) min_idx.push_back(i);
        }
        if (max_idx.size() < 2 || min_idx.size() < 2) { imf1 = x; break; }

        VectorXd upper = x, lower = x;
        for (size_t i = 0; i + 1 < max_idx.size(); ++i)
            upper.segment(max_idx[i], max_idx[i + 1] - max_idx[i] + 1)
            .setLinSpaced(max_idx[i + 1] - max_idx[i] + 1, x(max_idx[i]), x(max_idx[i + 1]));
        for (size_t i = 0; i + 1 < min_idx.size(); ++i)
            lower.segment(min_idx[i], min_idx[i + 1] - min_idx[i] + 1)
            .setLinSpaced(min_idx[i + 1] - min_idx[i] + 1, x(min_idx[i]), x(min_idx[i + 1]));

        VectorXd m = 0.5 * (upper + lower);
        imf1 = x - m;
        x = imf1;
    }
    return imf1;
}

// ================= Lambda to Sigma =================
void lambdaToSigma(const VectorXd& lambda_nm, const VectorXd& imf1, VectorXd& sigma, VectorXd& imf1_sigma) {
    sigma.resize(lambda_nm.size());
    imf1_sigma.resize(lambda_nm.size());
    for (Eigen::Index i = 0; i < lambda_nm.size(); ++i) {
        double lambda_um = lambda_nm(i) * 1e-3;
        sigma(i) = (lambda_um > 1e-6) ? (1.0 / lambda_um) : 0.0;
        imf1_sigma(i) = imf1(i);
    }
}

// ================= Calculate Power Spectrum =================
SpectrumData calculatePowerSpectrum(const VectorXd& sigma, const VectorXd& imf1_sigma, const AnalysisConfig& config) {
    int z_num = int((config.z_search_max - config.z_search_min) / config.z_search_step) + 1;

    SpectrumData result;
    result.z_axis.resize(z_num);
    result.power.resize(z_num);

    for (int i = 0; i < z_num; ++i) {
        double z = config.z_search_min + i * config.z_search_step;
        double sum_c = 0, sum_s = 0, sum_cc = 0, sum_ss = 0;

        for (int k = 0; k < sigma.size(); ++k) {
            double w = 2.0 * M_PI * z * sigma(k);
            double c = std::cos(w);
            double s = std::sin(w);
            double val = imf1_sigma(k);
            sum_c += val * c;
            sum_s += val * s;
            sum_cc += c * c;
            sum_ss += s * s;
        }

        double P = 0.0;
        if (sum_cc > 1e-9 && sum_ss > 1e-9)
            P = (sum_c * sum_c) / sum_cc + (sum_s * sum_s) / sum_ss;

        result.z_axis(i) = z;
        result.power(i) = P;
    }
    return result;
}

// ================= Find Peaks (Robust) =================
std::vector<PeakInfo> findPeaks(const SpectrumData& spectrum, const AnalysisConfig& config) {
    std::vector<PeakInfo> raw_peaks;
    double max_p = spectrum.power.maxCoeff();
    double threshold = max_p * config.peak_threshold;

    // 1. Local maxima detection
    for (int i = 1; i < spectrum.power.size() - 1; ++i) {
        if (spectrum.power(i) > spectrum.power(i - 1) &&
            spectrum.power(i) > spectrum.power(i + 1) &&
            spectrum.power(i) > threshold) {
            raw_peaks.push_back({ spectrum.z_axis(i), spectrum.power(i) });
        }
    }

    // 2. Sort by Power descending
    std::sort(raw_peaks.begin(), raw_peaks.end(), [](const PeakInfo& a, const PeakInfo& b) {
        return a.power > b.power;
        });

    // 3. Distance suppression (Non-Maximum Suppression)
    std::vector<PeakInfo> final_peaks;
    for (const auto& cand : raw_peaks) {
        bool keep = true;
        for (const auto& exist : final_peaks) {
            if (std::abs(cand.z_opd - exist.z_opd) < config.min_peak_dist) {
                keep = false; break;
            }
        }
        if (keep) final_peaks.push_back(cand);
    }

    // 4. Sort result by OPD for logic processing
    std::sort(final_peaks.begin(), final_peaks.end(), [](const PeakInfo& a, const PeakInfo& b) {
        return a.z_opd < b.z_opd;
        });
    return final_peaks;
}

// ================= Topology Solver (Smart Scoring) =================
LayerStructure solveLayerTopology(const std::vector<PeakInfo>& peaks, const AnalysisConfig& config) {
    LayerStructure best_result;
    best_result.is_double_layer = false;
    best_result.opd_layer_1 = 0;
    best_result.opd_layer_2 = 0;
    best_result.opd_total = 0;
    best_result.error = 99999.0;
    best_result.score = -1.0;

    int n = (int)peaks.size();
    if (n < 2) {
        if (n == 1) {
            best_result.opd_layer_1 = peaks[0].z_opd;
            best_result.opd_total = peaks[0].z_opd;
            best_result.error = 0.0;
        }
        return best_result;
    }

    // Iterate through all pairs to find A + B = C
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double sum = peaks[i].z_opd + peaks[j].z_opd;

            // Check if 'sum' exists in other peaks
            for (int k = 0; k < n; ++k) {
                // Constraints: 
                // 1. Total must be distinct from components
                if (k == i || k == j) continue;
                // 2. Total must be larger than components (Physics)
                if (peaks[k].z_opd <= peaks[i].z_opd || peaks[k].z_opd <= peaks[j].z_opd) continue;

                double diff = std::abs(peaks[k].z_opd - sum);

                if (diff < config.additivity_tol) {
                    // --- SCORING LOGIC ---
                    // We value HIGH INTENSITY and LOW ERROR.
                    // Total Power = P_A + P_B + P_Total
                    double total_power = peaks[i].power + peaks[j].power + peaks[k].power;

                    // Penalty Factor: 1.0 at 0 error, decaying linearly to 0.5 at max tolerance
                    double error_penalty = 1.0 - (diff / config.additivity_tol) * 0.5;
                    if (error_penalty < 0.1) error_penalty = 0.1; // Safety floor

                    double current_score = total_power * error_penalty;

                    // Update if this combination is better
                    if (current_score > best_result.score) {
                        best_result.is_double_layer = true;
                        best_result.opd_layer_1 = peaks[i].z_opd;
                        best_result.opd_layer_2 = peaks[j].z_opd;
                        best_result.opd_total = peaks[k].z_opd;
                        best_result.error = diff;
                        best_result.score = current_score;
                    }
                }
            }
        }
    }

    // Fallback: If no double layer found, return the strongest peak as single layer
    if (!best_result.is_double_layer && n > 0) {
        double max_p = -1;
        for (const auto& p : peaks) {
            if (p.power > max_p) {
                max_p = p.power;
                best_result.opd_layer_1 = p.z_opd;
                best_result.opd_total = p.z_opd;
            }
        }
        best_result.error = 0.0;
    }

    return best_result;
}