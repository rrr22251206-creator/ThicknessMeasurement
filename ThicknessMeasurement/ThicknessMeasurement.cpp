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
            double val = imf1_sigma(k); // 直接使用原始信号，不加窗

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

// ================= 物理约束拓扑求解器 =================
LayerStructure solveLayerTopology(const std::vector<PeakInfo>& peaks, const AnalysisConfig& config) {
    LayerStructure best_result;
    best_result.is_double_layer = false;
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

    // 1. 找到全局最强峰 (Dominant Peak)
    // 物理公理：在高折射率差的薄膜中，最强的干涉峰通常对应某个单层的物理厚度，
    // 而不是两层叠加的和频（和频能量衰减大）。
    double max_peak_power = -1.0;
    double dominant_opd = -1.0;
    for (const auto& p : peaks) {
        if (p.power > max_peak_power) {
            max_peak_power = p.power;
            dominant_opd = p.z_opd;
        }
    }

    // 2. 遍历寻找 A + B = C
    const double sigma = 15.0;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double sum = peaks[i].z_opd + peaks[j].z_opd;

            for (int k = 0; k < n; ++k) {
                // 基本约束
                if (k == i || k == j) continue;
                if (peaks[k].z_opd <= peaks[i].z_opd || peaks[k].z_opd <= peaks[j].z_opd) continue;

                double diff = std::abs(peaks[k].z_opd - sum);

                if (diff < config.additivity_tol) {

                    // --- 评分逻辑优化 ---
                    double total_power = peaks[i].power + peaks[j].power + peaks[k].power;
                    double gaussian_penalty = std::exp(-(diff * diff) / (2.0 * sigma * sigma));
                    double current_score = total_power * gaussian_penalty;

                    // --- 【关键逻辑】物理角色判定 ---
                    bool component_is_dominant =
                        (std::abs(peaks[i].z_opd - dominant_opd) < 1.0) ||
                        (std::abs(peaks[j].z_opd - dominant_opd) < 1.0);

                    bool sum_is_dominant =
                        (std::abs(peaks[k].z_opd - dominant_opd) < 1.0);

                    // 规则 1: 如果 Total Sum 是最强峰，这在物理上极不可能 (除非是特殊谐振腔)
                    // 我们对其进行严厉惩罚 (或者直接 continue 跳过)
                    if (sum_is_dominant) {
                        current_score *= 0.001; // 惩罚 1000 倍
                    }
                    // 规则 2: 如果 Layer 1 或 Layer 2 是最强峰，这非常合理，给予奖励
                    else if (component_is_dominant) {
                        current_score *= 10.0;  // 奖励 10 倍
                    }

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

    // 单层回退
    if (!best_result.is_double_layer && n > 0) {
        best_result.opd_layer_1 = dominant_opd;
        best_result.opd_total = dominant_opd;
        best_result.error = 0.0;
    }

    return best_result;
}