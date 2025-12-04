#include "ThicknessMeasurement.h"

int main() {
    // Generate wavelength array
    VectorXd lambda(PIXEL_NUM);
    for (Eigen::Index i = 0; i < PIXEL_NUM; ++i) {
        lambda(i) = LAMBDA_MIN + (LAMBDA_MAX - LAMBDA_MIN) * i / (PIXEL_NUM - 1);
    }

    // Simulation input: Assume true thickness d_true=9μm
    double d_true = 9.0;
    double refractive = 1.5;
    double z_true = 2 * d_true * refractive; // Round-trip optical path difference
    VectorXd I_raw = simulateInterferenceIntensity(0.5, 0.3, z_true, lambda);

    // Calculate some statistical information about I_raw
    cout << "I_raw mean: " << I_raw.mean() << endl;
    cout << "I_raw maximum: " << I_raw.maxCoeff() << endl;
    cout << "I_raw minimum: " << I_raw.minCoeff() << endl;

    // EMD filtering
    VectorXd imf1 = emdFilter(I_raw);
    
    // Debug: Check EMD filtered signal
    cout << "EMD filtered signal mean: " << imf1.mean() << endl;
    cout << "EMD filtered signal maximum: " << imf1.maxCoeff() << endl;
    cout << "EMD filtered signal minimum: " << imf1.minCoeff() << endl;
    
    // Wavenumber transformation
    VectorXd sigma, imf1_sigma;
    lambdaToSigma(lambda, imf1, sigma, imf1_sigma);
    
    // Debug: Check first few values of sigma and imf1_sigma
    cout << "First 5 sigma values: " << endl;
    for (int i = 0; i < 5; ++i) {
        cout << sigma(i) << " " << imf1_sigma(i) << endl;
    }

    // LSP analysis
    double z_peak = lombScargle(sigma, imf1_sigma);
    
    // Debug: Check LSP result
    cout << "LSP peak index: " << endl;

    // Calculate film thickness
    double d_calc = calculateThickness(z_peak, refractive);

    // Output results
    cout << "========================================" << endl;
    cout << "          Film Thickness Measurement Results              " << endl;
    cout << "========================================" << endl;
    cout << "True thickness:    " << fixed << setprecision(3) << d_true << " um" << endl;
    cout << "True optical path:  " << fixed << setprecision(3) << z_true << " um" << endl;
    cout << "Measured optical path:  " << fixed << setprecision(3) << z_peak << " um" << endl;
    cout << "Calculated thickness:    " << fixed << setprecision(3) << d_calc << " um" << endl;
    cout << "Measurement error:    " << fixed << setprecision(3) << abs(d_calc - d_true) << " um" << endl;
    cout << "Accurate measurement: " << (abs(d_calc - d_true) <= 0.1 ? "Yes (<=0.1 um)" : "No") << endl;
    cout << "========================================" << endl;

    return 0;
}
