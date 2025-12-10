#include <iostream>
#include <cmath>    // For math operations
#include <iomanip>  // For output formatting

using namespace std;

/*
 * ======================================================================================
 * Derivation of the Constraint Formula
 * ======================================================================================
 *
 * 1. Interference Principle:
 *    In Low Coherence Spectral Interferometry (LSP), reflections from the top and bottom
 *    surfaces of a thin film interfere. The condition for interference extrema is:
 *    OPD = 2 * n * d = m * lambda
 *    Where:
 *      n: Refractive index
 *      d: Physical thickness
 *      m: Interference order (integer)
 *      lambda: Wavelength
 *
 * 2. Spectral Periodicity:
 *    In the wavelength domain, the spacing between adjacent interference fringes
 *    (Delta_lambda) is approximated by:
 *    Delta_lambda_period approx lambda^2 / (2 * n * d)
 *    This implies: As thickness (d) increases, fringes become denser (Delta_lambda decreases).
 *
 * 3. Nyquist Sampling Theorem:
 *    To resolve the sinusoidal interference fringes, the spectrometer's sampling frequency
 *    must be at least twice the signal frequency.
 *    The sampling interval is the optical resolution (delta_lambda).
 *    According to Nyquist:
 *    delta_lambda <= Delta_lambda_period / 2
 *
 * 4. Final Calculation:
 *    Substituting Delta_lambda_period from step 2 into step 3:
 *    delta_lambda <= (lambda^2 / (2 * n * d)) / 2
 *    delta_lambda <= lambda^2 / (4 * n * d)
 *
 *    Rearranging to solve for maximum thickness d_max:
 *    4 * n * d * delta_lambda <= lambda^2
 *
 *    ==> d_max = lambda^2 / (4 * n * delta_lambda)
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