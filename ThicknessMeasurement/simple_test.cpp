#include <iostream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>

#include "ThicknessMeasurement.h"
#include "SpectrumLoader.h"

using namespace std;
using namespace Eigen;

int main()
{
    // ===== 1. 读取光谱文件 =====
    vector<double> lambda_nm_vec;
    vector<double> intensity_vec;

    if (!SpectrumLoader::load(
        "spectrum.txt",
        lambda_nm_vec,
        intensity_vec))
    {
        return -1;
    }

    // ===== 2. 转 Eigen =====
    int N = static_cast<int>(lambda_nm_vec.size());
    VectorXd lambda(N), I(N);

    for (int i = 0; i < N; ++i) {
        lambda(i) = lambda_nm_vec[i];   // nm
        I(i) = intensity_vec[i];
    }

    // ===== 3. λ → σ =====
    VectorXd sigma, I_sigma;
    lambdaToSigma(lambda, I, sigma, I_sigma);

    // ===== 4. LSP =====
    double z_peak = lombScargle(sigma, I_sigma);
    double d_calc = calculateThickness(z_peak);

    // ===== 5. 输出 =====
    cout << fixed << setprecision(4);
    cout << "Estimated thickness : "
        << d_calc << " um" << endl;

    system("pause"); // Keep console open
    return 0;
}
