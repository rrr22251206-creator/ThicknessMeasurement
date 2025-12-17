#include "SpectrumLoader.h"
#include <fstream>
#include <sstream>
#include <iostream>

bool SpectrumLoader::load(
    const std::string& filename,
    std::vector<double>& lambda,
    std::vector<double>& intensity,
    bool wavelengthInNm)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "[ERROR] Cannot open spectrum file: "
            << filename << std::endl;
        return false;
    }

    lambda.clear();
    intensity.clear();

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        // 支持 CSV / 空格
        for (char& c : line)
            if (c == ',') c = ' ';

        std::stringstream ss(line);
        double wl, I;
        ss >> wl >> I;

        if (ss.fail())
            continue;

        // 统一内部单位：nm
        if (!wavelengthInNm)
            wl *= 1e9;  // m → nm

        lambda.push_back(wl);
        intensity.push_back(I);
    }

    if (lambda.size() < 16) {
        std::cerr << "[ERROR] Spectrum data too short\n";
        return false;
    }

    std::cout << "[INFO] Loaded spectrum points: "
        << lambda.size() << std::endl;

    return true;
}
