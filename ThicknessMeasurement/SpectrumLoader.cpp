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
        if (line.empty())
            continue;

        // 跳过非数字开头的行（标题 / header）
        if (!(isdigit(line[0]) || line[0] == '.' || line[0] == '-'))
            continue;

        for (char& c : line)
            if (c == ',') c = ' ';

        std::stringstream ss(line);
        double wl, I;
        if (!(ss >> wl >> I))
            continue;

        // 统一：内部单位 = nm
        lambda.push_back(wl);
        intensity.push_back(I);
    }

    std::cout << "[INFO] Loaded spectrum points: "
        << lambda.size() << std::endl;

    if (lambda.size() < 16) {
        std::cerr << "[ERROR] Spectrum data too short\n";
        return false;
    }

    return true;
}
