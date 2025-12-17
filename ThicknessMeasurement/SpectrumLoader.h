#pragma once
#include <vector>
#include <string>

/**
 * @class SpectrumLoader
 * @brief 原始光谱数据加载模块（工程输入层）
 *
 * 支持 TXT / CSV 文件，格式示例：
 *   500.0  12345
 *   500.2, 12410
 *
 * 波长单位默认 nm
 */
class SpectrumLoader {
public:
    /**
     * @brief 加载光谱数据文件
     * @param filename 文件路径
     * @param lambda 输出波长数组（单位：nm）
     * @param intensity 输出强度数组
     * @param wavelengthInNm 输入是否为 nm（默认 true）
     * @return 是否成功
     */
    static bool load(
        const std::string& filename,
        std::vector<double>& lambda,
        std::vector<double>& intensity,
        bool wavelengthInNm = true
    );
};
