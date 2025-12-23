#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

// 定义数据点结构
struct Point3D {
    double x;
    double y;
    double z;
};

class SurfaceSmoother {
private:
    std::vector<Point3D> rawData;
    std::vector<Point3D> smoothedData;

public:
    // 1. 加载数据
    // 使用流式读取，忽略换行符的位置，每读取3个数字作为一个点
    // 这能有效解决OCR结果中换行混乱的问题
    bool loadData(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open file " << filename << std::endl;
            return false;
        }

        double val;
        std::vector<double> buffer;
        while (file >> val) {
            // 简单的过滤：忽略异常值
            // 在实际工程中，可以根据坐标范围过滤
            if (buffer.empty() && val > 200) continue; // 假设X坐标不会超过200

            buffer.push_back(val);
            if (buffer.size() == 3) {
                rawData.push_back({ buffer[0], buffer[1], buffer[2] });
                buffer.clear();
            }
        }
        std::cout << "Loaded " << rawData.size() << " points." << std::endl;
        return true;
    }

    // 2. 高斯平滑算法
    // sigma: 平滑力度。值越大越平滑，但可能丢失细节；值越小越接近原始数据。
    // searchRadius: 搜索邻居的半径，通常设为 3 * sigma。
    void applyGaussianSmoothing(double sigma) {
        if (rawData.empty()) return;

        smoothedData = rawData; // 复制坐标，Z值将被更新
        double searchRadius = 3.0 * sigma;
        double sigmaSq2 = 2.0 * sigma * sigma;

        for (size_t i = 0; i < rawData.size(); ++i) {
            double sumWeightedZ = 0.0;
            double sumWeight = 0.0;

            // 遍历所有点寻找邻居
            for (size_t j = 0; j < rawData.size(); ++j) {
                double dx = rawData[i].x - rawData[j].x;
                double dy = rawData[i].y - rawData[j].y;
                double distSq = dx * dx + dy * dy;

                // 如果在搜索半径内
                if (distSq < searchRadius * searchRadius) {
                    // 计算高斯权重
                    double weight = std::exp(-distSq / sigmaSq2);

                    sumWeightedZ += rawData[j].z * weight;
                    sumWeight += weight;
                }
            }

            if (sumWeight > 0.0) {
                smoothedData[i].z = sumWeightedZ / sumWeight;
            }
        }
        std::cout << "Smoothing completed with sigma = " << sigma << std::endl;
    }

    // 3. 导出数据 (CSV格式，方便用Excel查看曲线)
    void saveToCSV(const std::string& filename) {
        std::ofstream file(filename);
        file << "X,Y,Original_Z,Smoothed_Z\n";
        file << std::fixed << std::setprecision(3);

        for (size_t i = 0; i < rawData.size(); ++i) {
            file << rawData[i].x << ","
                << rawData[i].y << ","
                << rawData[i].z << ","
                << smoothedData[i].z << "\n";
        }
        std::cout << "Data saved to " << filename << std::endl;
    }

    // 简单的优化：剔除离群点 (Outlier Removal)
    // 如果某点与邻居的平均值差异过大，则将其拉回
    void removeOutliers(double threshold = 1.0) {
        // 此处逻辑与平滑类似，但在平滑前执行，用于去除极端的测量错误
        // 为保持代码简洁，这里主要依赖高斯平滑本身对噪声的抑制能力
    }
};

int main() {
    SurfaceSmoother smoother;

    // 假设OCR数据保存为 input.txt
    // 请确保 input.txt 中只包含数字，去除 "Start of PDF" 等文字
    if (smoother.loadData("input.txt")) {

        // 参数调整建议：
        // Sigma = 5.0 到 10.0 : 适用于保留较多局部特征
        // Sigma = 15.0 到 20.0 : 适用于生成非常平滑的曲面
        // 根据您的数据坐标间距（X间距约20，Y间距约5），建议 Sigma 设为 10.0 左右
        smoother.applyGaussianSmoothing(10.0);

        smoother.saveToCSV("output.csv");
    }
    else {
        // 如果没有文件，生成一些模拟数据测试
        std::cout << "Please provide 'input.txt' with your data." << std::endl;
    }

    return 0;
}