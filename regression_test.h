#pragma once
#include "cavity_common.h"
#include <string>
#include <vector>

// ======================== 回归测试数据结构 ========================

// 单个原始面在特征中的记录：faceId、是否被拆分、子区域数量
struct RegressionOriginalFaceRecord {
    int originalFaceId = 0;   // 原始面 ID
    bool isSplit = false;     // 该面在切片中是否被拆分为多个子区域
    int partCount = 1;        // 子区域数量（1=未拆分）
};

// 单个特征的回归记录：类型、Z范围、原始面集合、拆分信息、签名
struct RegressionFeatureRecord {
    std::string featureId;                        // 特征 ID
    std::string type;                             // 特征类型：OPEN / CLOSED / OTHER
    double topZ = 0.0;                            // 特征顶部 Z 坐标
    double bottomZ = 0.0;                         // 特征底部 Z 坐标
    std::vector<int> originalFaceIds;             // 参与该特征的所有原始面 ID 列表
    bool hasFaceSplit = false;                    // 是否有原始面被拆分
    std::vector<int> splitOriginalFaces;          // 被拆分的原始面 ID 列表
    std::vector<RegressionOriginalFaceRecord> originalFaces;  // 每个原始面的详细拆分记录
    std::string signature;                        // 特征签名字符串（用于快速比对）
};

// 单个模型的完整回归结果：模型名、运行 ID、所有特征、全局统计
struct RegressionModelResult {
    std::string model;                            // 模型文件名
    std::string runId;                            // 运行标识（如时间戳）
    std::vector<RegressionFeatureRecord> features;  // 所有识别到的特征
    int featureCount = 0;                         // 特征总数
    int splitOriginalFaceCount = 0;               // 被拆分的原始面总数（含重复计数）
    int sharedOriginalFaceCount = 0;              // 被多个特征共享的原始面数量
    std::string featureSignature;                 // 所有特征类型拼接的签名（如 OPEN|CLOSED）
};

// 对比结果枚举
enum class RegressionCompareResult {
    PASS,        // 与 baseline 完全一致
    CHANGE,      // 有变化但不算回归（如新增特征、拆分数量变化）
    REGRESSION   // 回归问题（如特征丢失、类型变化）
};

// 对比报告：总体结果 + 差异原因列表
struct RegressionDiffReport {
    RegressionCompareResult overallResult = RegressionCompareResult::PASS;
    std::vector<std::string> reasons;             // 导致判定的具体原因列表
};

// ======================== 回归测试函数 ========================

// 将本轮识别结果保存为 JSON 文件（包含每个特征的类型、Z范围、原始面集合、拆分信息）
void SaveModelRegressionResult(
    const std::string& outputPath,                // 输出 JSON 文件路径
    const std::string& modelName,                 // 模型文件名
    const std::string& runId,                     // 运行标识
    const std::vector<CavityFeature>& features,   // 识别到的特征列表
    const std::vector<std::vector<Face2D>>& allLayerFaces);  // 所有层的 Face2D（用于统计面拆分）

// 对比两个 JSON 结果文件，返回差异报告（PASS/CHANGE/REGRESSION）
RegressionDiffReport CompareModelRegressionResults(
    const std::string& baselineJsonPath,          // baseline JSON 文件路径
    const std::string& currentJsonPath,           // 当前结果 JSON 文件路径
    double zTol = 0.05);                          // Z 坐标容差
