#pragma once
#include "cavity_common.h"

// ======================== 调试/打印 ========================

// 打印所有 OneLine 线段的起终点坐标和 faceId（调试用）
void PrintAllLines(const std::vector<OneLine>& lines);

// 打印端点间隙诊断信息：统计开放端点数量、最大间隙等（调试用）
void PrintEndpointGapDiagnostics(const std::vector<OneLine>& lines, double snapTol);

// ======================== 线段去重/修复 ========================

// 对线段列表进行去重（基于起终点相同 + faceId 相同判定重复）
std::vector<OneLine> DeduplicateLines(const std::vector<OneLine>& lines);

// 修复开放端点间隙：将距离小于 healTol 的端点吸附到一起，返回修复的间隙数量
int HealOpenEndpointGaps(std::vector<OneLine>& lines, double healTol);

// ======================== 成环算法 ========================

// 将线段列表连接成闭合环（基于端点邻接关系构建），返回所有闭合环
std::vector<Loop> BuildLoops(const std::vector<OneLine>& inputLines);

// ======================== 拓扑构建 ========================

// 从线段列表完整构建拓扑：去重→修复→成环→嵌套分析→提取 Face2D
std::vector<Face2D> BuildTopologyAndExtractFaces(const std::vector<OneLine>& layerLines);

// 从嵌套环树中递归提取 Face2D（奇数层=外环，偶数层=内环）
void ExtractFacesFromTree(LoopNode* node, int depth, std::vector<Face2D>& faces);

// 将环列表导出为 BREP 文件（调试用，每个环单独编号）
void DebugExportLoops(const std::vector<Loop>& loops, const std::string& prefix, const std::string& savePath);
