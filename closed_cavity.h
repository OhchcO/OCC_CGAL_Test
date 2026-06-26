#pragma once
#include "cavity_common.h"

// ======================== 封闭型腔切分 ========================

// 递归切分：在指定 Z 高度处将一个 faceId 集合切分为上下两部分，递归处理直到所有 Z 层遍历完毕
void SplitCavityRecursive(const std::set<int>& currentFaceIds, TopTools_DataMapOfIntegerShape& idToFace, int& nextId, const std::vector<double>& splitZs, int zIndex, std::vector<std::set<int>>& result);

// 对一个 Compound 型腔执行完整 Z 轴切分，返回切分后的多个子 Compound
std::vector<TopoDS_Compound> SplitCavity(const TopoDS_Compound& cavity);

// 将一个 Compound 中互不相连的部分分离为独立的 Compound（基于面邻接关系判断连通性）
std::vector<TopoDS_Compound> SeparateDisconnectedCavities(const TopoDS_Compound& cavityCompound);

// ======================== 特征生成 ========================

// 将切分后的子 Compound 列表转换为 CavityFeature 列表（包含类型、Z范围、面ID等信息）
std::vector<CavityFeature> GenerateCavityFeatures(const std::vector<TopoDS_Compound>& isolatedCavities, const std::vector<std::vector<Face2D>>& allLayerFaces, const std::vector<double>& sortedSplitPoints, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// 从单个 CavityFeature 中提取真正的型腔面（壁面+盖面），用于可视化和导出
TopoDS_Compound ExtractTrueCavityFaces(const CavityFeature& feat, const TopoDS_Shape& mainShape, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// 将多个子 Compound 转换为独立的 CavityFeature（用于 SplitCavity 结果的后续处理）
std::vector<CavityFeature> ConvertPartsToFeatures(const std::vector<TopoDS_Compound>& parts, const CavityFeature& parentFeat, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// ======================== 主流程 ========================

// 封闭型腔处理主流程：从切层数据中提取封闭型腔、Z 轴切分、生成特征、输出结果
void ProcessAndSplitClosedCavityFeatures(const std::vector<std::vector<Face2D>>& allLayerFaces, const TopoDS_Shape& mainShape, const std::vector<double>& splitPoints, const TopTools_DataMapOfShapeInteger& faceToIdMap, const std::string& savePath, std::vector<TopoDS_Compound>& outTrueClosedCavities, std::vector<CavityFeature>& outFeatures);
