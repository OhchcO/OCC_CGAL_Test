#pragma once
#include "cavity_common.h"

// ======================== 开放型腔清洗 ========================

// 计算 Face2D 的紧凑度（面积/周长²），用于过滤过于狭长的噪声面
double CalculateFaceCompactness(const Face2D& face);

// 构建圆角半径映射：遍历模型所有面，识别圆角面并记录其 faceId 和半径
std::map<int, double> BuildFilletRadiusMap(const TopoDS_Shape& shape, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// 判断开放型腔面是否由小圆角产生（基于面中包含的圆角面半径判断），可选输出匹配的圆角信息
bool IsSmallFilletCavity(const Face2D& face, const std::map<int, double>& filletMap, double maxFilletRadius, std::vector<std::pair<int, double>>* matchedFillets = nullptr);

// 清洗并过滤开放型腔面：去除面积过小、紧凑度过低、刀具无法通过、圆角过小等不合格面
std::vector<Face2D> CleanAndFilterOpenCavities(const std::vector<Face2D>& rawOpenCavities, const std::map<int, double>& filletMap, const OpenCavityFilterParams& params = OpenCavityFilterParams{});

// ======================== 开放型腔提取 ========================

// 从所有层的开放型腔 Face2D 中收集所有出现过的 sourceFaceId
std::set<int> CollectOpenCavityFaceIds(const std::vector<std::vector<Face2D>>& allLayerFaces);

// 获取开放型腔的盖面（与型腔壁面相邻的水平面），用于可视化
TopoDS_Compound GetOpenCavityCapFaces(const TopoDS_Shape& solid, const std::set<int>& openFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, double angleTolerance = 0.001);

// 将开放型腔面导出为 BREP 文件（调试用）
void ExportOpenCavityFaces(const TopoDS_Shape& solid, const std::set<int>& openFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, const std::string& fileName);

// 获取开放型腔的完整 Compound（壁面），用于后续处理
TopoDS_Compound GetOpenCavityCompound(const TopoDS_Shape& solid, const std::set<int>& openFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// 将切分后的开放型腔子 Compound 列表转换为 CavityFeature 列表
std::vector<CavityFeature> GenerateOpenCavityFeatures(const std::vector<TopoDS_Compound>& isolatedCavities, const std::vector<std::vector<Face2D>>& allLayerFaces, const std::vector<double>& sortedSplitPoints, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// 从单个开放型腔 CavityFeature 中提取真正的型腔面（排除被封闭型腔占用的面）
TopoDS_Compound ExtractTrueOpenCavityFaces(const CavityFeature& openFeat, const TopoDS_Shape& mainShape, const TopTools_DataMapOfShapeInteger& faceToIdMap, const std::vector<std::vector<Face2D>>& allClosedLayerFaces, const std::vector<double>& splitPoints, const std::vector<TopoDS_Compound>& trueClosedCavities);

// ======================== 嵌套合并 ========================

// 判断 innerComp 是否完全被 outerComp 包含（精确几何包含检测）
bool IsCompoundInsideExact(const TopoDS_Compound& innerComp, const TopoDS_Compound& outerComp);

// 合并嵌套的岛状型腔：内层被外层完全包含时合并为一个 Compound
std::vector<TopoDS_Compound> MergeNestedIslands(const std::vector<TopoDS_Compound>& bfsParts);

// ======================== Compound Z 分割 ========================

// 在指定 Z 高度将 Compound 切分为上下两部分（使用 BRepAlgoAPI_Section）
std::pair<TopoDS_Compound, TopoDS_Compound> SplitCompoundAtZ(const TopoDS_Compound& cavity, double splitZ);

// 递归切分 Compound：在所有水平切分点处切分，直到每个子部分不再跨层
std::vector<TopoDS_Compound> RecursiveSplitCavity(const TopoDS_Compound& cavity);

// 合并距离过近的切分 Z 值（避免产生过薄的切片层）
std::vector<double> MergeCloseSplitZs(std::vector<double> splitZs, double minGap);

// ======================== 开放型腔 Z 切分 ========================

// 收集开放型腔 Compound 的切分 Z 值列表（基于型腔几何和盖面位置）
std::vector<double> CollectOpenCavitySplitZs(const TopoDS_Compound& cavity, const CavityFeature& feature, const OpenCavitySplitParams& params);

// 安全获取 Compound 的 Z 范围（即使 Compound 为空也不崩溃）
bool GetCompoundZRangeSafe(const TopoDS_Compound& compound, double& zmin, double& zmax);

// 对单个子部分收集切分 Z 值（用于递归切分场景）
std::vector<double> CollectOpenCavitySplitZsForPart(const TopoDS_Compound& part, const CavityFeature& feature, const OpenCavitySplitParams& params, const std::string& branchName);

// 将几何上相连但被 Z 切分打断的开放型腔子部分重新合并（避免过度碎片化）
std::vector<TopoDS_Compound> MergeGeometricallyConnectedOpenParts(const std::vector<TopoDS_Compound>& rawParts, double tolerance, int featureId, const std::string& branchName, int depth, bool exportDebugBreps);

// 在水平面位置对开放型腔进行 Z 切分
std::vector<TopoDS_Compound> SplitOpenCavityByZPlan(const TopoDS_Compound& cavity, const CavityFeature& feature, const OpenCavitySplitParams& params, int featureId);

// 保守策略切分开放型腔：先找切分点，再切分，最后合并连通部分
std::vector<TopoDS_Compound> SplitOpenCavityConservatively(const TopoDS_Compound& cavity, const CavityFeature& feature, const OpenCavitySplitParams& params);

// ======================== 主流程 ========================

// 开放型腔处理主流程：收集 ID → 构建 Compound → 合并嵌套 → Z 切分 → 生成特征
std::vector<CavityFeature> ProcessAndSplitOpenCavityFeatures(const std::vector<std::vector<Face2D>>& allOpenLayerFaces, const std::vector<std::vector<Face2D>>& allClosedLayerFaces, const std::vector<TopoDS_Compound>& trueClosedCavities, const TopoDS_Shape& mainShape, const std::vector<double>& splitPoints, const TopTools_DataMapOfShapeInteger& faceToIdMap, const std::string& savePath);
