#pragma once
#include "cavity_common.h"

// ======================== 中部开放型腔 ========================

// 获取 Face2D 的 Z 坐标（取外环第一个顶点的 Z 值）
double GetFace2DZForInteriorOpen(const Face2D& face);

// 判断是否存在中部开放型腔：检查开放型腔是否被封闭型腔在垂直方向上"夹住"
bool HasInteriorOpenCavity(const std::vector<std::vector<Face2D>>& allOpenLayerFaces, const std::vector<std::vector<Face2D>>& allClosedLayerFaces, double modelMaxZ, const std::vector<double>& splitPoints, double zTol = 0.5);

// 在切分点列表中查找紧邻当前 Z 下方的下一个切分点（用于确定中部开放区域的下边界）
double FindNextLowerSplitPoint(const std::vector<double>& splitPoints, double z, double fallbackZ);

// 将 Face2D 中所有线段的 faceId 添加到集合中（辅助函数）
void AddFaceIdsFromFace2D(const Face2D& face, std::set<int>& ids);

// 打印 faceId 集合（调试用）
void PrintIdSet(const std::string& label, const std::set<int>& ids);

// 收集中部开放区域周围的封闭型腔侧壁面 ID（用于判断是否被封闭型腔包围）
std::set<int> CollectClosedSideFaceIdsAroundOpenBand(const std::vector<std::vector<Face2D>>& allClosedLayerFaces, double openTopZ, double openBottomZ);

// 判断中部开放区域是否应该丢弃：如果与封闭型腔共享大部分侧壁面，则认为不需要单独处理
bool ShouldDiscardInteriorOpenRegion(const std::set<int>& openSideIds, const std::set<int>& closedSideIds, const TopoDS_Compound& topCapFaces, std::set<int>& sharedSideIds, double minSharedRatio = 0.8);

// 构建中部开放型腔的 CavityFeature：确定 Z 范围、收集面 ID、设置加工方向
CavityFeature BuildInteriorOpenCavityFeature(const std::vector<std::vector<Face2D>>& allOpenLayerFaces, const std::vector<std::vector<Face2D>>& allClosedLayerFaces, const std::vector<double>& splitPoints, double modelMinZ, const std::string& savePath);

// 根据 CavityFeature 的 sourceFaceIds 从原始模型中提取对应的 3D 面 Compound
TopoDS_Compound ExtractFacesByFeatureSourceIds(const CavityFeature& feature, const TopoDS_Shape& mainShape, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// 获取中部开放区域在 topZ 或 bottomZ 处的盖面（水平面）
TopoDS_Compound GetInteriorOpenBandCapFaces(const TopoDS_Shape& solid, const std::set<int>& sideFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, double topZ, double bottomZ, bool wantTopCap, double zTol = 0.1);

// 获取指定 Z 高度处与侧壁面相邻的水平盖面（用于确定型腔的顶部或底部封闭面）
TopoDS_Compound GetAdjacentHorizontalCapFacesAtZ(const TopoDS_Shape& solid, const std::set<int>& sideFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, double targetZ, bool normalUp, const std::string& debugLabel, double zTol = 0.1);

// 获取 Z 范围 [bandMinZ, bandMaxZ] 内与侧壁面相邻的所有水平面
TopoDS_Compound GetHorizontalFacesForBand(const TopoDS_Shape& solid, const std::set<int>& sideFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, double bandMinZ, double bandMaxZ, const std::string& debugLabel, double zTol = 0.1);

// 获取 Compound 的 XY 范围（最小/最大 X 和 Y）
void GetCompoundXYRange(const TopoDS_Compound& compound, double& xmin, double& ymin, double& xmax, double& ymax);

// 判断 OCC 面是否完全在指定 XY 矩形框内部
bool IsFaceInsideXYBox(const TopoDS_Face& face, double xmin, double ymin, double xmax, double ymax);

// 判断 OCC 面是否与指定 XY 矩形框交叉（部分在内）
bool IsFaceCrossingXYBox(const TopoDS_Face& face, double xmin, double ymin, double xmax, double ymax);

// 将面 Compound 按 XY 矩形框分为"完全在内"和"交叉"两组
void SplitFacesByXYBox(const TopoDS_Compound& inputFaces, double xmin, double ymin, double xmax, double ymax, TopoDS_Compound& insideFaces, TopoDS_Compound& crossingFaces);

// 从 Face2D 中收集所有线段的 faceId（去重后返回集合）
std::set<int> CollectFaceIdsFromFace2D(const Face2D& face);

// 从 Compound 中收集所有面的 faceId（通过 faceToIdMap 反查）
std::set<int> CollectFaceIdsFromCompound(const TopoDS_Compound& compound, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// 判断 3D 面是否与 2D Face2D 区域有 XY 重叠
bool IsFaceOverlappingFace2D(const TopoDS_Face& face3D, const Face2D& region2D);

// 按 2D Clipper 布尔区域将面 Compound 分为"可加工区域内的面"和"区域外的面"
void SplitFacesByClipperRegion(const TopoDS_Compound& inputFaces, const Face2D& machinableRegion, TopoDS_Compound& machinableFaces, TopoDS_Compound& nonMachinableFaces);

// 按 Z 范围将面 Compound 分为三部分：上部封闭区、中部开放区、下部封闭区
void SplitInteriorOpenSourceFacesByBand(const TopoDS_Compound& sourceFaces, double topZ, double bottomZ, TopoDS_Compound& upperClosed, TopoDS_Compound& middleOpen, TopoDS_Compound& lowerClosed);

// 中部开放型腔处理主流程：判断→分带→特征构建→上下分离→输出结果
void ProcessInteriorOpenCavityByBand(const std::vector<std::vector<Face2D>>& allClosedLayerFaces, const std::vector<std::vector<Face2D>>& allOpenLayerFaces, const TopoDS_Shape& mainShape, const std::vector<double>& splitPoints, const TopTools_DataMapOfShapeInteger& faceToIdMap, const std::string& savePath, std::vector<TopoDS_Compound>& outNewClosedParts, TopoDS_Compound& outNewOpenCavity, gp_Dir& outOpenToolDir);
