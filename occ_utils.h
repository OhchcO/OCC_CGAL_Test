#pragma once
#include "cavity_common.h"

// ======================== 平面/边分析 ========================

// 获取模型在指定 Z 高度处的所有水平面（法向量接近垂直的面），返回 Compound
TopoDS_Compound GetCoplanarFaces(const TopoDS_Shape& shape, double splitZ, double tol = 1e-4);

// 判断两个相邻面之间的公共边是否为凹边（用于识别型腔边界）
bool IsConcave(const TopoDS_Edge& E, const TopoDS_Face& F1, const TopoDS_Face& F2);

// 获取模型中所有凹边，返回 Compound（基于共面面列表分析）
TopoDS_Compound GetConcaveEdges(const TopoDS_Shape& solid, const TopoDS_Compound& coplanarFaces, const gp_Dir& toolDirection = gp_Dir(0, 0, 1));

// 判断两个相邻面之间的公共边是否为凸边
bool IsConvex(const TopoDS_Edge& E, const TopoDS_Face& F1, const TopoDS_Face& F2);

// 获取模型中所有来自对面的凸边，返回 Compound
TopoDS_Compound GetConvexEdgesFromOppositeFaces(const TopoDS_Shape& solid, const TopoDS_Compound& coplanarFaces, const gp_Dir& toolDirection = gp_Dir(0, 0, 1));

// 分类凹凸边并返回需要移除的边 Compound（用于切片线段提取）
TopoDS_Compound ClassifyAndGetEdgesToRemove(const TopoDS_Shape& solid, const TopoDS_Compound& coplanarFaces, const gp_Dir& toolDirection);

// ======================== 线段布尔运算 ========================

// 对两组线段做布尔差运算，并维护每条线段到原始面的映射关系
std::vector<OneEdge> SubtractLinesAndMapFaces(const TopoDS_Compound& linesA, const TopoDS_Compound& linesB, const TopTools_DataMapOfShapeShape& edgeToFaceMap);

// 将 OneEdge 列表（含边+源面信息）构建为 OCC Compound
TopoDS_Compound BuildCompoundFromOneEdgeVector(const std::vector<OneEdge>& edges);

// ======================== 全局 Compound ========================

// 初始化全局 Compound 容器（gAllFacesCompound、gAllLinesCompound），只执行一次
void InitGlobalCompounds();

// ======================== 面属性判断 ========================

// 判断 OCC 面是否为水平面（法向量 Z 分量接近 1 或 -1）
bool IsHorizontalFace(const TopoDS_Face& face, double tol = 1e-4);

// 获取 OCC 面在中心点处的法向量
gp_Dir GetFaceNormal(const TopoDS_Face& face);

// 获取水平面的 Z 坐标高度（非水平面返回 0）
double GetHorizontalFaceZ(const TopoDS_Face& face);

// 计算 OCC 面的面积（使用 GProp_GProps）
double CalculateFaceAreaOCC(const TopoDS_Face& face);

// 统计 Compound 中包含的面数量
int CountFacesInCompound(const TopoDS_Compound& compound);

// 将任意方向的法向量 Snap 到最近的轴对齐方向（±X/±Y/±Z）
gp_Dir SnapNormalToDirection(const gp_Dir& normal);

// 判断型腔是否可加工：统计垂直面与水平面比例，确定进刀方向
bool IsCavityMachinable(const TopoDS_Compound& cavity, gp_Dir& outToolDir);

// ======================== Face ID 提取 ========================

// 从所有层的 Face2D 中收集所有出现过的 sourceFaceId（去重）
std::set<int> CollectCavityFaceIds(const std::vector<std::vector<Face2D>>& allLayerFaces);

// 根据 faceId 集合从原始实体模型中提取对应的拓扑面，返回 Compound
TopoDS_Compound GetFacesByFaceIds(const TopoDS_Shape& solid, const std::set<int>& targetFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// 获取型腔的盖面（cap faces）：与型腔壁面相邻的水平面，用于封闭型腔
TopoDS_Compound GetCavityCapFaces(const TopoDS_Shape& solid, const std::set<int>& cavityFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, double angleTolerance = 0.001);

// 获取型腔的完整 Compound：壁面 + 盖面
TopoDS_Compound GetCavityCompound(const TopoDS_Shape& solid, const std::set<int>& cavityFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// ======================== Compound 工具 ========================

// 将 source Compound 中的所有面添加到 target Compound 中
void AddCompoundFaces(TopoDS_Compound& target, const TopoDS_Compound& source);

// 获取 Compound 中所有面的 Z 范围（最小值和最大值）
bool GetCompoundZRange(const TopoDS_Compound& compound, double& zmin, double& zmax);
