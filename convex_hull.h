#pragma once
#include "cavity_common.h"

// ======================== 凸包相关 ========================

// 将 Face2D 的外环顶点转换为 CGAL Exact_predicates_exact_constructions_kernel 的 Point_2 列表
std::vector<Point_2> ConvertFaceToPoints2D_EK(const Face2D& face);

// 将 CGAL 凸包计算结果（Point_2 列表）转换回 Face2D 闭合环，Z 坐标设为 zHeight
Face2D ConvertPointsToHullFace_EK(vector<Point_2>& hull_points, double zHeight);

// 对单个 Face2D 计算其 2D 凸包，返回凸包 Face2D
Face2D ComputeConvexHullFace(const Face2D& face);

// 将多个小凸包（由 memberIndices 指定）合并为一个大凸包 Face2D，Z 坐标设为 zHeight
Face2D ComputeMergedHullFace(const std::vector<HullItem>& hullItems, const std::vector<int>& memberIndices, double zHeight);

// 收集多个 Face2D 外环上的所有顶点，用于合并凸包计算
std::vector<Point_2> CollectFaceOuterPoints(const std::vector<Face2D>& faces);

// 根据距离阈值将相近的小凸包分组合并，返回合并后的凸包组
std::vector<HullGroup> BuildMergedHullGroups(const std::vector<HullItem>& hullItems, const OpenCavityFilterParams& params);
