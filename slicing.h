#pragma once
#include "cavity_common.h"

// ======================== Z 轴相关 ========================

// 获取单个 OCC 面在 Z 方向的最小值和最大值
void GetFaceZRange(const TopoDS_Face& face, double& zmin, double& zmax);

// 获取整个 OCC 形状在 Z 方向的最小值和最大值
void GetShapeZRange(const TopoDS_Shape& shape, double& zmin, double& zmax);

// 沿指定方向获取面的 Z 极值（用于非水平面的极值分析）
void GetExtremaZOfFace(const TopoDS_Face& face, const gp_Dir& direction, double& outMinZ, double& outMaxZ);

// 自适应计算模型的水平切分点列表（基于水平面检测和面积阈值筛选）
std::vector<double> GetSplitPointsAlongZ(const TopoDS_Shape& shape, double angleTolerance = 0.001, double mergeTol = 0.05, double areaThreshold = 1.0);

// ======================== 切片 ========================

// 在指定 Z 高度对模型进行水平切片，返回该层的所有 Face2D（含 SOLID 和 CAVITY 类型）
std::vector<Face2D> SliceModelAtZ(const TopoDS_Shape& shape, double splitZ, int sliceIndex, int maxIndex, const TopTools_DataMapOfShapeInteger& faceToIdMap);
