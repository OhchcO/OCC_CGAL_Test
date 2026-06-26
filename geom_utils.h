#pragma once
#include "cavity_common.h"

// ======================== 点/线段工具 ========================

// 判断两个 3D 点在容差范围内是否相等（只比较 x,y）
bool IsPointEqual(const Point3D& p1, const Point3D& p2, double tol = 1e-3);

// 计算向量 AB 与 AC 的 2D 叉积（Z 分量），用于判断方向/共线
double Cross2D(const Point3D& a, const Point3D& b, const Point3D& c);

// 判断 value 是否在 a 和 b 之间（含容差），不限定 a < b
bool IsValueBetween(double value, double a, double b, double tol = 1e-9);

// 判断点 p 是否在 2D 线段 ab 上（含容差）
bool IsPointOnSegment2D(const Point3D& p, const Point3D& a, const Point3D& b, double tol = 1e-9);

// 判断两条 2D 线段 ab 和 cd 是否相交（含端点接触）
bool SegmentsIntersect2D(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d);

// 计算 2D 点 p 到线段 ab 的最短距离
double PointSegmentDistance2D(const Point3D& p, const Point3D& a, const Point3D& b);

// 计算两条 2D 线段之间的最短距离
double SegmentSegmentDistance2D(const OneLine& a, const OneLine& b);

// 计算两个环之间的最短 2D 距离（逐线段对取最小值）
double LoopDistance2D(const Loop& a, const Loop& b);

// ======================== 环/面工具 ========================

// 使用鞋带公式计算 2D 环的有向面积（正值=逆时针）
double CalculateArea(const Loop& loop);

// 射线法判断 2D 点是否在环内部
bool IsPointInLoop(const Point3D& pt, const Loop& loop);

// 计算 Face2D 外环的几何中心点（各顶点均值）
Point3D GetFace2DCenter(const Face2D& face);

// 判断 innerFace 是否完全在 outerFace 内部（基于顶点包含检测）
bool IsFaceInsideFace(const Face2D& innerFace, const Face2D& outerFace);

// 计算环的 2D 轴对齐包围盒
LoopBBox2D GetLoopBBox(const Loop& loop);

// 计算两个包围盒之间的最短距离（0=重叠）
double BBoxDistance2D(const LoopBBox2D& a, const LoopBBox2D& b);

// ======================== Z 平移 ========================

// 将 Face2D 所有顶点的 Z 坐标平移到 newZ，返回新 Face2D
Face2D TranslateFace2DZ(const Face2D& face, double newZ);
