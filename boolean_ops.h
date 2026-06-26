#pragma once
#include "cavity_common.h"

// ======================== Clipper2 路径转换 ========================

using namespace Clipper2Lib;

// 将内部 Loop 结构转换为 Clipper2 的 Path64 格式（x,y 乘以 SCALE 取整）
Path64 LoopToPath(const Loop& loop);

// 将 Clipper2 的 Path64 结果转换回内部 Loop 结构（除以 SCALE 还原浮点坐标）
Loop PathToLoop(const Path64& path, double zHeight);

// 将多个 Face2D 转换为 Clipper2 的 Paths64 格式（仅外环）
Paths64 FacesToClipperPaths(const std::vector<Face2D>& faces);

// 从 Clipper2 的 PolyPath64 树中递归提取 Face2D（区分孔洞和外轮廓）
void ExtractFacesFromPolyNode(const PolyPath64* subjectPoly, const PolyPath64* clipPoly, ClipType clipType, double zHeight, std::vector<Face2D>& resultFaces);

// ======================== 布尔运算 ========================

// 对单对 Face2D 执行布尔运算（交/差/并/异或），返回结果 Face2D 列表
std::vector<Face2D> BooleanFacesSingle(const Face2D& subject, const Face2D& clip, ClipType clipType, double zHeight);

// 对多对 Face2D 执行批量布尔运算（subjects 中每个元素依次与 clips 中对应元素运算）
std::vector<Face2D> BooleanFaces(const std::vector<Face2D>& subjects, const std::vector<Face2D>& clips, ClipType clipType, double zHeight);

// 通过几何包含关系恢复布尔运算结果中各子区域的原始 faceId
std::vector<Face2D> RecoverOriginalFaceIdsByGeometry(const std::vector<Face2D>& pocketResults, const std::vector<Face2D>& currentSolidFaces);
