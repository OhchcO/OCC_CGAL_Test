#pragma once
#include "cavity_common.h"

// ======================== BREP 导出 ========================

// 将 OneLine 线段列表导出为 BREP 文件（用于调试查看切片线段）
void ExportOneLinesToBrep(const std::vector<OneLine>& lines, const std::string& fileName);

// 根据 faceId 集合从实体模型中提取型腔面，并导出为 BREP 文件
void ExportCavityFaces(const TopoDS_Shape& solid, const std::set<int>& cavityFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, const std::string& fileName);

// 将 Face2D 列表转换为 OCC Wire/Face 并导出为 BREP 文件（用于调试查看 2D 切片结果）
void ExportFace2DToBrep(const std::vector<Face2D>& faces, const std::string& fileName);

// 导出开放型腔面中的软边界线段（faceId == -1 的线段）为 BREP 文件
void ExportSoftEdgesToBrep(const std::vector<Face2D>& cavityFaces, const std::string& filePath);
