#include "open_cavity.h"
#include "geom_utils.h"
#include "occ_utils.h"
#include "export_utils.h"
#include "boolean_ops.h"
#include "slicing.h"
#include "closed_cavity.h"

// --- 开放型腔清洗 ---
double CalculateFaceCompactness(const Face2D& face)
{
    if (face.outerLoop.empty()) return 0.0;

    // 1. 利用鞋带公式计算面积（借用你原有的 CalculateArea）
    double area = CalculateArea(face.outerLoop);
    if (area < 1e-5) return 0.0;

    // 2. 累加外环总周长
    double perimeter = 0.0;
    for (const auto& line : face.outerLoop) {
        gp_Pnt p1(line.start.x, line.start.y, line.start.z);
        gp_Pnt p2(line.end.x, line.end.y, line.end.z);
        perimeter += p1.Distance(p2);
    }

    if (perimeter < 1e-5) return 0.0;

    // 3. 刚性等周商数学定义
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
    double compactness = (4.0 * M_PI * area) / (perimeter * perimeter);

    return compactness;
}

// (IsPureFilletSliverCavity_Advanced 已被 BuildFilletRadiusMap + IsSmallFilletCavity 替代)
std::map<int, double> BuildFilletRadiusMap(
    const TopoDS_Shape& shape,
    const TopTools_DataMapOfShapeInteger& faceToIdMap)
{
    std::map<int, double> filletMap;
    TopExp_Explorer exp(shape, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        if (!faceToIdMap.IsBound(face)) continue;
        int faceId = faceToIdMap.Find(face);

        BRepAdaptor_Surface surf(face);
        double radius = -1.0;

        if (surf.GetType() == GeomAbs_Cylinder) {
            radius = surf.Cylinder().Radius();
        }
        else if (surf.GetType() == GeomAbs_Torus) {
            radius = surf.Torus().MinorRadius();
        }

        filletMap[faceId] = radius;
    }
    return filletMap;
}


bool IsSmallFilletCavity(
    const Face2D& face,
    const std::map<int, double>& filletMap,
    double maxFilletRadius,
    std::vector<std::pair<int, double>>* matchedFillets)
{
    const auto& loop = face.outerLoop;
    if (loop.empty()) return false;

    bool hasVirtualEdge = false;
    for (const auto& line : loop) {
        if (line.faceId == -1) {
            hasVirtualEdge = true;
            break;
        }
    }
    if (!hasVirtualEdge) return false;

    std::set<int> edgeNeighborIds;
    for (size_t i = 0; i < loop.size(); i++) {
        if (loop[i].faceId != -1) continue;

        size_t prevIdx = (i == 0) ? loop.size() - 1 : i - 1;
        size_t nextIdx = (i + 1) % loop.size();

        if (loop[prevIdx].faceId >= 0) {
            edgeNeighborIds.insert(loop[prevIdx].faceId);
        }
        if (loop[nextIdx].faceId >= 0) {
            edgeNeighborIds.insert(loop[nextIdx].faceId);
        }
    }

    if (edgeNeighborIds.empty()) return false;

    std::vector<std::pair<int, double>> localMatchedFillets;
    for (int fid : edgeNeighborIds) {
        auto it = filletMap.find(fid);
        if (it == filletMap.end() || it->second <= 0 || it->second >= maxFilletRadius) {
            return false;
        }
        localMatchedFillets.push_back({ fid, it->second });
    }

    if (matchedFillets) {
        *matchedFillets = localMatchedFillets;
    }
    return true;
}


std::vector<Face2D> CleanAndFilterOpenCavities(
    const std::vector<Face2D>& rawOpenCavities,
    const std::map<int, double>& filletMap,
    const OpenCavityFilterParams& params)
{
    std::vector<Face2D> clearFeatures;

    int idx = 0;
    for (const auto& face : rawOpenCavities) {
        idx++;
        if (face.outerLoop.size() < 3) {
            std::cout << "  ❌ 型腔#" << idx << " 被过滤: 顶点数不足 ("
                << face.outerLoop.size() << " < 3)" << std::endl;
            continue;
        }

        double area = CalculateArea(face.outerLoop);
        if (area < params.minArea) {
            std::cout << "  ❌ 型腔#" << idx << " 被过滤: 面积过小 ("
                << area << " < " << params.minArea << ")" << std::endl;
            continue;
        }

        double compactness = CalculateFaceCompactness(face);
        if (compactness < params.minCompactness) {
            std::cout << "  ❌ 型腔#" << idx << " 被过滤: 等周商过小/细长 ("
                << compactness << " < " << params.minCompactness << ")" << std::endl;
            continue;
        }

        double xMin = 1e9, xMax = -1e9, yMin = 1e9, yMax = -1e9;
        for (const auto& line : face.outerLoop) {
            xMin = std::min(xMin, line.start.x); xMax = std::max(xMax, line.start.x);
            yMin = std::min(yMin, line.start.y); yMax = std::max(yMax, line.start.y);
        }
        if ((xMax - xMin) < params.minToolPassSpan && (yMax - yMin) < params.minToolPassSpan) {
            std::cout << "  ❌ 型腔#" << idx << " 被过滤: 包围盒跨度不足 ("
                << (xMax - xMin) << " x " << (yMax - yMin) << " < " << params.minToolPassSpan << ")" << std::endl;
            continue;
        }

        std::vector<std::pair<int, double>> matchedFillets;
        if (IsSmallFilletCavity(face, filletMap, params.maxFilletRadius, &matchedFillets)) {
            std::cout << "  ❌ 型腔#" << idx << " 被过滤: 小圆角特征 (半径 < " << params.maxFilletRadius << ")" << std::endl;
            std::cout << "      圆角面: ";
            for (const auto& [fid, radius] : matchedFillets) {
                std::cout << "faceId=" << fid << " radius=" << radius << " ";
            }
            std::cout << std::endl;
            continue;
        }

        std::cout << "  ✅ 型腔#" << idx << " 通过: 面积=" << area
            << " 等周商=" << compactness
            << " 跨度=" << (xMax - xMin) << "x" << (yMax - yMin) << std::endl;
        clearFeatures.push_back(face);
    }

    std::cout << "🏁 [工艺提纯报告] 原本共有型腔碎片: " << rawOpenCavities.size()
        << " 个 -> 过滤后最终存活型腔: " << clearFeatures.size() << " 个！" << std::endl;

    return clearFeatures;
}

//开放型腔提取流程
// 1. 【开放专用】提取属于开放型腔的表面 ID 集合

// --- 开放型腔提取 ---
std::set<int> CollectOpenCavityFaceIds(const std::vector<std::vector<Face2D>>& allLayerFaces) {
    std::set<int> cavityFaceIds;
    for (const auto& layerFaces : allLayerFaces) {
        for (const auto& face : layerFaces) {
            if (face.type != FaceType::OPENCAVITY) continue; // 🔓 严格限制为 OPENCAVITY
            for (const auto& line : face.outerLoop) {
                if (line.faceId >= 0) cavityFaceIds.insert(line.faceId);
            }
            for (const auto& innerLoop : face.innerLoops) {
                for (const auto& line : innerLoop) {
                    if (line.faceId >= 0) cavityFaceIds.insert(line.faceId);
                }
            }
        }
    }
    return cavityFaceIds;
}

// 2. 【开放专用】获取开放型腔的底面/开放槽平台面
// 2. 【开放专用】获取开放型腔的工艺底面（自适应有无底面，彻底杜绝顶面）
TopoDS_Compound GetOpenCavityCapFaces(const TopoDS_Shape& solid,
    const std::set<int>& openFaceIds,
    const TopTools_DataMapOfShapeInteger& faceToIdMap,
    double angleTolerance) {
    BRep_Builder builder;
    TopoDS_Compound capFaces;
    builder.MakeCompound(capFaces);

    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    double wallMinZ = 1e9;
    double wallMaxZ = -1e9;
    bool hasValidWalls = false;

    TopExp_Explorer wallExp(solid, TopAbs_FACE);
    for (; wallExp.More(); wallExp.Next()) {
        TopoDS_Face f = TopoDS::Face(wallExp.Current());
        if (faceToIdMap.IsBound(f)) {
            int id = faceToIdMap.Find(f);
            if (openFaceIds.count(id)) {
                double zmin, zmax;
                GetFaceZRange(f, zmin, zmax);
                wallMinZ = std::min(wallMinZ, zmin);
                wallMaxZ = std::max(wallMaxZ, zmax);
                hasValidWalls = true;
            }
        }
    }

    if (!hasValidWalls) return capFaces;

    const gp_Dir zAxis(0, 0, 1);
    const double zTol = 0.1;

    TopExp_Explorer exp(solid, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        BRepAdaptor_Surface surf(face);
        if (surf.GetType() != GeomAbs_Plane) continue;

        gp_Pln plane = surf.Plane();
        gp_Dir normal = plane.Axis().Direction();
        if (face.Orientation() == TopAbs_REVERSED) {
            normal.Reverse();
        }

        // 过滤1：只保留水平面（法线平行于Z轴）
        double angle = normal.Angle(zAxis);
        if (angle > angleTolerance && fabs(angle - M_PI) > angleTolerance) continue;

        int faceId = -1;
        if (faceToIdMap.IsBound(face)) faceId = faceToIdMap.Find(face);
        else continue;

        double faceZ = plane.Location().Z();

        // 过滤2：顶面（wallMaxZ附近）只保留法向朝下的面
        if (std::abs(faceZ - wallMaxZ) < zTol && normal.Z() > -0.99) {
            continue;
        }

        // 过滤3：底面（wallMinZ附近）只保留法向朝上的面
        if (std::abs(faceZ - wallMinZ) < zTol && normal.Z() < 0.99) {
            continue;
        }

        // 过滤4：邻接面必须属于型腔面集合
        bool hasOpenCavityNeighbor = false;
        TopExp_Explorer edgeExp(face, TopAbs_EDGE);
        for (; edgeExp.More(); edgeExp.Next()) {
            const TopTools_ListOfShape& adjacentFaces = edgeToFaces.FindFromKey(edgeExp.Current());
            TopTools_ListIteratorOfListOfShape it(adjacentFaces);
            for (; it.More(); it.Next()) {
                const TopoDS_Face& adjFace = TopoDS::Face(it.Value());
                if (adjFace.IsSame(face)) continue;

                if (faceToIdMap.IsBound(adjFace)) {
                    int adjId = faceToIdMap.Find(adjFace);
                    if (openFaceIds.count(adjId)) {
                        hasOpenCavityNeighbor = true;
                        break;
                    }
                }
            }
            if (hasOpenCavityNeighbor) break;
        }

        // 通过所有过滤的面加入结果集
        if (hasOpenCavityNeighbor && faceId >= 0) {
            builder.Add(capFaces, face);
        }
    }

    return capFaces;
}
// 3. 【开放专用】导出并保存完整开放型腔面
void ExportOpenCavityFaces(const TopoDS_Shape& solid,
    const std::set<int>& openFaceIds,
    const TopTools_DataMapOfShapeInteger& faceToIdMap,
    const std::string& fileName) {
    TopoDS_Compound wallFaces = GetFacesByFaceIds(solid, openFaceIds, faceToIdMap);
    TopoDS_Compound capFaces = GetOpenCavityCapFaces(solid, openFaceIds, faceToIdMap);

    BRep_Builder builder;
    TopoDS_Compound result;
    builder.MakeCompound(result);

    TopExp_Explorer exp(wallFaces, TopAbs_FACE);
    for (; exp.More(); exp.Next()) builder.Add(result, exp.Current());
    TopExp_Explorer expCap(capFaces, TopAbs_FACE);
    for (; expCap.More(); expCap.Next()) builder.Add(result, expCap.Current());

    BRepTools::Write(result, fileName.c_str());
    std::cout << "  🔓 完整开放型腔已保存至: " << fileName << std::endl;
}

// 4. 【开放专用】获取整个开放型腔的合并面 Compound
TopoDS_Compound GetOpenCavityCompound(const TopoDS_Shape& solid,
    const std::set<int>& openFaceIds,
    const TopTools_DataMapOfShapeInteger& faceToIdMap) {
    std::cout << "\n===== [DEBUG] GetOpenCavityCompound 开始 =====" << std::endl;
    std::cout << "  输入 openFaceIds 共 " << openFaceIds.size() << " 个: ";
    for (int fid : openFaceIds) std::cout << fid << " ";
    std::cout << std::endl;

    TopoDS_Compound wallFaces = GetFacesByFaceIds(solid, openFaceIds, faceToIdMap);

    int wallCount = 0;
    TopExp_Explorer wallDebug(wallFaces, TopAbs_FACE);
    for (; wallDebug.More(); wallDebug.Next()) {
        TopoDS_Face f = TopoDS::Face(wallDebug.Current());
        double zmin, zmax;
        GetFaceZRange(f, zmin, zmax);
        int fid = faceToIdMap.IsBound(f) ? faceToIdMap.Find(f) : -1;
        std::cout << "  [侧壁面] faceId=" << fid << "  Zrange=[" << zmin << ", " << zmax << "]" << std::endl;
        wallCount++;
    }
    std::cout << "  侧壁面总数: " << wallCount << std::endl;

    TopoDS_Compound capFaces = GetOpenCavityCapFaces(solid, openFaceIds, faceToIdMap);

    int capCount = 0;
    TopExp_Explorer capDebug(capFaces, TopAbs_FACE);
    for (; capDebug.More(); capDebug.Next()) {
        TopoDS_Face f = TopoDS::Face(capDebug.Current());
        double zmin, zmax;
        GetFaceZRange(f, zmin, zmax);
        int fid = faceToIdMap.IsBound(f) ? faceToIdMap.Find(f) : -1;
        std::cout << "  [盖板面] faceId=" << fid << "  Zrange=[" << zmin << ", " << zmax << "]" << std::endl;
        capCount++;
    }
    std::cout << "  盖板面总数: " << capCount << std::endl;

    BRep_Builder builder;
    TopoDS_Compound result;
    builder.MakeCompound(result);

    TopExp_Explorer exp(wallFaces, TopAbs_FACE);
    for (; exp.More(); exp.Next()) builder.Add(result, exp.Current());
    TopExp_Explorer expCap(capFaces, TopAbs_FACE);
    for (; expCap.More(); expCap.Next()) builder.Add(result, expCap.Current());

    std::cout << "  最终 Compound 总面数: " << wallCount + capCount << std::endl;
    std::cout << "===== [DEBUG] GetOpenCavityCompound 结束 =====\n" << std::endl;
    return result;
}

/**
 * @brief 【开放专用】将物理分离的开放型腔实体与所有切片的 2D 环进行匹配，重组为开放 CavityFeature 列表
 */
std::vector<CavityFeature> GenerateOpenCavityFeatures(
    const std::vector<TopoDS_Compound>& isolatedCavities,
    const std::vector<std::vector<Face2D>>& allLayerFaces,
    const std::vector<double>& sortedSplitPoints,
    const TopTools_DataMapOfShapeInteger& faceToIdMap)
{
    std::cout << "\n===== [DEBUG] GenerateOpenCavityFeatures 开始 =====" << std::endl;
    std::cout << "  isolatedCavities 数量: " << isolatedCavities.size() << std::endl;
    std::cout << "  allLayerFaces 层数: " << allLayerFaces.size() << std::endl;
    std::cout << "  sortedSplitPoints (" << sortedSplitPoints.size() << "个): ";
    for (size_t s = 0; s < sortedSplitPoints.size(); ++s) std::cout << sortedSplitPoints[s] << " ";
    std::cout << std::endl;

    std::vector<CavityFeature> finalFeatures;

    for (size_t i = 0; i < isolatedCavities.size(); ++i) {
        std::cout << "\n  --- [DEBUG] 处理 isolatedCavity [" << i << "] ---" << std::endl;
        CavityFeature feat;
        feat.featureId = (int)i;
        feat.type = CavityType::OPEN;
        feat.topZ = -1e9;
        feat.bottomZ = 1e9;
        double lastLoopZ = 1e9;

        TopExp_Explorer exp(isolatedCavities[i], TopAbs_FACE);
        int faceCountInCompound = 0;
        for (; exp.More(); exp.Next()) {
            if (faceToIdMap.IsBound(exp.Current())) {
                feat.sourceFaceIds.insert(faceToIdMap.Find(exp.Current()));
            }
            faceCountInCompound++;
        }
        std::cout << "  Compound 内总面数: " << faceCountInCompound << std::endl;
        std::cout << "  sourceFaceIds (" << feat.sourceFaceIds.size() << "个): ";
        for (int fid : feat.sourceFaceIds) std::cout << fid << " ";
        std::cout << std::endl;

        Bnd_Box box;
        BRepBndLib::Add(isolatedCavities[i], box);
        double x1, y1, zmin_box, x2, y2, zmax_box;
        box.Get(x1, y1, zmin_box, x2, y2, zmax_box);
        std::cout << "  3D包围盒: X=[" << x1 << ", " << x2 << "] Y=[" << y1 << ", " << y2 << "] Z=[" << zmin_box << ", " << zmax_box << "]" << std::endl;

        int layerIdx = 0;
        int skippedNonOpen = 0;
        int skippedIdMismatch = 0;
        int skippedBBoxMismatch = 0;
        int matchedLoops = 0;

        for (const auto& layer : allLayerFaces) {
            int faceIdxInLayer = 0;
            for (const auto& f2d : layer) {
                if (f2d.type != FaceType::OPENCAVITY || f2d.outerLoop.empty()) {
                    skippedNonOpen++;
                    faceIdxInLayer++;
                    continue;
                }

                // 投票法统计外环所有线段的 faceId，取出现频率最高的
                std::map<int, int> idVoteCount;
                for (const auto& l : f2d.outerLoop) {
                    if (l.faceId >= 0) {
                        idVoteCount[l.faceId]++;
                    }
                }
                int dominantFaceId = -1;
                int maxVotes = 0;
                for (const auto& [fid, cnt] : idVoteCount) {
                    if (cnt > maxVotes) {
                        maxVotes = cnt;
                        dominantFaceId = fid;
                    }
                }

                // ID 匹配检查
                if (dominantFaceId < 0 || !feat.sourceFaceIds.count(dominantFaceId)) {
                    if (dominantFaceId >= 0) {
                        std::cout << "    [层" << layerIdx << " 面" << faceIdxInLayer << "] OPENCAVITY dominantFaceId=" << dominantFaceId
                                  << " (投票=" << maxVotes << ") 不在 sourceFaceIds 中 → 跳过" << std::endl;
                    }
                    skippedIdMismatch++;
                    faceIdxInLayer++;
                    continue;
                }

                Point3D p = f2d.outerLoop[0].start;
                bool inBBox = (p.x >= x1 - 0.5 && p.x <= x2 + 0.5 && p.y >= y1 - 0.5 && p.y <= y2 + 0.5);
                if (!inBBox) {
                    std::cout << "    [层" << layerIdx << " 面" << faceIdxInLayer << "] OPENCAVITY dominantFaceId=" << dominantFaceId
                              << " ID匹配OK 但 XY不在包围盒内! 首点=(" << p.x << ", " << p.y << ", " << p.z << ") → 跳过" << std::endl;
                    skippedBBoxMismatch++;
                    faceIdxInLayer++;
                    continue;
                }

                double curZ = p.z;
                int innerCount = (int)f2d.innerLoops.size();

                CavityLoop out;
                out.lines = f2d.outerLoop; out.zHeight = curZ; out.isOuter = true;
                feat.stepLoops.push_back(out);

                for (const auto& inner : f2d.innerLoops) {
                    CavityLoop in;
                    in.lines = inner; in.zHeight = curZ; in.isOuter = false;
                    feat.stepLoops.push_back(in);
                }

                feat.topZ = std::max(feat.topZ, curZ);
                lastLoopZ = std::min(lastLoopZ, curZ);
                matchedLoops++;

                std::cout << "    [层" << layerIdx << " 面" << faceIdxInLayer << "] ✅ 匹配! dominantFaceId=" << dominantFaceId
                          << " curZ=" << curZ << " innerLoops=" << innerCount
                          << " outerLoop段数=" << f2d.outerLoop.size() << std::endl;

                faceIdxInLayer++;
            }
            layerIdx++;
        }

        std::cout << "  [DEBUG] 匹配统计: matchedLoops=" << matchedLoops
                  << " skippedNonOpen=" << skippedNonOpen
                  << " skippedIdMismatch=" << skippedIdMismatch
                  << " skippedBBoxMismatch=" << skippedBBoxMismatch << std::endl;
        std::cout << "  [DEBUG] 最终 stepLoops 数: " << feat.stepLoops.size()
                  << " topZ=" << feat.topZ << " lastLoopZ(最低环)=" << lastLoopZ << std::endl;

        if (!feat.stepLoops.empty()) {
            feat.bottomZ = lastLoopZ;
            bool foundSplitPoint = false;
            for (size_t k = 0; k < sortedSplitPoints.size(); ++k) {
                if (std::abs(sortedSplitPoints[k] - lastLoopZ) < 1e-4) {
                    foundSplitPoint = true;
                    if (k + 1 < sortedSplitPoints.size()) {
                        std::cout << "  [DEBUG] bottomZ修正: lastLoopZ=" << lastLoopZ
                                  << " 匹配到 splitPoints[" << k << "]=" << sortedSplitPoints[k]
                                  << " → 取下一个 splitPoints[" << k + 1 << "]=" << sortedSplitPoints[k + 1]
                                  << " 作为 bottomZ" << std::endl;
                        feat.bottomZ = sortedSplitPoints[k + 1];
                    } else {
                        std::cout << "  [DEBUG] bottomZ修正: lastLoopZ=" << lastLoopZ
                                  << " 匹配到 splitPoints[" << k << "]=" << sortedSplitPoints[k]
                                  << " 但已是最后一个点，无法再向下取 → bottomZ保持lastLoopZ=" << lastLoopZ << std::endl;
                    }
                    break;
                }
            }
            if (!foundSplitPoint) {
                std::cout << "  [DEBUG] ⚠️ bottomZ修正: lastLoopZ=" << lastLoopZ
                          << " 在 sortedSplitPoints 中没有找到匹配项(容差1e-4)! bottomZ保持lastLoopZ=" << lastLoopZ << std::endl;
            }
            feat.totalDepth = std::abs(feat.topZ - feat.bottomZ);
            std::sort(feat.stepLoops.begin(), feat.stepLoops.end(),
                [](const CavityLoop& a, const CavityLoop& b) { return a.zHeight > b.zHeight; });
            finalFeatures.push_back(feat);
        }
    }

    std::cout << "\n===== [DEBUG] GenerateOpenCavityFeatures 结果汇总 =====" << std::endl;
    for (const auto& feat : finalFeatures) {
        cout << ">>> [开放型腔 ID: " << feat.featureId << "]" << endl;
        cout << "    - 类型: " << (feat.type == CavityType::OPEN ? "OPEN" : feat.type == CavityType::CLOSED ? "CLOSED" : "OTHER") << endl;
        cout << "    - 进刀方向: (" << feat.toolDirection.X() << ", " << feat.toolDirection.Y() << ", " << feat.toolDirection.Z() << ")" << endl;
        cout << "    - 空间范围: Z_Top = " << feat.topZ << " / Z_Floor = " << feat.bottomZ << endl;
        cout << "    - 开粗深度: " << feat.totalDepth << " mm" << endl;
        cout << "    - 进刀方案: 外部侧向安全切入 (Side Entry)" << endl;

        cout << "    - 关联原始面 ID (" << feat.sourceFaceIds.size() << "个): ";
        for (int oid : feat.sourceFaceIds) cout << oid << " ";
        cout << endl;

        cout << "    - stepLoops 详情 (" << feat.stepLoops.size() << "个):" << endl;
        for (size_t sl = 0; sl < feat.stepLoops.size(); ++sl) {
            const auto& loop = feat.stepLoops[sl];
            cout << "      [" << sl << "] isOuter=" << (loop.isOuter ? "外环" : "内环")
                 << " zHeight=" << loop.zHeight
                 << " 线段数=" << loop.lines.size();
            if (!loop.lines.empty()) {
                Point3D s = loop.lines.front().start;
                Point3D e = loop.lines.back().end;
                cout << " 首点=(" << s.x << "," << s.y << "," << s.z << ")"
                     << " 末点=(" << e.x << "," << e.y << "," << e.z << ")";
            }
            cout << endl;
        }
    }
    std::cout << "===== [DEBUG] GenerateOpenCavityFeatures 结束 =====\n" << std::endl;

    return finalFeatures;
}


#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>

TopoDS_Compound ExtractTrueOpenCavityFaces(
    const CavityFeature& openFeat,
    const TopoDS_Shape& mainShape,
    const TopTools_DataMapOfShapeInteger& faceToIdMap,
    const std::vector<std::vector<Face2D>>& allClosedLayerFaces,
    const std::vector<double>& splitPoints,
    const std::vector<TopoDS_Compound>& trueClosedCavities)
{
    BRep_Builder builder;
    TopoDS_Compound refinedResult;
    builder.MakeCompound(refinedResult);

    // =========================================================================
    // 🧬 1. 建立确凿证据库：提取所有已知封闭型腔面的几何指纹（面积和重心）
    // =========================================================================
    struct FaceSignature {
        double area;
        gp_Pnt centroid;
    };
    std::vector<FaceSignature> closedSignatures;

    for (const auto& closedComp : trueClosedCavities) {
        TopExp_Explorer fExp(closedComp, TopAbs_FACE);
        for (; fExp.More(); fExp.Next()) {
            TopoDS_Face f = TopoDS::Face(fExp.Current());
            GProp_GProps props;
            BRepGProp::SurfaceProperties(f, props);
            if (props.Mass() > 1e-6) {
                closedSignatures.push_back({ props.Mass(), props.CentreOfMass() });
            }
        }
    }

    auto IsClosedCavityFace = [&](const TopoDS_Face& face) -> bool {
        GProp_GProps props;
        BRepGProp::SurfaceProperties(face, props);
        double area = props.Mass();
        if (area <= 1e-6) return false;
        gp_Pnt center = props.CentreOfMass();

        for (const auto& sig : closedSignatures) {
            // 容差设定：面积误差 < 1e-3，重心偏移 < 1e-3
            if (std::abs(sig.area - area) < 1e-3 && sig.centroid.Distance(center) < 1e-3) {
                return true;
            }
        }
        return false;
    };

    // =========================================================================
    // 📡 2. 雷达侦测：定位分水岭高度 cutZ（只为了切开跨越的侧壁）
    // =========================================================================
    double cutZ = openFeat.topZ;
    TopoDS_Wire closedTopWire;
    std::set<int> closedFaceIds = CollectCavityFaceIds(allClosedLayerFaces);

    if (!closedFaceIds.empty()) {
        TopoDS_Compound closedCompound = GetCavityCompound(mainShape, closedFaceIds, faceToIdMap);
        std::vector<TopoDS_Compound> isolatedClosedParts = SeparateDisconnectedCavities(closedCompound);
        std::vector<CavityFeature> allClosedFeatures = GenerateCavityFeatures(isolatedClosedParts, allClosedLayerFaces, splitPoints, faceToIdMap);

        double openXMin = 1e9, openXMax = -1e9, openYMin = 1e9, openYMax = -1e9;
        for (const auto& loop : openFeat.stepLoops) {
            for (const auto& line : loop.lines) {
                openXMin = std::min({ openXMin, line.start.x, line.end.x });
                openXMax = std::max({ openXMax, line.start.x, line.end.x });
                openYMin = std::min({ openYMin, line.start.y, line.end.y });
                openYMax = std::max({ openYMax, line.start.y, line.end.y });
            }
        }

        for (const auto& closedFeat : allClosedFeatures) {
            if (closedFeat.stepLoops.empty()) continue;
            Point3D closedTopPt = closedFeat.stepLoops.front().lines.front().start;

            if (closedTopPt.x >= openXMin - 0.5 && closedTopPt.x <= openXMax + 0.5 &&
                closedTopPt.y >= openYMin - 0.5 && closedTopPt.y <= openYMax + 0.5 &&
                closedFeat.topZ <= openFeat.topZ + 0.1 && closedFeat.topZ >= openFeat.bottomZ - 0.1)
            {
                cutZ = closedFeat.topZ;

                const auto& tLoop = closedFeat.stepLoops.front();
                BRepBuilderAPI_MakePolygon poly;
                for (const auto& l : tLoop.lines) poly.Add(gp_Pnt(l.start.x, l.start.y, l.start.z));
                if (!tLoop.lines.empty()) {
                    poly.Add(gp_Pnt(tLoop.lines.back().end.x, tLoop.lines.back().end.y, tLoop.lines.back().end.z));
                }
                if (poly.IsDone()) closedTopWire = poly.Wire();

                std::cout << "    [智能联动] 锁定内部嵌套封闭孔: 顶切分点 Z = " << cutZ << std::endl;
                break;
            }
        }
    }

    // =========================================================================
    // ⚔️ 3. 物理切割并使用指纹库比对剔除
    // =========================================================================
    gp_Pln cuttingPlane(gp_Pnt(0, 0, cutZ), gp_Dir(0, 0, 1));
    TopoDS_Face splitPlaneFace = BRepBuilderAPI_MakeFace(cuttingPlane);

    TopExp_Explorer exp(mainShape, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face originalFace = TopoDS::Face(exp.Current());
        if (!faceToIdMap.IsBound(originalFace)) continue;
        int oid = faceToIdMap.Find(originalFace);
        if (!openFeat.sourceFaceIds.count(oid)) continue;

        double zmin, zmax;
        GetFaceZRange(originalFace, zmin, zmax);

        // 如果面完全在 cutZ 上方
        if (zmin >= cutZ - 1e-4) {
            if (!IsClosedCavityFace(originalFace)) {
                builder.Add(refinedResult, originalFace);
            }
            continue;
        }

        // 如果面完全在 cutZ 下方（例如底面，如面 23）
        if (zmax <= cutZ + 1e-4) {
            if (!IsClosedCavityFace(originalFace)) {
                builder.Add(refinedResult, originalFace);
            }
            continue;
        }

        // 跨越 cutZ 的面，执行切割
        if (zmax > cutZ + 1e-4 && zmin < cutZ - 1e-4) {
            BRepAlgoAPI_Section section(originalFace, splitPlaneFace);
            section.Build();
            if (section.IsDone()) {
                BRepFeat_SplitShape splitter(originalFace);
                TopExp_Explorer eExp(section.Shape(), TopAbs_EDGE);
                for (; eExp.More(); eExp.Next())
                    splitter.Add(TopoDS::Edge(eExp.Current()), originalFace);
                splitter.Build();
                if (splitter.IsDone()) {
                    TopExp_Explorer subExp(splitter.Shape(), TopAbs_FACE);
                    for (; subExp.More(); subExp.Next()) {
                        TopoDS_Face subF = TopoDS::Face(subExp.Current());
                        if (!IsClosedCavityFace(subF)) {
                            builder.Add(refinedResult, subF);
                        }
                    }
                }
            }
        }
    }

    return refinedResult;
}


#include <BRepIntCurveSurface_Inter.hxx>
#include <gp_Lin.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
// 辅助函数：全向半山腰射线法判断 innerComp（底面/孤岛）是否在 outerComp（侧壁）内部

// --- 嵌套合并 ---
bool IsCompoundInsideExact(const TopoDS_Compound& innerComp, const TopoDS_Compound& outerComp) {

    // 1. 测算侧壁 (outerComp) 的绝对高度范围
    double outZmin = 1e9, outZmax = -1e9;
    TopExp_Explorer outZExp(outerComp, TopAbs_FACE);
    bool hasWall = false;
    for (; outZExp.More(); outZExp.Next()) {
        double z1, z2;
        GetFaceZRange(TopoDS::Face(outZExp.Current()), z1, z2);
        outZmin = std::min(outZmin, z1);
        outZmax = std::max(outZmax, z2);
        hasWall = true;
    }
    // 如果没有找到有效侧壁，直接返回 false
    if (!hasWall || outZmin > 1e8) return false;

    // 🎯 核心神技 1：算得侧壁的“半山腰”高度
    double safeMidZ = (outZmax + outZmin) / 2.0;

    // 2. 提取底面 (innerComp) 
    TopExp_Explorer fExp(innerComp, TopAbs_FACE);
    if (!fExp.More()) return false;
    TopoDS_Face firstFace = TopoDS::Face(fExp.Current());
    BRepAdaptor_Surface surf(firstFace);

    // 🎯 核心神技 2：使用 2D 分类器，确保提取的测试点 100% 踩在 U 型底面的真实材料上，而不是虚空里！
    BRepTopAdaptor_FClass2d classifier(firstFace, Precision::PConfusion());
    gp_Pnt P_test;
    bool foundValidPoint = false;

    // 在面上撒一张 10x10 的网，总能抓到一个真正内部的点
    for (int i = 1; i < 10; ++i) {
        for (int j = 1; j < 10; ++j) {
            Standard_Real u = surf.FirstUParameter() + i * (surf.LastUParameter() - surf.FirstUParameter()) / 10.0;
            Standard_Real v = surf.FirstVParameter() + j * (surf.LastVParameter() - surf.FirstVParameter()) / 10.0;

            if (classifier.Perform(gp_Pnt2d(u, v)) == TopAbs_IN) {
                P_test = surf.Value(u, v);
                foundValidPoint = true;
                break;
            }
        }
        if (foundValidPoint) break;
    }
    if (!foundValidPoint) return false;

    // 🎯 核心神技 3：把测试点直接拔高到侧壁的半山腰！
    // 彻底免疫底部圆角被删后留下的巨大物理缝隙！
    P_test.SetZ(safeMidZ);

    // 3. 全向雷达发射 (向 4 个方向射击)
    gp_Dir dirs[4] = { gp_Dir(1,0,0), gp_Dir(-1,0,0), gp_Dir(0,1,0), gp_Dir(0,-1,0) };
    int oddCount = 0;

    for (int d = 0; d < 4; ++d) {
        gp_Lin ray(P_test, dirs[d]);
        std::vector<gp_Pnt> hitPoints;

        TopExp_Explorer outerFaceExp(outerComp, TopAbs_FACE);
        for (; outerFaceExp.More(); outerFaceExp.Next()) {
            TopoDS_Face outFace = TopoDS::Face(outerFaceExp.Current());

            // 求交计算
            BRepIntCurveSurface_Inter inter;
            inter.Init(outFace, ray, 1e-4);

            while (inter.More()) {
                if (inter.W() > 1e-4) { // 只看正前方的交点
                    gp_Pnt hitP = inter.Pnt();

                    // 去重：防止刚好打在两面墙交界线上算成 2 次
                    bool isDuplicate = false;
                    for (const auto& p : hitPoints) {
                        if (p.Distance(hitP) < 1e-2) {
                            isDuplicate = true; break;
                        }
                    }
                    if (!isDuplicate) hitPoints.push_back(hitP);
                }
                inter.Next();
            }
        }

        // 如果在这个方向上打穿侧壁的次数为奇数，计一次有效包围
        if (hitPoints.size() % 2 != 0) {
            oddCount++;
        }
    }

    // 🎯 终极裁决：只要在至少 1 个方向上被确认为“在内部”（奇数穿透）
    // 就足以证明这个底面属于外面的侧壁！
    return oddCount >= 1;
}

// 核心工序：将错拆的“散件”重新合并
std::vector<TopoDS_Compound> MergeNestedIslands(const std::vector<TopoDS_Compound>& bfsParts) {
    if (bfsParts.size() <= 1) return bfsParts;

    std::vector<bool> isMerged(bfsParts.size(), false);
    std::vector<TopoDS_Compound> finalParts;
    BRep_Builder builder;

    for (size_t i = 0; i < bfsParts.size(); ++i) {
        if (isMerged[i]) continue;

        TopoDS_Compound currentMainPart = bfsParts[i];

        for (size_t j = i + 1; j < bfsParts.size(); ++j) {
            if (isMerged[j]) continue;

            // =======================================================
            // 保险 1: 若两者极端靠近（<0.5mm），直接合并
            // =======================================================
            BRepExtrema_DistShapeShape distCalc(bfsParts[j], currentMainPart);
            if (distCalc.IsDone() && distCalc.Value() < 0.5) {
                builder.Add(currentMainPart, bfsParts[j]);
                isMerged[j] = true;
                std::cout << "  [物理强吸] 发现微小裂缝缝合遗漏，已强行合并侧壁与底面！" << std::endl;
                continue;
            }

            // =======================================================
            // 保险 2: 半山腰全向射线法鉴定 (解决大缝隙悬空面/孤岛)
            // =======================================================
            if (IsCompoundInsideExact(bfsParts[j], currentMainPart)) {
                builder.Add(currentMainPart, bfsParts[j]);
                isMerged[j] = true;
                std::cout << "  [全向射线] 悬空底面/孤岛 (Sub_" << j << ") 属于主侧壁内部，已合并！" << std::endl;
            }
            else if (IsCompoundInsideExact(currentMainPart, bfsParts[j])) {
                TopoDS_Compound newMain;
                builder.MakeCompound(newMain);
                builder.Add(newMain, bfsParts[j]);
                builder.Add(newMain, currentMainPart);
                currentMainPart = newMain;
                isMerged[j] = true;
                std::cout << "  [全向射线] 主侧壁包含底面/孤岛，已合并！" << std::endl;
            }
        }
        finalParts.push_back(currentMainPart);
    }

    return finalParts;
}

/**
 * @brief 将物理切分后的零件块转换成独立的 CavityFeature 特征
 * @param parts 经过物理切分、缝合、孤岛合并后的独立零件列表
 * @param parentFeat 该零件所属的原始父型腔特征（包含完整的切片环信息）
 * @param faceToIdMap 原始模型 ID 映射
 * @return 转换后的、包含该零件专属特征信息的数组
 */

// --- Compound Z 分割 ---
std::pair<TopoDS_Compound, TopoDS_Compound> SplitCompoundAtZ(
    const TopoDS_Compound& cavity, double splitZ) {

    BRep_Builder builder;
    TopoDS_Compound upper, lower;
    builder.MakeCompound(upper);
    builder.MakeCompound(lower);

    TopTools_DataMapOfIntegerShape idToFace;
    std::set<int> aboveIds, belowIds, facesToSplit;
    int nextId = 1;

    TopExp_Explorer exp(cavity, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face f = TopoDS::Face(exp.Current());
        int id = nextId++;
        idToFace.Bind(id, f);

        if (IsHorizontalFace(f)) {
            double z = GetHorizontalFaceZ(f);
            if (std::abs(z - splitZ) <= 1e-4) {
                gp_Dir normal = GetFaceNormal(f);
                if (normal.Z() < -0.5) belowIds.insert(id);
                else aboveIds.insert(id);
            } else if (z > splitZ) {
                aboveIds.insert(id);
            } else {
                belowIds.insert(id);
            }
        } else {
            double zmin, zmax;
            GetFaceZRange(f, zmin, zmax);
            if (zmin >= splitZ - 1e-4) {
                aboveIds.insert(id);
            } else if (zmax <= splitZ + 1e-4) {
                belowIds.insert(id);
            } else {
                facesToSplit.insert(id);
            }
        }
    }

    gp_Pln cuttingPlane(gp_Pnt(0, 0, splitZ), gp_Dir(0, 0, 1));
    TopoDS_Face planeFace = BRepBuilderAPI_MakeFace(cuttingPlane);

    for (int id : facesToSplit) {
        TopoDS_Face faceToCut = TopoDS::Face(idToFace.Find(id));

        try {
            BRepAlgoAPI_Section section(faceToCut, planeFace);
            section.Build();

            if (!section.IsDone() || !TopExp_Explorer(section.Shape(), TopAbs_EDGE).More()) {
                belowIds.insert(id);
                continue;
            }

            BRepFeat_SplitShape splitter(faceToCut);
            TopExp_Explorer edgeExp(section.Shape(), TopAbs_EDGE);
            for (; edgeExp.More(); edgeExp.Next()) {
                splitter.Add(TopoDS::Edge(edgeExp.Current()), faceToCut);
            }
            splitter.Build();

            if (splitter.IsDone()) {
                TopExp_Explorer faceExp(splitter.Shape(), TopAbs_FACE);
                bool producedSubFace = false;
                for (; faceExp.More(); faceExp.Next()) {
                    TopoDS_Face subFace = TopoDS::Face(faceExp.Current());
                    double smin, smax;
                    GetFaceZRange(subFace, smin, smax);
                    double zMid = (smin + smax) / 2.0;

                    int newId = nextId++;
                    idToFace.Bind(newId, subFace);
                    producedSubFace = true;

                    if (zMid > splitZ) aboveIds.insert(newId);
                    else belowIds.insert(newId);
                }
                if (!producedSubFace) {
                    belowIds.insert(id);
                }
            } else {
                belowIds.insert(id);
            }
        }
        catch (const Standard_Failure& e) {
            std::cerr << "  [SplitCompoundAtZ] Z=" << splitZ
                << " 切分面ID=" << id
                << " OCCT异常: " << e.GetMessageString()
                << "，保留到下方" << std::endl;
            belowIds.insert(id);
        }
        catch (...) {
            std::cerr << "  [SplitCompoundAtZ] Z=" << splitZ
                << " 切分面ID=" << id
                << " 未知异常，保留到下方" << std::endl;
            belowIds.insert(id);
        }
    }

    for (int id : aboveIds) {
        if (idToFace.IsBound(id)) builder.Add(upper, idToFace.Find(id));
    }
    for (int id : belowIds) {
        if (idToFace.IsBound(id)) builder.Add(lower, idToFace.Find(id));
    }

    return { upper, lower };
}


std::vector<TopoDS_Compound> RecursiveSplitCavity(const TopoDS_Compound& cavity) {
    std::vector<TopoDS_Compound> result;

    int inputFaceCount = CountFacesInCompound(cavity);
    if (inputFaceCount == 0) {
        std::cout << "  [递推切割] 收到空型腔，跳过。" << std::endl;
        return result;
    }

    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(cavity, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    Bnd_Box cavityBox;
    BRepBndLib::Add(cavity, cavityBox);
    if (cavityBox.IsVoid()) {
        std::cout << "  [递推切割] 型腔包围盒为空，跳过。faces=" << inputFaceCount << std::endl;
        return result;
    }

    double xmin, ymin, zmin, xmax, ymax, topZ;
    try {
        cavityBox.Get(xmin, ymin, zmin, xmax, ymax, topZ);
    }
    catch (const Standard_Failure& e) {
        std::cerr << "  [递推切割] Bnd_Box::Get异常: " << e.GetMessageString()
            << " faces=" << inputFaceCount << std::endl;
        return result;
    }

    std::set<double> zSet;
    TopExp_Explorer exp(cavity, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face f = TopoDS::Face(exp.Current());
        if (!IsHorizontalFace(f)) continue;

        double zFace = GetHorizontalFaceZ(f);
        if (std::abs(zFace - topZ) < 1e-3) continue;

        bool isIntermediateStep = false;

        TopExp_Explorer edgeExp(f, TopAbs_EDGE);
        for (; edgeExp.More(); edgeExp.Next()) {
            const TopoDS_Edge& edge = TopoDS::Edge(edgeExp.Current());
            if (!edgeToFaces.Contains(edge)) continue;

            const TopTools_ListOfShape& neighbors = edgeToFaces.FindFromKey(edge);
            TopTools_ListIteratorOfListOfShape it(neighbors);
            for (; it.More(); it.Next()) {
                TopoDS_Face nFace = TopoDS::Face(it.Value());
                if (nFace.IsSame(f)) continue;

                double nZmin, nZmax;
                GetFaceZRange(nFace, nZmin, nZmax);
                if (nZmin < zFace - 1e-3) {
                    isIntermediateStep = true;
                    break;
                }
            }
            if (isIntermediateStep) break;
        }

        if (isIntermediateStep) zSet.insert(zFace);
    }

    if (zSet.empty()) {
        result.push_back(cavity);
        return result;
    }

    double splitZ = *zSet.rbegin();
    std::cout << "  [递推切割] Z = " << splitZ << std::endl;

    auto [upper, lower] = SplitCompoundAtZ(cavity, splitZ);
    int upperFaces = CountFacesInCompound(upper);
    int lowerFaces = CountFacesInCompound(lower);
    std::cout << "  [递推切割] 切后 upperFaces=" << upperFaces
        << " lowerFaces=" << lowerFaces << std::endl;

    if (upperFaces > 0) {
        result.push_back(upper);
    }

    if (lowerFaces == 0) {
        return result;
    }

    std::vector<TopoDS_Compound> lowerParts = SeparateDisconnectedCavities(lower);
    std::cout << "  [连通性] 下方拆分为 " << lowerParts.size() << " 个独立部分" << std::endl;

    if (lowerParts.size() > 1) {
        for (const auto& part : lowerParts) {
            if (CountFacesInCompound(part) == 0) continue;
            auto subResult = RecursiveSplitCavity(part);
            result.insert(result.end(), subResult.begin(), subResult.end());
        }
    } else {
        auto subResult = RecursiveSplitCavity(lower);
        result.insert(result.end(), subResult.begin(), subResult.end());
    }

    return result;
}


// --- 开放型腔 Z 切分 ---
std::vector<double> MergeCloseSplitZs(std::vector<double> splitZs, double minGap) {
    if (splitZs.empty()) return splitZs;

    std::sort(splitZs.begin(), splitZs.end(), std::greater<double>());
    std::vector<double> merged;
    for (double z : splitZs) {
        if (merged.empty() || std::abs(merged.back() - z) >= minGap) {
            merged.push_back(z);
        }
    }
    return merged;
}

std::vector<double> CollectOpenCavitySplitZs(
    const TopoDS_Compound& cavity,
    const CavityFeature& feature,
    const OpenCavitySplitParams& params)
{
    std::vector<double> splitZs;

    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(cavity, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    TopExp_Explorer exp(cavity, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        if (!IsHorizontalFace(face)) continue;

        double zFace = GetHorizontalFaceZ(face);
        if (zFace >= feature.topZ - params.zProtectionTol ||
            zFace <= feature.bottomZ + params.zProtectionTol) {
            std::cout << "  [开放切分候选] Z=" << zFace << " 被跳过: 顶/底保护区" << std::endl;
            continue;
        }

        double area = CalculateFaceAreaOCC(face);
        if (area < params.minStepFaceArea) {
            std::cout << "  [开放切分候选] Z=" << zFace << " 被跳过: 面积过小 area=" << area << std::endl;
            continue;
        }

        bool hasUpNeighbor = false;
        bool hasDownNeighbor = false;

        TopExp_Explorer edgeExp(face, TopAbs_EDGE);
        for (; edgeExp.More(); edgeExp.Next()) {
            const TopoDS_Edge& edge = TopoDS::Edge(edgeExp.Current());
            if (!edgeToFaces.Contains(edge)) continue;

            const TopTools_ListOfShape& neighbors = edgeToFaces.FindFromKey(edge);
            TopTools_ListIteratorOfListOfShape it(neighbors);
            for (; it.More(); it.Next()) {
                TopoDS_Face neighbor = TopoDS::Face(it.Value());
                if (neighbor.IsSame(face)) continue;

                double nZmin, nZmax;
                GetFaceZRange(neighbor, nZmin, nZmax);
                if (nZmax > zFace + params.neighborZTol) {
                    hasUpNeighbor = true;
                }
                if (nZmin < zFace - params.neighborZTol) {
                    hasDownNeighbor = true;
                }
            }
        }

        bool isMiddleStep = hasUpNeighbor && hasDownNeighbor;
        bool isDownOnlyStep = params.splitDownOnlyStepFaces && hasDownNeighbor && !hasUpNeighbor;

        if (!isMiddleStep && !isDownOnlyStep) {
            std::cout << "  [开放切分候选] Z=" << zFace << " 被跳过: 缺少上下连续邻接"
                << " up=" << hasUpNeighbor << " down=" << hasDownNeighbor
                << " area=" << area << std::endl;
            continue;
        }

        std::cout << "  [开放切分候选] Z=" << zFace << " 被采用: "
            << (isDownOnlyStep ? "仅下方邻接凸台/岛顶面" : "上下连续台阶")
            << " area=" << area << std::endl;
        splitZs.push_back(zFace);
    }

    return MergeCloseSplitZs(splitZs, params.minSplitGap);
}

bool GetCompoundZRangeSafe(const TopoDS_Compound& compound, double& zmin, double& zmax) {
    Bnd_Box box;
    BRepBndLib::Add(compound, box);
    if (box.IsVoid()) return false;

    double xmin, ymin, xmax, ymax;
    try {
        box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    }
    catch (...) {
        return false;
    }
    return true;
}

std::vector<double> CollectOpenCavitySplitZsForPart(
    const TopoDS_Compound& part,
    const CavityFeature& feature,
    const OpenCavitySplitParams& params,
    const std::string& branchName)
{
    double localZmin = 0.0;
    double localZmax = 0.0;
    if (!GetCompoundZRangeSafe(part, localZmin, localZmax)) {
        return {};
    }

    std::vector<double> raw = CollectOpenCavitySplitZs(part, feature, params);
    std::vector<double> local;
    for (double z : raw) {
        if (z >= localZmax - params.zProtectionTol ||
            z <= localZmin + params.zProtectionTol) {
            std::cout << "  [开放Z切分] branch=" << branchName
                << " 跳过局部顶/底切分点 Z=" << z
                << " localZ=[" << localZmin << "," << localZmax << "]" << std::endl;
            continue;
        }
        local.push_back(z);
    }

    local = MergeCloseSplitZs(local, params.minSplitGap);
    std::cout << "  [开放Z切分] branch=" << branchName
        << " localZ=[" << localZmin << "," << localZmax << "]"
        << " 局部候选数量=" << local.size() << ": ";
    for (double z : local) std::cout << z << " ";
    std::cout << std::endl;
    return local;
}

std::vector<TopoDS_Compound> MergeGeometricallyConnectedOpenParts(
    const std::vector<TopoDS_Compound>& rawParts,
    double tolerance,
    int featureId,
    const std::string& branchName,
    int depth,
    bool exportDebugBreps)
{
    if (rawParts.size() <= 1) return rawParts;

    UnionFind uf((int)rawParts.size());

    for (size_t i = 0; i < rawParts.size(); ++i) {
        for (size_t j = i + 1; j < rawParts.size(); ++j) {
            BRepExtrema_DistShapeShape distCalc(rawParts[i], rawParts[j]);
            if (!distCalc.IsDone()) {
                std::cout << "  [开放几何连通] feature=" << featureId
                    << " branch=" << branchName
                    << " depth=" << depth
                    << " part " << i << "-" << j
                    << " 距离计算失败" << std::endl;
                continue;
            }

            double distance = distCalc.Value();
            if (distance <= tolerance) {
                std::cout << "  [开放几何连通] feature=" << featureId
                    << " branch=" << branchName
                    << " depth=" << depth
                    << " 合并 part " << i << " + " << j
                    << " distance=" << distance
                    << " tol=" << tolerance << std::endl;
                uf.Unite((int)i, (int)j, distance);
            }
        }
    }

    std::map<int, std::vector<int>> groups;
    for (int i = 0; i < (int)rawParts.size(); ++i) {
        groups[uf.Find(i)].push_back(i);
    }

    std::vector<TopoDS_Compound> mergedParts;
    BRep_Builder builder;
    int mergedIndex = 0;
    for (const auto& entry : groups) {
        TopoDS_Compound merged;
        builder.MakeCompound(merged);

        for (int idx : entry.second) {
            TopExp_Explorer exp(rawParts[idx], TopAbs_FACE);
            for (; exp.More(); exp.Next()) {
                builder.Add(merged, exp.Current());
            }
        }

        if (CountFacesInCompound(merged) > 0) {
            if (exportDebugBreps) {
                std::string mergedPath = savePath + "OpenSplit_Feature_" + std::to_string(featureId)
                    + "_Branch_" + branchName
                    + "_Depth_" + std::to_string(depth)
                    + "_MergedLowerPart_" + std::to_string(mergedIndex)
                    + "_Members_" + std::to_string((int)entry.second.size()) + ".brep";
                BRepTools::Write(merged, mergedPath.c_str());
            }
            mergedParts.push_back(merged);
            mergedIndex++;
        }
    }

    std::cout << "  [开放几何连通] feature=" << featureId
        << " branch=" << branchName
        << " depth=" << depth
        << " rawParts=" << rawParts.size()
        << " mergedParts=" << mergedParts.size() << std::endl;

    return mergedParts;
}

std::vector<TopoDS_Compound> SplitOpenCavityByZPlan(
    const TopoDS_Compound& cavity,
    const CavityFeature& feature,
    const OpenCavitySplitParams& params,
    int featureId)
{
    std::vector<TopoDS_Compound> result;

    struct OpenSplitTask {
        TopoDS_Compound shape;
        int depth = 0;
        std::string branchName;
    };

    std::vector<OpenSplitTask> tasks;
    tasks.push_back({ cavity, 0, "root" });

    while (!tasks.empty()) {
        OpenSplitTask task = tasks.back();
        tasks.pop_back();

        int taskFaces = CountFacesInCompound(task.shape);
        if (taskFaces == 0) {
            continue;
        }

        if (task.depth >= params.maxRecursionDepth) {
            std::cout << "  [开放Z切分] branch=" << task.branchName
                << " 达到最大递归深度，作为最终特征保留。" << std::endl;
            result.push_back(task.shape);
            continue;
        }

        std::vector<double> localSplitZs =
            CollectOpenCavitySplitZsForPart(task.shape, feature, params, task.branchName);

        if (localSplitZs.empty()) {
            result.push_back(task.shape);
            continue;
        }

        double splitZ = localSplitZs.front();
        auto [upper, lower] = SplitCompoundAtZ(task.shape, splitZ);

        int upperFaces = CountFacesInCompound(upper);
        int lowerFaces = CountFacesInCompound(lower);
        std::cout << "  [开放Z切分] feature=" << featureId
            << " branch=" << task.branchName
            << " depth=" << task.depth
            << " Z=" << splitZ
            << " localCandidates=" << localSplitZs.size()
            << " inputFaces=" << taskFaces
            << " upperFaces=" << upperFaces
            << " lowerFaces=" << lowerFaces << std::endl;

        if (upperFaces == 0 || lowerFaces == 0) {
            std::cout << "  [开放Z切分] branch=" << task.branchName
                << " 切分未产生有效上下层，作为最终特征保留，避免重复递推。" << std::endl;
            result.push_back(task.shape);
            continue;
        }

        if (params.exportDebugBreps) {
            std::string upperPath = savePath + "OpenSplit_Feature_" + std::to_string(featureId)
                + "_Branch_" + task.branchName
                + "_Depth_" + std::to_string(task.depth) + "_Upper_Z" + std::to_string(splitZ) + ".brep";
            std::string lowerPath = savePath + "OpenSplit_Feature_" + std::to_string(featureId)
                + "_Branch_" + task.branchName
                + "_Depth_" + std::to_string(task.depth) + "_Lower_Z" + std::to_string(splitZ) + ".brep";
            BRepTools::Write(upper, upperPath.c_str());
            BRepTools::Write(lower, lowerPath.c_str());
        }

        if (upperFaces > 0) {
            result.push_back(upper);
        }

        if (lowerFaces == 0) {
            continue;
        }

        std::vector<TopoDS_Compound> rawLowerParts = SeparateDisconnectedCavities(lower);
        std::cout << "  [开放Z切分连通性] feature=" << featureId
            << " branch=" << task.branchName
            << " depth=" << task.depth
            << " rawParts=" << rawLowerParts.size() << std::endl;

        if (params.exportDebugBreps) {
            for (size_t p = 0; p < rawLowerParts.size(); ++p) {
                std::string rawPath = savePath + "OpenSplit_Feature_" + std::to_string(featureId)
                    + "_Branch_" + task.branchName
                    + "_Depth_" + std::to_string(task.depth)
                    + "_RawLowerPart_" + std::to_string(p) + ".brep";
                BRepTools::Write(rawLowerParts[p], rawPath.c_str());
            }
        }

        if (rawLowerParts.empty()) {
            tasks.push_back({ lower, task.depth + 1, task.branchName + "_L" });
            continue;
        }

        std::vector<TopoDS_Compound> lowerParts = MergeGeometricallyConnectedOpenParts(
            rawLowerParts,
            params.geometricConnectTol,
            featureId,
            task.branchName,
            task.depth,
            params.exportDebugBreps);

        if (lowerParts.empty()) {
            tasks.push_back({ lower, task.depth + 1, task.branchName + "_L" });
            continue;
        }

        for (size_t p = 0; p < lowerParts.size(); ++p) {
            int partFaces = CountFacesInCompound(lowerParts[p]);
            if (partFaces == 0) continue;

            std::string childBranch = task.branchName + "_L" + std::to_string(p);
            if (params.exportDebugBreps) {
                std::string partPath = savePath + "OpenSplit_Feature_" + std::to_string(featureId)
                    + "_Branch_" + childBranch
                    + "_AfterDepth_" + std::to_string(task.depth)
                    + "_LowerPart.brep";
                BRepTools::Write(lowerParts[p], partPath.c_str());
            }

            tasks.push_back({ lowerParts[p], task.depth + 1, childBranch });
        }
    }

    if (result.empty() && CountFacesInCompound(cavity) > 0) {
        result.push_back(cavity);
    }

    return result;
}

std::vector<TopoDS_Compound> SplitOpenCavityConservatively(
    const TopoDS_Compound& cavity,
    const CavityFeature& feature,
    const OpenCavitySplitParams& params)
{
    std::vector<double> splitZs = CollectOpenCavitySplitZs(cavity, feature, params);
    std::cout << "  [开放Z切分] feature=" << feature.featureId
        << " 候选切分点数量=" << splitZs.size() << ": ";
    for (double z : splitZs) {
        std::cout << z << " ";
    }
    std::cout << std::endl;

    if (splitZs.empty()) {
        return { cavity };
    }

    return SplitOpenCavityByZPlan(cavity, feature, params, feature.featureId);
}

/**
 * @brief 【开放型腔专属】全量特征提取与物理分割大管家
 * @param allOpenLayerFaces 开放型腔专用的切片数据容器
 * @param allClosedLayerFaces 封闭型腔专用的切片数据容器（用于嵌套检测）
 * @param mainShape 原始加载的完整 STEP 实体模型 (TopoDS_Shape)
 * @param splitPoints 自适应获取并排好序的所有 Z 轴切分高度
 * @param faceToIdMap 原始模型的面到唯一 ID 的绑定映射表
 * @param savePath 文件输出的绝对路径
 */
std::vector<CavityFeature> ProcessAndSplitOpenCavityFeatures(
    const std::vector<std::vector<Face2D>>& allOpenLayerFaces,
    const std::vector<std::vector<Face2D>>& allClosedLayerFaces,
    const std::vector<TopoDS_Compound>& trueClosedCavities,
    const TopoDS_Shape& mainShape,
    const std::vector<double>& splitPoints,
    const TopTools_DataMapOfShapeInteger& faceToIdMap,
    const std::string& savePath)
{
    // 1. 专属户口清查：只搜刮经过精密圆角穿透提纯后的开放型腔面 ID
    std::set<int> openFaceIds = CollectOpenCavityFaceIds(allOpenLayerFaces);
    if (openFaceIds.empty()) {
        return {};
    }

    // 2. 得到纯粹的开放特征壳体（包含专属的开放底面判定）
    std::string cavityFacesFileName = savePath + "OpenCavity_WallFaces.brep";
    ExportOpenCavityFaces(mainShape, openFaceIds, faceToIdMap, cavityFacesFileName);

    // 获取对应的开放合并面 Compound
    TopoDS_Compound cavityCompound = GetOpenCavityCompound(mainShape, openFaceIds, faceToIdMap);


    // 保存为 brep 文件
    std::string cavityCompoundFileName = savePath + "OpenCavityCompound.brep";
    BRepTools::Write(cavityCompound, cavityCompoundFileName.c_str());
    std::cout << "  开放型腔合并面已保存至: " << cavityCompoundFileName << std::endl;

    // 3. BFS 拓扑连通体打散
    std::vector<TopoDS_Compound> isolatedCavities = SeparateDisconnectedCavities(cavityCompound);
    cout << "    [拓扑分析] 成功将开放合并面打散为 " << isolatedCavities.size() << " 个空间独立的开放区域。" << endl;

    // 4. 生成高级特征表达与智能剪裁（三维壳体精密还原）
    std::vector<CavityFeature> cavityFeatures = GenerateOpenCavityFeatures(isolatedCavities, allOpenLayerFaces, splitPoints, faceToIdMap);
    
    vector<TopoDS_Compound> trueOpenCavities;

    for (const auto& feat : cavityFeatures) {
        TopoDS_Compound trueCavity = ExtractTrueOpenCavityFaces(feat, mainShape, faceToIdMap, allClosedLayerFaces, splitPoints, trueClosedCavities);
		trueOpenCavities.push_back(trueCavity);
        std::string filePath = savePath + "Final_True_OpenCavity_" + std::to_string(feat.featureId) + ".brep";
        BRepTools::Write(trueCavity, filePath.c_str());
    }

    std::vector<TopoDS_Compound> filteredOpenCavities;
    std::vector<CavityFeature> filteredFeatures;
    for (size_t i = 0; i < trueOpenCavities.size(); ++i) {
        gp_Dir toolDir;
        bool machinable = IsCavityMachinable(trueOpenCavities[i], toolDir);
        std::cout << "  [加工可行性赛选] OpenCavity[" << i << "] 进刀方向=("
            << toolDir.X() << "," << toolDir.Y() << "," << toolDir.Z()
            << ") -> " << (machinable ? "可加工，保留" : "不可加工，丢弃") << std::endl;
        if (machinable) {
            filteredOpenCavities.push_back(trueOpenCavities[i]);
            CavityFeature updatedFeat = cavityFeatures[i];
            updatedFeat.type = CavityType::OPEN;
            updatedFeat.toolDirection = toolDir;
            filteredFeatures.push_back(updatedFeat);
        }
    }
    trueOpenCavities = filteredOpenCavities;
    cavityFeatures = filteredFeatures;
    std::cout << "  [加工可行性赛选] 赛选后保留 " << trueOpenCavities.size() << " 个开放型腔。" << std::endl;
 
    // 5. 沿 Z 轴自适应多级物理切割
    int totalPartCount = 0;
    std::vector<CavityFeature> finalMachinableOpenFeatures;
    std::vector<TopoDS_Compound> allOpenPartCompounds;
    OpenCavitySplitParams openSplitParams;

    for (size_t i = 0; i < trueOpenCavities.size(); ++i) {
        bool isZDirection = (cavityFeatures[i].toolDirection.IsEqual(gp_Dir(0, 0, -1), 1e-6));

        if (!isZDirection) {
            std::cout << "  [开放分割跳过] feature=" << cavityFeatures[i].featureId
                << " 进刀方向=(" << cavityFeatures[i].toolDirection.X()
                << "," << cavityFeatures[i].toolDirection.Y()
                << "," << cavityFeatures[i].toolDirection.Z()
                << ") 非(0,0,-1)，跳过Z轴分割，整体保留。" << std::endl;

            BRepBuilderAPI_Sewing sewer(1e-2);
            sewer.Add(trueOpenCavities[i]);
            sewer.Perform();
            TopoDS_Shape sewedShape = sewer.SewedShape();
            TopoDS_Compound sewedCompound;
            BRep_Builder compBuilder;
            compBuilder.MakeCompound(sewedCompound);
            TopExp_Explorer faceExp(sewedShape, TopAbs_FACE);
            for (; faceExp.More(); faceExp.Next()) {
                compBuilder.Add(sewedCompound, faceExp.Current());
            }
            std::string fileName = savePath + "Final_CAM_OpenPart_" + std::to_string(totalPartCount) + ".brep";
            BRepTools::Write(sewedCompound, fileName.c_str());

            CavityFeature convertedFeat = cavityFeatures[i];
            convertedFeat.featureId = totalPartCount;
            finalMachinableOpenFeatures.push_back(convertedFeat);
            totalPartCount++;
            continue;
        }

        // 获取当前父型腔沿 Z 轴保守切分后的所有层
        std::vector<TopoDS_Compound> zLayers =
            SplitOpenCavityConservatively(trueOpenCavities[i], cavityFeatures[i], openSplitParams);

        for (size_t j = 0; j < zLayers.size(); ++j) {

            // 🎯 【核心拓扑修复】：缝合
            BRepBuilderAPI_Sewing sewer(1e-2);
            sewer.Add(zLayers[j]);
            sewer.Perform();
            TopoDS_Shape sewedShape = sewer.SewedShape();

            TopoDS_Compound sewedCompound;
            BRep_Builder compBuilder;
            compBuilder.MakeCompound(sewedCompound);
            TopExp_Explorer faceExp(sewedShape, TopAbs_FACE);
            for (; faceExp.More(); faceExp.Next()) {
                compBuilder.Add(sewedCompound, faceExp.Current());
            }
            std::string fileName = savePath + "temp_sew" + std::to_string(totalPartCount) + ".brep";
            BRepTools::Write(sewedCompound, fileName.c_str());


            std::vector<TopoDS_Compound> isolatedSubParts;

            // 开放型腔的同一Z层可能在拓扑上断开，但加工语义上仍属于同一特征。
            isolatedSubParts.push_back(sewedCompound);
            std::cout << "  [开放层保持合并] feature=" << cavityFeatures[i].featureId
                << " layer=" << j
                << " faces=" << CountFacesInCompound(sewedCompound) << std::endl;
        

            std::vector<CavityFeature> converted = ConvertPartsToFeatures(isolatedSubParts, cavityFeatures[i], faceToIdMap);
           
            for (size_t k = 0; k < converted.size(); ++k) {
                // 1. 设置正确的 ID
                converted[k].featureId = totalPartCount;

                // 2. 保存几何体 (这里用 k 访问 isolatedSubParts)
                std::string fileName = savePath + "Final_CAM_OpenPart_" + std::to_string(totalPartCount) + ".brep";
                BRepTools::Write(isolatedSubParts[k], fileName.c_str());

                // 3. 存入总表
                finalMachinableOpenFeatures.push_back(converted[k]);

                totalPartCount++;
            }
        }
    }

    for (const auto& feat : finalMachinableOpenFeatures) {
        std::cout << "\n>>> [独立可加工开放型腔特征 ID: " << feat.featureId << "]" << std::endl;
        std::cout << "    - 类型: " << (feat.type == CavityType::OPEN ? "OPEN" : feat.type == CavityType::CLOSED ? "CLOSED" : "OTHER") << std::endl;
        std::cout << "    - 进刀方向: (" << feat.toolDirection.X() << ", " << feat.toolDirection.Y() << ", " << feat.toolDirection.Z() << ")" << std::endl;
        std::cout << "    - 空间范围: Z_Top = " << feat.topZ << " / Z_Floor = " << feat.bottomZ << std::endl;
        std::cout << "    - 深度: " << feat.totalDepth << " mm" << std::endl;
        std::cout << "    - 包含 2D 层级切片数量: " << feat.stepLoops.size() << " 圈" << std::endl;
    }

    return finalMachinableOpenFeatures;
}

// 封闭型腔提取流程
/**
 * @brief 【型腔全自动特征提取与物理分割引擎】
 * 输入全层剖分面数据，自动抽离型腔，提取三维面，按连通体打散，并执行 Z 轴物理切割输出组件零件。
 * * @param allLayerFaces 包含了每层所有 Face2D（SOLID, CAVITY, HULL）的全局大容器
 * @param mainShape 原始加载的完整 STEP 实体模型 (TopoDS_Shape)
 * @param splitPoints 之前自适应获取并排好序的所有 Z 轴切分高度
 * @param faceToIdMap 原始模型的面到唯一 ID 的绑定映射表
 * @param savePath 文件输出的绝对路径
 */
