#include "interior_open.h"
#include "geom_utils.h"
#include "occ_utils.h"
#include "export_utils.h"
#include "boolean_ops.h"
#include "slicing.h"
#include "closed_cavity.h"
#include "open_cavity.h"

// --- 中部开放型腔 ---
double GetFace2DZForInteriorOpen(const Face2D& face) {
    if (!face.outerLoop.empty()) return face.outerLoop.front().start.z;
    for (const auto& inner : face.innerLoops) {
        if (!inner.empty()) return inner.front().start.z;
    }
    return 0.0;
}

bool HasInteriorOpenCavity(
    const std::vector<std::vector<Face2D>>& allOpenLayerFaces,
    const std::vector<std::vector<Face2D>>& allClosedLayerFaces,
    double modelMaxZ,
    const std::vector<double>& splitPoints,
    double zTol)
{
    std::cout << "\n===== [中部开放旁路检测] modelMaxZ=" << modelMaxZ
        << " zTol=" << zTol << " =====" << std::endl;

    bool hasOpenFace = false;
    double topOpenZ = -std::numeric_limits<double>::max();
    for (const auto& layer : allOpenLayerFaces) {
        for (const auto& face : layer) {
            if (face.type != FaceType::OPENCAVITY || face.outerLoop.empty()) continue;
            double openZ = GetFace2DZForInteriorOpen(face);
            if (!hasOpenFace || openZ > topOpenZ) { hasOpenFace = true; topOpenZ = openZ; }
        }
    }

    if (!hasOpenFace) {
        std::cout << "  未检测到 OPENCAVITY 切片，保持原有分支。" << std::endl;
        std::cout << "===== [中部开放旁路检测] 保持原有分支 =====\n" << std::endl;
        return false;
    }

    double gapToTop = modelMaxZ - topOpenZ;
    bool hasInterior = gapToTop > zTol;
    std::cout << "  topOpenZ=" << topOpenZ
        << " modelMaxZ-topOpenZ=" << gapToTop
        << " threshold=" << zTol
        << " -> " << (hasInterior ? "中部开放" : "顶部开放") << std::endl;

    if (!hasInterior) {
        std::cout << "===== [中部开放旁路检测] 保持原有分支 =====\n" << std::endl;
        return false;
    }

    // 找到紧邻上方的切分点，判断该切片是否是封闭型腔
    double nextSplitZ = modelMaxZ;
    for (double sp : splitPoints) {
        if (sp > topOpenZ && sp < nextSplitZ) nextSplitZ = sp;
    }

    bool foundAdjacentClosed = false;
    for (const auto& layer : allClosedLayerFaces) {
        for (const auto& face : layer) {
            if (face.type != FaceType::CAVITY || face.outerLoop.empty()) continue;
            double z = face.outerLoop.front().start.z;
            if (std::abs(z - nextSplitZ) <= zTol) {
                foundAdjacentClosed = true;
                std::cout << "  紧邻上方切分点 Z=" << nextSplitZ << " -> 封闭型腔切片，可进入中部开放分支。" << std::endl;
                break;
            }
        }
        if (foundAdjacentClosed) break;
    }

    if (!foundAdjacentClosed) {
        std::cout << "  紧邻上方切分点 Z=" << nextSplitZ << " -> 非封闭型腔切片（实体或其他），保持原有分支。" << std::endl;
        std::cout << "===== [中部开放旁路检测] 保持原有分支 =====\n" << std::endl;
        return false;
    }

    std::cout << "===== [中部开放旁路检测] 进入新旁路 =====\n" << std::endl;
    return true;
}

double FindNextLowerSplitPoint(
    const std::vector<double>& splitPoints,
    double z,
    double fallbackZ)
{
    double best = -std::numeric_limits<double>::max();
    for (double sp : splitPoints) {
        if (sp < z - 1e-4 && sp > best) {
            best = sp;
        }
    }
    if (best > -std::numeric_limits<double>::max() / 2.0) return best;
    return fallbackZ;
}

void AddFaceIdsFromFace2D(const Face2D& face, std::set<int>& ids) {
    for (const auto& line : face.outerLoop) {
        if (line.faceId >= 0) ids.insert(line.faceId);
    }
    for (const auto& inner : face.innerLoops) {
        for (const auto& line : inner) {
            if (line.faceId >= 0) ids.insert(line.faceId);
        }
    }
}

void PrintIdSet(const std::string& label, const std::set<int>& ids) {
    std::cout << label << " (" << ids.size() << "): ";
    for (int id : ids) std::cout << id << " ";
    std::cout << std::endl;
}

std::set<int> CollectClosedSideFaceIdsAroundOpenBand(
    const std::vector<std::vector<Face2D>>& allClosedLayerFaces,
    double openTopZ,
    double openBottomZ)
{
    std::set<int> closedSideIds;
    int closedFaceCount = 0;

    for (const auto& layer : allClosedLayerFaces) {
        for (const auto& face : layer) {
            if (face.type != FaceType::CAVITY || face.outerLoop.empty()) continue;

            double z = GetFace2DZForInteriorOpen(face);
            if (z <= openTopZ + 1e-4 && z >= openBottomZ - 1e-4) {
                continue;
            }

            AddFaceIdsFromFace2D(face, closedSideIds);
            closedFaceCount++;
        }
    }

    std::cout << "  [同源验证] 开放band外侧封闭CAVITY切片数=" << closedFaceCount << std::endl;
    return closedSideIds;
}

bool ShouldDiscardInteriorOpenRegion(
    const std::set<int>& openSideIds,
    const std::set<int>& closedSideIds,
    const TopoDS_Compound& topCapFaces,
    std::set<int>& sharedSideIds,
    double minSharedRatio)
{
    sharedSideIds.clear();
    for (int id : openSideIds) {
        if (closedSideIds.count(id)) sharedSideIds.insert(id);
    }

    double openCoverage = openSideIds.empty()
        ? 0.0
        : (double)sharedSideIds.size() / (double)openSideIds.size();
    double closedCoverage = closedSideIds.empty()
        ? 0.0
        : (double)sharedSideIds.size() / (double)closedSideIds.size();
    bool sameOriginalWall = closedCoverage >= minSharedRatio;
    bool blockedByTopCap = CountFacesInCompound(topCapFaces) > 0;
    bool discardMiddleOpen = sameOriginalWall && blockedByTopCap;

    PrintIdSet("  [同源验证] openSideIds", openSideIds);
    PrintIdSet("  [同源验证] closedSideIds", closedSideIds);
    PrintIdSet("  [同源验证] sharedSideIds", sharedSideIds);
    std::cout << "  [同源验证] openCoverage(shared/open)=" << openCoverage << std::endl;
    std::cout << "  [同源验证] closedCoverage(shared/closed)=" << closedCoverage
        << " threshold=" << minSharedRatio
        << " sameOriginalWall=" << (sameOriginalWall ? "true" : "false") << std::endl;
    std::cout << "  [可加工性判断] blockedByTopCap="
        << (blockedByTopCap ? "true" : "false")
        << " topCapFaces=" << CountFacesInCompound(topCapFaces) << std::endl;
    std::cout << "  [重分类] discardMiddleOpen="
        << (discardMiddleOpen ? "true" : "false") << std::endl;

    return discardMiddleOpen;
}

CavityFeature BuildInteriorOpenCavityFeature(
    const std::vector<std::vector<Face2D>>& allOpenLayerFaces,
    const std::vector<std::vector<Face2D>>& allClosedLayerFaces,
    const std::vector<double>& splitPoints,
    double modelMinZ,
    const std::string& savePath)
{
    CavityFeature feat;
    feat.featureId = 0;
    feat.type = CavityType::OTHER;
    feat.topZ = -1e9;
    feat.bottomZ = 1e9;

    std::vector<Face2D> allOpenFaces;
    double minOpenSliceZ = 1e9;

    for (const auto& layer : allOpenLayerFaces) {
        for (const auto& face : layer) {
            if (face.type != FaceType::OPENCAVITY || face.outerLoop.empty()) continue;

            double z = GetFace2DZForInteriorOpen(face);
            CavityLoop outer;
            outer.lines = face.outerLoop;
            outer.zHeight = z;
            outer.isOuter = true;
            feat.stepLoops.push_back(outer);

            for (const auto& innerLoop : face.innerLoops) {
                CavityLoop inner;
                inner.lines = innerLoop;
                inner.zHeight = z;
                inner.isOuter = false;
                feat.stepLoops.push_back(inner);
            }

            feat.topZ = std::max(feat.topZ, z);
            minOpenSliceZ = std::min(minOpenSliceZ, z);
            AddFaceIdsFromFace2D(face, feat.sourceFaceIds);
            allOpenFaces.push_back(face);
        }
    }

    for (const auto& layer : allClosedLayerFaces) {
        for (const auto& face : layer) {
            if (face.type != FaceType::CAVITY || face.outerLoop.empty()) continue;
            AddFaceIdsFromFace2D(face, feat.sourceFaceIds);
        }
    }

    if (!feat.stepLoops.empty()) {
        feat.bottomZ = FindNextLowerSplitPoint(splitPoints, minOpenSliceZ, modelMinZ);
        feat.totalDepth = std::abs(feat.topZ - feat.bottomZ);
        std::sort(feat.stepLoops.begin(), feat.stepLoops.end(),
            [](const CavityLoop& a, const CavityLoop& b) {
                return a.zHeight > b.zHeight;
            });
    }

    std::string allOpenPath = savePath + "InteriorOpen_AllOpenSlices.brep";
    ExportFace2DToBrep(allOpenFaces, allOpenPath);
    std::cout << "  [中部开放旁路] 全部开放切片已保存: " << allOpenPath
        << " 数量=" << allOpenFaces.size() << std::endl;

    return feat;
}


TopoDS_Compound ExtractFacesByFeatureSourceIds(
    const CavityFeature& feature,
    const TopoDS_Shape& mainShape,
    const TopTools_DataMapOfShapeInteger& faceToIdMap)
{
    return GetFacesByFaceIds(mainShape, feature.sourceFaceIds, faceToIdMap);
}

TopoDS_Compound GetInteriorOpenBandCapFaces(
    const TopoDS_Shape& solid,
    const std::set<int>& sideFaceIds,
    const TopTools_DataMapOfShapeInteger& faceToIdMap,
    double topZ,
    double bottomZ,
    bool wantTopCap,
    double zTol)
{
    BRep_Builder builder;
    TopoDS_Compound capFaces;
    builder.MakeCompound(capFaces);

    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    const double targetZ = wantTopCap ? topZ : bottomZ;
    const gp_Dir toolDirection(0, 0, -1);
    int capCount = 0;

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

        if (std::abs(std::abs(normal.Z()) - 1.0) > 1e-3) continue;

        double faceZ = plane.Location().Z();
        if (std::abs(faceZ - targetZ) > zTol) continue;

        // 顶面：法向和进刀方向相同，即向下；底面：法向和进刀方向相反，即向上。
        if (wantTopCap) {
            if (normal.Dot(toolDirection) < 0.99) continue;
        } else {
            if (normal.Dot(toolDirection) > -0.99) continue;
        }

        bool hasSideNeighbor = false;
        TopExp_Explorer edgeExp(face, TopAbs_EDGE);
        for (; edgeExp.More(); edgeExp.Next()) {
            const TopoDS_Edge& edge = TopoDS::Edge(edgeExp.Current());
            if (!edgeToFaces.Contains(edge)) continue;

            const TopTools_ListOfShape& adjacentFaces = edgeToFaces.FindFromKey(edge);
            TopTools_ListIteratorOfListOfShape it(adjacentFaces);
            for (; it.More(); it.Next()) {
                TopoDS_Face adjFace = TopoDS::Face(it.Value());
                if (adjFace.IsSame(face)) continue;
                if (!faceToIdMap.IsBound(adjFace)) continue;

                int adjId = faceToIdMap.Find(adjFace);
                if (sideFaceIds.count(adjId)) {
                    hasSideNeighbor = true;
                    break;
                }
            }
            if (hasSideNeighbor) break;
        }

        if (hasSideNeighbor) {
            builder.Add(capFaces, face);
            capCount++;
        }
    }

    std::cout << "  [中部开放旁路] "
        << (wantTopCap ? "TopCap(法向向下)" : "BottomCap(法向向上)")
        << " targetZ=" << targetZ
        << " count=" << capCount << std::endl;
    return capFaces;
}

TopoDS_Compound GetAdjacentHorizontalCapFacesAtZ(
    const TopoDS_Shape& solid,
    const std::set<int>& sideFaceIds,
    const TopTools_DataMapOfShapeInteger& faceToIdMap,
    double targetZ,
    bool normalUp,
    const std::string& debugLabel,
    double zTol)
{
    BRep_Builder builder;
    TopoDS_Compound capFaces;
    builder.MakeCompound(capFaces);

    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    int capCount = 0;
    TopExp_Explorer exp(solid, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        BRepAdaptor_Surface surf(face);
        if (surf.GetType() != GeomAbs_Plane) continue;

        gp_Pln plane = surf.Plane();
        gp_Dir normal = plane.Axis().Direction();
        if (face.Orientation() == TopAbs_REVERSED) normal.Reverse();
        if (std::abs(std::abs(normal.Z()) - 1.0) > 1e-3) continue;

        double faceZ = plane.Location().Z();
        if (std::abs(faceZ - targetZ) > zTol) continue;
        if (normalUp && normal.Z() < 0.99) continue;
        if (!normalUp && normal.Z() > -0.99) continue;

        bool hasSideNeighbor = false;
        TopExp_Explorer edgeExp(face, TopAbs_EDGE);
        for (; edgeExp.More(); edgeExp.Next()) {
            const TopoDS_Edge& edge = TopoDS::Edge(edgeExp.Current());
            if (!edgeToFaces.Contains(edge)) continue;

            const TopTools_ListOfShape& adjacentFaces = edgeToFaces.FindFromKey(edge);
            TopTools_ListIteratorOfListOfShape it(adjacentFaces);
            for (; it.More(); it.Next()) {
                TopoDS_Face adjFace = TopoDS::Face(it.Value());
                if (adjFace.IsSame(face)) continue;
                if (!faceToIdMap.IsBound(adjFace)) continue;
                int adjId = faceToIdMap.Find(adjFace);
                if (sideFaceIds.count(adjId)) {
                    hasSideNeighbor = true;
                    break;
                }
            }
            if (hasSideNeighbor) break;
        }

        if (hasSideNeighbor) {
            builder.Add(capFaces, face);
            capCount++;
        }
    }

    std::cout << "  [中部开放旁路] " << debugLabel
        << " targetZ=" << targetZ
        << " normal=" << (normalUp ? "up" : "down")
        << " count=" << capCount << std::endl;
    return capFaces;
}

TopoDS_Compound GetHorizontalFacesForBand(
    const TopoDS_Shape& solid,
    const std::set<int>& sideFaceIds,
    const TopTools_DataMapOfShapeInteger& faceToIdMap,
    double bandMinZ,
    double bandMaxZ,
    const std::string& debugLabel,
    double zTol)
{
    BRep_Builder builder;
    TopoDS_Compound horizontalFaces;
    builder.MakeCompound(horizontalFaces);

    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    int candidateCount = 0;
    int keptCount = 0;
    int removedTopUp = 0;
    int removedBottomDown = 0;
    int skippedNoSideNeighbor = 0;

    TopExp_Explorer exp(solid, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        BRepAdaptor_Surface surf(face);
        if (surf.GetType() != GeomAbs_Plane) continue;

        gp_Pln plane = surf.Plane();
        gp_Dir normal = plane.Axis().Direction();
        if (face.Orientation() == TopAbs_REVERSED) normal.Reverse();
        if (std::abs(std::abs(normal.Z()) - 1.0) > 1e-3) continue;

        double faceZ = plane.Location().Z();
        if (faceZ < bandMinZ - zTol || faceZ > bandMaxZ + zTol) continue;
        candidateCount++;

        bool hasSideNeighbor = false;
        TopExp_Explorer edgeExp(face, TopAbs_EDGE);
        for (; edgeExp.More(); edgeExp.Next()) {
            const TopoDS_Edge& edge = TopoDS::Edge(edgeExp.Current());
            if (!edgeToFaces.Contains(edge)) continue;

            const TopTools_ListOfShape& adjacentFaces = edgeToFaces.FindFromKey(edge);
            TopTools_ListIteratorOfListOfShape it(adjacentFaces);
            for (; it.More(); it.Next()) {
                TopoDS_Face adjFace = TopoDS::Face(it.Value());
                if (adjFace.IsSame(face)) continue;
                if (!faceToIdMap.IsBound(adjFace)) continue;
                int adjId = faceToIdMap.Find(adjFace);
                if (sideFaceIds.count(adjId)) {
                    hasSideNeighbor = true;
                    break;
                }
            }
            if (hasSideNeighbor) break;
        }

        if (!hasSideNeighbor) {
            skippedNoSideNeighbor++;
            continue;
        }

        bool isBandTop = std::abs(faceZ - bandMaxZ) <= zTol;
        bool isBandBottom = std::abs(faceZ - bandMinZ) <= zTol;
        bool normalUp = normal.Z() > 0.99;
        bool normalDown = normal.Z() < -0.99;

        if (isBandTop && !normalDown) {
            removedTopUp++;
            continue;
        }
        if (isBandBottom && !normalUp) {
            removedBottomDown++;
            continue;
        }

        builder.Add(horizontalFaces, face);
        keptCount++;
    }

    std::cout << "  [中部开放旁路] " << debugLabel
        << " band=[" << bandMinZ << ", " << bandMaxZ << "]"
        << " candidates=" << candidateCount
        << " kept=" << keptCount
        << " removedTopUp=" << removedTopUp
        << " removedBottomDown=" << removedBottomDown
        << " skippedNoSideNeighbor=" << skippedNoSideNeighbor
        << std::endl;
    return horizontalFaces;
}


void GetCompoundXYRange(const TopoDS_Compound& compound, double& xmin, double& ymin, double& xmax, double& ymax) {
    Bnd_Box box;
    BRepBndLib::Add(compound, box);
    if (box.IsVoid()) {
        xmin = ymin = xmax = ymax = 0.0;
        return;
    }
    double zmin, zmax;
    box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
}

bool IsFaceInsideXYBox(const TopoDS_Face& face, double xmin, double ymin, double xmax, double ymax) {
    Bnd_Box faceBox;
    BRepBndLib::Add(face, faceBox);
    if (faceBox.IsVoid()) {
        return false;
    }
    double fxmin, fymin, fzmin, fxmax, fymax, fzmax;
    faceBox.Get(fxmin, fymin, fzmin, fxmax, fymax, fzmax);
    return (fxmin >= xmin && fxmax <= xmax && fymin >= ymin && fymax <= ymax);
}

bool IsFaceCrossingXYBox(const TopoDS_Face& face, double xmin, double ymin, double xmax, double ymax) {
    Bnd_Box faceBox;
    BRepBndLib::Add(face, faceBox);
    if (faceBox.IsVoid()) {
        return false;
    }
    double fxmin, fymin, fzmin, fxmax, fymax, fzmax;
    faceBox.Get(fxmin, fymin, fzmin, fxmax, fymax, fzmax);
    bool xOverlap = (fxmin < xmax && fxmax > xmin);
    bool yOverlap = (fymin < ymax && fymax > ymin);
    if (!xOverlap || !yOverlap) {
        return false;
    }
    bool xInside = (fxmin >= xmin && fxmax <= xmax);
    bool yInside = (fymin >= ymin && fymax <= ymax);
    return !(xInside && yInside);
}

void SplitFacesByXYBox(
    const TopoDS_Compound& inputFaces,
    double xmin, double ymin, double xmax, double ymax,
    TopoDS_Compound& insideFaces,
    TopoDS_Compound& crossingFaces) {
    BRep_Builder builder;
    builder.MakeCompound(insideFaces);
    builder.MakeCompound(crossingFaces);

    TopExp_Explorer explorer(inputFaces, TopAbs_FACE);
    for (; explorer.More(); explorer.Next()) {
        TopoDS_Face face = TopoDS::Face(explorer.Current());
        if (IsFaceInsideXYBox(face, xmin, ymin, xmax, ymax)) {
            builder.Add(insideFaces, face);
        } else if (IsFaceCrossingXYBox(face, xmin, ymin, xmax, ymax)) {
            builder.Add(crossingFaces, face);
        } else {
            builder.Add(crossingFaces, face);
        }
    }
}


std::set<int> CollectFaceIdsFromFace2D(const Face2D& face) {
    std::set<int> ids;
    for (const auto& line : face.outerLoop) {
        if (line.faceId >= 0) ids.insert(line.faceId);
    }
    for (const auto& inner : face.innerLoops) {
        for (const auto& line : inner) {
            if (line.faceId >= 0) ids.insert(line.faceId);
        }
    }
    return ids;
}

bool IsFaceOverlappingFace2D(const TopoDS_Face& face3D, const Face2D& region2D) {
    Bnd_Box faceBox;
    BRepBndLib::Add(face3D, faceBox);
    if (faceBox.IsVoid()) return false;

    double fxmin, fymin, fzmin, fxmax, fymax, fzmax;
    faceBox.Get(fxmin, fymin, fzmin, fxmax, fymax, fzmax);

    double cx = (fxmin + fxmax) / 2.0;
    double cy = (fymin + fymax) / 2.0;
    double z = (fzmin + fzmax) / 2.0;

    double xSpan = fxmax - fxmin;
    double ySpan = fymax - fymin;

    double halfW = xSpan / 2.0;
    double halfH = ySpan / 2.0;
    double minHalf = 0.05;
    if (halfW < minHalf) halfW = minHalf;
    if (halfH < minHalf) halfH = minHalf;

    double bxMin = cx - halfW;
    double bxMax = cx + halfW;
    double byMin = cy - halfH;
    double byMax = cy + halfH;

    Face2D face2DBox;
    face2DBox.type = FaceType::SOLID;
    face2DBox.outerLoop = {
        {{bxMin, byMin, z}, {bxMax, byMin, z}, -1},
        {{bxMax, byMin, z}, {bxMax, byMax, z}, -1},
        {{bxMax, byMax, z}, {bxMin, byMax, z}, -1},
        {{bxMin, byMax, z}, {bxMin, byMin, z}, -1}
    };

    std::vector<Face2D> intersection = BooleanFacesSingle(face2DBox, region2D, ClipType::Intersection, z);
    return !intersection.empty();
}

void SplitFacesByClipperRegion(
    const TopoDS_Compound& inputFaces,
    const Face2D& machinableRegion,
    TopoDS_Compound& machinableFaces,
    TopoDS_Compound& nonMachinableFaces) {
    BRep_Builder builder;
    builder.MakeCompound(machinableFaces);
    builder.MakeCompound(nonMachinableFaces);

    TopExp_Explorer explorer(inputFaces, TopAbs_FACE);
    for (; explorer.More(); explorer.Next()) {
        TopoDS_Face face = TopoDS::Face(explorer.Current());
        if (IsFaceOverlappingFace2D(face, machinableRegion)) {
            builder.Add(machinableFaces, face);
        } else {
            builder.Add(nonMachinableFaces, face);
        }
    }
}

std::set<int> CollectFaceIdsFromCompound(
    const TopoDS_Compound& compound,
    const TopTools_DataMapOfShapeInteger& faceToIdMap) {
    std::set<int> ids;
    TopExp_Explorer exp(compound, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        if (faceToIdMap.IsBound(face)) {
            ids.insert(faceToIdMap.Find(face));
        }
    }
    return ids;
}

void SplitInteriorOpenSourceFacesByBand(
    const TopoDS_Compound& sourceFaces,
    double topZ,
    double bottomZ,
    TopoDS_Compound& upperClosed,
    TopoDS_Compound& middleOpen,
    TopoDS_Compound& lowerClosed)
{
    auto topSplit = SplitCompoundAtZ(sourceFaces, topZ);
    upperClosed = topSplit.first;

    auto bottomSplit = SplitCompoundAtZ(topSplit.second, bottomZ);
    middleOpen = bottomSplit.first;
    lowerClosed = bottomSplit.second;
}

void ProcessInteriorOpenCavityByBand(
    const std::vector<std::vector<Face2D>>& allClosedLayerFaces,
    const std::vector<std::vector<Face2D>>& allOpenLayerFaces,
    const TopoDS_Shape& mainShape,
    const std::vector<double>& splitPoints,
    const TopTools_DataMapOfShapeInteger& faceToIdMap,
    const std::string& savePath,
    std::vector<TopoDS_Compound>& outNewClosedParts,
    TopoDS_Compound& outNewOpenCavity,
    gp_Dir& outOpenToolDir)
{
    std::cout << "\n===== [中部开放型腔旁路] 开始 =====" << std::endl;

    double modelMinZ = 0.0, modelMaxZ = 0.0;
    GetShapeZRange(mainShape, modelMinZ, modelMaxZ);

    CavityFeature openFeature =
        BuildInteriorOpenCavityFeature(allOpenLayerFaces, allClosedLayerFaces, splitPoints, modelMinZ, savePath);

    std::cout << "  [中部开放旁路] openFeature.topZ=" << openFeature.topZ
        << " bottomZ=" << openFeature.bottomZ
        << " sourceFaceIds=" << openFeature.sourceFaceIds.size() << std::endl;
    std::cout << "  [中部开放旁路] sourceFaceIds: ";
    for (int id : openFeature.sourceFaceIds) std::cout << id << " ";
    std::cout << std::endl;

    if (openFeature.stepLoops.empty() || openFeature.sourceFaceIds.empty()) {
        std::cout << "  [中部开放旁路] 开放feature为空，跳过。" << std::endl;
        return;
    }

    TopoDS_Compound sourceFaces =
        ExtractFacesByFeatureSourceIds(openFeature, mainShape, faceToIdMap);
    std::string sourcePath = savePath + "InteriorOpen_SourceSideFaces.brep";
    BRepTools::Write(sourceFaces, sourcePath.c_str());

    TopoDS_Compound upperClosed;
    TopoDS_Compound middleOpen;
    TopoDS_Compound lowerClosed;
    SplitInteriorOpenSourceFacesByBand(
        sourceFaces,
        openFeature.topZ,
        openFeature.bottomZ,
        upperClosed,
        middleOpen,
        lowerClosed);

    double upperZMin = 0.0, upperZMax = 0.0;
    double middleZMin = 0.0, middleZMax = 0.0;
    double lowerZMin = 0.0, lowerZMax = 0.0;
    bool hasUpperRange = GetCompoundZRange(upperClosed, upperZMin, upperZMax);
    bool hasMiddleRange = GetCompoundZRange(middleOpen, middleZMin, middleZMax);
    bool hasLowerRange = GetCompoundZRange(lowerClosed, lowerZMin, lowerZMax);

    std::set<int> upperSideIds = CollectFaceIdsFromCompound(upperClosed, faceToIdMap);
    std::set<int> middleSideIds = CollectFaceIdsFromCompound(middleOpen, faceToIdMap);
    std::set<int> lowerSideIds = CollectFaceIdsFromCompound(lowerClosed, faceToIdMap);

    std::cout << "  [DEBUG-水平面] upperSideIds(size=" << upperSideIds.size() << "): ";
    for (int id : upperSideIds) std::cout << id << " ";
    std::cout << std::endl;
    std::cout << "  [DEBUG-水平面] middleSideIds(size=" << middleSideIds.size() << "): ";
    for (int id : middleSideIds) std::cout << id << " ";
    std::cout << std::endl;
    std::cout << "  [DEBUG-水平面] lowerSideIds(size=" << lowerSideIds.size() << "): ";
    for (int id : lowerSideIds) std::cout << id << " ";
    std::cout << std::endl;

    auto debugHorizontalFaces = [&](const TopoDS_Compound& compound, const std::string& label) {
        TopExp_Explorer exp(compound, TopAbs_FACE);
        for (; exp.More(); exp.Next()) {
            TopoDS_Face f = TopoDS::Face(exp.Current());
            BRepAdaptor_Surface surf(f);
            if (surf.GetType() != GeomAbs_Plane) continue;
            gp_Pln plane = surf.Plane();
            gp_Dir normal = plane.Axis().Direction();
            if (f.Orientation() == TopAbs_REVERSED) normal.Reverse();
            if (std::abs(std::abs(normal.Z()) - 1.0) > 1e-3) continue;
            double faceZ = plane.Location().Z();
            double area = CalculateFaceAreaOCC(f);
            int faceId = faceToIdMap.IsBound(f) ? faceToIdMap.Find(f) : -1;
            std::cout << "  [DEBUG-" << label << "] 水平面 ID=" << faceId
                << " Z=" << faceZ
                << " 法向=(" << normal.X() << "," << normal.Y() << "," << normal.Z() << ")"
                << " 面积=" << area << std::endl;
        }
    };

    debugHorizontalFaces(upperClosed, "上区水平面");
    debugHorizontalFaces(middleOpen, "中区水平面");
    debugHorizontalFaces(lowerClosed, "下区水平面");

    TopoDS_Compound upperHorizontalFaces = GetHorizontalFacesForBand(
        mainShape,
        upperSideIds.empty() ? openFeature.sourceFaceIds : upperSideIds,
        faceToIdMap,
        hasUpperRange ? upperZMin : openFeature.topZ,
        hasUpperRange ? upperZMax : openFeature.topZ,
        "UpperClosedHorizontalFaces");
    TopoDS_Compound middleHorizontalFaces = GetHorizontalFacesForBand(
        mainShape,
        middleSideIds.empty() ? openFeature.sourceFaceIds : middleSideIds,
        faceToIdMap,
        hasMiddleRange ? middleZMin : openFeature.bottomZ,
        hasMiddleRange ? middleZMax : openFeature.topZ,
        "MiddleOpenHorizontalFaces");
    TopoDS_Compound lowerHorizontalFaces = GetHorizontalFacesForBand(
        mainShape,
        lowerSideIds.empty() ? openFeature.sourceFaceIds : lowerSideIds,
        faceToIdMap,
        hasLowerRange ? lowerZMin : openFeature.bottomZ,
        hasLowerRange ? lowerZMax : openFeature.bottomZ,
        "LowerClosedHorizontalFaces");

    AddCompoundFaces(upperClosed, upperHorizontalFaces);
    AddCompoundFaces(middleOpen, middleHorizontalFaces);
    AddCompoundFaces(lowerClosed, lowerHorizontalFaces);

    // ================= 中间开放型腔：封闭型腔轮廓作为可加工区域边界 =================
    double topOpenZ = -1e9;
    for (const auto& layer : allOpenLayerFaces) {
        for (const auto& face : layer) {
            if (face.type != FaceType::OPENCAVITY || face.outerLoop.empty()) continue;
            double z = face.outerLoop.front().start.z;
            if (z > topOpenZ) topOpenZ = z;
        }
    }
    std::cout << "  [DEBUG-分离] 最高开放型腔Z=" << topOpenZ << std::endl;

    Face2D closedClipFace;
    bool foundClosedClip = false;
    int closedClipLayerIdx = -1;
    for (size_t layerIdx = 0; layerIdx < allClosedLayerFaces.size(); ++layerIdx) {
        for (const auto& face : allClosedLayerFaces[layerIdx]) {
            if (face.type != FaceType::CAVITY || face.outerLoop.empty()) continue;
            double z = face.outerLoop.front().start.z;
            if (z > topOpenZ) {
                closedClipFace = face;
                foundClosedClip = true;
                closedClipLayerIdx = (int)layerIdx;
                break;
            }
        }
        if (foundClosedClip) break;
    }

    if (foundClosedClip) {
        double cxMin = 1e9, cyMin = 1e9, cxMax = -1e9, cyMax = -1e9;
        for (const auto& line : closedClipFace.outerLoop) {
            cxMin = std::min(cxMin, std::min(line.start.x, line.end.x));
            cyMin = std::min(cyMin, std::min(line.start.y, line.end.y));
            cxMax = std::max(cxMax, std::max(line.start.x, line.end.x));
            cyMax = std::max(cyMax, std::max(line.start.y, line.end.y));
        }
        std::cout << "  [DEBUG-分离] 封闭型腔边界切片 layerIdx=" << closedClipLayerIdx
            << " Z=" << closedClipFace.outerLoop.front().start.z
            << " outerLoop边数=" << closedClipFace.outerLoop.size()
            << " innerLoops数=" << closedClipFace.innerLoops.size() << std::endl;
        std::cout << "  [DEBUG-分离] 封闭型腔边界 XY范围: ("
            << cxMin << "," << cyMin << ") -> (" << cxMax << "," << cyMax << ")" << std::endl;
        std::string closedClipFacePath = savePath + "DEBUG_closedClipFace.brep";
        ExportFace2DToBrep({closedClipFace}, closedClipFacePath);
    } else {
        std::cout << "  [DEBUG-分离] 未找到Z > topOpenZ的封闭型腔切片！" << std::endl;
    }

    TopoDS_Compound machinableFaces;
    TopoDS_Compound nonMachinableFaces;
    BRep_Builder splitBuilder;
    splitBuilder.MakeCompound(machinableFaces);
    splitBuilder.MakeCompound(nonMachinableFaces);

    if (foundClosedClip) {
        int faceIdx = 0;
        TopExp_Explorer explorer(middleOpen, TopAbs_FACE);
        for (; explorer.More(); explorer.Next(), ++faceIdx) {
            TopoDS_Face face = TopoDS::Face(explorer.Current());
            Bnd_Box faceBox;
            BRepBndLib::Add(face, faceBox);
            double fxmin, fymin, fzmin, fxmax, fymax, fzmax;
            faceBox.Get(fxmin, fymin, fzmin, fxmax, fymax, fzmax);
            bool isMachinable = IsFaceOverlappingFace2D(face, closedClipFace);
            std::cout << "    [DEBUG-分离] middleOpen面[" << faceIdx
                << "] Z=(" << fzmin << "," << fzmax
                << ") XY=(" << fxmin << "," << fymin << ")->(" << fxmax << "," << fymax
                << ") -> " << (isMachinable ? "可加工" : "不可加工") << std::endl;
            if (isMachinable) {
                splitBuilder.Add(machinableFaces, face);
            } else {
                splitBuilder.Add(nonMachinableFaces, face);
            }
        }
    } else {
        std::cout << "  [DEBUG-分离] 无封闭型腔边界，全部归为不可加工。" << std::endl;
        AddCompoundFaces(nonMachinableFaces, middleOpen);
    }

    std::string machinablePath = savePath + "InteriorOpen_MachinableFaces.brep";
    std::string nonMachinablePath = savePath + "InteriorOpen_NonMachinableFaces.brep";
    BRepTools::Write(machinableFaces, machinablePath.c_str());
    BRepTools::Write(nonMachinableFaces, nonMachinablePath.c_str());

    std::cout << "  [可加工区域分离] MiddleOpen总面数=" << CountFacesInCompound(middleOpen) << std::endl;
    std::cout << "  [可加工区域分离] 可加工区域面数=" << CountFacesInCompound(machinableFaces)
        << " -> " << machinablePath << std::endl;
    std::cout << "  [可加工区域分离] 不可加工区域面数=" << CountFacesInCompound(nonMachinableFaces)
        << " -> " << nonMachinablePath << std::endl;

    BRep_Builder reclassifiedBuilder;
    TopoDS_Compound newClosedCavity;
    reclassifiedBuilder.MakeCompound(newClosedCavity);
    AddCompoundFaces(newClosedCavity, upperClosed);
    AddCompoundFaces(newClosedCavity, upperHorizontalFaces);
    AddCompoundFaces(newClosedCavity, lowerClosed);
    AddCompoundFaces(newClosedCavity, lowerHorizontalFaces);
    AddCompoundFaces(newClosedCavity, machinableFaces);

    TopoDS_Compound newOpenCavity;
    reclassifiedBuilder.MakeCompound(newOpenCavity);
    AddCompoundFaces(newOpenCavity, nonMachinableFaces);

    std::string newClosedPath = savePath + "InteriorOpen_NewClosedCavity.brep";
    std::string newOpenPath = savePath + "InteriorOpen_NewOpenCavity.brep";
    BRepTools::Write(newClosedCavity, newClosedPath.c_str());
    BRepTools::Write(newOpenCavity, newOpenPath.c_str());

    std::cout << "  [重分类结果] 新封闭型腔面数=" << CountFacesInCompound(newClosedCavity)
        << " -> " << newClosedPath << std::endl;
    std::cout << "  [重分类结果] 新开放型腔面数=" << CountFacesInCompound(newOpenCavity)
        << " -> " << newOpenPath << std::endl;

    gp_Dir openToolDir;
    bool openMachinable = IsCavityMachinable(newOpenCavity, openToolDir);
    std::cout << "  [重分类结果] 不可加工区域加工可行性: -> "
        << (openMachinable ? "可行" : "不可行")
        << " 进刀方向=(" << openToolDir.X() << "," << openToolDir.Y() << "," << openToolDir.Z() << ")" << std::endl;
    outOpenToolDir = openToolDir;

    // ================= 新封闭型腔：沿Z轴分割得到简单型腔 =================
    std::cout << "\n  [新封闭型腔Z轴分割] 开始对新封闭型腔进行Z轴递推切割..." << std::endl;
    std::vector<TopoDS_Compound> newClosedParts = RecursiveSplitCavity(newClosedCavity);
    std::cout << "  [新封闭型腔Z轴分割] 切割完成，共得到 " << newClosedParts.size() << " 个简单型腔。" << std::endl;

    for (size_t partIdx = 0; partIdx < newClosedParts.size(); ++partIdx) {
        std::string partPath = savePath + "InteriorOpen_NewClosedCavity_Part_" + std::to_string(partIdx) + ".brep";
        BRepTools::Write(newClosedParts[partIdx], partPath.c_str());
        std::cout << "    - Part[" << partIdx << "] faces=" << CountFacesInCompound(newClosedParts[partIdx])
            << " -> " << partPath << std::endl;
    }

    std::string upperPath = savePath + "InteriorOpen_UpperClosedSideFaces.brep";
    std::string middlePath = savePath + "InteriorOpen_MiddleOpenSideFaces.brep";
    std::string lowerPath = savePath + "InteriorOpen_LowerClosedSideFaces.brep";
    std::string upperHorizontalPath = savePath + "InteriorOpen_UpperClosedHorizontalFaces.brep";
    std::string middleHorizontalPath = savePath + "InteriorOpen_MiddleOpenHorizontalFaces.brep";
    std::string lowerHorizontalPath = savePath + "InteriorOpen_LowerClosedHorizontalFaces.brep";
    BRepTools::Write(upperHorizontalFaces, upperHorizontalPath.c_str());
    BRepTools::Write(middleHorizontalFaces, middleHorizontalPath.c_str());
    BRepTools::Write(lowerHorizontalFaces, lowerHorizontalPath.c_str());
    BRepTools::Write(upperClosed, upperPath.c_str());
    BRepTools::Write(middleOpen, middlePath.c_str());
    BRepTools::Write(lowerClosed, lowerPath.c_str());

    std::cout << "  [中部开放旁路] Source faces=" << CountFacesInCompound(sourceFaces) << std::endl;
    std::cout << "  [中部开放旁路] UpperClosed faces=" << CountFacesInCompound(upperClosed)
        << " -> " << upperPath << std::endl;
    std::cout << "  [中部开放旁路] MiddleOpen faces=" << CountFacesInCompound(middleOpen)
        << " -> " << middlePath << std::endl;
    std::cout << "  [中部开放旁路] LowerClosed faces=" << CountFacesInCompound(lowerClosed)
        << " -> " << lowerPath << std::endl;

    outNewClosedParts = newClosedParts;
    outNewOpenCavity = newOpenCavity;

    std::cout << "===== [中部开放型腔旁路] 结束 =====\n" << std::endl;
}
