#include "occ_utils.h"
#include "geom_utils.h"

// --- 平面/边分析 ---
TopoDS_Compound GetCoplanarFaces(const TopoDS_Shape& shape, double splitZ, double tol) {
    BRep_Builder builder;
    TopoDS_Compound comp;
    builder.MakeCompound(comp);
    
    const gp_Dir zAxis(0, 0, 1);
    
    TopExp_Explorer exp(shape, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        BRepAdaptor_Surface surf(face);
        
        // 必须是平面
        if (surf.GetType() != GeomAbs_Plane) continue;
        
        gp_Pln plane = surf.Plane();
        gp_Dir normal = plane.Axis().Direction();
        
        // 判断法向量是否平行于 Z 轴（即水平面）
        double angle = normal.Angle(zAxis);
        if (angle < tol || fabs(angle - M_PI) < tol) {
            // 判断平面的精确 Z 高度是否等于 splitZ
            double z = plane.Location().Z();
            if (std::abs(z - splitZ) < tol) {
                builder.Add(comp, face);
            }
        }
    }
    return comp;
}

// 判断边是否为凹边（仅用于非光顺边/硬棱边）
// 用 crossN·T 几何判断：crossN = N1×N2
bool IsConcave(const TopoDS_Edge& E, const TopoDS_Face& F1, const TopoDS_Face& F2) {
    Standard_Real f, l;
    Handle(Geom_Curve) C = BRep_Tool::Curve(E, f, l);
    if (C.IsNull()) return false;

    gp_Pnt P;
    gp_Vec T;
    C->D1((f + l) / 2.0, P, T);

    // 找到边在 F1 中的方向
    TopExp_Explorer expE(F1, TopAbs_EDGE);
    TopoDS_Edge E_in_F1;
    bool found = false;
    for (; expE.More(); expE.Next()) {
        if (expE.Current().IsSame(E)) {
            E_in_F1 = TopoDS::Edge(expE.Current());
            found = true;
            break;
        }
    }
    if (!found) return false;

    if (E_in_F1.Orientation() == TopAbs_REVERSED) T.Reverse();

    // 获取 F1 的法向
    BRepAdaptor_Surface AS1(F1, Standard_True);
    Standard_Real u1, v1;
    GeomAPI_ProjectPointOnSurf proj1(P, AS1.Surface().Surface());
    if (!proj1.IsDone()) return false;
    proj1.LowerDistanceParameters(u1, v1);
    gp_Pnt p1; gp_Vec d1u, d1v;
    AS1.D1(u1, v1, p1, d1u, d1v);
    gp_Vec N1 = d1u.Crossed(d1v);
    if (F1.Orientation() == TopAbs_REVERSED) N1.Reverse();
    if (N1.Magnitude() > 1e-7) N1.Normalize();

    // 获取 F2 的法向
    BRepAdaptor_Surface AS2(F2, Standard_True);
    Standard_Real u2, v2;
    GeomAPI_ProjectPointOnSurf proj2(P, AS2.Surface().Surface());
    if (!proj2.IsDone()) return false;
    proj2.LowerDistanceParameters(u2, v2);
    gp_Pnt p2; gp_Vec d2u, d2v;
    AS2.D1(u2, v2, p2, d2u, d2v);
    gp_Vec N2 = d2u.Crossed(d2v);
    if (F2.Orientation() == TopAbs_REVERSED) N2.Reverse();
    if (N2.Magnitude() > 1e-7) N2.Normalize();

    // 计算法向叉乘
    gp_Vec crossN = N1.Crossed(N2);
    double dot = crossN.Dot(T);

    // 如果点乘小于0，说明是凹边（材料内部角度 > 180度）
    return dot < -1e-5;
}

// 获得法向量与进刀方向相同的平面的凹边
TopoDS_Compound GetConcaveEdges(const TopoDS_Shape& solid, const TopoDS_Compound& coplanarFaces, const gp_Dir& toolDirection) {
    BRep_Builder builder;
    TopoDS_Compound concaveEdges;
    builder.MakeCompound(concaveEdges);

    // 建立边到面的映射关系
    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    // 先建立面方向的快速查找：对于 coplanarFaces 中的每个面，
    // 如果法向量与进刀方向相同，才检查它的凹边
    TopTools_IndexedMapOfShape targetFaces;
    TopExp_Explorer faceExp(coplanarFaces, TopAbs_FACE);
    for (; faceExp.More(); faceExp.Next()) {
        TopoDS_Face face = TopoDS::Face(faceExp.Current());
        BRepAdaptor_Surface surf(face);
        if (surf.GetType() != GeomAbs_Plane) continue;
        
        gp_Pln plane = surf.Plane();
        gp_Dir normal = plane.Axis().Direction();
        if (face.Orientation() == TopAbs_REVERSED) normal.Reverse();
        
        // 法向量与进刀方向点积 > 0 → 法向量指向进刀方向 → 处理这个面的凹边
        if (normal.Dot(toolDirection) > 0) {
            targetFaces.Add(face);
        }
    }

    // 遍历 targetFaces 的所有边，收集凹边
    for (int i = 1; i <= targetFaces.Extent(); i++) {
        const TopoDS_Shape& shapeRef = targetFaces(i);
        TopoDS_Face face = TopoDS::Face(shapeRef);
        TopExp_Explorer edgeExp(face, TopAbs_EDGE);
        for (; edgeExp.More(); edgeExp.Next()) {
            TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());
            
            const TopTools_ListOfShape& faceList = edgeToFaces.FindFromKey(edge);
            if (faceList.Extent() == 2) {
                TopoDS_Face f1 = TopoDS::Face(faceList.First());
                TopoDS_Face f2 = TopoDS::Face(faceList.Last());
                
                if (IsConcave(edge, f1, f2)) {
                    builder.Add(concaveEdges, edge);
                }
            }
        }
    }
    
    return concaveEdges;
}

// 判断边是否为凸边
bool IsConvex(const TopoDS_Edge& E, const TopoDS_Face& F1, const TopoDS_Face& F2) {
    Standard_Real f, l;
    Handle(Geom_Curve) C = BRep_Tool::Curve(E, f, l);
    if (C.IsNull()) return false;

    gp_Pnt P;
    gp_Vec T;
    C->D1((f + l) / 2.0, P, T);

    // 找到边在 F1 中的方向
    TopExp_Explorer expE(F1, TopAbs_EDGE);
    TopoDS_Edge E_in_F1;
    bool found = false;
    for (; expE.More(); expE.Next()) {
        if (expE.Current().IsSame(E)) {
            E_in_F1 = TopoDS::Edge(expE.Current());
            found = true;
            break;
        }
    }
    if (!found) return false;

    if (E_in_F1.Orientation() == TopAbs_REVERSED) T.Reverse();

    // 获取 F1 的法向
    BRepAdaptor_Surface AS1(F1, Standard_True);
    Standard_Real u1, v1;
    GeomAPI_ProjectPointOnSurf proj1(P, AS1.Surface().Surface());
    if (!proj1.IsDone()) return false;
    proj1.LowerDistanceParameters(u1, v1);
    gp_Pnt p1; gp_Vec d1u, d1v;
    AS1.D1(u1, v1, p1, d1u, d1v);
    gp_Vec N1 = d1u.Crossed(d1v);
    if (F1.Orientation() == TopAbs_REVERSED) N1.Reverse();
    if (N1.Magnitude() > 1e-7) N1.Normalize();

    // 获取 F2 的法向
    BRepAdaptor_Surface AS2(F2, Standard_True);
    Standard_Real u2, v2;
    GeomAPI_ProjectPointOnSurf proj2(P, AS2.Surface().Surface());
    if (!proj2.IsDone()) return false;
    proj2.LowerDistanceParameters(u2, v2);
    gp_Pnt p2; gp_Vec d2u, d2v;
    AS2.D1(u2, v2, p2, d2u, d2v);
    gp_Vec N2 = d2u.Crossed(d2v);
    if (F2.Orientation() == TopAbs_REVERSED) N2.Reverse();
    if (N2.Magnitude() > 1e-7) N2.Normalize();

    // 计算法向叉乘
    gp_Vec crossN = N1.Crossed(N2);
    double dot = crossN.Dot(T);

    // 凸边判断：点乘 > 0（与凹边相反）
    return dot > 1e-5;
}

// 获得法向量与进刀方向相反的平面的凸边
TopoDS_Compound GetConvexEdgesFromOppositeFaces(
    const TopoDS_Shape& solid,
    const TopoDS_Compound& coplanarFaces,
    const gp_Dir& toolDirection)
{
    BRep_Builder builder;
    TopoDS_Compound convexEdges;
    builder.MakeCompound(convexEdges);

    // 边 -> 面映射（共用边）
    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    // 找出：法向量 与 进刀方向 相反 的平面
    TopTools_IndexedMapOfShape oppositeFaces;
    TopExp_Explorer faceExp(coplanarFaces, TopAbs_FACE);
    for (; faceExp.More(); faceExp.Next()) {
        TopoDS_Face face = TopoDS::Face(faceExp.Current());
        BRepAdaptor_Surface surf(face);
        if (surf.GetType() != GeomAbs_Plane) continue;

        gp_Pln plane = surf.Plane();
        gp_Dir normal = plane.Axis().Direction();
        if (face.Orientation() == TopAbs_REVERSED)
            normal.Reverse();

        // 关键：法向量与进刀方向 相反（点积 < 0）
        if (normal.Dot(toolDirection) < 0) {
            oppositeFaces.Add(face);
        }
    }

    // 遍历这些面，提取 凸边
    for (int i = 1; i <= oppositeFaces.Extent(); ++i) {
        const TopoDS_Shape& shapeRef = oppositeFaces(i);
        TopoDS_Face face = TopoDS::Face(shapeRef);
        TopExp_Explorer edgeExp(face, TopAbs_EDGE);

        for (; edgeExp.More(); edgeExp.Next()) {
            TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());

            const TopTools_ListOfShape& faceList = edgeToFaces.FindFromKey(edge);
            if (faceList.Extent() == 2) {
                TopoDS_Face f1 = TopoDS::Face(faceList.First());
                TopoDS_Face f2 = TopoDS::Face(faceList.Last());

                // 关键：判断为 凸边 才加入
                if (IsConvex(edge, f1, f2)) {
                    builder.Add(convexEdges, edge);
                }
            }
        }
    }

    return convexEdges;
}

// 统一清洗函数：对 coplanarFaces 中的所有边进行分类，
/**
 * @brief 识别并获取需要移除的边缘
 * 逻辑：
 * 1. 进刀方向 toolDirection (通常为 0,0,-1)
 * 2. 朝向刀具的面 (isFacingTool, Dot < 0): 移除凹边 (isConcave)
 * 3. 背向刀具的面 (isFacingTool == false): 移除凸边 (!isConcave)
 */
TopoDS_Compound ClassifyAndGetEdgesToRemove(
    const TopoDS_Shape& solid,
    const TopoDS_Compound& coplanarFaces,
    const gp_Dir& toolDirection)
{
    BRep_Builder builder;
    TopoDS_Compound result;
    builder.MakeCompound(result);

    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    TopExp_Explorer faceExp(coplanarFaces, TopAbs_FACE);
    for (; faceExp.More(); faceExp.Next()) {
        TopoDS_Face face = TopoDS::Face(faceExp.Current());
        BRepAdaptor_Surface surf(face);
        if (surf.GetType() != GeomAbs_Plane) continue;

        gp_Pln plane = surf.Plane();
        double splitZ = plane.Location().Z();
        gp_Dir planeNormal = plane.Axis().Direction();
        if (face.Orientation() == TopAbs_REVERSED) planeNormal.Reverse();

        // 进刀(0,0,-1) 与 槽底法向(0,0,1) 点积为负 -> 朝向刀具的面
        bool isFacingTool = (planeNormal.Dot(toolDirection) < 0);

        TopExp_Explorer edgeExp(face, TopAbs_EDGE);
        for (; edgeExp.More(); edgeExp.Next()) {
            TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());

            if (!edgeToFaces.Contains(edge)) continue;
            const TopTools_ListOfShape& faceList = edgeToFaces.FindFromKey(edge);
            if (faceList.Extent() != 2) continue;

            TopoDS_Face f1 = TopoDS::Face(faceList.First());
            TopoDS_Face f2 = TopoDS::Face(faceList.Last());
            TopoDS_Face neighborFace = f1.IsSame(face) ? f2 : f1;

            bool isConcave = false;
            GeomAbs_Shape continuity = BRep_Tool::Continuity(edge, face, neighborFace);

            if (continuity >= GeomAbs_G1) {
                // --- 【策略 A：光滑边使用中点高度判定法】 ---
                // 光滑边（如圆角面）的 BndBox 中点能很好地反映其相对于基准面的偏移方向
                Bnd_Box bbox;
                BRepBndLib::Add(neighborFace, bbox);
                Standard_Real xmin, ymin, zmin, xmax, ymax, zmax;
                bbox.Get(xmin, ymin, zmin, xmax, ymax, zmax);

                double zCenter = (zmin + zmax) * 0.5;
                bool isUpwards = (zCenter > splitZ);
                double Nz = planeNormal.Z();

                if (isUpwards) {
                    isConcave = (Nz > 0); // 向上弯曲且面向上 -> 凹圆角
                }
                else {
                    isConcave = (Nz < 0); // 向下弯曲且面向下 -> 凹圆角
                }
            }
            else {
                // --- 【策略 B：硬棱边使用传统向量判断】 ---
                // 使用你提供的 IsConcave 逻辑进行精确计算
                isConcave = IsConcave(edge, f1, f2);
            }

            // --- 策略：底面删凹边，顶面删凸边 ---
            // isFacingTool (true)  -> shouldRemove = isConcave
            // isFacingTool (false) -> shouldRemove = !isConcave
            bool shouldRemove = isFacingTool ? isConcave : !isConcave;

            if (shouldRemove) {
                builder.Add(result, edge);
            }
        }
    }

    return result;
}
// 3. 线段布尔减 (A - B) 并保留来源面映射
std::vector<OneEdge> SubtractLinesAndMapFaces(
    const TopoDS_Compound& linesA, 
    const TopoDS_Compound& linesB,
    const TopTools_DataMapOfShapeShape& edgeToFaceMap) 
{
    std::vector<OneEdge> result;
    BRepAlgoAPI_Cut cutAlgo(linesA, linesB);
    cutAlgo.Build();
    
    if (!cutAlgo.IsDone()) return result;
    
    TopoDS_Shape cutShape = cutAlgo.Shape();
    
    // 收集布尔减后最终存在的所有边
    TopTools_IndexedMapOfShape resultEdgesMap;
    TopExp::MapShapes(cutShape, TopAbs_EDGE, resultEdgesMap);
    
    // 遍历原始交线，追踪它们在布尔减之后的变化
    TopExp_Explorer expOrig(linesA, TopAbs_EDGE);
    for (; expOrig.More(); expOrig.Next()) {
        TopoDS_Edge origE = TopoDS::Edge(expOrig.Current());
        
        // 找到这条原始交线对应的来源面
        if (!edgeToFaceMap.IsBound(origE)) continue;
        TopoDS_Face origF = TopoDS::Face(edgeToFaceMap.Find(origE));
        
        // 如果这条边被完全删除了，跳过
        if (cutAlgo.IsDeleted(origE)) continue;
        
        // 检查这条边是否被修改（被切成了多段）
        const TopTools_ListOfShape& modified = cutAlgo.Modified(origE);
        if (!modified.IsEmpty()) {
            TopTools_ListIteratorOfListOfShape it(modified);
            for (; it.More(); it.Next()) {
                TopoDS_Edge newE = TopoDS::Edge(it.Value());
                // 确保新边确实在最终结果中
                if (resultEdgesMap.Contains(newE)) {
                    OneEdge oe;
                    oe.edge = newE;
                    oe.sourceFace = origF;
                    result.push_back(oe);
                }
            }
        } else {
            // 如果没有被修改，且在最终结果中，直接添加
            if (resultEdgesMap.Contains(origE)) {
                OneEdge oe;
                oe.edge = origE;
                oe.sourceFace = origF;
                result.push_back(oe);
            }
        }
    }
    return result;
}

// 辅助函数：把 vector<OneEdge> 转成 TopoDS_Compound
TopoDS_Compound BuildCompoundFromOneEdgeVector(const std::vector<OneEdge>& edges)
{
    BRep_Builder builder;
    TopoDS_Compound comp;
    builder.MakeCompound(comp);
    for (const OneEdge& e : edges) {
        builder.Add(comp, e.edge);
    }
    return comp;
}

// ======================== Global Compound Functions ========================

void InitGlobalCompounds() {
    if (!gCompoundsInitialized) {
        gBuilder.MakeCompound(gAllFacesCompound);
        gBuilder.MakeCompound(gAllLinesCompound);
        gCompoundsInitialized = true;
    }
}

// 判断两个点是否重合（带容差）
//bool IsPointEqual(const Point3D& p1, const Point3D& p2, double tol) {
//    return std::abs(p1.x - p2.x) < tol && std::abs(p1.y - p2.y) < tol;
//}

// --- Face ID 提取 (CollectCavityFaceIds) ---
std::set<int> CollectCavityFaceIds(const std::vector<std::vector<Face2D>>& allLayerFaces) {
    std::set<int> cavityFaceIds;

    // 遍历每一层
    for (int layerIdx = 0; layerIdx < allLayerFaces.size(); ++layerIdx)
    {
        const auto& layerFaces = allLayerFaces[layerIdx];
        std::set<int> currentLayerFaceIds; // 只存当前层的faceId

        // 遍历当前层所有面
        for (const auto& face : layerFaces) {
            if (face.type != FaceType::CAVITY)
                continue;

            // 外环线
            for (const auto& line : face.outerLoop) {
                if (line.faceId >= 0) {
                    currentLayerFaceIds.insert(line.faceId);
                    cavityFaceIds.insert(line.faceId);
                }
            }
            // 内环线
            for (const auto& innerLoop : face.innerLoops) {
                for (const auto& line : innerLoop) {
                    if (line.faceId >= 0) {
                        currentLayerFaceIds.insert(line.faceId);
                        cavityFaceIds.insert(line.faceId);
                    }
                }
            }
        }

        // 输出：当前层得到了哪些 faceId
        std::cout << "Layer " << layerIdx << " faces: ";
        for (int id : currentLayerFaceIds) {
            std::cout << id << " ";
        }
        std::cout << std::endl;
    }

    // 最终总结果
    std::cout << "All cavity face IDs: ";
    for (int id : cavityFaceIds) {
        std::cout << id << " ";
    }
    std::cout << std::endl;

    return cavityFaceIds;
}

// 根据 faceID 集合，从原始实体模型中提取对应的拓扑面

// --- Face ID 提取 (GetFacesByFaceIds, GetCavityCapFaces) ---
TopoDS_Compound GetFacesByFaceIds(const TopoDS_Shape& solid,
                                   const std::set<int>& targetFaceIds,
                                   const TopTools_DataMapOfShapeInteger& faceToIdMap) {
    BRep_Builder builder;
    TopoDS_Compound comp;
    builder.MakeCompound(comp);
    
    TopExp_Explorer exp(solid, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        const TopoDS_Face& face = TopoDS::Face(exp.Current());
        if (faceToIdMap.IsBound(face)) {
            int id = faceToIdMap.Find(face);
            if (targetFaceIds.find(id) != targetFaceIds.end()) {
                builder.Add(comp, face);
            }
        }
    }
    return comp;
}

// 获取封闭型腔的顶面/底面/中间平台面
// 条件：水平面（法线平行于Z轴），且它的所有邻接面都属于 cavityFaceIds 集合
TopoDS_Compound GetCavityCapFaces(const TopoDS_Shape& solid,
                                   const std::set<int>& cavityFaceIds,
                                   const TopTools_DataMapOfShapeInteger& faceToIdMap,
                                   double angleTolerance) {
    BRep_Builder builder;
    TopoDS_Compound capFaces;
    builder.MakeCompound(capFaces);
    
    // 建立边 -> 面的映射（用于查找邻接关系）
    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);
    
    const gp_Dir zAxis(0, 0, 1);
    
    TopExp_Explorer exp(solid, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        
        // 检查是否是水平面
        BRepAdaptor_Surface surf(face);
        if (surf.GetType() != GeomAbs_Plane) continue;
        
        gp_Pln plane = surf.Plane();
        gp_Dir normal = plane.Axis().Direction();
        double angle = normal.Angle(zAxis);
        if (angle > angleTolerance && fabs(angle - M_PI) > angleTolerance) continue;
        
        // 获取当前面的 ID
        int faceId = -1;
        if (faceToIdMap.IsBound(face)) {
            faceId = faceToIdMap.Find(face);
        } else {
            continue;
        }
        
        // 检查所有邻接面是否都属于 cavityFaceIds
        bool allAdjacentAreCavity = true;
        TopExp_Explorer edgeExp(face, TopAbs_EDGE);
        for (; edgeExp.More(); edgeExp.Next()) {
            const TopTools_ListOfShape& adjacentFaces = edgeToFaces.FindFromKey(edgeExp.Current());
            TopTools_ListIteratorOfListOfShape it(adjacentFaces);
            for (; it.More(); it.Next()) {
                const TopoDS_Face& adjFace = TopoDS::Face(it.Value());
                if (adjFace.IsSame(face)) continue; // 跳过自己
                
                if (faceToIdMap.IsBound(adjFace)) {
                    int adjId = faceToIdMap.Find(adjFace);
                    // 如果邻接面不在 cavityFaceIds 中，说明这个水平面不是型腔的封闭面
                    if (cavityFaceIds.find(adjId) == cavityFaceIds.end()) {
                        allAdjacentAreCavity = false;
                        break;
                    }
                } else {
                    allAdjacentAreCavity = false;
                    break;
                }
            }
            if (!allAdjacentAreCavity) break;
        }
        
        if (allAdjacentAreCavity && faceId >= 0) {
            builder.Add(capFaces, face);
        }
    }
    
    return capFaces;
}

// 将纯数学线段集合保存为 BREP 文件，用于可视化验证

// --- GetCavityCompound ---
TopoDS_Compound GetCavityCompound(const TopoDS_Shape& solid,
                                   const std::set<int>& cavityFaceIds,
                                   const TopTools_DataMapOfShapeInteger& faceToIdMap) {
    TopoDS_Compound wallFaces = GetFacesByFaceIds(solid, cavityFaceIds, faceToIdMap);
    TopoDS_Compound capFaces = GetCavityCapFaces(solid, cavityFaceIds, faceToIdMap);

    BRep_Builder builder;
    TopoDS_Compound result;
    builder.MakeCompound(result);

    TopExp_Explorer exp(wallFaces, TopAbs_FACE);
    for (; exp.More(); exp.Next()) builder.Add(result, exp.Current());

    TopExp_Explorer expCap(capFaces, TopAbs_FACE);
    for (; expCap.More(); expCap.Next()) builder.Add(result, expCap.Current());

    return result;
}

// ============================================================

// 步骤 1：将打散的线段拼接成封闭环
// 去重：删除完全重复的线段（关键！解决假环）
vector<OneLine> DeduplicateLines(const vector<OneLine>& lines) {
    vector<OneLine> res;
    int duplicateCount = 0;
    int reverseDuplicateCount = 0;

    auto MergePoint = [](Point3D& target, const Point3D& incoming) {
        target.x = (target.x + incoming.x) * 0.5;
        target.y = (target.y + incoming.y) * 0.5;
        target.z = (target.z + incoming.z) * 0.5;
    };

    for (const auto& line : lines) {
        bool dup = false;
        for (auto& r : res) {
            if (IsPointEqual(line.start, r.start) && IsPointEqual(line.end, r.end)) {
                MergePoint(r.start, line.start);
                MergePoint(r.end, line.end);
                if (r.faceId < 0 && line.faceId >= 0) r.faceId = line.faceId;
                dup = true;
                duplicateCount++;
                break;
            }
            if (IsPointEqual(line.start, r.end) && IsPointEqual(line.end, r.start)) {
                MergePoint(r.start, line.end);
                MergePoint(r.end, line.start);
                if (r.faceId < 0 && line.faceId >= 0) r.faceId = line.faceId;
                dup = true;
                reverseDuplicateCount++;
                break;
            }
        }
        if (!dup) res.push_back(line);
    }
    if (duplicateCount > 0 || reverseDuplicateCount > 0) {
        cout << "  [DeduplicateLines] 输入=" << lines.size()
            << " 输出=" << res.size()
            << " 同向重复=" << duplicateCount
            << " 反向重复=" << reverseDuplicateCount << endl;
    }
    return res;
}


// --- 面属性判断 ---
bool IsHorizontalFace(const TopoDS_Face& face, double tol) {
    BRepAdaptor_Surface surf(face);
    if (surf.GetType() != GeomAbs_Plane) return false;
    gp_Dir normal = surf.Plane().Axis().Direction();
    if (face.Orientation() == TopAbs_REVERSED) normal.Reverse();
    return normal.IsParallel(gp_Dir(0, 0, 1), tol);
}

// 获取面在中心点处的法线方向
gp_Dir GetFaceNormal(const TopoDS_Face& face) {
    BRepAdaptor_Surface surf(face);
    Standard_Real umin, umax, vmin, vmax;
    BRepTools::UVBounds(face, umin, umax, vmin, vmax);

    gp_Pnt p;
    gp_Vec du, dv;
    surf.D1((umin + umax) / 2.0, (vmin + vmax) / 2.0, p, du, dv);

    gp_Vec normal = du.Crossed(dv);
    if (face.Orientation() == TopAbs_REVERSED) {
        normal.Reverse();
    }

    if (normal.Magnitude() > 1e-6) {
        return gp_Dir(normal);
    }
    return gp_Dir(0, 0, 1); // 默认向上
}


double GetHorizontalFaceZ(const TopoDS_Face& face) {
    BRepAdaptor_Surface surf(face);
    return surf.Plane().Location().Z();
}

// ============================================================
// 核心逻辑：递归切分型腔
// ============================================================

// --- OCC 面属性工具 ---
double CalculateFaceAreaOCC(const TopoDS_Face& face) {
    GProp_GProps props;
    BRepGProp::SurfaceProperties(face, props);
    return props.Mass();
}

int CountFacesInCompound(const TopoDS_Compound& compound) {
    int count = 0;
    TopExp_Explorer exp(compound, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        count++;
    }
    return count;
}

gp_Dir SnapNormalToDirection(const gp_Dir& normal) {
    double nx = std::abs(normal.X());
    double ny = std::abs(normal.Y());
    double nz = std::abs(normal.Z());
    if (nx >= ny && nx >= nz) return gp_Dir(normal.X() > 0.0 ? 1.0 : -1.0, 0.0, 0.0);
    if (ny >= nx && ny >= nz) return gp_Dir(0.0, normal.Y() > 0.0 ? 1.0 : -1.0, 0.0);
    return gp_Dir(0.0, 0.0, normal.Z() > 0.0 ? 1.0 : -1.0);
}

bool IsCavityMachinable(const TopoDS_Compound& cavity, gp_Dir& outToolDir) {
    struct FaceInfo {
        TopoDS_Face face;
        double area;
        gp_Dir dir;
    };
    std::vector<FaceInfo> faces;

    TopExp_Explorer exp(cavity, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face f = TopoDS::Face(exp.Current());
        BRepAdaptor_Surface surf(f);
        if (surf.GetType() != GeomAbs_Plane) continue;
        gp_Pln plane = surf.Plane();
        gp_Dir normal = plane.Axis().Direction();
        if (f.Orientation() == TopAbs_REVERSED) normal.Reverse();
        faces.push_back({f, CalculateFaceAreaOCC(f), SnapNormalToDirection(normal)});
    }

    if (faces.empty()) { outToolDir = gp_Dir(0, 0, -1); return false; }

    std::sort(faces.begin(), faces.end(),
        [](const FaceInfo& a, const FaceInfo& b) { return a.area > b.area; });

    for (const auto& fi : faces) {
        gp_Dir feedDir = fi.dir.Reversed();
        bool blocked = false;
        for (const auto& other : faces) {
            if (other.face.IsSame(fi.face)) continue;
            if (other.dir.IsEqual(feedDir, 1e-6)) {
                blocked = true;
                break;
            }
        }
        std::cout << "    [加工可行性] 面法向=(" << fi.dir.X() << "," << fi.dir.Y() << "," << fi.dir.Z()
            << ") 进刀方向=(" << feedDir.X() << "," << feedDir.Y() << "," << feedDir.Z()
            << ") 面积=" << fi.area
            << " -> " << (blocked ? "不可行" : "可行") << std::endl;
        if (!blocked) {
            outToolDir = feedDir;
            return true;
        }
    }
    outToolDir = gp_Dir(0, 0, -1);
    return false;
}


// --- Compound 工具 ---
void AddCompoundFaces(TopoDS_Compound& target, const TopoDS_Compound& source) {
    BRep_Builder builder;
    TopExp_Explorer exp(source, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        builder.Add(target, exp.Current());
    }
}

bool GetCompoundZRange(const TopoDS_Compound& compound, double& zmin, double& zmax) {
    Bnd_Box box;
    BRepBndLib::Add(compound, box);
    if (box.IsVoid()) {
        zmin = zmax = 0.0;
        return false;
    }

    double xmin, ymin, xmax, ymax;
    box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    return true;
}

