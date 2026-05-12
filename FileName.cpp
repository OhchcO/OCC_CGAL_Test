#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <cmath>

//CGAL 库用于计算凸包
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>

// OCC 核心基础
#include <gp_Pnt.hxx>
#include <gp_Pln.hxx>
#include <gp_Dir.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#pragma comment(lib, "TKFeat.lib")    // 包含 BRepFeat_SplitShape
#pragma comment(lib, "TKBool.lib")    // 包含布尔运算
#pragma comment(lib, "TKTopAlgo.lib") // 包含基础拓扑算法
#pragma comment(lib, "TKGeomAlgo.lib")// 包含几何算法
// STEP 读取
#include <STEPControl_Reader.hxx>
#include <IFSelect_ReturnStatus.hxx>

// 拓扑遍历与计算
#include <TopExp_Explorer.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>

// 求交算法头文件
#include <BRepAlgoAPI_Section.hxx>

// 几何分析
#include <BRepAdaptor_Surface.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_DataMapOfShapeShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopExp.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GCPnts_QuasiUniformDeflection.hxx>
#include <Geom2d_Curve.hxx>
#include <gp_Vec2d.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <TopTools_DataMapOfIntegerShape.hxx>
#include <algorithm>
#include <BRepFeat_SplitShape.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <TopoDS_Iterator.hxx>
using namespace std;

// 分割算法
#include <TopTools_MapOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <BRepExtrema_DistShapeShape.hxx>

#include <GCPnts_UniformAbscissa.hxx>
#include <CPnts_AbscissaPoint.hxx>

struct Point3D {
    double x, y, z;
};

struct OneEdge {
    TopoDS_Edge edge;
    TopoDS_Face sourceFace;
};

struct OneLine {
    Point3D start;
    Point3D end;
    int faceId;
};

// ---------------- 数据结构定义 ----------------
// 一个封闭环就是一组首尾相连的线段
typedef std::vector<OneLine> Loop;

// 树节点：用于表示环的嵌套关系
struct LoopNode {
    Loop loopLines;             // 当前环的线段
    double area;                // 环的面积（用于排序）
    std::vector<LoopNode*> children; // 包含在它内部的子环
    
    LoopNode() : area(0.0) {}
};

enum class FaceType {
    SOLID,  // 实体面（有材料）
    CAVITY,  // 空腔面（无材料/孔洞）
    HULL     // 凸包面
};

// 最终提取出的 2D 面（全空间剖分）
struct Face2D {
    FaceType type;              // 面类型
    Loop outerLoop;             // 外环（唯一）
    std::vector<Loop> innerLoops; // 内环/孔洞（可能有多个）
};
// ----------------------------------------------

// 型腔环，型腔特征定义
struct CavityLoop {
    std::vector<OneLine> lines;
    double zHeight;
    bool isOuter;
};

struct CavityFeature {
    int featureId;
    std::vector<CavityLoop> stepLoops;
    double topZ = -1e9;
    double bottomZ = 1e9;
    double totalDepth = 0.0;
    gp_Dir toolDirection = gp_Dir(0, 0, -1);
    std::set<int> sourceFaceIds; // 关联的原始面 ID 集合
};

#include <fstream>
#include <cstdlib>

#include <fstream>
#include <cstdlib>
#include <string>
#include <direct.h> // 用于获取当前路径

std::string savePath = "E:\\soft\\code\\Project1\\output\\";
std::string inputPath = "E:\\soft\\code\\cMake_test\\input\\";
//std::string inputPath = "E:\\soft\\code\\Project1\\input\\";

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;

/**
 * @brief 对 Face2D 的外环计算二维凸包，并将结果封装回一个新的 Face2D 结构
 * @param face 输入的原始面（SOLID 类型）
 * @return Face2D 类型为 HULL 的新面，包含封闭的凸包外环
 */
Face2D ComputeConvexHullFace(const Face2D& face) {
    Face2D resultFace;
    resultFace.type = FaceType::HULL;      // 设置面类型为 HULL
    resultFace.innerLoops.clear();         // 明确内环为空

    if (face.outerLoop.empty()) {
        return resultFace;
    }

    std::vector<Point_2> points_2d;
    double currentZ = face.outerLoop[0].start.z;

    // 1. 收集所有端点投影到 2D
    for (const auto& line : face.outerLoop) {
        points_2d.emplace_back(line.start.x, line.start.y);
        points_2d.emplace_back(line.end.x, line.end.y);
    }

    // 2. 计算 2D 凸包
    std::vector<Point_2> hull_points;
    CGAL::convex_hull_2(points_2d.begin(), points_2d.end(), std::back_inserter(hull_points));

    // 3. 还原为 OneLine 闭合环并存入 outerLoop
    for (size_t i = 0; i < hull_points.size(); ++i) {
        OneLine line;
        // 起点
        line.start = {
            CGAL::to_double(hull_points[i].x()),
            CGAL::to_double(hull_points[i].y()),
            currentZ
        };
        // 终点 (取模连回起点)
        size_t nextIdx = (i + 1) % hull_points.size();
        line.end = {
            CGAL::to_double(hull_points[nextIdx].x()),
            CGAL::to_double(hull_points[nextIdx].y()),
            currentZ
        };

        line.faceId = -1; // 标记为凸包边缘
        resultFace.outerLoop.push_back(line);
    }

    return resultFace;
}

//打印lines
void PrintAllLines(const std::vector<OneLine>& lines)
{
    printf("========== ALL LINES (%d) ==========\n", (int)lines.size());
    for (int i = 0; i < lines.size(); i++) {
        const auto& l = lines[i];
        printf("[%4d] (%.6f, %.6f) -> (%.6f, %.6f)\n",
            i,
            l.start.x, l.start.y,
            l.end.x, l.end.y);
    }
    printf("=====================================\n");
}
// 获取面的 Z 范围
void GetFaceZRange(const TopoDS_Face& face, double& zmin, double& zmax) {
    Bnd_Box box;
    BRepBndLib::Add(face, box);
    if (!box.IsVoid()) {
        double xmin, ymin, xmax, ymax;
        box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    }
    else {
        zmin = zmax = 0.0;
    }
}
// 自适应划分，获取所有切分点
// 1. 所有水平面（任意高度，包括台阶、内台、顶面、底面）
// 2. 所有非平面的最低点 Z（侧壁底部、圆角底部、斜面底部等）
// 最后去重合并，按从大到小排序
std::vector<double> GetSplitPointsAlongZ(const TopoDS_Shape& shape,
    double angleTolerance = 0.001,
    double mergeTol = 1e-3)
{
    std::set<double> zPoints;
    const gp_Dir zAxis(0, 0, 1);

    TopExp_Explorer exp(shape, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        BRepAdaptor_Surface surf(face);

        // ==============================
        // 第一步：判断是不是水平面
        // ==============================
        bool isHorizontal = false;
        double zPlane = 0.0;

        if (surf.GetType() == GeomAbs_Plane) {
            gp_Pln plane = surf.Plane();
            gp_Dir normal = plane.Axis().Direction();

            if (face.Orientation() == TopAbs_REVERSED)
                normal.Reverse();

            double angle = normal.Angle(zAxis);
            if (angle < angleTolerance || fabs(angle - M_PI) < angleTolerance) {
                isHorizontal = true;
                zPlane = plane.Location().Z();
            }
        }

        // ==============================
        // 规则 1：水平面 → 加 Z
        // ==============================
        if (isHorizontal) {
            zPoints.insert(zPlane);
        }
        // ==============================
        // 规则 2：所有其他面（竖直平面、斜面、圆角、曲面）→ 加 最低点 Z
        // ==============================
        else {
            Bnd_Box box;
            BRepBndLib::Add(face, box);
            if (!box.IsVoid()) {
                double xmin, ymin, zmin, xmax, ymax, zmax;
                box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
                zPoints.insert(zmin);
            }
        }
    }

    // 从大到小排序
    std::vector<double> sortedZ(zPoints.begin(), zPoints.end());
    std::sort(sortedZ.begin(), sortedZ.end(), std::greater<double>());

    // 按容差去重
    std::vector<double> finalPoints;
    for (double z : sortedZ) {
        if (finalPoints.empty()) {
            finalPoints.push_back(z);
        }
        else {
            if (finalPoints.back() - z > mergeTol) {
                finalPoints.push_back(z);
            }
        }
    }

    return finalPoints;
}

// 1. 获得平面集合（基于精确几何方程，忽略包围盒容差）
TopoDS_Compound GetCoplanarFaces(const TopoDS_Shape& shape, double splitZ, double tol = 1e-4) {
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
TopoDS_Compound GetConcaveEdges(const TopoDS_Shape& solid, const TopoDS_Compound& coplanarFaces, const gp_Dir& toolDirection = gp_Dir(0, 0, 1)) {
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
    const gp_Dir& toolDirection = gp_Dir(0, 0, 1))
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

BRep_Builder gBuilder;
TopoDS_Compound gAllFacesCompound;
TopoDS_Compound gAllLinesCompound;
bool gCompoundsInitialized = false;

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
bool IsPointEqual(const Point3D& p1, const Point3D& p2, double tol = 1e-3) {
    return fabs(p1.x - p2.x) < tol &&
        fabs(p1.y - p2.y) < tol;
}
// ================== 封闭型腔特征提取函数 ==================

// 从所有分层切片中提取出 CAVITY 面的 faceId 集合
std::set<int> CollectCavityFaceIds(const std::vector<std::vector<Face2D>>& allLayerFaces) {
    std::set<int> cavityFaceIds;
    for (const auto& layerFaces : allLayerFaces) {
        for (const auto& face : layerFaces) {
            if (face.type != FaceType::CAVITY) continue;
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

// 根据 faceID 集合，从原始实体模型中提取对应的拓扑面
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
                                   double angleTolerance = 0.001) {
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

// 保存整个封闭型腔（侧壁 + 顶/底面）为 BREP
void ExportCavityFaces(const TopoDS_Shape& solid,
                        const std::set<int>& cavityFaceIds,
                        const TopTools_DataMapOfShapeInteger& faceToIdMap,
                        const std::string& fileName) {
    // 1. 获取侧壁面
    TopoDS_Compound wallFaces = GetFacesByFaceIds(solid, cavityFaceIds, faceToIdMap);
    
    // 2. 获取顶/底面/中间平台面
    TopoDS_Compound capFaces = GetCavityCapFaces(solid, cavityFaceIds, faceToIdMap);
    
    // 3. 合并
    BRep_Builder builder;
    TopoDS_Compound result;
    builder.MakeCompound(result);
    
    TopExp_Explorer exp(wallFaces, TopAbs_FACE);
    for (; exp.More(); exp.Next()) builder.Add(result, exp.Current());
    
    TopExp_Explorer expCap(capFaces, TopAbs_FACE);
    for (; expCap.More(); expCap.Next()) builder.Add(result, expCap.Current());
    
    BRepTools::Write(result, fileName.c_str());
    std::cout << "  完整封闭型腔已保存至: " << fileName << std::endl;
}

// 获取整个封闭型腔的合并面（侧壁 + 顶/底面/中间平台），返回 Compound
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
    for (const auto& line : lines) {
        bool dup = false;
        for (const auto& r : res) {
            if (IsPointEqual(line.start, r.start) && IsPointEqual(line.end, r.end)) {
                dup = true;
                break;
            }
            if (IsPointEqual(line.start, r.end) && IsPointEqual(line.end, r.start)) {
                dup = true;
                break;
            }
        }
        if (!dup) res.push_back(line);
    }
    return res;
}

// 核心成环算法：只找大环，自动忽略碎线
vector<Loop> BuildLoops(const vector<OneLine>& inputLines) {
    vector<Loop> loops;
    auto lines = DeduplicateLines(inputLines); // 先去重！
    vector<bool> visited(lines.size(), false);
    const double tol = 0.1;

    for (size_t i = 0; i < lines.size(); ++i) {
        if (visited[i]) continue;

        Loop currentLoop;
        currentLoop.push_back(lines[i]);
        visited[i] = true;
        Point3D start = currentLoop.front().start;
        Point3D currEnd = currentLoop.back().end;

        while (true) {
            bool found = false;
            for (size_t j = 0; j < lines.size(); ++j) {
                if (visited[j]) continue;
                const auto& l = lines[j];

                if (IsPointEqual(currEnd, l.start, tol)) {
                    currentLoop.push_back(l);
                    currEnd = l.end;
                    visited[j] = true;
                    found = true;
                    break;
                }
                if (IsPointEqual(currEnd, l.end, tol)) {
                    OneLine rev = l;
                    swap(rev.start, rev.end);
                    currentLoop.push_back(rev);
                    currEnd = rev.end;
                    visited[j] = true;
                    found = true;
                    break;
                }
            }
            if (!found) break;
            if (IsPointEqual(currEnd, start, tol)) break;
        }

        // ===================== 关键过滤 =====================
        // 只保留线段数 >=3 的大环，自动过滤假环、碎环
        if (IsPointEqual(currentLoop.back().end, currentLoop.front().start, tol)
            && currentLoop.size() >= 3)
        {
            loops.push_back(currentLoop);
        }
    }

    // 输出日志
    cout << "\n================ 最终结果 =================" << endl;
    cout << " 成功生成环数量：" << loops.size()  << endl;
    for (int i = 0; i < loops.size(); i++) {
        cout << "  环 " << i + 1 << "：" << loops[i].size() << " 条线段" << endl;
    }
    cout << "==========================================\n" << endl;

    return loops;
}



// 计算多边形面积（鞋带公式 Shoelace Formula）
double CalculateArea(const Loop& loop) {
    double area = 0.0;
    for (const auto& line : loop) {
        area += (line.start.x * line.end.y - line.end.x * line.start.y);
    }
    return std::abs(area) / 2.0;
}

// 步骤 2：射线法判断点是否在多边形（环）内部
bool IsPointInLoop(const Point3D& pt, const Loop& loop) {
    bool inside = false;
    // 为了避免射线刚好打在多边形顶点上的经典 bug，我们给测试点的 Y 坐标加一个极小的偏移量
    double testY = pt.y + 1e-7; 
    
    for (const auto& line : loop) {
        double x1 = line.start.x, y1 = line.start.y;
        double x2 = line.end.x, y2 = line.end.y;
        
        // 射线向 +X 方向发射，统计交点个数
        if (((y1 > testY) != (y2 > testY)) &&
            (pt.x < (x2 - x1) * (testY - y1) / (y2 - y1) + x1)) {
            inside = !inside;
        }
    }
    return inside;
}

// 递归提取面（全空间剖分：无论奇偶层都提取为面）
void ExtractFacesFromTree(LoopNode* node, int depth, std::vector<Face2D>& faces) {
    // 根节点（depth == 0）是虚拟节点，没有自身的线段，直接跳过其自身的构造
    if (depth > 0) {
        Face2D newFace;
        // 奇数层代表实体，偶数层代表空腔
        newFace.type = (depth % 2 != 0) ? FaceType::SOLID : FaceType::CAVITY;
        newFace.outerLoop = node->loopLines;
        
        // 它的直接子节点就是它的内环
        for (LoopNode* child : node->children) {
            newFace.innerLoops.push_back(child->loopLines);
        }
        faces.push_back(newFace);
    }
    
    // 继续递归遍历子节点
    for (LoopNode* child : node->children) {
        ExtractFacesFromTree(child, depth + 1, faces);
    }
}

// 对每个环可视化导出为独立的 BREP 文件
void DebugExportLoops(const std::vector<Loop>& loops, const std::string& prefix, const std::string& savePath) {
    for (size_t i = 0; i < loops.size(); ++i) {
        BRep_Builder builder;
        TopoDS_Compound comp;
        builder.MakeCompound(comp);
        
        for (const auto& line : loops[i]) {
            gp_Pnt p1(line.start.x, line.start.y, line.start.z);
            gp_Pnt p2(line.end.x, line.end.y, line.end.z);
            if (!p1.IsEqual(p2, 1e-7)) {
                TopoDS_Edge anEdge = BRepBuilderAPI_MakeEdge(p1, p2);
                builder.Add(comp, anEdge);
            }
        }
        
        double area = CalculateArea(loops[i]);
        std::string fileName = savePath + prefix + "_Loop_" + std::to_string(i) + "_Area" + std::to_string((int)area) + ".brep";
        BRepTools::Write(comp, fileName.c_str());
        std::cout << "  环 " << i << " 面积=" << area << " -> " << fileName << std::endl;
    }
}

// 构建嵌套树并提取 Face2D
std::vector<Face2D> BuildTopologyAndExtractFaces(const std::vector<OneLine>& layerLines) {
    // 1. 组装成环
    PrintAllLines(layerLines);
    std::vector<Loop> loops = BuildLoops(layerLines);
    
    // 调试：导出每个环
    //DebugExportLoops(loops, "DebugLayer", savePath);
    
    // 2. 创建树节点并计算面积
    std::vector<LoopNode*> nodes;
    for (const auto& loop : loops) {
        LoopNode* node = new LoopNode();
        node->loopLines = loop;
        node->area = CalculateArea(loop);
        nodes.push_back(node);
    }
    
    // 3. 按面积从大到小排序
    std::sort(nodes.begin(), nodes.end(), [](LoopNode* a, LoopNode* b) {
        return a->area > b->area;
    });
    
    // 4. 构建嵌套树
    LoopNode root; // 虚拟根节点（第 0 层）
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        bool foundParent = false;
        // 从比它大一点的环开始往上找（即从 i-1 倒序遍历到 0）
        // 找到的第一个包含它的环，就是它的直接父节点（面积最小的包围者）
        for (int j = (int)i - 1; j >= 0; --j) {
            // 取当前环的第一个点进行测试
            if (IsPointInLoop(nodes[i]->loopLines.front().start, nodes[j]->loopLines)) {
                nodes[j]->children.push_back(nodes[i]);
                foundParent = true;
                break;
            }
        }
        // 如果没有环包含它，它就是最外层，挂在虚拟根节点下
        if (!foundParent) {
            root.children.push_back(nodes[i]);
        }
    }
    
    // 5. 提取 Face2D (全空间剖分)
    std::vector<Face2D> finalFaces;
    // 从虚拟根节点开始提取，根节点深度为0
    ExtractFacesFromTree(&root, 0, finalFaces); 
    
    // 清理内存
    for (LoopNode* node : nodes) {
        delete node;
    }
    
    return finalFaces;
}

// 将提取出的 Face2D 转换回 OCCT 的 TopoDS_Face 并保存为 BREP
void ExportFace2DToBrep(const std::vector<Face2D>& faces, const std::string& fileName) {
    if (faces.empty()) return;

    BRep_Builder builder;
    TopoDS_Compound comp;
    builder.MakeCompound(comp);

    for (const auto& face2d : faces) {
        try {
            if (face2d.outerLoop.size() < 3) continue;

            // 1. 构建外环实心面
            BRepBuilderAPI_MakePolygon outerPoly;
            for (const auto& line : face2d.outerLoop) {
                outerPoly.Add(gp_Pnt(line.start.x, line.start.y, line.start.z));
            }
            outerPoly.Close(); 
            if (!outerPoly.IsDone()) continue;
            
            TopoDS_Wire outerWire = outerPoly.Wire();
            if (outerWire.IsNull() || !outerWire.Closed()) continue;

            BRepBuilderAPI_MakeFace outerFaceMaker(outerWire);
            if (!outerFaceMaker.IsDone()) continue;
            TopoDS_Shape currentShape = outerFaceMaker.Face();

            // 2. 用布尔减法挖去所有内环（孔洞）
            for (const auto& innerLoop : face2d.innerLoops) {
                if (innerLoop.size() < 3) continue;
                
                BRepBuilderAPI_MakePolygon innerPoly;
                for (const auto& line : innerLoop) {
                    innerPoly.Add(gp_Pnt(line.start.x, line.start.y, line.start.z));
                }
                innerPoly.Close();
                
                if (innerPoly.IsDone()) {
                    TopoDS_Wire innerWire = innerPoly.Wire();
                    if (!innerWire.IsNull() && innerWire.Closed()) {
                        BRepBuilderAPI_MakeFace innerFaceMaker(innerWire);
                        if (innerFaceMaker.IsDone()) {
                            // 核心修复：使用极其稳定的布尔减法，大面减小面
                            BRepAlgoAPI_Cut cutAlgo(currentShape, innerFaceMaker.Face());
                            cutAlgo.Build();
                            if (cutAlgo.IsDone()) {
                                currentShape = cutAlgo.Shape();
                            }
                        }
                    }
                }
            }

            // 3. 将最终生成好的面加入 Compound
            builder.Add(comp, currentShape);
        }
        catch (const Standard_Failure& e) {
            std::cerr << "  [警告] 生成面时发生 OCCT 异常: " << e.GetMessageString() << std::endl;
        }
        catch (...) {
            std::cerr << "  [警告] 生成面时发生未知异常" << std::endl;
        }
    }

    BRepTools::Write(comp, fileName.c_str());
    std::cout << "  重建的面已保存至: " << fileName << std::endl;
}

// 将纯数学线段集合保存为 BREP 文件，用于可视化验证
void ExportOneLinesToBrep(const std::vector<OneLine>& lines, const std::string& fileName) {
    if (lines.empty()) return;
    
    BRep_Builder builder;
    TopoDS_Compound comp;
    builder.MakeCompound(comp);
    
    for (const auto& line : lines) {
        gp_Pnt p1(line.start.x, line.start.y, line.start.z);
        gp_Pnt p2(line.end.x, line.end.y, line.end.z);
        
        // 防止起点和终点重合导致 MakeEdge 失败
        if (!p1.IsEqual(p2, 1e-7)) {
            TopoDS_Edge anEdge = BRepBuilderAPI_MakeEdge(p1, p2);
            builder.Add(comp, anEdge);
        }
    }
    
    BRepTools::Write(comp, fileName.c_str());
    cout << "  打散后的线段已保存至: " << fileName << endl;
}

//获取切分结果：切分面 + 切分线（减去凹边）
std::vector<Face2D> SliceModelAtZ(const TopoDS_Shape& shape, double splitZ, int sliceIndex, int maxIndex, const TopTools_DataMapOfShapeInteger& faceToIdMap) {
    InitGlobalCompounds();
    
    gp_Pln cuttingPlane(gp_Pnt(0, 0, splitZ), gp_Dir(0, 0, 1));
    TopoDS_Face algoPlane = BRepBuilderAPI_MakeFace(cuttingPlane);
    
    BRep_Builder builder;
    TopoDS_Compound selectedFacesCompound;
    builder.MakeCompound(selectedFacesCompound);
    
    TopoDS_Compound intersectionLinesCompound;
    builder.MakeCompound(intersectionLinesCompound);
    
    int faceCount = 0;
    int lineCount = 0;
    
    TopTools_DataMapOfShapeShape edgeToFaceMap;
    
    // 1. 获取表面（共面的水平面）
    TopoDS_Compound coplanarFaces = GetCoplanarFaces(shape, splitZ);
    std::string coplanarFacesFileName = savePath + "coplanarFaces.brep";
    BRepTools::Write(coplanarFaces, coplanarFacesFileName.c_str());
    
    // 为了快速判断一个面是否是共面水平面，我们把它存到一个集合中
    TopTools_IndexedMapOfShape coplanarFacesMap;
    TopExp::MapShapes(coplanarFaces, TopAbs_FACE, coplanarFacesMap);

    TopExp_Explorer exp(shape, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        
        // 【关键修复】：如果这个面本身就是我们要切的水平面（顶面/底面/平台），直接跳过求交！
        // 这样可以避免共面求交生成重复的边界线。
        if (coplanarFacesMap.Contains(face)) {
            continue;
        }

        double zmin, zmax;
        GetFaceZRange(face, zmin, zmax);
        
        
        if (splitZ > zmin  && splitZ < zmax ) {
            builder.Add(selectedFacesCompound, face);
            gBuilder.Add(gAllFacesCompound, face);
            faceCount++;
            
            BRepAlgoAPI_Section section(face, algoPlane, Standard_True);
            section.Build();
            
            if (section.IsDone()) {
                TopoDS_Shape intersectionShape = section.Shape();
                TopExp_Explorer edgeExp(intersectionShape, TopAbs_EDGE);
                for (; edgeExp.More(); edgeExp.Next()) {
                    TopoDS_Edge E = TopoDS::Edge(edgeExp.Current());
                    builder.Add(intersectionLinesCompound, E);
                    edgeToFaceMap.Bind(E, face); // 记录交线来源的面
                    lineCount++;
                }
            }
        }
    }
    
    // 2. 统一分类清洗：对 coplanarFaces 中的边做全量分类，
    //    光顺边用平面法向·邻面法向等效判断；非光顺边用 S₁·S₂ 几何判断
    //    只保留需要的边（朝下保留凹边，朝上保留凸边）
    gp_Dir toolDirection(0, 0, -1);
    TopoDS_Compound edgesToRemove;
    BRep_Builder builderToRemove; 
    builderToRemove.MakeCompound(edgesToRemove); // 默认初始化为一个空的容器
    if (sliceIndex != maxIndex)
    {
        edgesToRemove = ClassifyAndGetEdgesToRemove(shape, coplanarFaces, toolDirection);
    }
    std::string edgesToRemoveFileName = savePath + "edgesToRemove.brep";
    BRepTools::Write(edgesToRemove, edgesToRemoveFileName.c_str());
    
    // 3. 从交线中布尔减掉要删除的边，得到最终干净的边缘集合
    std::vector<OneEdge> finalEdges = SubtractLinesAndMapFaces(intersectionLinesCompound, edgesToRemove, edgeToFaceMap);

    // 重新构建 finalLinesCompound 用于保存
    TopoDS_Compound finalLinesCompound;
    builder.MakeCompound(finalLinesCompound);
    
    std::vector<OneLine> layerLines; // 存储打散后的纯数学线段
    
    for (const auto& oneEdge : finalEdges) {
        builder.Add(finalLinesCompound, oneEdge.edge);
        gBuilder.Add(gAllLinesCompound, oneEdge.edge);
        
        // 获取来源面的 ID
        int faceId = -1;
        if (faceToIdMap.IsBound(oneEdge.sourceFace)) {
            faceId = faceToIdMap.Find(oneEdge.sourceFace);
        }
        
        // 将复杂边打散为简单线段
        BRepAdaptor_Curve bac(oneEdge.edge);

        //// 0.01 是离散化容差，值越小曲线被切分得越细
        //GCPnts_QuasiUniformDeflection discretizer(bac, 0.1); 
        //if (discretizer.IsDone()) {
        //    int nbPoints = discretizer.NbPoints();
        //    for (int i = 1; i < nbPoints; ++i) {
        //        gp_Pnt p1 = discretizer.Value(i);
        //        gp_Pnt p2 = discretizer.Value(i + 1);
        //        
        //        OneLine line;
        //        line.start.x = p1.X();
        //        line.start.y = p1.Y();
        //        line.start.z = p1.Z();
        //        line.end.x = p2.X();
        //        line.end.y = p2.Y();
        //        line.end.z = p2.Z();
        //        line.faceId = faceId;
        //        
        //        layerLines.push_back(line);
        //    }
        //}

        // 1. 获取总弧长
        Standard_Real totalLen = CPnts_AbscissaPoint::Length(bac);

        // 2. 设定固定步长（建议 0.5mm，解决你 8000 条碎线太密的问题）
        Standard_Real myStep = 0.5;

        // 3. 执行离散化
        GCPnts_UniformAbscissa discretizer;
        discretizer.Initialize(bac, myStep);

        if (discretizer.IsDone()) {
            int nbPoints = discretizer.NbPoints();

            // 关键点：UniformAbscissa 存储的是参数 u
            // 我们需要用曲线对象 bac.Value(u) 来换取物理坐标
            for (int i = 1; i < nbPoints; ++i) {
                // 获取第 i 个点和第 i+1 个点的曲线参数
                Standard_Real u1 = discretizer.Parameter(i);
                Standard_Real u2 = discretizer.Parameter(i + 1);

                // 通过参数计算真实的 3D 点
                gp_Pnt p1 = bac.Value(u1);
                gp_Pnt p2 = bac.Value(u2);

                // 过滤零长度线段
                if (p1.Distance(p2) < 1e-6) continue;

                OneLine line;
                line.start = { p1.X(), p1.Y(), p1.Z() };
                line.end = { p2.X(), p2.Y(), p2.Z() };
                line.faceId = faceId;

                layerLines.push_back(line);
            }
        }

    }
    


    cout << "  打散为简单线段数量: " << layerLines.size() << endl;

    // 执行拓扑重建
    std::vector<Face2D> layerFaces = BuildTopologyAndExtractFaces(layerLines);
    cout << "  成功重建拓扑，提取出 " << layerFaces.size() << " 个带孔洞的 2D 面。" << endl;

    // 将重建的 2D 面保存为 BREP 文件（分别保存实体和空腔）
    std::vector<Face2D> layerSolid, layerCavity;
    for (const auto& f : layerFaces) {
        if (f.type == FaceType::SOLID) layerSolid.push_back(f);
        else layerCavity.push_back(f);
    }

    std::string solidFacesFileName = savePath + "Slice_" + std::to_string(sliceIndex) + "_SolidFaces_Z" + std::to_string(splitZ) + ".brep";
    std::string cavityFacesFileName = savePath + "Slice_" + std::to_string(sliceIndex) + "_CavityFaces_Z" + std::to_string(splitZ) + ".brep";
    ExportFace2DToBrep(layerSolid, solidFacesFileName);
    ExportFace2DToBrep(layerCavity, cavityFacesFileName);

    // 调用新封装的函数保存打散后的线段
    std::string discretizedFileName = savePath + "Slice_" + std::to_string(sliceIndex) + "_Discretized_Z" + std::to_string(splitZ) + ".brep";
    ExportOneLinesToBrep(layerLines, discretizedFileName);

    std::string faceFileName = savePath + "Slice_" + std::to_string(sliceIndex) + "_Faces_Z" + std::to_string(splitZ) + ".brep";
    std::string templineFileName = savePath + "Slice_" + std::to_string(sliceIndex) + "_Temp_Lines_Z" + std::to_string(splitZ) + ".brep";
    std::string lineFileName = savePath + "Slice_" + std::to_string(sliceIndex) + "_Lines_Z" + std::to_string(splitZ) + ".brep";
    
    if (faceCount > 0) {
        BRepTools::Write(selectedFacesCompound, faceFileName.c_str());
        cout << "  筛选出面数: " << faceCount << " -> " << faceFileName << endl;
    }
    


    if (lineCount > 0) {
        BRepTools::Write( intersectionLinesCompound, templineFileName.c_str());
        cout << "  交线数: " << lineCount << " -> " << templineFileName << endl;
        BRepTools::Write(finalLinesCompound, lineFileName.c_str());
        cout << "  交线已修正" << endl;
    }
    else {
        cout << "  未生成交线" << endl;
    }
    int pa = 0;
    
    return layerFaces;
}

// 判断面是否为水平面
bool IsHorizontalFace(const TopoDS_Face& face, double tol = 1e-4) {
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
void SplitCavityRecursive(const std::set<int>& currentFaceIds,
    TopTools_DataMapOfIntegerShape& idToFace,
    int& nextId,
    const std::vector<double>& splitZs,
    int zIndex,
    std::vector<std::set<int>>& result) {

    // 终止条件：没有切分点或没有待处理的面
    if (zIndex >= splitZs.size() || currentFaceIds.empty()) {
        if (!currentFaceIds.empty()) result.push_back(currentFaceIds);
        return;
    }

    double curSplitZ = splitZs[zIndex];
    std::set<int> aboveIds, belowIds;
    std::set<int> facesToSplit;

    // --- 1. 初步筛选：区分哪些面需要被“切一刀” ---
    for (int id : currentFaceIds) {
        TopoDS_Face face = TopoDS::Face(idToFace.Find(id));

        // 水平面分类
        if (IsHorizontalFace(face)) {
            double z = GetHorizontalFaceZ(face);

            // 如果面不在当前切分高度，按常规高度判断
            if (std::abs(z - curSplitZ) > 1e-4) {
                if (z > curSplitZ) aboveIds.insert(id);
                else belowIds.insert(id);
                continue;
            }

            // --- 关键修改：处理刚好在切分点处的水平面 ---
            gp_Dir normal = GetFaceNormal(face);
            double dot = normal.Z(); // 法线与 [0,0,1] 的点积

            if (dot < -0.5) {
                // 1. 法线向下 [0,0,-1]，与进刀方向相同
                // 它是下一层的“盖子”，传给下方
                belowIds.insert(id);
                cout << "  [水平面判定] ID:" << id << " 法线向下 -> 归入下方 (盖子)" << endl;
            }
            else {
                // 2. 法线向上 [0,0,1]，与进刀方向相反
                // 它是这一层的“底面”，归入上方
                aboveIds.insert(id);
                cout << "  [水平面判定] ID:" << id << " 法线向上 -> 归入上方 (底面)" << endl;
            }
            continue;
        }

        // 侧面根据包围盒分类
        double zmin, zmax;
        GetFaceZRange(face, zmin, zmax);

        if (zmin >= curSplitZ - 1e-4) {
            aboveIds.insert(id);
        }
        else if (zmax <= curSplitZ + 1e-4) {
            belowIds.insert(id);
        }
        else {
            facesToSplit.insert(id); // 跨越切分点的侧面
        }
    }

    // --- 2. 构造切分工具：无限大平面 ---
    gp_Pln cuttingPlane(gp_Pnt(0, 0, curSplitZ), gp_Dir(0, 0, 1));
    TopoDS_Face planeFace = BRepBuilderAPI_MakeFace(cuttingPlane);

    // --- 3. 核心切割循环 ---
    for (int id : facesToSplit) {
        TopoDS_Face faceToCut = TopoDS::Face(idToFace.Find(id));

        // A. 求取交线 (Section)
        BRepAlgoAPI_Section section(faceToCut, planeFace);
        section.Build();

        if (!section.IsDone() || !TopExp_Explorer(section.Shape(), TopAbs_EDGE).More()) {
            belowIds.insert(id); // 异常情况：归入下方
            continue;
        }

        // B. 使用 SplitShape 物理分割面
        BRepFeat_SplitShape splitter(faceToCut);
        TopExp_Explorer edgeExp(section.Shape(), TopAbs_EDGE);
        for (; edgeExp.More(); edgeExp.Next()) {
            splitter.Add(TopoDS::Edge(edgeExp.Current()), faceToCut);
        }
        splitter.Build();

        if (splitter.IsDone()) {
            // C. 收集并分类新产生的所有“子面”
            BRep_Builder debugBuilder;
            TopoDS_Compound splitResultComp;
            debugBuilder.MakeCompound(splitResultComp);

            TopExp_Explorer faceExp(splitter.Shape(), TopAbs_FACE);
            for (; faceExp.More(); faceExp.Next()) {
                TopoDS_Face subFace = TopoDS::Face(faceExp.Current());

                double smin, smax;
                GetFaceZRange(subFace, smin, smax);
                double zMid = (smin + smax) / 2.0; // 用中点判定归属最稳健

                int newId = nextId++;
                idToFace.Bind(newId, subFace);
                debugBuilder.Add(splitResultComp, subFace);

                if (zMid > curSplitZ) {
                    aboveIds.insert(newId);
                }
                else {
                    belowIds.insert(newId);
                }
            }

            // D. 可视化导出：这一层当前面被切开后的结果
            std::string debugPath = savePath + "Split_Z" + std::to_string((int)curSplitZ) + "_ID" + std::to_string(id) + ".brep";
            BRepTools::Write(splitResultComp, debugPath.c_str());
        }
        else {
            belowIds.insert(id); // 切割失败保留原始面
        }
    }

    // --- 4. 分类结果可视化 (Above vs Below) ---
    if (!aboveIds.empty() || !belowIds.empty()) {
        BRep_Builder resultBuilder;

        // 导出“上方集合” (当前层的结果)
        TopoDS_Compound compAbove;
        resultBuilder.MakeCompound(compAbove);
        for (int id : aboveIds) {
            resultBuilder.Add(compAbove, idToFace.Find(id));
        }
        std::string abovePath = savePath + "Layer_Z" + std::to_string((int)curSplitZ) + "_ABOVE_FINAL.brep";
        BRepTools::Write(compAbove, abovePath.c_str());

        // 导出“下方集合” (留给下一层递归的形状)
        TopoDS_Compound compBelow;
        resultBuilder.MakeCompound(compBelow);
        for (int id : belowIds) {
            resultBuilder.Add(compBelow, idToFace.Find(id));
        }
        std::string belowPath = savePath + "Layer_Z" + std::to_string((int)curSplitZ) + "_BELOW_REMAINING.brep";
        BRepTools::Write(compBelow, belowPath.c_str());

        cout << "  [Visualize] Z=" << curSplitZ << " 层提取完成: " << endl;
        cout << "    - 上方导出: " << abovePath << " (" << aboveIds.size() << " 个面)" << endl;
        cout << "    - 下方导出: " << belowPath << " (" << belowIds.size() << " 个面)" << endl;
    }

    // --- 4. 结果收集与递归 ---
    // 将当前层（切分线以上的所有面）存入结果
    if (!aboveIds.empty()) {
        result.push_back(aboveIds);
    }

    // 将切分线以下的所有面交给下一个切分点递归处理
    SplitCavityRecursive(belowIds, idToFace, nextId, splitZs, zIndex + 1, result);
}

// 入口函数
// 将返回值改为 Compound 数组
std::vector<TopoDS_Compound> SplitCavity(const TopoDS_Compound& cavity) {
    TopTools_DataMapOfIntegerShape idToFace;
    std::set<int> allIds;
    int nextId = 1;

    // 1. 初始化映射
    TopExp_Explorer exp(cavity, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        int id = nextId++;
        idToFace.Bind(id, TopoDS::Face(exp.Current()));
        allIds.insert(id);
    }

    // 2. 获取并排序切分点
    std::set<double> zSet;
    for (int id : allIds) {
        TopoDS_Face f = TopoDS::Face(idToFace.Find(id));
        if (IsHorizontalFace(f)) zSet.insert(GetHorizontalFaceZ(f));
    }

    std::vector<double> splitZs(zSet.begin(), zSet.end());
    std::sort(splitZs.begin(), splitZs.end(), std::greater<double>());
    if (!splitZs.empty()) splitZs.pop_back();

    // 3. 执行递归切分
    std::vector<std::set<int>> idGroups; // 存储各组的 ID
    SplitCavityRecursive(allIds, idToFace, nextId, splitZs, 0, idGroups);

    // 4. 【关键步骤】在 idToFace 销毁前，将 ID 转换为真正的几何体
    std::vector<TopoDS_Compound> finalShapes;
    BRep_Builder builder;

    for (const auto& group : idGroups) {
        TopoDS_Compound comp;
        builder.MakeCompound(comp);
        for (int id : group) {
            if (idToFace.IsBound(id)) {
                builder.Add(comp, idToFace.Find(id));
            }
        }
        finalShapes.push_back(comp);
    }

    return finalShapes; // 现在返回的是实实在在的形状
}

// 将一个包含多个独立区域的 Compound 拆分为多个独立的 Compound
std::vector<TopoDS_Compound> SeparateDisconnectedCavities(const TopoDS_Compound& cavityCompound) {
    std::vector<TopoDS_Compound> individualCavities;

    // 1. 建立 Edge -> Face 的映射，用于寻找邻接面
    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopExp::MapShapesAndAncestors(cavityCompound, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

    TopTools_MapOfShape visitedFaces;
    TopExp_Explorer exp(cavityCompound, TopAbs_FACE);

    for (; exp.More(); exp.Next()) {
        TopoDS_Face startFace = TopoDS::Face(exp.Current());
        if (visitedFaces.Contains(startFace)) continue;

        // 2. 发现一个新的连通区域，使用广度优先搜索 (BFS)
        BRep_Builder builder;
        TopoDS_Compound singleCavity;
        builder.MakeCompound(singleCavity);

        TopTools_ListOfShape queue;
        queue.Append(startFace);
        visitedFaces.Add(startFace);

        while (!queue.IsEmpty()) {
            TopoDS_Face currentFace = TopoDS::Face(queue.First());
            queue.RemoveFirst();
            builder.Add(singleCavity, currentFace);

            // 遍历当前面的所有边，找到邻接面
            TopExp_Explorer edgeExp(currentFace, TopAbs_EDGE);
            for (; edgeExp.More(); edgeExp.Next()) {
                const TopoDS_Shape& edge = edgeExp.Current();
                if (edgeFaceMap.Contains(edge)) {
                    const TopTools_ListOfShape& adjFaces = edgeFaceMap.FindFromKey(edge);
                    for (TopTools_ListIteratorOfListOfShape it(adjFaces); it.More(); it.Next()) {
                        const TopoDS_Shape& adjFace = it.Value();
                        if (!visitedFaces.Contains(adjFace)) {
                            visitedFaces.Add(adjFace);
                            queue.Append(adjFace);
                        }
                    }
                }
            }
        }
        individualCavities.push_back(singleCavity);
    }
    return individualCavities;
}


/**
 * 功能：将物理分离的型腔实体与所有切片的 2D 环进行匹配，重组为 CavityFeature 列表
 * 修正：bottomZ 自动向匹配下方最近的一个总切分点对齐
 */
std::vector<CavityFeature> GenerateCavityFeatures(
    const std::vector<TopoDS_Compound>& isolatedCavities,
    const std::vector<std::vector<Face2D>>& allLayerFaces,
    const std::vector<double>& sortedSplitPoints, // 传入降序排列的总切分点
    const TopTools_DataMapOfShapeInteger& faceToIdMap)
{
    std::vector<CavityFeature> finalFeatures;

    for (size_t i = 0; i < isolatedCavities.size(); ++i) {
        CavityFeature feat;
        feat.featureId = (int)i;
        feat.topZ = -1e9;   // 初始化极小值
        feat.bottomZ = 1e9; // 初始化极大值，稍后修正为下一个切分点
        double lastLoopZ = 1e9; // 记录探测到的物理环的最低 Z

        // --- 步骤 A: 收集当前型腔实体的面 ID 集合 ---
        TopExp_Explorer exp(isolatedCavities[i], TopAbs_FACE);
        for (; exp.More(); exp.Next()) {
            if (faceToIdMap.IsBound(exp.Current())) {
                feat.sourceFaceIds.insert(faceToIdMap.Find(exp.Current()));
            }
        }

        // --- 步骤 B: 空间 XY 包围盒计算 ---
        Bnd_Box box;
        BRepBndLib::Add(isolatedCavities[i], box);
        double x1, y1, zmin_box, x2, y2, zmax_box;
        box.Get(x1, y1, zmin_box, x2, y2, zmax_box);

        // --- 步骤 C: 全层扫描与匹配 ---
        for (const auto& layer : allLayerFaces) {
            for (const auto& f2d : layer) {
                // 1. 类型过滤：只看空腔类型的切片
                if (f2d.type != FaceType::CAVITY || f2d.outerLoop.empty()) continue;

                // 2. ID 匹配：检查该切片环的来源面 ID 是否属于本型腔
                int lineFaceId = f2d.outerLoop[0].faceId;
                if (feat.sourceFaceIds.count(lineFaceId)) {

                    // 3. 空间匹配：检查切片环的点是否在型腔实体的 XY 范围内
                    Point3D p = f2d.outerLoop[0].start;
                    if (p.x >= x1 - 0.5 && p.x <= x2 + 0.5 && p.y >= y1 - 0.5 && p.y <= y2 + 0.5) {

                        double curZ = p.z;

                        // 转换外环
                        CavityLoop out;
                        out.lines = f2d.outerLoop;
                        out.zHeight = curZ;
                        out.isOuter = true;
                        feat.stepLoops.push_back(out);

                        // 转换内环（孔洞/孤岛）
                        for (const auto& inner : f2d.innerLoops) {
                            CavityLoop in;
                            in.lines = inner;
                            in.zHeight = curZ;
                            in.isOuter = false;
                            feat.stepLoops.push_back(in);
                        }

                        // 更新探测到的 Z 极限
                        feat.topZ = std::max(feat.topZ, curZ);
                        lastLoopZ = std::min(lastLoopZ, curZ);
                    }
                }
            }
        }

        // --- 步骤 D: 特征收尾与底面修正 ---
        if (!feat.stepLoops.empty()) {
            // 修正 bottomZ：寻找检测到的最低环 Z 对应的下一个切分点
            feat.bottomZ = lastLoopZ; // 默认值
            for (size_t k = 0; k < sortedSplitPoints.size(); ++k) {
                // 如果当前点与最低环高度重合（带容差）
                if (std::abs(sortedSplitPoints[k] - lastLoopZ) < 1e-4) {
                    // 如果存在更低的一个切分点，则那个点才是型腔的物理底面
                    if (k + 1 < sortedSplitPoints.size()) {
                        feat.bottomZ = sortedSplitPoints[k + 1];
                    }
                    break;
                }
            }

            feat.totalDepth = std::abs(feat.topZ - feat.bottomZ);

            // 按高度从高到低排序，模拟加工顺序
            std::sort(feat.stepLoops.begin(), feat.stepLoops.end(),
                [](const CavityLoop& a, const CavityLoop& b) { return a.zHeight > b.zHeight; });

            finalFeatures.push_back(feat);
        }
    }

    // ==================== 详细信息打印开始 ====================
    cout << "\n============================================================" << endl;
    cout << "                型腔特征识别报告 (底面修正版)                 " << endl;
    cout << "============================================================" << endl;

    for (const auto& feat : finalFeatures) {
        cout << ">>> [型腔 ID: " << feat.featureId << "]" << endl;
        cout << "    - 空间范围: Z_Top = " << feat.topZ << " / Z_Floor = " << feat.bottomZ << endl;
        cout << "    - 加工深度: " << feat.totalDepth << endl;
        cout << "    - 进刀方向: [0, 0, -1]" << endl;

        // 打印原始面 ID 清单
        cout << "    - 关联原始面 ID (" << feat.sourceFaceIds.size() << "个): ";
        for (int oid : feat.sourceFaceIds) cout << oid << " ";
        cout << endl;

        // 打印层级环信息
        cout << "    - 层级切片环清单 (按 Z 降序):" << endl;

        double lastZ = 999999.9;
        for (size_t k = 0; k < feat.stepLoops.size(); ++k) {
            const auto& loop = feat.stepLoops[k];

            if (fabs(loop.zHeight - lastZ) > 1e-4) {
                cout << "      [高度 Z = " << loop.zHeight << "]" << endl;
                lastZ = loop.zHeight;
            }

            string typeStr = loop.isOuter ? "外边界 (Outer)" : "孤岛/内孔 (Inner)";
            cout << "        * " << typeStr << ": 线段数 " << loop.lines.size();

            if (!loop.lines.empty()) {
                const auto& startPt = loop.lines.front().start;
                const auto& endPt = loop.lines.back().end;
                cout << " | 闭合性检校: 起(" << startPt.x << "," << startPt.y
                    << ") -> 终(" << endPt.x << "," << endPt.y << ")";
            }
            cout << endl;
        }
        cout << "------------------------------------------------------------" << endl;
    }
    cout << "============================================================\n" << endl;
    // ==================== 详细信息打印结束 ====================

    return finalFeatures;
}


// 根据识别到的型腔特征，从原始模型中提取真正属于型腔的面，形成新的 Compound
TopoDS_Compound ExtractTrueCavityFaces(
    const CavityFeature& feat,
    const TopoDS_Shape& mainShape,
    const TopTools_DataMapOfShapeInteger& faceToIdMap)
{
    BRep_Builder builder;
    TopoDS_Compound refinedResult;
    builder.MakeCompound(refinedResult);

    // 1. 准备“开口环（最上方环）”用于子面的水平位置校验
    TopoDS_Wire topWire;
    if (!feat.stepLoops.empty()) {
        const auto& tLoop = feat.stepLoops.front();
        BRepBuilderAPI_MakePolygon poly;
        for (const auto& l : tLoop.lines) poly.Add(gp_Pnt(l.start.x, l.start.y, l.start.z));
        const auto& lastL = tLoop.lines.back();
        poly.Add(gp_Pnt(lastL.end.x, lastL.end.y, lastL.end.z));
        if (poly.IsDone()) topWire = poly.Wire();
    }

    if (topWire.IsNull()) return refinedResult;

    TopExp_Explorer exp(mainShape, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face originalFace = TopoDS::Face(exp.Current());
        if (!faceToIdMap.IsBound(originalFace)) continue;

        int oid = faceToIdMap.Find(originalFace);
        // 只处理属于该型腔特征关联的面
        if (feat.sourceFaceIds.count(oid)) {

            double zmin, zmax;
            GetFaceZRange(originalFace, zmin, zmax);

            // --- 策略 A：处理跨越切分点(topZ)的面 ---
            if (zmax > feat.topZ + 1e-4 && zmin < feat.topZ - 1e-4) {

                gp_Pln splitPln(gp_Pnt(0, 0, feat.topZ), gp_Dir(0, 0, 1));
                TopoDS_Face algoPlane = BRepBuilderAPI_MakeFace(splitPln);
                BRepAlgoAPI_Section section(originalFace, algoPlane);
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

                            double sMin, sMax;
                            GetFaceZRange(subF, sMin, sMax);

                            // 校验 1: 子面必须位于 topZ 以下
                            if (sMax <= feat.topZ + 1e-4) {
                                // 校验 2: 子面必须在水平方向上靠近型腔开口环
                                BRepExtrema_DistShapeShape distCalc(subF, topWire);
                                if (distCalc.IsDone() && distCalc.Value() < 1e-2) {
                                    builder.Add(refinedResult, subF);
                                }
                            }
                        }
                    }
                }
            }
            // --- 策略 B：处理完全在 topZ 以下的面 (不在待切分列表中的面) ---
            // 只要面在 topZ 以下且属于 sourceFaceIds，则全部直接加入
            else if (zmax <= feat.topZ + 1e-4) {
                builder.Add(refinedResult, originalFace);
                // 打印 Debug 信息，确认这包含水平面和其他非跨层侧面
                // cout << ">>>> [直接加入] 面 ID: " << oid << " 完全在高度范围内。" << endl;
            }
        }
    }

    return refinedResult;
}




#if 1
int main() {
    std::string stepfile = "8_stp.stp";

    STEPControl_Reader reader;
    std::string inputFileName = inputPath + stepfile;
    if (reader.ReadFile(inputFileName.c_str()) != IFSelect_RetDone) {
        cout << "无法读取 STEP 文件！" << endl;
        return 1;
    }
    reader.TransferRoots();
    TopoDS_Shape mainShape = reader.OneShape();
    
    cout << "模型已加载，开始自适应获取切分点..." << endl;
    
    // 为模型中的每个面分配一个唯一的 ID
    TopTools_DataMapOfShapeInteger faceToIdMap;
    int currentFaceId = 1;
    TopExp_Explorer expFace(mainShape, TopAbs_FACE);
    for (; expFace.More(); expFace.Next()) {
        if (!faceToIdMap.IsBound(expFace.Current())) {
            faceToIdMap.Bind(expFace.Current(), currentFaceId++);
        }
    }
    cout << "已为 " << currentFaceId - 1 << " 个面分配了 ID。" << endl;
    
    std::vector<double> splitPoints = GetSplitPointsAlongZ(mainShape);
    
    cout << "检测到 " << splitPoints.size() << " 个切分点:" << endl;
    for (int i = 0; i < splitPoints.size(); i++) {
        cout << "  [" << i << "] Z = " << splitPoints[i] << endl;
    }
    
    // ================= 新增：建立最大毛坯包围体并切分 =================
    cout << "\n=== 建立最大毛坯包围体并分层切割 ===" << endl;
    Bnd_Box totalBox;
    BRepBndLib::Add(mainShape, totalBox);
    if (!totalBox.IsVoid()) {
        double xMin, yMin, zMin, xMax, yMax, zMax;
        totalBox.Get(xMin, yMin, zMin, xMax, yMax, zMax);

        // 给 X/Y 方向增加一点加工余量 (可选，这里默认为 0)
        double offset = 2.0;
        xMin -= offset; yMin -= offset;
        xMax += offset; yMax += offset;

        // 创建长方体毛坯
        TopoDS_Shape stockShape = BRepPrimAPI_MakeBox(gp_Pnt(xMin, yMin, zMin), gp_Pnt(xMax, yMax, zMax)).Shape();
        std::string StockFileName = savePath + "Stock_Body.brep";
        BRepTools::Write(stockShape, StockFileName.c_str());
        cout << "  已生成整体毛坯实体: Stock_Body.brep" << endl;

        // 对毛坯进行分层切割
        for (int i = (int)splitPoints.size() - 1; i >= 0; i--) {
            double splitZ = splitPoints[i];
            int reversedIndex = splitPoints.size() - 1 - i;
            
            gp_Pln cuttingPlane(gp_Pnt(0, 0, splitZ), gp_Dir(0, 0, 1));
            TopoDS_Face algoPlane = BRepBuilderAPI_MakeFace(cuttingPlane);

            BRepAlgoAPI_Section section(stockShape, algoPlane, Standard_True);
            section.Build();

            std::vector<OneLine> stockLines;
            if (section.IsDone()) {
                TopExp_Explorer expEdge(section.Shape(), TopAbs_EDGE);
                for (; expEdge.More(); expEdge.Next()) {
                    TopoDS_Edge E = TopoDS::Edge(expEdge.Current());
                    BRepAdaptor_Curve bac(E);
                    GCPnts_QuasiUniformDeflection discretizer(bac, 0.01);
                    if (discretizer.IsDone()) {
                        for (int k = 1; k < discretizer.NbPoints(); ++k) {
                            gp_Pnt p1 = discretizer.Value(k);
                            gp_Pnt p2 = discretizer.Value(k + 1);
                            OneLine line;
                            line.start.x = p1.X(); line.start.y = p1.Y(); line.start.z = p1.Z();
                            line.end.x = p2.X();   line.end.y = p2.Y();   line.end.z = p2.Z();
                            line.faceId = -888; 
                            stockLines.push_back(line);
                        }
                    }
                }
            }

            std::vector<Face2D> stockFaces = BuildTopologyAndExtractFaces(stockLines);
            std::vector<Face2D> stockSolids;
            for (const auto& f : stockFaces) {
                if (f.type == FaceType::SOLID) stockSolids.push_back(f);
            }

            std::string stockName = savePath + "Stock_Slice_" + std::to_string(reversedIndex) + ".brep";
            ExportFace2DToBrep(stockSolids, stockName);
        }
        cout << "  毛坯分层切割完成！" << endl;
    }
    // =================================================================

    cout << "\n开始逐层切分模型 (从上往下)..." << endl;
    
    // 用于累加最大轮廓面
    //std::vector<Face2D> maxSilhouetteFaces;
    
    // 收集所有层的所有数据用于后续特征识别
    std::vector<std::vector<Face2D>> allLayerFaces;
    
    for (int i = (int)splitPoints.size() - 1; i >= 0; i--) {
        double splitZ = splitPoints[i];
        int reversedIndex = splitPoints.size() - 1 - i;
        cout << "\n切分 (倒序 #" << reversedIndex << "/" << splitPoints.size() << ") (Z = " << splitZ << ")" << endl;
        
        if (i == 4)
        {
            // 提取本层的面
            std::vector<Face2D> currentLayerFaces = SliceModelAtZ(mainShape, splitZ, i, splitPoints.size()-1, faceToIdMap);

            // 将本层面数据保存到全局列表
            allLayerFaces.push_back(currentLayerFaces);

            // 1. 分离出实体面参与运算
            std::vector<Face2D> currentSolidFaces;
            std::vector<Face2D> currentConvexHullFaces;
            for (const auto& f : currentLayerFaces) {
                if (f.type == FaceType::SOLID) {
                    currentSolidFaces.push_back(f);
                    Face2D HullFace = ComputeConvexHullFace(f);
                    currentConvexHullFaces.push_back(HullFace);
                }
            }

            std::vector<Face2D> currentCavityFaces;
            for (const auto& f : currentLayerFaces) {
                if (f.type == FaceType::CAVITY) {
                    currentCavityFaces.push_back(f);
                }
            }

			//DEBUG: 导出当前层的实体面和腔面，检查切分结果
            std::string currentHullFaceName = savePath + "Slice_" + std::to_string(i) + "currentHullFace.brep";
            ExportFace2DToBrep(currentConvexHullFaces, currentHullFaceName);
            std::string currentSolidFaceName = savePath + "Slice_" + std::to_string(i) + "currentSolidFace.brep";
            ExportFace2DToBrep(currentSolidFaces, currentSolidFaceName);
            std::string currentCavityFaceName = savePath + "Slice_" + std::to_string(i) + "currentCavityFacee.brep";
            ExportFace2DToBrep(currentCavityFaces, currentCavityFaceName);
            int aaa = 0;
        }

    }
    
    // ================= 提取封闭型腔特征面 =================
    cout << "\n=== 开始提取封闭型腔特征面 ===" << endl;
    std::set<int> cavityFaceIds = CollectCavityFaceIds(allLayerFaces);
    cout << "  检测到 " << cavityFaceIds.size() << " 个型腔面 ID" << endl;
    
    // 从原始实体模型中提取完整封闭型腔（侧壁 + 顶/底面/中间平台）
    std::string cavityFacesFileName = savePath + "CavityWallFaces.brep";
    ExportCavityFaces(mainShape, cavityFaceIds, faceToIdMap, cavityFacesFileName);
    cout << "  完整封闭型腔提取完成！" << endl;

    // 获取型腔合并面 Compound，供后续切割使用
    TopoDS_Compound cavityCompound = GetCavityCompound(mainShape, cavityFaceIds, faceToIdMap);
    
    // 拆分为独立的型腔区域 (独立连通体)
    std::vector<TopoDS_Compound> isolatedCavities = SeparateDisconnectedCavities(cavityCompound);
    cout << "识别到 " << isolatedCavities.size() << " 个独立的型腔区域。" << endl;

	// 得到型腔特征的表示（包含每个型腔对应的切片环信息）
    std::vector<CavityFeature> cavityFeatures = GenerateCavityFeatures(isolatedCavities, allLayerFaces, splitPoints, faceToIdMap);

    for (const auto& feat : cavityFeatures) {
        // 提取真正的型腔几何面
        TopoDS_Compound trueCavity = ExtractTrueCavityFaces(feat, mainShape, faceToIdMap);

        // 保存可视化
        string filePath = savePath + "Final_TrueCavity_" + to_string(feat.featureId) + ".brep";
        BRepTools::Write(trueCavity, filePath.c_str());
    }

    // 3. 对每个独立型腔分别进行 Z 轴切分
    int totalPartCount = 0;
    for (size_t i = 0; i < isolatedCavities.size(); ++i) {
        cout << "正在处理第 " << i << " 个型腔..." << endl;

        // 对当前独立区域调用你的 SplitCavity
        std::vector<TopoDS_Compound> parts = SplitCavity(isolatedCavities[i]);

        // 4. 保存结果
        for (size_t j = 0; j < parts.size(); ++j) {
            std::string fileName = savePath + "Cavity_" + std::to_string(i) + "_Part_" + std::to_string(j) + ".brep";
            BRepTools::Write(parts[j], fileName.c_str());
            totalPartCount++;
        }
    }
    cout << "全部处理完成，共生成 " << totalPartCount << " 个型腔零件。" << endl;

    
    cout << "\n所有切分完成！" << endl;
    cout << "生成的 BREP 文件保存在当前目录。" << endl;
    

    std::string AllFacesFileName = savePath + "AllFaces.brep";
    BRepTools::Write(gAllFacesCompound, AllFacesFileName.c_str());
    std::string AllLinesFileName = savePath + "AllLines.brep";
    BRepTools::Write(gAllLinesCompound, AllLinesFileName.c_str());

    cout << "\n合并文件已生成:" << endl;
    cout << "  AllFaces.brep - 所有切分面的集合" << endl;
    cout << "  AllLines.brep - 所有切分线的集合" << endl;
    cout << "\n使用 BrepViewer.exe 打开当前目录查看结果。" << endl;
    
    system("pause");
    return 0;
}
#endif