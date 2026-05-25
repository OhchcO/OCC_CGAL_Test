#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
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
#include <BRepBuilderAPI_Sewing.hxx>
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
    CAVITY,  // 封闭空腔面（无材料/孔洞）
	OPENCAVITY, // 开口空腔面（边界开口的孔洞）
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

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/number_utils.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
typedef CGAL::Point_2<EK> Point_2;

// 将 Face2D 转换为 2D 点集合
std::vector<Point_2> ConvertFaceToPoints2D_EK(const Face2D& face) {
    std::vector<Point_2> points_2d;

    for (const auto& line : face.outerLoop) {
        points_2d.emplace_back(line.start.x, line.start.y);
        points_2d.emplace_back(line.end.x, line.end.y);
    }
    return points_2d;
}

/**
 * @brief 【核心还原函数】将 CGAL 计算出的凸包高精点集，重新包装为你原生的 Face2D 凸包结构
 * @param hull_points CGAL 计算出的凸包轮廓顶点序列
 * @param zHeight 当前切层的绝对高度 Z 值
 * @return Face2D 组装好外环（outerLoop）的 HULL 类型面
 */
Face2D ConvertPointsToHullFace_EK(
    vector<Point_2>& hull_points,
    double zHeight)
{
    Face2D resultFace;
    resultFace.type = FaceType::HULL;      // 明确面类型为 HULL（凸包毛坯面）
    resultFace.innerLoops.clear();         // 明确凸包没有嵌套内孔

    if (hull_points.empty()) {
        return resultFace;
    }

    // 顺着凸包的顶点序列，首尾相连重建你原生的 OneLine 闭合线段环
    for (size_t i = 0; i < hull_points.size(); ++i) {
        OneLine line;

        // 1. 将当前折点提取出来作为线段的起点
        // 使用 CGAL::to_double 将高精度的有理数/代数数坐标，安全转换为你结构体需要的 double 物理坐标
        line.start = {
            CGAL::to_double(hull_points[i].x()),
            CGAL::to_double(hull_points[i].y()),
            zHeight
        };

        // 2. 将下一个折点作为线段的终点
        // 当循环到最后一个点时，(i + 1) % size 会自动变成 0，从而全自动连回起点，完美闭合！
        size_t nextIdx = (i + 1) % hull_points.size();
        line.end = {
            CGAL::to_double(hull_points[nextIdx].x()),
            CGAL::to_double(hull_points[nextIdx].y()),
            zHeight
        };

        line.faceId = -1; // 标记为凸包逻辑算法生成的虚拟边缘

        // 塞入外环数组
        resultFace.outerLoop.push_back(line);
    }

    return resultFace;
}

/**
 * @brief 对 Face2D 的外环计算二维凸包，并将结果封装回一个新的 Face2D 结构
 * @param face 输入的原始面（SOLID 类型）
 * @return Face2D 类型为 HULL 的新面，包含封闭的凸包外环
 */
Face2D ComputeConvexHullFace(const Face2D& face) {

    double currentZ = face.outerLoop[0].start.z;

    // 1. 转成Point_2 的格式，准备计算二维凸包
    vector<Point_2> points_2d = ConvertFaceToPoints2D_EK(face);

    // 2. 计算 2D 凸包
    std::vector<Point_2> hull_points;
    CGAL::convex_hull_2(points_2d.begin(), points_2d.end(), std::back_inserter(hull_points));

    // 3. 还原为 OneLine 闭合环并存入 outerLoop
    Face2D resultFace = ConvertPointsToHullFace_EK(hull_points, currentZ);

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
void GetExtremaZOfFace(const TopoDS_Face& face, const gp_Dir& direction, double& outMinZ, double& outMaxZ)
{
    double minProjectedValue = std::numeric_limits<double>::max();
    double maxProjectedValue = -std::numeric_limits<double>::max();
    gp_Pnt lowestPoint(0, 0, 0), highestPoint(0, 0, 0);
    bool found = false;

    // 遍历当前面的所有边界（ACIS 风格：边缘曲线投影）
    TopExp_Explorer edgeExp(face, TopAbs_EDGE);
    for (; edgeExp.More(); edgeExp.Next()) {
        TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());
        BRepAdaptor_Curve curve(edge);

        double firstParam = curve.FirstParameter();
        double lastParam = curve.LastParameter();
        gp_Pnt pFirst = curve.Value(firstParam);
        gp_Pnt pLast = curve.Value(lastParam);

        // 利用 .Dot() 计算点在 Z 轴上的纯解析投影距离
        double distFirst = gp_Vec(pFirst.XYZ()).Dot(gp_Vec(direction.XYZ()));
        double distLast = gp_Vec(pLast.XYZ()).Dot(gp_Vec(direction.XYZ()));

        if (distFirst < minProjectedValue) { minProjectedValue = distFirst; lowestPoint = pFirst; found = true; }
        if (distFirst > maxProjectedValue) { maxProjectedValue = distFirst; highestPoint = pFirst; found = true; }
        if (distLast < minProjectedValue) { minProjectedValue = distLast;  lowestPoint = pLast;  found = true; }
        if (distLast > maxProjectedValue) { maxProjectedValue = distLast;  highestPoint = pLast;  found = true; }

        // 如果边界是圆弧、样条曲线，启动一维黄金分割精密搜索切点导数零点
        if (curve.GetType() != GeomAbs_Line) {
            double uMin = firstParam, uMax = lastParam;
            const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
            const double resphi = 2.0 - phi;
            double u1 = uMin + resphi * (uMax - uMin), u2 = uMax - resphi * (uMax - uMin);
            gp_Pnt p1 = curve.Value(u1), p2 = curve.Value(u2);
            double f1 = gp_Vec(p1.XYZ()).Dot(gp_Vec(direction.XYZ())), f2 = gp_Vec(p2.XYZ()).Dot(gp_Vec(direction.XYZ()));

            for (int i = 0; i < 15; ++i) {
                if (f1 < f2) { uMax = u2; u2 = u1; p2 = p1; f2 = f1; u1 = uMin + resphi * (uMax - uMin); p1 = curve.Value(u1); f1 = gp_Vec(p1.XYZ()).Dot(gp_Vec(direction.XYZ())); }
                else { uMin = u1; u1 = u2; p1 = p2; f1 = f2; u2 = uMax - resphi * (uMax - uMin); p2 = curve.Value(u2); f2 = gp_Vec(p2.XYZ()).Dot(gp_Vec(direction.XYZ())); }
            }
            gp_Pnt pExtrema = curve.Value((uMin + uMax) * 0.5);
            double distExtrema = gp_Vec(pExtrema.XYZ()).Dot(gp_Vec(direction.XYZ()));
            if (distExtrema < minProjectedValue) { minProjectedValue = distExtrema; lowestPoint = pExtrema; found = true; }
            if (distExtrema > maxProjectedValue) { maxProjectedValue = distExtrema; highestPoint = pExtrema; found = true; }
        }
    }
    if (!found) { BRepAdaptor_Surface surf(face); lowestPoint = surf.Value(surf.FirstUParameter(), surf.FirstVParameter()); highestPoint = lowestPoint; }
    outMinZ = lowestPoint.Z(); outMaxZ = highestPoint.Z();
}

// =================================================================
// 主功能函数：通用的自适应 Z 轴切分点生成引擎
// =================================================================
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
std::vector<double> GetSplitPointsAlongZ(const TopoDS_Shape& shape,
    double angleTolerance = 0.001,
    double mergeTol = 0.05,
    double areaThreshold = 1.0)
{
    if (shape.IsNull()) return {};

    std::set<double> zPoints;
    std::set<double> masterPlanes; // 水平大平面骨架白名单
    const gp_Dir zAxis(0, 0, 1);

    // -------------------------------------------------------------
    // 【第一轮遍历】：精确提取所有水平大平面基准
    // -------------------------------------------------------------
    TopExp_Explorer exp(shape, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        BRepAdaptor_Surface surf(face);

        if (surf.GetType() == GeomAbs_Plane) {
            gp_Pln plane = surf.Plane();
            gp_Dir normal = plane.Axis().Direction();
            if (face.Orientation() == TopAbs_REVERSED) normal.Reverse();

            double angle = normal.Angle(zAxis);
            if (angle < angleTolerance || std::fabs(angle - M_PI) < angleTolerance) {
                double zVal = plane.Location().Z();
                zPoints.insert(zVal);
                masterPlanes.insert(zVal); // 登入白名单
            }
        }
    }

    // -------------------------------------------------------------
    // 【第二轮遍历】：处理非水平特征（应用 ACIS 风格的面积与边界解析过滤）
    // -------------------------------------------------------------
    const double planeNoiseTol = 0.4; // 贴近大平面的网格过渡噪声盲区隔绝带

    exp.ReInit();
    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current()); // 如果你的Explorer叫exp，这里保持一致
        BRepAdaptor_Surface surf(face);

        // 如果是第一步处理过的水平面，直接跳过
        bool isHorizontal = false;
        if (surf.GetType() == GeomAbs_Plane) {
            gp_Pln plane = surf.Plane(); gp_Dir normal = plane.Axis().Direction();
            if (face.Orientation() == TopAbs_REVERSED) normal.Reverse();
            double angle = normal.Angle(zAxis);
            if (angle < angleTolerance || std::fabs(angle - M_PI) < angleTolerance) isHorizontal = true;
        }
        if (isHorizontal) continue;

        // 🔘 工业级核心滤镜一：基于通用表面积的碎面拦截
        // 彻底根治 51.242 这种由于建模、缝合产生的不具备切层价值的局部碎面噪声
        GProp_GProps gprops;
        BRepGProp::SurfaceProperties(face, gprops);
        if (gprops.Mass() < areaThreshold) {
            continue;
        }

        // 🔘 工业级核心滤镜二：基于边界曲线的纯解析极值引擎（100%代替Bnd_Box）
        // 彻底消灭 51.58、51.52 等由于网格离散没能掉到谷底的弦高误差
        double realZmin = 0.0;
        double realZmax = 0.0;
        GetExtremaZOfFace(face, zAxis, realZmin, realZmax);

        // 🔘 工业级核心滤镜三：主平面保护带智能关联
        bool isMinNoise = false;
        bool isMaxNoise = false;
        for (double mp : masterPlanes) {
            if (std::abs(realZmin - mp) > 1e-5 && std::abs(realZmin - mp) < planeNoiseTol) isMinNoise = true;
            if (std::abs(realZmax - mp) > 1e-5 && std::abs(realZmax - mp) < planeNoiseTol) isMaxNoise = true;
        }

        if (realZmin > 0.0 && !isMinNoise) zPoints.insert(realZmin);
        if (realZmax > 0.0 && !isMaxNoise) zPoints.insert(realZmax);
    }

    if (zPoints.empty()) return {};

    // -------------------------------------------------------------
    // 【第三阶段】：由高到低，多特征并排智能特征去重
    // -------------------------------------------------------------
    std::vector<double> sortedZ(zPoints.begin(), zPoints.end());
    std::sort(sortedZ.begin(), sortedZ.end(), std::greater<double>()); // 从大到小降序

    std::vector<double> finalPoints;

    // 核心工艺间距：因为你的碗状模型内外壁 51.1737 和 51.1219 相差 0.0518mm。
    // 我们把去重分辨率降到 0.02mm，这样间距大于 0.02mm 的真实物理特征都会并排独立保留！
    const double precisionResolution = 0.02;

    for (double z : sortedZ) {
        // 白名单大平面优先精准校准对齐（将 52.001 强制纠偏成标准的 52.0）
        for (double mp : masterPlanes) {
            if (std::abs(z - mp) < mergeTol) {
                z = mp;
                break;
            }
        }

        if (finalPoints.empty()) {
            finalPoints.push_back(z);
        }
        else {
            // 各回各家条件：只有落差跨越了工艺临界值，才被允许登记为新切层高度
            if (finalPoints.back() - z > precisionResolution) {
                finalPoints.push_back(z);
            }
        }
    }

    // 最后一轮清理：擦除由于校准产生的重复大平面连续项
    finalPoints.erase(std::unique(finalPoints.begin(), finalPoints.end(),
        [](double a, double b) { return std::abs(a - b) < 0.05; }), finalPoints.end());

    return finalPoints;
}

//std::vector<double> GetSplitPointsAlongZ(const TopoDS_Shape& shape,
//    double angleTolerance = 0.001,
//    double mergeTol = 1e-3)
//{
//    std::set<double> zPoints;
//    const gp_Dir zAxis(0, 0, 1);
//
//    TopExp_Explorer exp(shape, TopAbs_FACE);
//    for (; exp.More(); exp.Next()) {
//        TopoDS_Face face = TopoDS::Face(exp.Current());
//        BRepAdaptor_Surface surf(face);
//
//        // ==============================
//        // 第一步：判断是不是水平面
//        // ==============================
//        bool isHorizontal = false;
//        double zPlane = 0.0;
//
//        if (surf.GetType() == GeomAbs_Plane) {
//            gp_Pln plane = surf.Plane();
//            gp_Dir normal = plane.Axis().Direction();
//
//            if (face.Orientation() == TopAbs_REVERSED)
//                normal.Reverse();
//
//            double angle = normal.Angle(zAxis);
//            if (angle < angleTolerance || fabs(angle - M_PI) < angleTolerance) {
//                isHorizontal = true;
//                zPlane = plane.Location().Z();
//            }
//        }
//
//        // ==============================
//        // 规则 1：水平面 → 加 Z
//        // ==============================
//        if (isHorizontal) {
//            zPoints.insert(zPlane);
//        }
//        // ==============================
//        // 规则 2：所有其他面（竖直平面、斜面、圆角、曲面）→ 加 最低点 Z
//        // ==============================
//        else {
//            Bnd_Box box;
//            BRepBndLib::Add(face, box);
//            if (!box.IsVoid()) {
//                double xmin, ymin, zmin, xmax, ymax, zmax;
//                box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
//                zPoints.insert(zmin);
//            }
//        }
//    }
//
//    // 从大到小排序
//    std::vector<double> sortedZ(zPoints.begin(), zPoints.end());
//    std::sort(sortedZ.begin(), sortedZ.end(), std::greater<double>());
//
//    // 按容差去重
//    std::vector<double> finalPoints;
//    for (double z : sortedZ) {
//        if (finalPoints.empty()) {
//            finalPoints.push_back(z);
//        }
//        else {
//            if (finalPoints.back() - z > mergeTol) {
//                finalPoints.push_back(z);
//            }
//        }
//    }
//
//    return finalPoints;
//}

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
    auto lines = DeduplicateLines(inputLines);
    if (lines.empty()) return loops;

    vector<bool> visited(lines.size(), false);
    const double snapTol = 1e-3;

    auto Snap = [](const Point3D& p, double tol) -> std::pair<int64_t, int64_t> {
        return { (int64_t)round(p.x / tol), (int64_t)round(p.y / tol) };
    };

    std::unordered_map<int64_t, std::vector<size_t>> endMap;
    for (size_t i = 0; i < lines.size(); ++i) {
        auto sk = Snap(lines[i].start, snapTol);
        auto ek = Snap(lines[i].end, snapTol);
        int64_t startKey = sk.first * 1000000007LL + sk.second;
        int64_t endKey = ek.first * 1000000007LL + ek.second;
        endMap[startKey].push_back(i);
        if (startKey != endKey) {
            endMap[endKey].push_back(i);
        }
    }

    for (size_t i = 0; i < lines.size(); ++i) {
        if (visited[i]) continue;

        Loop currentLoop;
        currentLoop.push_back(lines[i]);
        visited[i] = true;
        Point3D startPt = lines[i].start;
        Point3D currEnd = lines[i].end;

        while (true) {
            auto key = Snap(currEnd, snapTol);
            int64_t mapKey = key.first * 1000000007LL + key.second;

            size_t bestIdx = SIZE_MAX;
            double bestDist = 1e9;
            bool bestReverse = false;

            auto it = endMap.find(mapKey);
            if (it != endMap.end()) {
                for (size_t idx : it->second) {
                    if (visited[idx]) continue;
                    const auto& l = lines[idx];
                    double dS = fabs(currEnd.x - l.start.x) + fabs(currEnd.y - l.start.y);
                    double dE = fabs(currEnd.x - l.end.x) + fabs(currEnd.y - l.end.y);
                    if (dS < snapTol && dS < bestDist) {
                        bestDist = dS; bestIdx = idx; bestReverse = false;
                    }
                    if (dE < snapTol && dE < bestDist) {
                        bestDist = dE; bestIdx = idx; bestReverse = true;
                    }
                }
            }

            if (bestIdx == SIZE_MAX) break;

            visited[bestIdx] = true;
            if (bestReverse) {
                OneLine rev = lines[bestIdx];
                swap(rev.start, rev.end);
                currentLoop.push_back(rev);
                currEnd = rev.end;
            } else {
                currentLoop.push_back(lines[bestIdx]);
                currEnd = lines[bestIdx].end;
            }

            if (IsPointEqual(currEnd, startPt, snapTol)) break;
        }

        if (IsPointEqual(currentLoop.back().end, currentLoop.front().start, snapTol)
            && currentLoop.size() >= 3)
        {
            loops.push_back(currentLoop);
        }
    }

    cout << "\n================ 最终结果 =================" << endl;
    cout << " 成功生成环数量：" << loops.size() << endl;
    for (size_t i = 0; i < loops.size(); i++) {
        cout << "  环 " << i + 1 << "：" << loops[i].size() << " 条线段" << endl;
    }
    int usedCount = 0;
    for (bool v : visited) if (v) usedCount++;
    cout << "  输入线段: " << lines.size() << " 已用: " << usedCount << " 未用: " << lines.size() - usedCount << endl;
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

// 获取 Face2D 的重心（基于外环顶点）
Point3D GetFace2DCenter(const Face2D& face) {
    Point3D center = {0.0, 0.0, 0.0};
    int ptCount = 0;
    for (const auto& line : face.outerLoop) {
        center.x += line.start.x;
        center.y += line.start.y;
        center.z += line.start.z;
        ptCount++;
    }
    if (ptCount > 0) {
        center.x /= ptCount;
        center.y /= ptCount;
        center.z /= ptCount;
    }
    return center;
}

// 判断一个 Face2D 是否完全在另一个 Face2D 的外环内部
bool IsFaceInsideFace(const Face2D& innerFace, const Face2D& outerFace) {
    // 检查 innerFace 的所有顶点是否都在 outerFace 的外环内部
    for (const auto& line : innerFace.outerLoop) {
        if (!IsPointInLoop(line.start, outerFace.outerLoop)) {
            return false; // 有顶点在外环外部
        }
    }
    return true; // 所有顶点都在内部
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
    //PrintAllLines(layerLines);
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
        // 取当前环的重心作为测试点（比只取第一个点更稳健）
        Point3D center = {0.0, 0.0, 0.0};
        int ptCount = 0;
        for (const auto& l : nodes[i]->loopLines) {
            center.x += l.start.x;
            center.y += l.start.y;
            center.z += l.start.z;
            ptCount++;
        }
        if (ptCount > 0) {
            center.x /= ptCount;
            center.y /= ptCount;
            center.z /= ptCount;
        }

        bool foundParent = false;
        // 从比它大一点的环开始往上找（即从 i-1 倒序遍历到 0）
        // 找到的第一个包含它的环，就是它的直接父节点（面积最小的包围者）
        for (int j = (int)i - 1; j >= 0; --j) {
            // 用重心测试：重心在多边形内部的概率远高于任意一个顶点
            if (IsPointInLoop(center, nodes[j]->loopLines)) {
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
                    if (edgeToFaceMap.IsBound(E)) {
                        // 两条边重叠：选择更低的面（更贴近 splitZ 的面）
                        TopoDS_Face existingFace = TopoDS::Face(edgeToFaceMap.Find(E));
                        double eZmin, eZmax, fZmin, fZmax;
                        GetFaceZRange(existingFace, eZmin, eZmax);
                        GetFaceZRange(face, fZmin, fZmax);
                        if (fZmax < eZmax) {
                            edgeToFaceMap.Bind(E, face); // 新面更低，替换
                        }
                    } else {
                        edgeToFaceMap.Bind(E, face);
                    }
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
        //std::string abovePath = savePath + "Layer_Z" + std::to_string((int)curSplitZ) + "_ABOVE_FINAL.brep";
        //BRepTools::Write(compAbove, abovePath.c_str());

        // 导出“下方集合” (留给下一层递归的形状)
        TopoDS_Compound compBelow;
        resultBuilder.MakeCompound(compBelow);
        for (int id : belowIds) {
            resultBuilder.Add(compBelow, idToFace.Find(id));
        }
        //std::string belowPath = savePath + "Layer_Z" + std::to_string((int)curSplitZ) + "_BELOW_REMAINING.brep";
        //BRepTools::Write(compBelow, belowPath.c_str());

        //cout << "  [Visualize] Z=" << curSplitZ << " 层提取完成: " << endl;
        //cout << "    - 上方导出: " << abovePath << " (" << aboveIds.size() << " 个面)" << endl;
        //cout << "    - 下方导出: " << belowPath << " (" << belowIds.size() << " 个面)" << endl;
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
// ============================================================
// 入口函数：针对单个独立型腔的自适应切分（智能台阶过滤版）
// ============================================================
std::vector<TopoDS_Compound> SplitCavity(const TopoDS_Compound& cavity) {
    TopTools_DataMapOfIntegerShape idToFace;
    std::set<int> allIds;
    int nextId = 1;

    // 1. 初始化 ID 到面的映射
    TopExp_Explorer exp(cavity, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        int id = nextId++;
        idToFace.Bind(id, TopoDS::Face(exp.Current()));
        allIds.insert(id);
    }

    // ==============================================================
    // 🌟 新增：建立 Edge -> Face 的拓扑邻接图，用于打探邻居面的情报
    // ==============================================================
    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(cavity, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    // 2. 获取并排序切分点
    std::set<double> zSet;
    for (int id : allIds) {
        TopoDS_Face f = TopoDS::Face(idToFace.Find(id));

        // 我们只考察水平平面是否能成为切分点
        if (IsHorizontalFace(f)) {
            double zFace = GetHorizontalFaceZ(f);

            // 默认假设它是真正的底面（不论是封闭底还是开放悬空底），默认不切！
            bool isIntermediateStep = false;

            // 扫荡这个平面的所有边界边，看看邻居是往上长还是往下走
            TopExp_Explorer edgeExp(f, TopAbs_EDGE);
            for (; edgeExp.More(); edgeExp.Next()) {
                const TopoDS_Edge& edge = TopoDS::Edge(edgeExp.Current());

                // 如果这条边有相邻的面
                if (edgeToFaces.Contains(edge)) {
                    const TopTools_ListOfShape& neighbors = edgeToFaces.FindFromKey(edge);
                    TopTools_ListIteratorOfListOfShape it(neighbors);
                    for (; it.More(); it.Next()) {
                        TopoDS_Face nFace = TopoDS::Face(it.Value());

                        // 排除自己，只看真正的邻居面
                        if (!nFace.IsSame(f)) {
                            // 探测邻居面在 Z 轴上的极限范围
                            double nZmin, nZmax;
                            GetFaceZRange(nFace, nZmin, nZmax);

                            // 🎯 核心判据：如果发现任何一个邻居面比我还低（容差 1e-3 防止浮点抖动）
                            // 铁证如山：我下方还有深坑（比如小腔体侧壁）或者外围有更低的侧壁！
                            // 我绝不是底面，我是一个“中间台阶面 (Step)”，必须在这里切一刀！
                            if (nZmin < zFace - 1e-3) {
                                isIntermediateStep = true;
                                break;
                            }
                        }
                    }
                }
                // 只要发现一个向下的邻居，就已经实锤是台阶了，立刻停止排查这层面的其他边，提升效率
                if (isIntermediateStep) break;
            }

            // 判决时刻：只有被判定为“中间台阶”的平面，才允许登记为切分点
            if (isIntermediateStep) {
                zSet.insert(zFace);
            }
            else {
                std::cout << "  [智能拦截] 侦测到水平面 (Z=" << zFace
                    << ") 所有邻面均朝上或为悬空边界，确认为极限大底面，已免于切分！" << std::endl;
            }
        }
    }

    std::vector<double> splitZs(zSet.begin(), zSet.end());
    // 从大到小降序排列（从型腔最顶端往下排）
    std::sort(splitZs.begin(), splitZs.end(), std::greater<double>());

    // 【已被移除】：if (!splitZs.empty()) splitZs.pop_back(); 
    // 原因：上面的智能拦截网已经把真正的闭合底面和开放底面全部保护起来了，
    // zSet 里剩下的【全都是必须切分的中间台阶】，如果硬弹出一个，反而会漏切！

    // 3. 带着极其精准的切分高度列表，执行递归切分
    std::vector<std::set<int>> idGroups; // 存储各组切分后的面 ID
    SplitCavityRecursive(allIds, idToFace, nextId, splitZs, 0, idGroups);

    // 4. 【关键步骤】在 idToFace 销毁前，将打散的 ID 转换回真正的 3D 几何壳体
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

    // 返回由物理切分产生的一组实实在在的加工层级壳体
    return finalShapes;
}
// 将一个包含多个独立区域的 Compound 拆分为多个独立的 Compound
std::vector<TopoDS_Compound> SeparateDisconnectedCavities(const TopoDS_Compound& cavityCompound) {
    std::vector<TopoDS_Compound> individualCavities;

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopExp::MapShapesAndAncestors(cavityCompound, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

    // 额外建立 Vertex -> Face 的映射，处理"只共点不共边"的情况
    TopTools_IndexedDataMapOfShapeListOfShape vertexFaceMap;
    TopExp::MapShapesAndAncestors(cavityCompound, TopAbs_VERTEX, TopAbs_FACE, vertexFaceMap);

    TopTools_MapOfShape visitedFaces;
    TopExp_Explorer exp(cavityCompound, TopAbs_FACE);

    for (; exp.More(); exp.Next()) {
        TopoDS_Face startFace = TopoDS::Face(exp.Current());
        if (visitedFaces.Contains(startFace)) continue;

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

            // 通过共享边查找邻接面
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

            // 通过共享顶点查找邻接面（处理切分后只共点不共边的情况）
            TopExp_Explorer vertExp(currentFace, TopAbs_VERTEX);
            for (; vertExp.More(); vertExp.Next()) {
                const TopoDS_Shape& vertex = vertExp.Current();
                if (vertexFaceMap.Contains(vertex)) {
                    const TopTools_ListOfShape& adjFaces = vertexFaceMap.FindFromKey(vertex);
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

                // 2. ID 匹配：投票法统计外环所有线段的 faceId，取出现频率最高的
                //    比只看第一个线段更稳健，避免被虚拟软边界(faceId==-1)干扰
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

                if (dominantFaceId >= 0 && feat.sourceFaceIds.count(dominantFaceId)) {

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
        }
    }

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



// 用于面与面布尔运算的工业级引擎，完美支持交集、并集、差集，输出全空间剖分的 Face2D 结果
#include "clipper2/clipper.h"

using namespace Clipper2Lib;
// ================== 2. 数据转换工具 ==================

// 放大倍数，保留 6 位小数的精度（满足机械加工需求）
const double SCALE = 1000000.0;

// 将自定义的 Loop 转为 Clipper 的 Path64
Path64 LoopToPath(const Loop& loop) {
    Path64 path;
    for (const auto& line : loop) {
        path.push_back(Point64(line.start.x * SCALE, line.start.y * SCALE));
    }
    return path;
}

// 将 Clipper 的 Path64 还原为 Loop
Loop PathToLoop(const Path64& path, double zHeight) {
    Loop loop;
    if (path.empty()) return loop;

    for (size_t i = 0; i < path.size(); ++i) {
        size_t next_i = (i + 1) % path.size();
        OneLine line;
        line.start.x = static_cast<double>(path[i].x) / SCALE;
        line.start.y = static_cast<double>(path[i].y) / SCALE;
        line.start.z = zHeight;

        line.end.x = static_cast<double>(path[next_i].x) / SCALE;
        line.end.y = static_cast<double>(path[next_i].y) / SCALE;
        line.end.z = zHeight;

        line.faceId = -1; // 布尔运算后的新线段丢失了原始面ID
        loop.push_back(line);
    }
    return loop;
}

// 将 Face2D 转换为 Clipper 的输入路径集合
Paths64 FacesToClipperPaths(const std::vector<Face2D>& faces) {
    Paths64 result;
    for (const auto& face : faces) {
        // 只提取实体面参与运算
        if (face.outerLoop.empty()) continue;

        Path64 outer = LoopToPath(face.outerLoop);
        if (!IsPositive(outer)) std::reverse(outer.begin(), outer.end()); // 外环必须正向
        result.push_back(outer);

        for (const auto& inner : face.innerLoops) {
            Path64 innerPath = LoopToPath(inner);
            if (IsPositive(innerPath)) std::reverse(innerPath.begin(), innerPath.end()); // 内环必须反向
            result.push_back(innerPath);
        }
    }
    return result;
}

// ================== 3. 结果还原：从 PolyTree 重建 Face2D ==================

// 递归遍历 PolyTree，还原出全空间剖分的 Face2D（和之前的手写树逻辑完全一致！）
void ExtractFacesFromPolyNode(const PolyPath64* node, int depth, double zHeight, std::vector<Face2D>& faces) {
    if (depth > 0) {
        // 在布尔运算的结果中，奇数层(depth 1, 3, 5...)代表正向实体区域
        // 偶数层(depth 2, 4, 6...)代表孔洞区域
        // 因为 Face2D 结构已经包含了 innerLoops 来存储孔洞
        // 所以我们只需要在奇数层创建 Face2D，并将其直接子节点作为孔洞加入
        if (depth % 2 != 0) {
            Face2D newFace;
            newFace.type = FaceType::SOLID; // 提取出来的都是正向实体
            newFace.outerLoop = PathToLoop(node->Polygon(), zHeight);

            // 它的直接子节点就是它的孔洞 (depth + 1, 为偶数层)
            for (size_t i = 0; i < node->Count(); ++i) {
                newFace.innerLoops.push_back(PathToLoop(node->Child(i)->Polygon(), zHeight));
            }
            faces.push_back(newFace);
        }
    }

    // 递归处理子节点
    for (size_t i = 0; i < node->Count(); ++i) {
        ExtractFacesFromPolyNode(node->Child(i), depth + 1, zHeight, faces);
    }
}

// ================== 4. 终极封装：面面布尔运算 ==================

/**
 * 执行单个面与单个面的布尔运算（避免多个 clips 路径方向相互抵消）
 * @param subject 目标面 (A)
 * @param clip 裁剪面 (B)
 * @param clipType 运算类型：ClipType::Intersection(交集), Union(并集), Difference(A-B)
 * @param zHeight 输出面的 Z 高度
 */
std::vector<Face2D> BooleanFacesSingle(
    const Face2D& subject,
    const Face2D& clip,
    ClipType clipType,
    double zHeight)
{
    Clipper64 clipper;

    // 1. 添加单个运算主体 A 和单个裁剪面 B
    std::vector<Face2D> subjects = { subject };
    std::vector<Face2D> clips = { clip };
    clipper.AddSubject(FacesToClipperPaths(subjects));
    clipper.AddClip(FacesToClipperPaths(clips));

    // 2. 执行布尔运算，要求输出保留树状结构的 PolyTree
    PolyTree64 solutionTree;
    // 使用 NonZero 填充规则，完美识别由于方向导致的实体与孔洞
    clipper.Execute(clipType, FillRule::NonZero, solutionTree);

    // 3. 从 PolyTree 完美重建所有的 Face2D
    std::vector<Face2D> resultFaces;
    ExtractFacesFromPolyNode(&solutionTree, 0, zHeight, resultFaces);

    return resultFaces;
}

/**
 * 执行面与面的布尔运算（支持多个 subjects 和多个 clips）
 * @param subjects 目标面集合 (A)
 * @param clips 裁剪面集合 (B)
 * @param clipType 运算类型：ClipType::Intersection(交集), Union(并集), Difference(A-B)
 * @param zHeight 输出面的 Z 高度
 */
std::vector<Face2D> BooleanFaces(
    const std::vector<Face2D>& subjects,
    const std::vector<Face2D>& clips,
    ClipType clipType,
    double zHeight)
{
    // 如果 clips 只有一个 Face2D，直接使用原始方法
    if (clips.size() <= 1) {
        Clipper64 clipper;
        clipper.AddSubject(FacesToClipperPaths(subjects));
        clipper.AddClip(FacesToClipperPaths(clips));
        
        PolyTree64 solutionTree;
        clipper.Execute(clipType, FillRule::NonZero, solutionTree);
        
        std::vector<Face2D> resultFaces;
        ExtractFacesFromPolyNode(&solutionTree, 0, zHeight, resultFaces);
        return resultFaces;
    }
    
    // 如果 clips 有多个 Face2D，逐个处理，避免路径方向相互抵消
    std::vector<Face2D> result = subjects;
    for (const auto& clip : clips) {
        std::vector<Face2D> newResult;
        for (const auto& subject : result) {
            std::vector<Face2D> diffResult = BooleanFacesSingle(subject, clip, clipType, zHeight);
            newResult.insert(newResult.end(), diffResult.begin(), diffResult.end());
        }
        result = newResult;
    }
    return result;
}

/**
 * @brief 【特征追溯引擎】空间比对碰撞法：找回被 Clipper2 弄丢的原始面 ID
 * @param pocketResults Clipper2 布尔运算直接吐出来的原始型腔面集合（此时 faceId 均为 -1）
 * @param currentSolidFaces 这一层最原始的工件实体面集合（内部保留了 STEP 模型的真实 faceId）
 * @return 恢复了原始 faceId（硬边界）且将空气边严格标记为 -1（软边界）的全新 Face2D 型腔集合
 */
std::vector<Face2D> RecoverOriginalFaceIdsByGeometry(
    const std::vector<Face2D>& pocketResults,
    const std::vector<Face2D>& currentSolidFaces)
{
    std::vector<Face2D> recoveredCavities = pocketResults; // 拷贝一份准备改写
    const double match_tol = 0.05; // 50微米空间碰撞重合容差

    // 1. 遍历每一个型腔区域
    for (auto& cavityFace : recoveredCavities) {

        // 2. 盘查当前型腔区域的外环线段
        for (auto& cLine : cavityFace.outerLoop) {
            gp_Pnt cStart(cLine.start.x, cLine.start.y, cLine.start.z);
            gp_Pnt cEnd(cLine.end.x, cLine.end.y, cLine.end.z);

            bool isMatched = false;

            // 3. 去最原始的实体面阵营里进行空间高精碰撞
            for (const auto& sFace : currentSolidFaces) {
                for (const auto& sLine : sFace.outerLoop) {
                    gp_Pnt sStart(sLine.start.x, sLine.start.y, sLine.start.z);
                    gp_Pnt sEnd(sLine.end.x, sLine.end.y, sLine.end.z);

                    // 🎯 空间几何对碰成功（支持正向重合或反向首尾重合）
                    if ((cStart.Distance(sStart) < match_tol && cEnd.Distance(sEnd) < match_tol) ||
                        (cStart.Distance(sEnd) < match_tol && cEnd.Distance(sStart) < match_tol))
                    {
                        // 🟢 【黄金继承】：将型腔这根线的 ID，完美恢复成它亲生父母在 STEP 里的原始面 ID！
                        cLine.faceId = sLine.faceId;
                        isMatched = true;
                        break;
                    }
                }
                if (isMatched) break;
            }

            // 4. 🔴 如果遍历了所有的实体边界都碰不上，铁证如山：它就是悬空的【虚拟软边界】
            if (!isMatched) {
                cLine.faceId = -1; // 强制给它盖章为软边界标记
            }
        }

        // 5. 同样的逻辑，顺手盘查可能存在的内孔边界（如果有的话）
        for (auto& innerLoop : cavityFace.innerLoops) {
            for (auto& cLine : innerLoop) {
                gp_Pnt cStart(cLine.start.x, cLine.start.y, cLine.start.z);
                gp_Pnt cEnd(cLine.end.x, cLine.end.y, cLine.end.z);
                bool isMatched = false;

                for (const auto& sFace : currentSolidFaces) {
                    for (const auto& sLine : sFace.outerLoop) {
                        gp_Pnt sStart(sLine.start.x, sLine.start.y, sLine.start.z);
                        gp_Pnt sEnd(sLine.end.x, sLine.end.y, sLine.end.z);

                        if ((cStart.Distance(sStart) < match_tol && cEnd.Distance(sEnd) < match_tol) ||
                            (cStart.Distance(sEnd) < match_tol && cEnd.Distance(sStart) < match_tol))
                        {
                            cLine.faceId = sLine.faceId;
                            isMatched = true;
                            break;
                        }
                    }
                    if (isMatched) break;
                }
                if (!isMatched) {
                    cLine.faceId = -1;
                }
            }
        }
    }

    return recoveredCavities;
}

/**
 * @brief 【调试专用】将型腔集合中所有 faceId == -1 的虚拟软边界单独导出为 BREP 线框
 * @param cavityFaces 已经过 ID 恢复（RecoverOriginalFaceIdsByGeometry）的型腔面集合
 * @param filePath 导出的文件路径（例如 "D:/debug_soft_edges.brep"）
 */
void ExportSoftEdgesToBrep(const std::vector<Face2D>& cavityFaces, const std::string& filePath)
{
    std::vector<OneLine> softLines;

    // 1. 搜刮所有型腔面的外环
    for (const auto& face : cavityFaces) {
        for (const auto& line : face.outerLoop) {
            if (line.faceId == -1) {
                softLines.push_back(line);
            }
        }

        // 2. 顺手搜刮可能存在的内孔（理论上凸包求差的内孔一般贴着实体，但安全起见也扫一遍）
        for (const auto& innerLoop : face.innerLoops) {
            for (const auto& line : innerLoop) {
                if (line.faceId == -1) {
                    softLines.push_back(line);
                }
            }
        }
    }

    if (softLines.empty()) {
        std::cout << "⚠️ [软边界可视化] Z 轴该层未检测到任何 faceId == -1 的软边界线段！" << std::endl;
        return;
    }

    // 3. 借用你原有的 ExportOneLinesToBrep 刚性写入磁盘
    ExportOneLinesToBrep(softLines, filePath);
    std::cout << "🚀 [软边界可视化成功] 共提取出 " << softLines.size() << " 条软边界线段 -> " << filePath << std::endl;
}
//开放型腔清洗流程
/**
 * @brief 【特征量化】已知 Face2D 的外环计算其几何紧实度 (Isoperimetric Quotient)
 * @param face 输入的 2D 面结构
 * @return 紧实度值，范围 (0, 1]。越接近 1 越胖（趋近于圆/正方形），越接近 0 越骨感（趋近于长毛刺）
 */
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

/**
 * @brief 【精密细化版】找准软边界宏观两端，通过硬壁趋势刺探圆角垃圾
 * @param face 已经通过 RecoverOriginalFaceIdsByGeometry 恢复了硬壁 ID、标记了 -1 软边的型腔面
 * @return true 代表是纯外壁圆角夹缝（需清洗）；false 代表是有效型腔（保留）
 */
bool IsPureFilletSliverCavity_Advanced(const Face2D& face)
{
    size_t nSize = face.outerLoop.size();
    if (nSize < 8) return false; // 边数太少直接放行，不构成复杂圆角条件

    // ----------------────────────────────────────────=================
    // 步骤 1：第一轮盘查，找出所有连续为 -1 的“软边界线段链”的宏观起点和终点下标
    // ----------------------------------------------------------------=
    struct SoftChain {
        size_t macroStartIdx = 0; // 整个软边界链的第一根线段在 outerLoop 中的下标
        size_t macroEndIdx = 0;   // 整个软边界链的最后一根线段在 outerLoop 中的下标
    };
    std::vector<SoftChain> chains;

    bool inChain = false;
    SoftChain currentChain;

    for (size_t i = 0; i < nSize; ++i) {
        if (face.outerLoop[i].faceId == -1) {
            if (!inChain) {
                currentChain.macroStartIdx = i;
                inChain = true;
            }
            currentChain.macroEndIdx = i;
        }
        else {
            if (inChain) {
                chains.push_back(currentChain);
                inChain = false;
            }
        }
    }
    // 闭环处理：如果首尾在 -1 上连起来了
    if (inChain) {
        chains.push_back(currentChain);
    }
    // 如果首尾正好跨越了 0 点粘连，做一次安全合并
    if (chains.size() > 1 && chains.front().macroStartIdx == 0 && chains.back().macroEndIdx == nSize - 1) {
        chains.front().macroStartIdx = chains.back().macroStartIdx;
        chains.pop_back();
    }

    // 如果这个圈压根没有软边界（全都是硬墙），说明是封闭内腔，绝对不能删
    if (chains.empty()) return false;

    // ----------------────────────────────────────────=================
    // 步骤 2：针对每一条宏观软边界链，看它的“两端趋势”是否都被圆角包夹
    // ----------------------------------------------------------------=
    for (const auto& chain : chains) {

        // =============================================================
        // 🔎 探测【最初的起点端】：逆着拓扑环，从 macroStartIdx 往前穿透看 3 根硬线
        // =============================================================
        bool startSideIsFillet = false;

        size_t sIdx = chain.macroStartIdx;
        size_t prev1 = (sIdx == 0) ? nSize - 1 : sIdx - 1;
        size_t prev2 = (prev1 == 0) ? nSize - 1 : prev1 - 1;
        size_t prev3 = (prev2 == 0) ? nSize - 1 : prev2 - 1;

        // 确保探测到的不是其他软边界线
        if (face.outerLoop[prev1].faceId >= 0 && face.outerLoop[prev2].faceId >= 0 && face.outerLoop[prev3].faceId >= 0)
        {
            gp_XY pDir1(face.outerLoop[prev1].end.x - face.outerLoop[prev1].start.x, face.outerLoop[prev1].end.y - face.outerLoop[prev1].start.y);
            gp_XY pDir2(face.outerLoop[prev2].end.x - face.outerLoop[prev2].start.x, face.outerLoop[prev2].end.y - face.outerLoop[prev2].start.y);
            gp_XY pDir3(face.outerLoop[prev3].end.x - face.outerLoop[prev3].start.x, face.outerLoop[prev3].end.y - face.outerLoop[prev3].start.y);

            if (pDir1.SquareModulus() > 1e-6 && pDir2.SquareModulus() > 1e-6 && pDir3.SquareModulus() > 1e-6) {
                pDir1.Normalize(); pDir2.Normalize(); pDir3.Normalize();

                double dotPrev1 = pDir1.Dot(pDir2);
                double dotPrev2 = pDir2.Dot(pDir3);
                double totalPrevDot = pDir1.Dot(pDir3); // 宏观大方向改变趋势

                // ⚖️ 趋势判据：每一步都是平缓连接（>165°），但走完3步后，总方向发生了明确偏转（totalDot < 0.985），说明踩中圆角趋势
                if (dotPrev1 > 0.98 && dotPrev2 > 0.98 && (totalPrevDot > 0.88 && totalPrevDot < 0.985)) {
                    startSideIsFillet = true;
                }
            }
        }

        // =============================================================
        // 🔎 探测【最初的终点端】：顺着拓扑环，从 macroEndIdx 往后穿透看 3 根硬线
        // =============================================================
        bool endSideIsFillet = false;

        size_t eIdx = chain.macroEndIdx;
        size_t next1 = (eIdx + 1) % nSize;
        size_t next2 = (next1 + 1) % nSize;
        size_t next3 = (next2 + 1) % nSize;

        if (face.outerLoop[next1].faceId >= 0 && face.outerLoop[next2].faceId >= 0 && face.outerLoop[next3].faceId >= 0)
        {
            gp_XY nDir1(face.outerLoop[next1].end.x - face.outerLoop[next1].start.x, face.outerLoop[next1].end.y - face.outerLoop[next1].start.y);
            gp_XY nDir2(face.outerLoop[next2].end.x - face.outerLoop[next2].start.x, face.outerLoop[next2].end.y - face.outerLoop[next2].start.y);
            gp_XY nDir3(face.outerLoop[next3].end.x - face.outerLoop[next3].start.x, face.outerLoop[next3].end.y - face.outerLoop[next3].start.y);

            if (nDir1.SquareModulus() > 1e-6 && nDir2.SquareModulus() > 1e-6 && nDir3.SquareModulus() > 1e-6) {
                nDir1.Normalize(); nDir2.Normalize(); nDir3.Normalize();

                double dotNext1 = nDir1.Dot(nDir2);
                double dotNext2 = nDir2.Dot(nDir3);
                double totalNextDot = nDir1.Dot(nDir3); // 宏观大方向改变趋势

                if (dotNext1 > 0.98 && dotNext2 > 0.98 && (totalNextDot > 0.88 && totalNextDot < 0.985)) {
                    endSideIsFillet = true;
                }
            }
        }

        // =============================================================
        // ⚖️ 【工艺终审裁决】
        // 如果当前这条完整的宏观软边界链，它的最初起点和最初终点，同时咬死在圆角的转弯趋势里
        // 证实它是由于工件外壁圆角大擦边产生的梭形伪腔，判定为垃圾面，立刻返回 true 下死手！
        // =============================================================
        if (startSideIsFillet && endSideIsFillet) {
            return true;
        }
    }

    return false; // 通过了考核，是个含有直壁的良民型腔
}
/**
 * @brief 【工艺总接口】全量开放型腔多边形综合清洗大管家
 * @param rawOpenCavities Clipper2 计算出来并排除了重叠内腔的原始开口型腔面总集
 * @return 过滤清洗干净、可以直接下刀生成刀轨的真正高质量加工型腔数组
 */
std::vector<Face2D> CleanAndFilterOpenCavities(const std::vector<Face2D>& rawOpenCavities)
{
    std::vector<Face2D> clearFeatures;

    // --- 🔑 工艺大闸门参数定义 ---
    const double MIN_AREA = 2.5;         // 门槛 1：面积小于 2.5 平方毫米的纳米残渣直接蒸发
    const double MIN_COMPACTNESS = 0.015; // 门槛 2：等周商小于 0.015 的极其骨感的超级长尾巴单线直接蒸发
    const double MIN_TOOL_PASS_SPAN = 1.0; // 门槛 3：包围盒跨度小于 1mm 的微观伪胖正方形直接剔除

    for (const auto& face : rawOpenCavities) {
        if (face.outerLoop.size() < 3) continue;

        // 🟢 1. 拦截极小面积碎屑
        double area = CalculateArea(face.outerLoop);
        if (area < MIN_AREA) continue;

        // 🟢 2. 拦截骨感僵尸长尾巴
        double compactness = CalculateFaceCompactness(face);
        if (compactness < MIN_COMPACTNESS) continue;

        // 🟢 3. 拦截微观过小面积盲区（应用快速外接框跨度校验）
        double xMin = 1e9, xMax = -1e9, yMin = 1e9, yMax = -1e9;
        for (const auto& line : face.outerLoop) {
            xMin = std::min(xMin, line.start.x); xMax = std::max(xMax, line.start.x);
            yMin = std::min(yMin, line.start.y); yMax = std::max(yMax, line.start.y);
        }
        if ((xMax - xMin) < MIN_TOOL_PASS_SPAN && (yMax - yMin) < MIN_TOOL_PASS_SPAN) {
            continue;
        }

        // 🟢 4. 🚀 🔥【终极绝杀】：执行你提出来的软边界双端连通圆角穿透检查
        // 如果证明它只是外壁大圆角的顺滑擦边伪区域，直接就地物理毁灭！
        if (IsPureFilletSliverCavity_Advanced(face)) {
            continue;
        }

        // 🏆 恭喜它！闯过四道鬼门关，证明是绝对需要开粗的良民多边形，准予放行！
        clearFeatures.push_back(face);
    }

    std::cout << "🏁 [工艺提纯报告] 原本共有型腔碎片: " << rawOpenCavities.size()
        << " 个 -> 经过四大滤网刚性清算后，最终存活黄金型腔: " << clearFeatures.size() << " 个！" << std::endl;

    return clearFeatures;
}

//开放型腔提取流程
// 1. 【开放专用】提取属于开放型腔的表面 ID 集合
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
    double angleTolerance = 0.001,
    gp_Dir toolDirection = gp_Dir(0, 0, -1)) {
    BRep_Builder builder;
    TopoDS_Compound capFaces;
    builder.MakeCompound(capFaces);

    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    // =========================================================================
    // ⚖️ 步骤 A：宏观清算 —— 算出当前开放侧壁在 3D 空间中的真正 Z 轴极限
    // =========================================================================
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
                GetFaceZRange(f, zmin, zmax); // 借用你原有的 Z 范围函数
                wallMinZ = std::min(wallMinZ, zmin);
                wallMaxZ = std::max(wallMaxZ, zmax);
                hasValidWalls = true;
            }
        }
    }

    // 如果这一层连开放侧壁都丢了，直接退出
    if (!hasValidWalls) return capFaces;

    // =========================================================================
    // ⚖️ 步骤 B：微观盘查 —— 带着法向大闸门筛选真正的工艺底面
    // =========================================================================
    const gp_Dir zAxis(0, 0, 1);
    TopExp_Explorer exp(solid, TopAbs_FACE);

    for (; exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        BRepAdaptor_Surface surf(face);
        if (surf.GetType() != GeomAbs_Plane) continue;

        gp_Pln plane = surf.Plane();
        gp_Dir normal = plane.Axis().Direction();

        // 关键：获取面在中心点处的准确法线方向（处理好 REVERSED 取反）
        if (face.Orientation() == TopAbs_REVERSED) {
            normal.Reverse();
        }

        // 1. 🔍 水平大闸门：法线必须平行于 Z 轴
        double angle = normal.Angle(zAxis);
        if (angle > angleTolerance && fabs(angle - M_PI) > angleTolerance) continue;

        int faceId = -1;
        if (faceToIdMap.IsBound(face)) faceId = faceToIdMap.Find(face);
        else continue;

        // 2. 🔍 顶面判断：用法向与进刀方向的点积关系来决定
        //    dot(normal, toolDirection) > 0 → 法向与进刀同向 → 天花板/悬垂面 → 保留
        //    dot(normal, toolDirection) < 0 → 法向与进刀反向 → 天空开口面 → 丢弃
        double faceZ = plane.Location().Z();
        double dotWithTool = normal.X() * toolDirection.X() + normal.Y() * toolDirection.Y() + normal.Z() * toolDirection.Z();
        if (std::abs(faceZ - wallMaxZ) < 0.1 && dotWithTool < 0) {
            continue;
        }

        // 3. 🔍 底面法向判据
        // 模型的实体法向【必须严格朝上】(normal.Z() > 0.99)！
        // 如果通槽下方悬空导致邻接到了工件大底面，大底面法向朝下(Z = -1)，在这里会被直接过滤掉，完美防错！
        if (normal.Z() < 0.99) {
            continue;
        }

        // 4. 🔍 拓扑比对：检查这个法向朝上的底面是否真的贴着开放侧壁
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
                    // 只要有一侧邻接了我们认可的有效开放型腔侧壁
                    if (openFaceIds.count(adjId)) {
                        hasOpenCavityNeighbor = true;
                        break;
                    }
                }
            }
            if (hasOpenCavityNeighbor) break;
        }

        // 5. 🔍 极限区间校准：底面的高度，必须落在侧壁的合理跨度区间内
        if (hasOpenCavityNeighbor && faceId >= 0) {
            if (faceZ >= wallMinZ - 0.1 && faceZ <= wallMaxZ + 0.1) {
                builder.Add(capFaces, face);
            }
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
std::vector<CavityFeature> ConvertPartsToFeatures(
    const std::vector<TopoDS_Compound>& parts,
    const CavityFeature& parentFeat,
    const TopTools_DataMapOfShapeInteger& faceToIdMap)
{
    std::vector<CavityFeature> subFeatures;

    for (const auto& part : parts) {
        CavityFeature subFeat;
        subFeat.featureId = -1; // 建议在外部循环中赋值
        subFeat.toolDirection = parentFeat.toolDirection;

        // 1. 获取物理边界 & ID 映射
        double subZmin = 1e9, subZmax = -1e9;
        TopExp_Explorer exp(part, TopAbs_FACE);
        for (; exp.More(); exp.Next()) {
            TopoDS_Face f = TopoDS::Face(exp.Current());
            double fzmin, fzmax;
            GetFaceZRange(f, fzmin, fzmax);
            subZmin = std::min(subZmin, fzmin);
            subZmax = std::max(subZmax, fzmax);

            if (faceToIdMap.IsBound(f)) {
                subFeat.sourceFaceIds.insert(faceToIdMap.Find(f));
            }
        }
        subFeat.topZ = subZmax;
        subFeat.bottomZ = subZmin;
        subFeat.totalDepth = std::abs(subFeat.topZ - subFeat.bottomZ);

        // 2. 几何嵌套过滤：从父特征中继承属于本零件的轮廓线段
        Bnd_Box subBox;
        BRepBndLib::Add(part, subBox);
        double bx1, by1, bz1, bx2, by2, bz2;
        subBox.Get(bx1, by1, bz1, bx2, by2, bz2);

        for (const auto& loop : parentFeat.stepLoops) {
            // 高度重合检测
            if (loop.zHeight <= subFeat.topZ + 1e-3 && loop.zHeight >= subFeat.bottomZ - 1e-3) {
                if (!loop.lines.empty()) {
                    Point3D pt = loop.lines.front().start;
                    // XY 包围盒包含检测
                    if (pt.x >= bx1 - 1e-3 && pt.x <= bx2 + 1e-3 &&
                        pt.y >= by1 - 1e-3 && pt.y <= by2 + 1e-3) {
                        subFeat.stepLoops.push_back(loop);
                    }
                }
            }
        }
        subFeatures.push_back(subFeat);
    }
    return subFeatures;
}



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
            for (; faceExp.More(); faceExp.Next()) {
                TopoDS_Face subFace = TopoDS::Face(faceExp.Current());
                double smin, smax;
                GetFaceZRange(subFace, smin, smax);
                double zMid = (smin + smax) / 2.0;

                int newId = nextId++;
                idToFace.Bind(newId, subFace);

                if (zMid > splitZ) aboveIds.insert(newId);
                else belowIds.insert(newId);
            }
        } else {
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

std::vector<TopoDS_Compound> RecursiveSplitOpenCavity(const TopoDS_Compound& cavity) {
    std::vector<TopoDS_Compound> result;

    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(cavity, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    std::set<double> zSet;
    TopExp_Explorer exp(cavity, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face f = TopoDS::Face(exp.Current());
        if (!IsHorizontalFace(f)) continue;

        double zFace = GetHorizontalFaceZ(f);
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
    result.push_back(upper);

    std::vector<TopoDS_Compound> lowerParts = SeparateDisconnectedCavities(lower);
    std::cout << "  [连通性] 下方拆分为 " << lowerParts.size() << " 个独立部分" << std::endl;

    if (lowerParts.size() > 1) {
        for (const auto& part : lowerParts) {
            auto subResult = RecursiveSplitOpenCavity(part);
            result.insert(result.end(), subResult.begin(), subResult.end());
        }
    } else {
        auto subResult = RecursiveSplitOpenCavity(lower);
        result.insert(result.end(), subResult.begin(), subResult.end());
    }

    return result;
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
void ProcessAndSplitOpenCavityFeatures(
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
        return;
    }

    // 2. 得到纯粹的开放特征壳体（包含专属的开放底面判定）
    std::string cavityFacesFileName = savePath + "OpenCavity_WallFaces.brep";
    ExportOpenCavityFaces(mainShape, openFaceIds, faceToIdMap, cavityFacesFileName);

    // 获取对应的开放合并面 Compound
    TopoDS_Compound cavityCompound = GetOpenCavityCompound(mainShape, openFaceIds, faceToIdMap);


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
 
    // 5. 沿 Z 轴自适应多级物理切割
    int totalPartCount = 0;
    std::vector<CavityFeature> finalMachinableOpenFeatures;

    for (size_t i = 0; i < trueOpenCavities.size(); ++i) {

        // 获取当前父型腔沿 Z 轴递推切分后的所有层
        std::vector<TopoDS_Compound> zLayers = RecursiveSplitOpenCavity(trueOpenCavities[i]);

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


            // 🎯 【连通性检测与孤岛修复】
            std::vector<TopoDS_Compound> rawSubParts = SeparateDisconnectedCavities(sewedCompound);
            std::vector<TopoDS_Compound> isolatedSubParts = MergeNestedIslands(rawSubParts);

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
        std::cout << "    - 空间范围: Z_Top = " << feat.topZ << " / Z_Floor = " << feat.bottomZ << std::endl;
        std::cout << "    - 深度: " << feat.totalDepth << " mm" << std::endl;
        std::cout << "    - 包含 2D 层级切片数量: " << feat.stepLoops.size() << " 圈" << std::endl;
    }
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
std::vector<TopoDS_Compound> ProcessAndSplitClosedCavityFeatures(
    const std::vector<std::vector<Face2D>>& allLayerFaces,
    const TopoDS_Shape& mainShape,
    const std::vector<double>& splitPoints,
    const TopTools_DataMapOfShapeInteger& faceToIdMap,
    const std::string& savePath)
{
    std::vector<TopoDS_Compound> trueClosedCavities;

    // 1. 提取所有腔体面 ID 集合
    std::set<int> cavityFaceIds = CollectCavityFaceIds(allLayerFaces);
    if (cavityFaceIds.empty()) {
        return trueClosedCavities;
    }

    // 2. 得到纯粹的“型腔壳体”
    std::string cavityFacesFileName = savePath + "CavityWallFaces.brep";
    ExportCavityFaces(mainShape, cavityFaceIds, faceToIdMap, cavityFacesFileName);

    // 获取型腔合并面 Compound，供后续所有三维切割和特征工程使用
    TopoDS_Compound cavityCompound = GetCavityCompound(mainShape, cavityFaceIds, faceToIdMap);

    // 3. 利用 BFS（广度优先搜索）拓扑连通算法，把“拼在一起”的壳体按物理区域拆开
    std::vector<TopoDS_Compound> isolatedCavities = SeparateDisconnectedCavities(cavityCompound);
    cout << "    [拓扑分析] 成功将合并壳体拆分为 " << isolatedCavities.size() << " 个空间独立的型腔连通区域。" << endl;

    // 4. 生成高级特征表达
    std::vector<CavityFeature> cavityFeatures = GenerateCavityFeatures(isolatedCavities, allLayerFaces, splitPoints, faceToIdMap);
    
    for (const auto& feat : cavityFeatures) {
        // 智能剪裁：利用 topZ 解析切割，剔除超出开粗顶部的侧壁面（还原真正的工艺加工面）
        TopoDS_Compound trueCavity = ExtractTrueCavityFaces(feat, mainShape, faceToIdMap);
		trueClosedCavities.push_back(trueCavity);
        string filePath = savePath + "Final_True_ClosedCavity_" + to_string(feat.featureId) + ".brep";
        BRepTools::Write(trueCavity, filePath.c_str());
    }

    // 5. 沿z轴分割
    int totalPartCount = 0;
    std::vector<CavityFeature> finalMachinableClosedFeatures;
    for (size_t i = 0; i < trueClosedCavities.size(); ++i) {
        cout << "在处理第 [" << i << "] 个独立型腔区域..." << endl;

        // 获取当前父型腔沿 Z 轴递推切分后的所有层
        std::vector<TopoDS_Compound> parts = RecursiveSplitOpenCavity(trueClosedCavities[i]);

        std::vector<CavityFeature> converted = ConvertPartsToFeatures(parts, cavityFeatures[i], faceToIdMap);

        for (size_t k = 0; k < converted.size(); ++k) {
            converted[k].featureId = totalPartCount;

            std::string fileName = savePath + "Final_CAM_ClosedPart_" + std::to_string(totalPartCount) + ".brep";
            BRepTools::Write(parts[k], fileName.c_str());

            finalMachinableClosedFeatures.push_back(converted[k]);

            totalPartCount++;
        }

    }

    for (const auto& feat : finalMachinableClosedFeatures) {
        std::cout << "\n>>> [独立可加工闭合型腔特征 ID: " << feat.featureId << "]" << std::endl;
        std::cout << "    - 空间范围: Z_Top = " << feat.topZ << " / Z_Floor = " << feat.bottomZ << std::endl;
        std::cout << "    - 深度: " << feat.totalDepth << " mm" << std::endl;
        std::cout << "    - 包含 2D 层级切片数量: " << feat.stepLoops.size() << " 圈" << std::endl;
    }

    return trueClosedCavities;
}


#if 1
int main() {
    std::string stepfile = "7_stp_stp.stp";

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

    //// ================= 新增：建立最大毛坯包围体并切分 =================
    //cout << "\n=== 建立最大毛坯包围体并分层切割 ===" << endl;
    //Bnd_Box totalBox;
    //BRepBndLib::Add(mainShape, totalBox);
    //if (!totalBox.IsVoid()) {
    //    double xMin, yMin, zMin, xMax, yMax, zMax;
    //    totalBox.Get(xMin, yMin, zMin, xMax, yMax, zMax);

    //    // 给 X/Y 方向增加一点加工余量 (可选，这里默认为 0)
    //    double offset = 2.0;
    //    xMin -= offset; yMin -= offset;
    //    xMax += offset; yMax += offset;

    //    // 创建长方体毛坯
    //    TopoDS_Shape stockShape = BRepPrimAPI_MakeBox(gp_Pnt(xMin, yMin, zMin), gp_Pnt(xMax, yMax, zMax)).Shape();
    //    std::string StockFileName = savePath + "Stock_Body.brep";
    //    BRepTools::Write(stockShape, StockFileName.c_str());
    //    cout << "  已生成整体毛坯实体: Stock_Body.brep" << endl;

    //    // 对毛坯进行分层切割
    //    for (int i = (int)splitPoints.size() - 1; i >= 0; i--) {
    //        double splitZ = splitPoints[i];
    //        int reversedIndex = splitPoints.size() - 1 - i;
    //        
    //        gp_Pln cuttingPlane(gp_Pnt(0, 0, splitZ), gp_Dir(0, 0, 1));
    //        TopoDS_Face algoPlane = BRepBuilderAPI_MakeFace(cuttingPlane);

    //        BRepAlgoAPI_Section section(stockShape, algoPlane, Standard_True);
    //        section.Build();

    //        std::vector<OneLine> stockLines;
    //        if (section.IsDone()) {
    //            TopExp_Explorer expEdge(section.Shape(), TopAbs_EDGE);
    //            for (; expEdge.More(); expEdge.Next()) {
    //                TopoDS_Edge E = TopoDS::Edge(expEdge.Current());
    //                BRepAdaptor_Curve bac(E);
    //                GCPnts_QuasiUniformDeflection discretizer(bac, 0.01);
    //                if (discretizer.IsDone()) {
    //                    for (int k = 1; k < discretizer.NbPoints(); ++k) {
    //                        gp_Pnt p1 = discretizer.Value(k);
    //                        gp_Pnt p2 = discretizer.Value(k + 1);
    //                        OneLine line;
    //                        line.start.x = p1.X(); line.start.y = p1.Y(); line.start.z = p1.Z();
    //                        line.end.x = p2.X();   line.end.y = p2.Y();   line.end.z = p2.Z();
    //                        line.faceId = -888; 
    //                        stockLines.push_back(line);
    //                    }
    //                }
    //            }
    //        }

    //        std::vector<Face2D> stockFaces = BuildTopologyAndExtractFaces(stockLines);
    //        std::vector<Face2D> stockSolids;
    //        for (const auto& f : stockFaces) {
    //            if (f.type == FaceType::SOLID) stockSolids.push_back(f);
    //        }

    //        std::string stockName = savePath + "Stock_Slice_" + std::to_string(reversedIndex) + ".brep";
    //        ExportFace2DToBrep(stockSolids, stockName);
    //    }
    //    cout << "  毛坯分层切割完成！" << endl;
    //}
    //// =================================================================

    cout << "\n开始逐层切分模型 (从上往下)..." << endl;
    
    
    // 分离存储：封闭型腔和开放型腔的数据使用独立容器，避免混合污染
    std::vector<std::vector<Face2D>> allClosedLayerFaces;  // 封闭型腔专用
    std::vector<std::vector<Face2D>> allOpenLayerFaces;    // 开放型腔专用

    for (int i = (int)splitPoints.size() - 1; i >= 0; i--) {
        double splitZ = splitPoints[i];
        int reversedIndex = splitPoints.size() - 1 - i;
        cout << "\n切分 (倒序 #" << reversedIndex << "/" << splitPoints.size() << ") (Z = " << splitZ << ")" << endl;
        
        if (1)
        {
            // 提取本层的面（包含 SOLID + CAVITY）
            std::vector<Face2D> currentLayerFaces = SliceModelAtZ(mainShape, splitZ, i, splitPoints.size()-1, faceToIdMap);

            // 将封闭型腔相关的切片数据保存到专用容器
            allClosedLayerFaces.push_back(currentLayerFaces);

            // 1. 分离出实体面参与运算
            std::vector<Face2D> currentSolidFaces;
            std::vector<Face2D> currentConvexHullFaces;
            for (const auto& f : currentLayerFaces) {
                if (f.type == FaceType::SOLID) {
                    currentSolidFaces.push_back(f);

                    // 1. 转成Point_2 的格式，准备计算二维凸包
                    vector<Point_2> points_2d = ConvertFaceToPoints2D_EK(f);

                    // 2. 计算 2D 凸包
                    std::vector<Point_2> hull_points;
                    CGAL::convex_hull_2(points_2d.begin(), points_2d.end(), std::back_inserter(hull_points));

                    // 3. 还原为 OneLine 闭合环并存入 outerLoop
                    Face2D HullFace = ConvertPointsToHullFace_EK(hull_points, splitZ);
                    currentConvexHullFaces.push_back(HullFace);
                }
            }

            std::vector<Face2D> currentCavityFaces;
            for (const auto& f : currentLayerFaces) {
                if (f.type == FaceType::CAVITY) {
                    currentCavityFaces.push_back(f);
                }
            }

			// 分离出开口型腔区域
            // 新方案：找出最外层的凸包（不被其他凸包包围的凸包）
            // 只对最外层的凸包进行布尔运算，小凸包的结果自然包含在大凸包中
            std::vector<Face2D> outerHulls;
            for (const auto& hullFace : currentConvexHullFaces) {
                bool isInsideOtherHull = false;
                for (const auto& otherHull : currentConvexHullFaces) {
                    if (&otherHull == &hullFace) continue; // 跳过自身
                    if (IsFaceInsideFace(hullFace, otherHull)) {
                        isInsideOtherHull = true;
                        break;
                    }
                }
                // 只有不被其他凸包包围的才是最外层凸包
                if (!isInsideOtherHull) {
                    outerHulls.push_back(hullFace);
                }
            }
            
            std::vector<Face2D> pocketResults;
            
            for (const auto& hullFace : outerHulls) {
                // 找出在当前凸包内部的所有实体面
                std::vector<Face2D> solidsInsideHull;
                for (const auto& solidFace : currentSolidFaces) {
                    solidsInsideHull.push_back(solidFace);
                }
                
                // 找出在当前凸包内部的所有空腔面
                std::vector<Face2D> cavitiesInsideHull;
                for (const auto& cavityFace : currentCavityFaces) {
                    cavitiesInsideHull.push_back(cavityFace);
                }
                
                // 用当前凸包依次减去它内部的实体面
                std::vector<Face2D> hullResult = { hullFace };
                for (const auto& solidFace : solidsInsideHull) {
                    std::vector<Face2D> newResult;
                    for (const auto& resultFace : hullResult) {
                        std::vector<Face2D> diffResult = BooleanFacesSingle(resultFace, solidFace, ClipType::Difference, splitZ);
                        newResult.insert(newResult.end(), diffResult.begin(), diffResult.end());
                    }
                    hullResult = newResult;
                }
                
                // 继续减去它内部的空腔面
                for (const auto& cavityFace : cavitiesInsideHull) {
                    std::vector<Face2D> newResult;
                    for (const auto& resultFace : hullResult) {
                        std::vector<Face2D> diffResult = BooleanFacesSingle(resultFace, cavityFace, ClipType::Difference, splitZ);
                        newResult.insert(newResult.end(), diffResult.begin(), diffResult.end());
                    }
                    hullResult = newResult;
                }
                
                // 合并当前凸包的处理结果
                pocketResults.insert(pocketResults.end(), hullResult.begin(), hullResult.end());
            }
            
            std::vector<Face2D> tempOpenCavityFaces = RecoverOriginalFaceIdsByGeometry(pocketResults, currentSolidFaces);
            for (auto& f : tempOpenCavityFaces) {
                f.type = FaceType::OPENCAVITY;
            }

			//到处软边界线段，检查 ID 恢复和软边界标记是否正确（faceId == -1 的线段应该就是软边界）
            std::string softEdgesName = savePath + "Slice_" + std::to_string(i) + "_ONLY_SoftEdges_Z" + std::to_string(splitZ) + ".brep";
            ExportSoftEdgesToBrep(tempOpenCavityFaces, softEdgesName);

            vector<Face2D> currentOpenCavityFaces = CleanAndFilterOpenCavities(tempOpenCavityFaces);

            //debug 
            for(auto& f : currentOpenCavityFaces) {
				vector<Face2D> singleFaceVec = { f };
				std::string singleFaceName = savePath + "Slice_" + "SingleOpenCavityFace_ID" + ".brep";
				ExportFace2DToBrep(singleFaceVec, singleFaceName);
                int cc = 0;
            }

            // 将开放型腔数据保存到专用容器
            allOpenLayerFaces.push_back(currentOpenCavityFaces);



			//DEBUG: 导出当前层的实体面和腔面，检查切分结果
            std::string currentLayerFacesName = savePath + "Slice_" + std::to_string(i) + "currentLayerFaces.brep";
            ExportFace2DToBrep(currentLayerFaces, currentLayerFacesName);
            std::string currentHullFaceName = savePath + "Slice_" + std::to_string(i) + "currentHullFace.brep";
            ExportFace2DToBrep(currentConvexHullFaces, currentHullFaceName);
            std::string currentSolidFaceName = savePath + "Slice_" + std::to_string(i) + "currentSolidFace.brep";
            ExportFace2DToBrep(currentSolidFaces, currentSolidFaceName);
            std::string currentCavityFaceName = savePath + "Slice_" + std::to_string(i) + "currentCavityFacee.brep";
            ExportFace2DToBrep(currentCavityFaces, currentCavityFaceName);
            std::string tempOpenCavityFaceName = savePath + "Slice_" + std::to_string(i) + "temp_OpenCavityFaces_Z" + std::to_string(splitZ) + ".brep";
            ExportFace2DToBrep(tempOpenCavityFaces, tempOpenCavityFaceName);
            std::string currentOpenCavityFaceName = savePath + "Slice_" + std::to_string(i) + "CGAL_OpenCavityFaces_Z" + std::to_string(splitZ) + ".brep";
            ExportFace2DToBrep(currentOpenCavityFaces, currentOpenCavityFaceName);
            int aaa = 0;
        }

    }
    // ================= 提取封闭型腔特征面 =================
    auto trueClosedCavities = ProcessAndSplitClosedCavityFeatures(allClosedLayerFaces, mainShape, splitPoints, faceToIdMap, savePath);

    // ================= 提取开放型腔特征面 =================
    ProcessAndSplitOpenCavityFeatures(allOpenLayerFaces, allClosedLayerFaces, trueClosedCavities, mainShape, splitPoints, faceToIdMap, savePath);
    
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