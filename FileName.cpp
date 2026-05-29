#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <cmath>
#include <limits>

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

enum class CavityType {
    OPEN,    // 开放型腔
    CLOSED,  // 封闭型腔
    OTHER    // 其他
};

struct CavityFeature {
    int featureId;
    CavityType type = CavityType::CLOSED;
    std::vector<CavityLoop> stepLoops;
    double topZ = -1e9;
    double bottomZ = 1e9;
    double totalDepth = 0.0;
    gp_Dir toolDirection = gp_Dir(0, 0, -1);
    std::set<int> sourceFaceIds; // 关联的原始面 ID 集合
};

//开放型腔筛选条件
struct OpenCavityFilterParams {
    double minArea          = 2.5; //面积
    double minCompactness   = 0.015; //等周商，删去细长区域
    double minToolPassSpan  = 1.0; //包围盒跨度，删去扁平的区域
    double maxFilletRadius  = 0.5; //圆角半径
    double hullMergeDistance = 360.0; //独立凸包链式合并距离阈值
    double bboxPrecheckMargin = 0.0; //凸包距离粗筛额外余量
};

struct OpenCavitySplitParams {
    double minStepFaceArea = 5.0; //开放型腔切分台阶面的最小面积
    double zProtectionTol = 0.05; //保护顶面/底面，避免贴边切碎
    double neighborZTol = 1e-3; //判断邻接面上下延伸的 Z 容差
    double minSplitGap = 0.2; //过近切分高度合并
    double geometricConnectTol = 0.05; //开放型腔拓扑断开后的几何连通合并容差
    int maxRecursionDepth = 32; //防止异常模型无限递推
    bool splitDownOnlyStepFaces = true; //凸台顶面/岛顶面只有下方邻接时也作为切分面
    bool exportDebugBreps = true;
};

struct HullItem {
    Face2D hull;
    Face2D solid;
    int index = -1;
};

struct HullGroup {
    std::vector<int> memberIndices;
    double triggerDistance = 0.0;
};

struct LoopBBox2D {
    double xMin = 0.0;
    double yMin = 0.0;
    double xMax = 0.0;
    double yMax = 0.0;
    bool valid = false;
};

struct UnionFind {
    std::vector<int> parent;
    std::vector<int> rank;
    std::vector<double> triggerDistance;

    UnionFind(int n = 0) {
        Reset(n);
    }

    void Reset(int n) {
        parent.resize(n);
        rank.assign(n, 0);
        triggerDistance.assign(n, std::numeric_limits<double>::max());
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    int Find(int x) {
        if (parent[x] != x) {
            parent[x] = Find(parent[x]);
        }
        return parent[x];
    }

    void Unite(int a, int b, double distance) {
        int rootA = Find(a);
        int rootB = Find(b);
        if (rootA == rootB) {
            triggerDistance[rootA] = std::min(triggerDistance[rootA], distance);
            return;
        }
        if (rank[rootA] < rank[rootB]) {
            std::swap(rootA, rootB);
        }
        parent[rootB] = rootA;
        triggerDistance[rootA] = std::min(triggerDistance[rootA], triggerDistance[rootB]);
        triggerDistance[rootA] = std::min(triggerDistance[rootA], distance);
        if (rank[rootA] == rank[rootB]) {
            rank[rootA]++;
        }
    }
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

void GetShapeZRange(const TopoDS_Shape& shape, double& zmin, double& zmax) {
    Bnd_Box box;
    BRepBndLib::Add(shape, box);
    if (!box.IsVoid()) {
        double xmin, ymin, xmax, ymax;
        box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    } else {
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

void PrintEndpointGapDiagnostics(const std::vector<OneLine>& lines, double snapTol) {
    struct EndpointInfo {
        Point3D p;
        int lineIndex = -1;
        bool isStart = true;
        int degree = 0;
    };

    std::vector<EndpointInfo> endpoints;
    endpoints.reserve(lines.size() * 2);
    for (int i = 0; i < (int)lines.size(); ++i) {
        endpoints.push_back({ lines[i].start, i, true, 0 });
        endpoints.push_back({ lines[i].end, i, false, 0 });
    }

    for (size_t i = 0; i < endpoints.size(); ++i) {
        for (size_t j = i + 1; j < endpoints.size(); ++j) {
            double dx = endpoints[i].p.x - endpoints[j].p.x;
            double dy = endpoints[i].p.y - endpoints[j].p.y;
            double dist = std::sqrt(dx * dx + dy * dy);
            if (dist < snapTol) {
                endpoints[i].degree++;
                endpoints[j].degree++;
            }
        }
    }

    int openEndpointCount = 0;
    double maxNearestGap = 0.0;
    std::cout << "  [BuildLoops诊断] 未用线段端点断口:" << std::endl;
    for (size_t i = 0; i < endpoints.size(); ++i) {
        if (endpoints[i].degree > 0) continue;
        openEndpointCount++;

        double bestDist = std::numeric_limits<double>::max();
        int bestIndex = -1;
        for (size_t j = 0; j < endpoints.size(); ++j) {
            if (i == j) continue;
            double dx = endpoints[i].p.x - endpoints[j].p.x;
            double dy = endpoints[i].p.y - endpoints[j].p.y;
            double dist = std::sqrt(dx * dx + dy * dy);
            if (dist < bestDist) {
                bestDist = dist;
                bestIndex = (int)j;
            }
        }

        maxNearestGap = std::max(maxNearestGap, bestDist);
        std::cout << "    line=" << endpoints[i].lineIndex
            << (endpoints[i].isStart ? ".start" : ".end")
            << " p=(" << endpoints[i].p.x << "," << endpoints[i].p.y << "," << endpoints[i].p.z << ")"
            << " nearest=" << bestDist;
        if (bestIndex >= 0) {
            std::cout << " -> line=" << endpoints[bestIndex].lineIndex
                << (endpoints[bestIndex].isStart ? ".start" : ".end")
                << " p=(" << endpoints[bestIndex].p.x << "," << endpoints[bestIndex].p.y << "," << endpoints[bestIndex].p.z << ")";
        }
        std::cout << std::endl;
    }

    std::cout << "  [BuildLoops诊断] openEndpointCount=" << openEndpointCount
        << " maxNearestGap=" << maxNearestGap
        << " snapTol=" << snapTol << std::endl;
}

int HealOpenEndpointGaps(std::vector<OneLine>& lines, double healTol) {
    struct EndpointRef {
        int lineIndex = -1;
        bool isStart = true;
        Point3D p;
        int degree = 0;
    };

    auto Distance2DLocal = [](const Point3D& a, const Point3D& b) {
        double dx = a.x - b.x;
        double dy = a.y - b.y;
        return std::sqrt(dx * dx + dy * dy);
    };

    std::vector<EndpointRef> endpoints;
    endpoints.reserve(lines.size() * 2);
    for (int i = 0; i < (int)lines.size(); ++i) {
        endpoints.push_back({ i, true, lines[i].start, 0 });
        endpoints.push_back({ i, false, lines[i].end, 0 });
    }

    for (size_t i = 0; i < endpoints.size(); ++i) {
        for (size_t j = i + 1; j < endpoints.size(); ++j) {
            if (Distance2DLocal(endpoints[i].p, endpoints[j].p) < 1e-3) {
                endpoints[i].degree++;
                endpoints[j].degree++;
            }
        }
    }

    int healedCount = 0;
    std::vector<bool> endpointHealed(endpoints.size(), false);
    while (true) {
        int bestA = -1;
        int bestB = -1;
        double bestDist = healTol;

        for (int a = 0; a < (int)endpoints.size(); ++a) {
            if (endpointHealed[a] || endpoints[a].degree > 0) continue;
            for (int b = a + 1; b < (int)endpoints.size(); ++b) {
                if (endpointHealed[b] || endpoints[b].degree > 0) continue;
                if (endpoints[a].lineIndex == endpoints[b].lineIndex) continue;

                double dist = Distance2DLocal(endpoints[a].p, endpoints[b].p);
                if (dist <= bestDist) {
                    bestDist = dist;
                    bestA = a;
                    bestB = b;
                }
            }
        }

        if (bestA < 0 || bestB < 0) break;

        Point3D mergedPoint;
        mergedPoint.x = (endpoints[bestA].p.x + endpoints[bestB].p.x) * 0.5;
        mergedPoint.y = (endpoints[bestA].p.y + endpoints[bestB].p.y) * 0.5;
        mergedPoint.z = (endpoints[bestA].p.z + endpoints[bestB].p.z) * 0.5;

        auto ApplyEndpoint = [&](const EndpointRef& endpoint) {
            if (endpoint.isStart) {
                lines[endpoint.lineIndex].start = mergedPoint;
            }
            else {
                lines[endpoint.lineIndex].end = mergedPoint;
            }
        };

        ApplyEndpoint(endpoints[bestA]);
        ApplyEndpoint(endpoints[bestB]);
        endpointHealed[bestA] = true;
        endpointHealed[bestB] = true;
        healedCount++;
    }

    return healedCount;
}

// 核心成环算法：只找大环，自动忽略碎线
vector<Loop> BuildLoops(const vector<OneLine>& inputLines) {
    vector<Loop> loops;
    vector<OneLine> lines = DeduplicateLines(inputLines);
    int healedCount = HealOpenEndpointGaps(lines, 0.02);
    if (healedCount > 0) {
        cout << "  [BuildLoops] 端点断口愈合数量: " << healedCount
            << " (healTol=0.02mm)" << endl;
    }
    //debug 
    //std::string Name = savePath + "Slice_debug_line.brep";
    //ExportOneLinesToBrep(lines, Name);

    if (lines.empty()) return loops;

    vector<bool> visited(lines.size(), false);
    const double snapTol = 0.001;

    auto Snap = [](const Point3D& p, double tol) -> std::pair<int64_t, int64_t> {
        return { (int64_t)round(p.x / tol), (int64_t)round(p.y / tol) };
    };

    auto MakeSnapKey = [](int64_t x, int64_t y) -> int64_t {
        return x * 1000000007LL + y;
    };

    std::unordered_map<int64_t, std::vector<size_t>> endMap;
    for (size_t i = 0; i < lines.size(); ++i) {
        auto sk = Snap(lines[i].start, snapTol);
        auto ek = Snap(lines[i].end, snapTol);
        int64_t startKey = MakeSnapKey(sk.first, sk.second);
        int64_t endKey = MakeSnapKey(ek.first, ek.second);
        endMap[startKey].push_back(i);
        if (startKey != endKey) {
            endMap[endKey].push_back(i);
        }
    }

    for (size_t i = 0; i < lines.size(); ++i) {
        if (visited[i]) continue;

        Loop currentLoop;
        std::vector<size_t> trialIndices;
        std::vector<bool> trialUsed(lines.size(), false);

        currentLoop.push_back(lines[i]);
        trialIndices.push_back(i);
        trialUsed[i] = true;
        Point3D startPt = lines[i].start;
        Point3D currEnd = lines[i].end;

        while (true) {
            auto key = Snap(currEnd, snapTol);

            size_t bestIdx = SIZE_MAX;
            double bestDist = 1e9;
            bool bestReverse = false;

            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    int64_t mapKey = MakeSnapKey(key.first + dx, key.second + dy);
                    auto it = endMap.find(mapKey);
                    if (it == endMap.end()) continue;

                    for (size_t idx : it->second) {
                        if (visited[idx] || trialUsed[idx]) continue;
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
            }

            if (bestIdx == SIZE_MAX) break;

            trialUsed[bestIdx] = true;
            trialIndices.push_back(bestIdx);
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
            for (size_t idx : trialIndices) {
                visited[idx] = true;
            }
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
    if (usedCount < (int)lines.size()) {
        std::vector<OneLine> unusedLines;
        for (size_t i = 0; i < lines.size(); ++i) {
            if (!visited[i]) {
                unusedLines.push_back(lines[i]);
            }
        }
        std::string unusedName = savePath + "BuildLoops_UnusedLines_" + std::to_string((int)lines.size())
            + "_unused_" + std::to_string((int)unusedLines.size()) + ".brep";
        ExportOneLinesToBrep(unusedLines, unusedName);
        cout << "  [BuildLoops] 未用线段已导出: " << unusedName << endl;
        PrintEndpointGapDiagnostics(unusedLines, snapTol);
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

LoopBBox2D GetLoopBBox(const Loop& loop) {
    LoopBBox2D box;
    if (loop.empty()) return box;

    box.xMin = box.xMax = loop.front().start.x;
    box.yMin = box.yMax = loop.front().start.y;
    box.valid = true;

    for (const auto& line : loop) {
        box.xMin = std::min(box.xMin, line.start.x);
        box.xMax = std::max(box.xMax, line.start.x);
        box.yMin = std::min(box.yMin, line.start.y);
        box.yMax = std::max(box.yMax, line.start.y);

        box.xMin = std::min(box.xMin, line.end.x);
        box.xMax = std::max(box.xMax, line.end.x);
        box.yMin = std::min(box.yMin, line.end.y);
        box.yMax = std::max(box.yMax, line.end.y);
    }

    return box;
}

double BBoxDistance2D(const LoopBBox2D& a, const LoopBBox2D& b) {
    if (!a.valid || !b.valid) return std::numeric_limits<double>::max();

    double dx = 0.0;
    if (a.xMax < b.xMin) dx = b.xMin - a.xMax;
    else if (b.xMax < a.xMin) dx = a.xMin - b.xMax;

    double dy = 0.0;
    if (a.yMax < b.yMin) dy = b.yMin - a.yMax;
    else if (b.yMax < a.yMin) dy = a.yMin - b.yMax;

    return std::sqrt(dx * dx + dy * dy);
}

double Cross2D(const Point3D& a, const Point3D& b, const Point3D& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool IsValueBetween(double value, double a, double b, double tol = 1e-9) {
    return value >= std::min(a, b) - tol && value <= std::max(a, b) + tol;
}

bool IsPointOnSegment2D(const Point3D& p, const Point3D& a, const Point3D& b, double tol = 1e-9) {
    return std::abs(Cross2D(a, b, p)) <= tol
        && IsValueBetween(p.x, a.x, b.x, tol)
        && IsValueBetween(p.y, a.y, b.y, tol);
}

bool SegmentsIntersect2D(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) {
    const double tol = 1e-9;
    double c1 = Cross2D(a, b, c);
    double c2 = Cross2D(a, b, d);
    double c3 = Cross2D(c, d, a);
    double c4 = Cross2D(c, d, b);

    if (((c1 > tol && c2 < -tol) || (c1 < -tol && c2 > tol)) &&
        ((c3 > tol && c4 < -tol) || (c3 < -tol && c4 > tol))) {
        return true;
    }

    return IsPointOnSegment2D(c, a, b, tol)
        || IsPointOnSegment2D(d, a, b, tol)
        || IsPointOnSegment2D(a, c, d, tol)
        || IsPointOnSegment2D(b, c, d, tol);
}

double PointSegmentDistance2D(const Point3D& p, const Point3D& a, const Point3D& b) {
    double vx = b.x - a.x;
    double vy = b.y - a.y;
    double wx = p.x - a.x;
    double wy = p.y - a.y;
    double lenSq = vx * vx + vy * vy;

    if (lenSq < 1e-18) {
        double dx = p.x - a.x;
        double dy = p.y - a.y;
        return std::sqrt(dx * dx + dy * dy);
    }

    double t = (wx * vx + wy * vy) / lenSq;
    t = std::max(0.0, std::min(1.0, t));
    double projX = a.x + t * vx;
    double projY = a.y + t * vy;
    double dx = p.x - projX;
    double dy = p.y - projY;
    return std::sqrt(dx * dx + dy * dy);
}

double SegmentSegmentDistance2D(const OneLine& a, const OneLine& b) {
    if (SegmentsIntersect2D(a.start, a.end, b.start, b.end)) {
        return 0.0;
    }

    double d1 = PointSegmentDistance2D(a.start, b.start, b.end);
    double d2 = PointSegmentDistance2D(a.end, b.start, b.end);
    double d3 = PointSegmentDistance2D(b.start, a.start, a.end);
    double d4 = PointSegmentDistance2D(b.end, a.start, a.end);
    return std::min(std::min(d1, d2), std::min(d3, d4));
}

double LoopDistance2D(const Loop& a, const Loop& b) {
    if (a.empty() || b.empty()) return std::numeric_limits<double>::max();

    double best = std::numeric_limits<double>::max();
    for (const auto& lineA : a) {
        for (const auto& lineB : b) {
            best = std::min(best, SegmentSegmentDistance2D(lineA, lineB));
            if (best <= 1e-9) return 0.0;
        }
    }
    return best;
}

std::vector<Point_2> CollectFaceOuterPoints(const std::vector<Face2D>& faces) {
    std::vector<Point_2> points;
    for (const auto& face : faces) {
        for (const auto& line : face.outerLoop) {
            points.emplace_back(line.start.x, line.start.y);
            points.emplace_back(line.end.x, line.end.y);
        }
    }
    return points;
}

Face2D ComputeMergedHullFace(
    const std::vector<HullItem>& hullItems,
    const std::vector<int>& memberIndices,
    double zHeight)
{
    std::vector<Face2D> memberSolids;
    for (int memberIndex : memberIndices) {
        if (memberIndex >= 0 && memberIndex < (int)hullItems.size()) {
            memberSolids.push_back(hullItems[memberIndex].solid);
        }
    }

    std::vector<Point_2> points_2d = CollectFaceOuterPoints(memberSolids);
    std::vector<Point_2> hull_points;
    if (!points_2d.empty()) {
        CGAL::convex_hull_2(points_2d.begin(), points_2d.end(), std::back_inserter(hull_points));
    }

    return ConvertPointsToHullFace_EK(hull_points, zHeight);
}

std::vector<HullGroup> BuildMergedHullGroups(
    const std::vector<HullItem>& hullItems,
    const OpenCavityFilterParams& params)
{
    std::vector<HullGroup> groups;
    const int n = (int)hullItems.size();
    if (n == 0) return groups;

    UnionFind uf(n);
    std::vector<LoopBBox2D> boxes(n);
    for (int i = 0; i < n; ++i) {
        boxes[i] = GetLoopBBox(hullItems[i].hull.outerLoop);
    }

    const double mergeLimit = params.hullMergeDistance;
    const double precheckLimit = mergeLimit + params.bboxPrecheckMargin;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double bboxDistance = BBoxDistance2D(boxes[i], boxes[j]);
            if (bboxDistance > precheckLimit) {
                continue;
            }

            double hullDistance = std::numeric_limits<double>::max();
            if (IsFaceInsideFace(hullItems[i].hull, hullItems[j].hull) ||
                IsFaceInsideFace(hullItems[j].hull, hullItems[i].hull)) {
                hullDistance = 0.0;
            }
            else {
                hullDistance = LoopDistance2D(hullItems[i].hull.outerLoop, hullItems[j].hull.outerLoop);
            }

            if (hullDistance <= mergeLimit) {
                uf.Unite(i, j, hullDistance);
            }
        }
    }

    std::map<int, HullGroup> groupedByRoot;
    for (int i = 0; i < n; ++i) {
        int root = uf.Find(i);
        groupedByRoot[root].memberIndices.push_back(i);
    }

    for (auto& entry : groupedByRoot) {
        int root = uf.Find(entry.first);
        entry.second.triggerDistance = uf.triggerDistance[root];
        if (entry.second.triggerDistance == std::numeric_limits<double>::max()) {
            entry.second.triggerDistance = 0.0;
        }
        groups.push_back(entry.second);
    }

    return groups;
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
        feat.type = CavityType::CLOSED;
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
    std::vector<std::pair<int, double>>* matchedFillets = nullptr)
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
    const OpenCavityFilterParams& params = OpenCavityFilterParams{})
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
    double angleTolerance = 0.001) {
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
        subFeat.type = parentFeat.type;
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

int CountFacesInCompound(const TopoDS_Compound& compound);

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
void ProcessAndSplitClosedCavityFeatures(
    const std::vector<std::vector<Face2D>>& allLayerFaces,
    const TopoDS_Shape& mainShape,
    const std::vector<double>& splitPoints,
    const TopTools_DataMapOfShapeInteger& faceToIdMap,
    const std::string& savePath,
    std::vector<TopoDS_Compound>& outTrueClosedCavities,
    std::vector<CavityFeature>& outFeatures)
{
    std::vector<TopoDS_Compound> trueClosedCavities;

    // 1. 提取所有腔体面 ID 集合
    std::set<int> cavityFaceIds = CollectCavityFaceIds(allLayerFaces);
    if (cavityFaceIds.empty()) {
        return;
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
        std::vector<TopoDS_Compound> parts = RecursiveSplitCavity(trueClosedCavities[i]);

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

    outTrueClosedCavities = trueClosedCavities;
    outFeatures = finalMachinableClosedFeatures;
}

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
    double zTol = 0.5)
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
    double minSharedRatio = 0.8)
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
    double zTol = 0.1)
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
    double zTol = 0.1)
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
    double zTol = 0.1)
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

Face2D TranslateFace2DZ(const Face2D& face, double newZ) {
    Face2D result = face;
    for (auto& line : result.outerLoop) {
        line.start.z = newZ;
        line.end.z = newZ;
    }
    for (auto& inner : result.innerLoops) {
        for (auto& line : inner) {
            line.start.z = newZ;
            line.end.z = newZ;
        }
    }
    return result;
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

    std::map<int, double> filletMap = BuildFilletRadiusMap(mainShape, faceToIdMap);
    cout << "已建立圆角半径映射，共 " << filletMap.size() << " 个面。" << endl;

    OpenCavityFilterParams openCavityParams;

    std::vector<double> splitPoints = GetSplitPointsAlongZ(mainShape);
    
    cout << "检测到 " << splitPoints.size() << " 个切分点:" << endl;
    for (int i = 0; i < splitPoints.size(); i++) {
        cout << "  [" << i << "] Z = " << splitPoints[i] << endl;
    }

    cout << "\n开始逐层切分模型 (从上往下)..." << endl;
    
    
    // 分离存储：封闭型腔和开放型腔的数据使用独立容器，避免混合污染
    std::vector<std::vector<Face2D>> allClosedLayerFaces;  // 封闭型腔专用
    std::vector<std::vector<Face2D>> allOpenLayerFaces;    // 开放型腔专用

    for (int i = (int)splitPoints.size() - 1; i >= 0; i--) {
        double splitZ = splitPoints[i];
        int reversedIndex = splitPoints.size() - 1 - i;
        cout << "\n切分 (倒序 #" << reversedIndex << "/" << splitPoints.size() << ") (Z = " << splitZ << ")" << endl;
        
        if ( 2 )
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

            std::vector<HullItem> hullItems;
            for (size_t h = 0; h < currentConvexHullFaces.size() && h < currentSolidFaces.size(); ++h) {
                HullItem item;
                item.hull = currentConvexHullFaces[h];
                item.solid = currentSolidFaces[h];
                item.index = (int)h;
                hullItems.push_back(item);
            }

            std::vector<HullGroup> mergedHullGroups = BuildMergedHullGroups(hullItems, openCavityParams);
            std::vector<Face2D> currentMergedHullFaces;

            cout << "  [开放型腔凸包合并] 原始凸包数量: " << currentConvexHullFaces.size()
                << " -> 合并后大凸包数量: " << mergedHullGroups.size()
                << " (阈值=" << openCavityParams.hullMergeDistance << "mm)" << endl;

            std::vector<Face2D> pocketResults;

            for (size_t groupIdx = 0; groupIdx < mergedHullGroups.size(); ++groupIdx) {
                const HullGroup& group = mergedHullGroups[groupIdx];
                Face2D mergedHullFace = ComputeMergedHullFace(hullItems, group.memberIndices, splitZ);
                if (mergedHullFace.outerLoop.empty()) {
                    continue;
                }
                currentMergedHullFaces.push_back(mergedHullFace);

                cout << "    - 大凸包组[" << groupIdx << "] 成员数: " << group.memberIndices.size()
                    << " 触发距离: " << group.triggerDistance << endl;

                std::vector<Face2D> solidsInsideHull;
                for (int memberIndex : group.memberIndices) {
                    if (memberIndex >= 0 && memberIndex < (int)hullItems.size()) {
                        solidsInsideHull.push_back(hullItems[memberIndex].solid);
                    }
                }

                std::vector<Face2D> cavitiesInsideHull;
                for (const auto& cavityFace : currentCavityFaces) {
                    bool shouldSubtract = IsFaceInsideFace(cavityFace, mergedHullFace);
                    if (!shouldSubtract) {
                        std::vector<Face2D> intersection =
                            BooleanFacesSingle(cavityFace, mergedHullFace, ClipType::Intersection, splitZ);
                        shouldSubtract = !intersection.empty();
                    }
                    if (shouldSubtract) {
                        cavitiesInsideHull.push_back(cavityFace);
                    }
                }

                std::vector<Face2D> hullResult = { mergedHullFace };
                for (const auto& solidFace : solidsInsideHull) {
                    std::vector<Face2D> newResult;
                    for (const auto& resultFace : hullResult) {
                        std::vector<Face2D> diffResult = BooleanFacesSingle(resultFace, solidFace, ClipType::Difference, splitZ);
                        newResult.insert(newResult.end(), diffResult.begin(), diffResult.end());
                    }
                    hullResult = newResult;
                }

                for (const auto& cavityFace : cavitiesInsideHull) {
                    std::vector<Face2D> newResult;
                    for (const auto& resultFace : hullResult) {
                        std::vector<Face2D> diffResult = BooleanFacesSingle(resultFace, cavityFace, ClipType::Difference, splitZ);
                        newResult.insert(newResult.end(), diffResult.begin(), diffResult.end());
                    }
                    hullResult = newResult;
                }

                pocketResults.insert(pocketResults.end(), hullResult.begin(), hullResult.end());
            }
            
            std::vector<Face2D> tempOpenCavityFaces = RecoverOriginalFaceIdsByGeometry(pocketResults, currentSolidFaces);
            for (auto& f : tempOpenCavityFaces) {
                f.type = FaceType::OPENCAVITY;
            }

			//到处软边界线段，检查 ID 恢复和软边界标记是否正确（faceId == -1 的线段应该就是软边界）
            std::string softEdgesName = savePath + "Slice_" + std::to_string(i) + "_ONLY_SoftEdges_Z" + std::to_string(splitZ) + ".brep";
            ExportSoftEdgesToBrep(tempOpenCavityFaces, softEdgesName);

            vector<Face2D> currentOpenCavityFaces = CleanAndFilterOpenCavities(tempOpenCavityFaces, filletMap, openCavityParams);

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
            std::string currentSmallHullFaceName = savePath + "Slice_" + std::to_string(i) + "_IndependentSmallHulls_Z" + std::to_string(splitZ) + ".brep";
            ExportFace2DToBrep(currentConvexHullFaces, currentSmallHullFaceName);
            std::cout << "  [DEBUG] 本层独立小凸包已保存: " << currentSmallHullFaceName
                << " 数量=" << currentConvexHullFaces.size() << std::endl;

            std::string currentMergedHullFaceName = savePath + "Slice_" + std::to_string(i) + "_MergedBigHulls_Z" + std::to_string(splitZ) + ".brep";
            ExportFace2DToBrep(currentMergedHullFaces, currentMergedHullFaceName);
            std::cout << "  [DEBUG] 本层合并大凸包已保存: " << currentMergedHullFaceName
                << " 数量=" << currentMergedHullFaces.size() << std::endl;
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
    double modelMinZ = 0.0, modelMaxZ = 0.0;
    GetShapeZRange(mainShape, modelMinZ, modelMaxZ);
    bool isInteriorOpen = HasInteriorOpenCavity(allOpenLayerFaces, allClosedLayerFaces, modelMaxZ, splitPoints, 0.5);
    std::vector<TopoDS_Compound> trueClosedCavities;
    std::vector<CavityFeature> openCavityFeatures;
    std::vector<CavityFeature> closedCavityFeatures;
    std::vector<TopoDS_Compound> interiorNewClosedParts;
    TopoDS_Compound interiorNewOpenCavity;
    gp_Dir interiorOpenToolDir;

    if (isInteriorOpen) {
        // ================= 中部开放型腔独立旁路 =================
        ProcessInteriorOpenCavityByBand(allClosedLayerFaces, allOpenLayerFaces, mainShape, splitPoints, faceToIdMap, savePath, interiorNewClosedParts, interiorNewOpenCavity, interiorOpenToolDir);
    } else {
        // ================= 提取封闭型腔特征面 =================
        ProcessAndSplitClosedCavityFeatures(allClosedLayerFaces, mainShape, splitPoints, faceToIdMap, savePath, trueClosedCavities, closedCavityFeatures);

        // ================= 提取开放型腔特征面 =================
        openCavityFeatures = ProcessAndSplitOpenCavityFeatures(allOpenLayerFaces, allClosedLayerFaces, trueClosedCavities, mainShape, splitPoints, faceToIdMap, savePath);
    }

    cout << "\n";
    cout << "======================================================================\n";
    cout << "                      特征提取结果汇总\n";
    cout << "======================================================================\n\n";

    if (isInteriorOpen) {
        cout << "【分支类型】中部开放型腔分支 (InteriorOpen)\n\n";

        cout << "--- 重组封闭型腔 (Z轴分割后) ---\n";
        cout << "  数量: " << interiorNewClosedParts.size() << " 个\n";
        for (size_t i = 0; i < interiorNewClosedParts.size(); ++i) {
            Bnd_Box box;
            BRepBndLib::Add(interiorNewClosedParts[i], box);
            double cxmin, cymin, czmin, cxmax, cymax, czmax;
            box.Get(cxmin, cymin, czmin, cxmax, cymax, czmax);
            cout << "  Part_" << i
                << "  类型=OTHER"
                << "  进刀方向=(0,0,-1)"
                << "  Z范围=[" << czmax << "," << czmin << "]"
                << "  深度=" << (czmax - czmin) << "mm\n";
        }

        cout << "\n--- 不可加工区域 (侧向加工) ---\n";
        Bnd_Box openBox;
        BRepBndLib::Add(interiorNewOpenCavity, openBox);
        double oxmin, oymin, ozmin, oxmax, oymax, ozmax;
        openBox.Get(oxmin, oymin, ozmin, oxmax, oymax, ozmax);
        cout << "  类型=OTHER"
            << "  进刀方向=(" << interiorOpenToolDir.X() << "," << interiorOpenToolDir.Y() << "," << interiorOpenToolDir.Z() << ")"
            << "  Z范围=[" << ozmax << "," << ozmin << "]"
            << "  深度=" << (ozmax - ozmin) << "mm\n";
    } else {
        cout << "【分支类型】正常分支\n\n";

        cout << "--- 封闭型腔 (Z轴分割后) ---\n";
        cout << "  数量: " << closedCavityFeatures.size() << " 个\n";
        for (const auto& feat : closedCavityFeatures) {
            cout << "  ID=" << feat.featureId
                << "  类型=CLOSED  进刀方向=(0,0,-1)"
                << "  Z范围=[" << feat.topZ << "," << feat.bottomZ << "]"
                << "  深度=" << feat.totalDepth << "mm"
                << "  切片层数=" << feat.stepLoops.size() << "\n";
        }

        cout << "\n--- 开放型腔 (加工可行性赛选后) ---\n";
        cout << "  保留: " << openCavityFeatures.size() << " 个\n";
        for (const auto& feat : openCavityFeatures) {
            cout << "  ID=" << feat.featureId
                << "  类型=" << (feat.type == CavityType::OPEN ? "OPEN" : feat.type == CavityType::CLOSED ? "CLOSED" : "OTHER")
                << "  进刀方向=(" << feat.toolDirection.X() << "," << feat.toolDirection.Y() << "," << feat.toolDirection.Z() << ")"
                << "  Z范围=[" << feat.topZ << "," << feat.bottomZ << "]"
                << "  深度=" << feat.totalDepth << "mm";
            if (feat.toolDirection.IsEqual(gp_Dir(0, 0, -1), 1e-6)) {
                cout << "  Z轴分割=" << feat.stepLoops.size() << "层";
            } else {
                cout << "  跳过Z轴分割(侧向)";
            }
            cout << "\n";
        }
    }

    cout << "\n======================================================================\n";
    cout << "生成的 BREP 文件保存在 " << savePath << " 目录。\n";
    cout << "======================================================================\n";

    system("pause");
    return 0;
}
#endif
