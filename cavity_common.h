#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#define NOMINMAX
#include <cmath>
#include <limits>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <direct.h>
#include <windows.h>

// CGAL 库
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/number_utils.h>

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
#pragma comment(lib, "TKFeat.lib")
#pragma comment(lib, "TKBool.lib")
#pragma comment(lib, "TKTopAlgo.lib")
#pragma comment(lib, "TKGeomAlgo.lib")

// STEP 读取
#include <STEPControl_Reader.hxx>
#include <IFSelect_ReturnStatus.hxx>

// 拓扑遍历与计算
#include <TopExp_Explorer.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>

// 求交算法
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
#include <BRepFeat_SplitShape.hxx>
#include <TopoDS_Iterator.hxx>

// 分割算法
#include <TopTools_MapOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <BRepExtrema_DistShapeShape.hxx>

#include <GCPnts_UniformAbscissa.hxx>
#include <CPnts_AbscissaPoint.hxx>

#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>

#include <BRepIntCurveSurface_Inter.hxx>
#include <gp_Lin.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>

// Clipper2 库
#include "clipper2/clipper.h"

using namespace std;

// ======================== 数据结构定义 ========================

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

typedef std::vector<OneLine> Loop;

struct LoopNode {
    Loop loopLines;
    double area;
    std::vector<LoopNode*> children;
    LoopNode() : area(0.0) {}
};

enum class FaceType {
    SOLID,      // 实体面（有材料）
    CAVITY,     // 封闭空腔面（无材料/孔洞）
    OPENCAVITY, // 开口空腔面（边界开口的孔洞）
    HULL        // 凸包面
};

struct Face2D {
    FaceType type;
    Loop outerLoop;
    std::vector<Loop> innerLoops;
};

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
    std::set<int> sourceFaceIds;
};

struct OpenCavityFilterParams {
    double minArea          = 2.5;
    double minCompactness   = 0.015;
    double minToolPassSpan  = 1.0;
    double maxFilletRadius  = 0.5;
    double hullMergeDistance = 360.0;
    double bboxPrecheckMargin = 0.0;
};

struct OpenCavitySplitParams {
    double minStepFaceArea = 5.0;
    double zProtectionTol = 0.05;
    double neighborZTol = 1e-3;
    double minSplitGap = 0.2;
    double geometricConnectTol = 0.05;
    int maxRecursionDepth = 32;
    bool splitDownOnlyStepFaces = true;
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

    UnionFind(int n = 0) { Reset(n); }

    void Reset(int n) {
        parent.resize(n);
        rank.assign(n, 0);
        triggerDistance.assign(n, std::numeric_limits<double>::max());
        for (int i = 0; i < n; ++i) parent[i] = i;
    }

    int Find(int x) {
        if (parent[x] != x) parent[x] = Find(parent[x]);
        return parent[x];
    }

    void Unite(int a, int b, double distance) {
        int rootA = Find(a);
        int rootB = Find(b);
        if (rootA == rootB) {
            triggerDistance[rootA] = std::min(triggerDistance[rootA], distance);
            return;
        }
        if (rank[rootA] < rank[rootB]) std::swap(rootA, rootB);
        parent[rootB] = rootA;
        triggerDistance[rootA] = std::min(triggerDistance[rootA], triggerDistance[rootB]);
        triggerDistance[rootA] = std::min(triggerDistance[rootA], distance);
        if (rank[rootA] == rank[rootB]) rank[rootA]++;
    }
};

// ======================== CGAL 类型 ========================

typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
typedef CGAL::Point_2<EK> Point_2;

// ======================== 全局变量 (extern) ========================

extern std::string inputPath;
extern std::string savePath;

extern const double SCALE;

extern BRep_Builder gBuilder;
extern TopoDS_Compound gAllFacesCompound;
extern TopoDS_Compound gAllLinesCompound;
extern bool gCompoundsInitialized;
