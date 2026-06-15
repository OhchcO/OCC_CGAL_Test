#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <direct.h>

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

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/number_utils.h>
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

// ======================== 函数声明 ========================

// --- 凸包相关 ---
std::vector<Point_2> ConvertFaceToPoints2D_EK(const Face2D& face);
Face2D ConvertPointsToHullFace_EK(vector<Point_2>& hull_points, double zHeight);
Face2D ComputeConvexHullFace(const Face2D& face);
Face2D ComputeMergedHullFace(const std::vector<HullItem>& hullItems, const std::vector<int>& memberIndices, double zHeight);
std::vector<Point_2> CollectFaceOuterPoints(const std::vector<Face2D>& faces);
std::vector<HullGroup> BuildMergedHullGroups(const std::vector<HullItem>& hullItems, const OpenCavityFilterParams& params);

// --- 调试/打印 ---
void PrintAllLines(const std::vector<OneLine>& lines);
void PrintEndpointGapDiagnostics(const std::vector<OneLine>& lines, double snapTol);

// --- Z 轴相关 ---
void GetFaceZRange(const TopoDS_Face& face, double& zmin, double& zmax);
void GetShapeZRange(const TopoDS_Shape& shape, double& zmin, double& zmax);
void GetExtremaZOfFace(const TopoDS_Face& face, const gp_Dir& direction, double& outMinZ, double& outMaxZ);
std::vector<double> GetSplitPointsAlongZ(const TopoDS_Shape& shape, double angleTolerance = 0.001, double mergeTol = 0.05, double areaThreshold = 1.0);

// --- 平面/边相关 ---
TopoDS_Compound GetCoplanarFaces(const TopoDS_Shape& shape, double splitZ, double tol = 1e-4);
bool IsConcave(const TopoDS_Edge& E, const TopoDS_Face& F1, const TopoDS_Face& F2);
TopoDS_Compound GetConcaveEdges(const TopoDS_Shape& solid, const TopoDS_Compound& coplanarFaces, const gp_Dir& toolDirection = gp_Dir(0, 0, 1));
bool IsConvex(const TopoDS_Edge& E, const TopoDS_Face& F1, const TopoDS_Face& F2);
TopoDS_Compound GetConvexEdgesFromOppositeFaces(const TopoDS_Shape& solid, const TopoDS_Compound& coplanarFaces, const gp_Dir& toolDirection = gp_Dir(0, 0, 1));
TopoDS_Compound ClassifyAndGetEdgesToRemove(const TopoDS_Shape& solid, const TopoDS_Compound& coplanarFaces, const gp_Dir& toolDirection);

// --- 线段布尔运算 ---
std::vector<OneEdge> SubtractLinesAndMapFaces(const TopoDS_Compound& linesA, const TopoDS_Compound& linesB, const TopTools_DataMapOfShapeShape& edgeToFaceMap);
TopoDS_Compound BuildCompoundFromOneEdgeVector(const std::vector<OneEdge>& edges);

// --- 全局 Compound ---
void InitGlobalCompounds();

// --- 点/线段工具 ---
bool IsPointEqual(const Point3D& p1, const Point3D& p2, double tol = 1e-3);
double Cross2D(const Point3D& a, const Point3D& b, const Point3D& c);
bool IsValueBetween(double value, double a, double b, double tol = 1e-9);
bool IsPointOnSegment2D(const Point3D& p, const Point3D& a, const Point3D& b, double tol = 1e-9);
bool SegmentsIntersect2D(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d);
double PointSegmentDistance2D(const Point3D& p, const Point3D& a, const Point3D& b);
double SegmentSegmentDistance2D(const OneLine& a, const OneLine& b);
double LoopDistance2D(const Loop& a, const Loop& b);

// --- 环/面工具 ---
double CalculateArea(const Loop& loop);
bool IsPointInLoop(const Point3D& pt, const Loop& loop);
Point3D GetFace2DCenter(const Face2D& face);
bool IsFaceInsideFace(const Face2D& innerFace, const Face2D& outerFace);
LoopBBox2D GetLoopBBox(const Loop& loop);
double BBoxDistance2D(const LoopBBox2D& a, const LoopBBox2D& b);

// --- 封闭型腔面 ID ---
std::set<int> CollectCavityFaceIds(const std::vector<std::vector<Face2D>>& allLayerFaces);
TopoDS_Compound GetFacesByFaceIds(const TopoDS_Shape& solid, const std::set<int>& targetFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap);
TopoDS_Compound GetCavityCapFaces(const TopoDS_Shape& solid, const std::set<int>& cavityFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, double angleTolerance = 0.001);
TopoDS_Compound GetCavityCompound(const TopoDS_Shape& solid, const std::set<int>& cavityFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// --- BREP 导出 ---
void ExportOneLinesToBrep(const std::vector<OneLine>& lines, const std::string& fileName);
void ExportCavityFaces(const TopoDS_Shape& solid, const std::set<int>& cavityFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, const std::string& fileName);
void ExportFace2DToBrep(const std::vector<Face2D>& faces, const std::string& fileName);
void ExportSoftEdgesToBrep(const std::vector<Face2D>& cavityFaces, const std::string& filePath);

// --- 线段去重/修复 ---
std::vector<OneLine> DeduplicateLines(const std::vector<OneLine>& lines);
int HealOpenEndpointGaps(std::vector<OneLine>& lines, double healTol);

// --- 成环算法 ---
std::vector<Loop> BuildLoops(const std::vector<OneLine>& inputLines);

// --- 拓扑构建 ---
std::vector<Face2D> BuildTopologyAndExtractFaces(const std::vector<OneLine>& layerLines);
void ExtractFacesFromTree(LoopNode* node, int depth, std::vector<Face2D>& faces);
void DebugExportLoops(const std::vector<Loop>& loops, const std::string& prefix, const std::string& savePath);

// --- 切片 ---
std::vector<Face2D> SliceModelAtZ(const TopoDS_Shape& shape, double splitZ, int sliceIndex, int maxIndex, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// --- 面属性判断 ---
bool IsHorizontalFace(const TopoDS_Face& face, double tol = 1e-4);
gp_Dir GetFaceNormal(const TopoDS_Face& face);
double GetHorizontalFaceZ(const TopoDS_Face& face);
double CalculateFaceAreaOCC(const TopoDS_Face& face);
int CountFacesInCompound(const TopoDS_Compound& compound);
gp_Dir SnapNormalToDirection(const gp_Dir& normal);
bool IsCavityMachinable(const TopoDS_Compound& cavity, gp_Dir& outToolDir);

// --- 封闭型腔切分 ---
void SplitCavityRecursive(const std::set<int>& currentFaceIds, TopTools_DataMapOfIntegerShape& idToFace, int& nextId, const std::vector<double>& splitZs, int zIndex, std::vector<std::set<int>>& result);
std::vector<TopoDS_Compound> SplitCavity(const TopoDS_Compound& cavity);
std::vector<TopoDS_Compound> SeparateDisconnectedCavities(const TopoDS_Compound& cavityCompound);

// --- 特征生成 ---
std::vector<CavityFeature> GenerateCavityFeatures(const std::vector<TopoDS_Compound>& isolatedCavities, const std::vector<std::vector<Face2D>>& allLayerFaces, const std::vector<double>& sortedSplitPoints, const TopTools_DataMapOfShapeInteger& faceToIdMap);
TopoDS_Compound ExtractTrueCavityFaces(const CavityFeature& feat, const TopoDS_Shape& mainShape, const TopTools_DataMapOfShapeInteger& faceToIdMap);
std::vector<CavityFeature> ConvertPartsToFeatures(const std::vector<TopoDS_Compound>& parts, const CavityFeature& parentFeat, const TopTools_DataMapOfShapeInteger& faceToIdMap);

// --- Clipper2 布尔运算 ---
using namespace Clipper2Lib;
Path64 LoopToPath(const Loop& loop);
Loop PathToLoop(const Path64& path, double zHeight);
Paths64 FacesToClipperPaths(const std::vector<Face2D>& faces);
void ExtractFacesFromPolyNode(const PolyPath64* node, int depth, double zHeight, std::vector<Face2D>& faces);
std::vector<Face2D> BooleanFacesSingle(const Face2D& subject, const Face2D& clip, ClipType clipType, double zHeight);
std::vector<Face2D> BooleanFaces(const std::vector<Face2D>& subjects, const std::vector<Face2D>& clips, ClipType clipType, double zHeight);
std::vector<Face2D> RecoverOriginalFaceIdsByGeometry(const std::vector<Face2D>& pocketResults, const std::vector<Face2D>& currentSolidFaces);

// --- 开放型腔清洗 ---
double CalculateFaceCompactness(const Face2D& face);
std::map<int, double> BuildFilletRadiusMap(const TopoDS_Shape& shape, const TopTools_DataMapOfShapeInteger& faceToIdMap);
bool IsSmallFilletCavity(const Face2D& face, const std::map<int, double>& filletMap, double maxFilletRadius, std::vector<std::pair<int, double>>* matchedFillets = nullptr);
std::vector<Face2D> CleanAndFilterOpenCavities(const std::vector<Face2D>& rawOpenCavities, const std::map<int, double>& filletMap, const OpenCavityFilterParams& params = OpenCavityFilterParams{});

// --- 开放型腔提取 ---
std::set<int> CollectOpenCavityFaceIds(const std::vector<std::vector<Face2D>>& allLayerFaces);
TopoDS_Compound GetOpenCavityCapFaces(const TopoDS_Shape& solid, const std::set<int>& openFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, double angleTolerance = 0.001);
void ExportOpenCavityFaces(const TopoDS_Shape& solid, const std::set<int>& openFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, const std::string& fileName);
TopoDS_Compound GetOpenCavityCompound(const TopoDS_Shape& solid, const std::set<int>& openFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap);
std::vector<CavityFeature> GenerateOpenCavityFeatures(const std::vector<TopoDS_Compound>& isolatedCavities, const std::vector<std::vector<Face2D>>& allLayerFaces, const std::vector<double>& sortedSplitPoints, const TopTools_DataMapOfShapeInteger& faceToIdMap);
TopoDS_Compound ExtractTrueOpenCavityFaces(const CavityFeature& openFeat, const TopoDS_Shape& mainShape, const TopTools_DataMapOfShapeInteger& faceToIdMap, const std::vector<std::vector<Face2D>>& allClosedLayerFaces, const std::vector<double>& splitPoints, const std::vector<TopoDS_Compound>& trueClosedCavities);

// --- 嵌套合并 ---
bool IsCompoundInsideExact(const TopoDS_Compound& innerComp, const TopoDS_Compound& outerComp);
std::vector<TopoDS_Compound> MergeNestedIslands(const std::vector<TopoDS_Compound>& bfsParts);

// --- Compound Z 分割 ---
std::pair<TopoDS_Compound, TopoDS_Compound> SplitCompoundAtZ(const TopoDS_Compound& cavity, double splitZ);
std::vector<TopoDS_Compound> RecursiveSplitCavity(const TopoDS_Compound& cavity);
std::vector<double> MergeCloseSplitZs(std::vector<double> splitZs, double minGap);

// --- 开放型腔 Z 切分 ---
std::vector<double> CollectOpenCavitySplitZs(const TopoDS_Compound& cavity, const CavityFeature& feature, const OpenCavitySplitParams& params);
bool GetCompoundZRangeSafe(const TopoDS_Compound& compound, double& zmin, double& zmax);
std::vector<double> CollectOpenCavitySplitZsForPart(const TopoDS_Compound& part, const CavityFeature& feature, const OpenCavitySplitParams& params, const std::string& branchName);
std::vector<TopoDS_Compound> MergeGeometricallyConnectedOpenParts(const std::vector<TopoDS_Compound>& rawParts, double tolerance, int featureId, const std::string& branchName, int depth, bool exportDebugBreps);
std::vector<TopoDS_Compound> SplitOpenCavityByZPlan(const TopoDS_Compound& cavity, const CavityFeature& feature, const OpenCavitySplitParams& params, int featureId);
std::vector<TopoDS_Compound> SplitOpenCavityConservatively(const TopoDS_Compound& cavity, const CavityFeature& feature, const OpenCavitySplitParams& params);

// --- 主流程 ---
std::vector<CavityFeature> ProcessAndSplitOpenCavityFeatures(const std::vector<std::vector<Face2D>>& allOpenLayerFaces, const std::vector<std::vector<Face2D>>& allClosedLayerFaces, const std::vector<TopoDS_Compound>& trueClosedCavities, const TopoDS_Shape& mainShape, const std::vector<double>& splitPoints, const TopTools_DataMapOfShapeInteger& faceToIdMap, const std::string& savePath);
void ProcessAndSplitClosedCavityFeatures(const std::vector<std::vector<Face2D>>& allLayerFaces, const TopoDS_Shape& mainShape, const std::vector<double>& splitPoints, const TopTools_DataMapOfShapeInteger& faceToIdMap, const std::string& savePath, std::vector<TopoDS_Compound>& outTrueClosedCavities, std::vector<CavityFeature>& outFeatures);

// --- 中部开放型腔 ---
double GetFace2DZForInteriorOpen(const Face2D& face);
bool HasInteriorOpenCavity(const std::vector<std::vector<Face2D>>& allOpenLayerFaces, const std::vector<std::vector<Face2D>>& allClosedLayerFaces, double modelMaxZ, const std::vector<double>& splitPoints, double zTol = 0.5);
double FindNextLowerSplitPoint(const std::vector<double>& splitPoints, double z, double fallbackZ);
void AddFaceIdsFromFace2D(const Face2D& face, std::set<int>& ids);
void PrintIdSet(const std::string& label, const std::set<int>& ids);
std::set<int> CollectClosedSideFaceIdsAroundOpenBand(const std::vector<std::vector<Face2D>>& allClosedLayerFaces, double openTopZ, double openBottomZ);
bool ShouldDiscardInteriorOpenRegion(const std::set<int>& openSideIds, const std::set<int>& closedSideIds, const TopoDS_Compound& topCapFaces, std::set<int>& sharedSideIds, double minSharedRatio = 0.8);
CavityFeature BuildInteriorOpenCavityFeature(const std::vector<std::vector<Face2D>>& allOpenLayerFaces, const std::vector<std::vector<Face2D>>& allClosedLayerFaces, const std::vector<double>& splitPoints, double modelMinZ, const std::string& savePath);
TopoDS_Compound ExtractFacesByFeatureSourceIds(const CavityFeature& feature, const TopoDS_Shape& mainShape, const TopTools_DataMapOfShapeInteger& faceToIdMap);
TopoDS_Compound GetInteriorOpenBandCapFaces(const TopoDS_Shape& solid, const std::set<int>& sideFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, double topZ, double bottomZ, bool wantTopCap, double zTol = 0.1);
TopoDS_Compound GetAdjacentHorizontalCapFacesAtZ(const TopoDS_Shape& solid, const std::set<int>& sideFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, double targetZ, bool normalUp, const std::string& debugLabel, double zTol = 0.1);
TopoDS_Compound GetHorizontalFacesForBand(const TopoDS_Shape& solid, const std::set<int>& sideFaceIds, const TopTools_DataMapOfShapeInteger& faceToIdMap, double bandMinZ, double bandMaxZ, const std::string& debugLabel, double zTol = 0.1);
void AddCompoundFaces(TopoDS_Compound& target, const TopoDS_Compound& source);
bool GetCompoundZRange(const TopoDS_Compound& compound, double& zmin, double& zmax);
void GetCompoundXYRange(const TopoDS_Compound& compound, double& xmin, double& ymin, double& xmax, double& ymax);
bool IsFaceInsideXYBox(const TopoDS_Face& face, double xmin, double ymin, double xmax, double ymax);
bool IsFaceCrossingXYBox(const TopoDS_Face& face, double xmin, double ymin, double xmax, double ymax);
void SplitFacesByXYBox(const TopoDS_Compound& inputFaces, double xmin, double ymin, double xmax, double ymax, TopoDS_Compound& insideFaces, TopoDS_Compound& crossingFaces);
Face2D TranslateFace2DZ(const Face2D& face, double newZ);
std::set<int> CollectFaceIdsFromFace2D(const Face2D& face);
std::set<int> CollectFaceIdsFromCompound(const TopoDS_Compound& compound, const TopTools_DataMapOfShapeInteger& faceToIdMap);
bool IsFaceOverlappingFace2D(const TopoDS_Face& face3D, const Face2D& region2D);
void SplitFacesByClipperRegion(const TopoDS_Compound& inputFaces, const Face2D& machinableRegion, TopoDS_Compound& machinableFaces, TopoDS_Compound& nonMachinableFaces);
void SplitInteriorOpenSourceFacesByBand(const TopoDS_Compound& sourceFaces, double topZ, double bottomZ, TopoDS_Compound& upperClosed, TopoDS_Compound& middleOpen, TopoDS_Compound& lowerClosed);
void ProcessInteriorOpenCavityByBand(const std::vector<std::vector<Face2D>>& allClosedLayerFaces, const std::vector<std::vector<Face2D>>& allOpenLayerFaces, const TopoDS_Shape& mainShape, const std::vector<double>& splitPoints, const TopTools_DataMapOfShapeInteger& faceToIdMap, const std::string& savePath, std::vector<TopoDS_Compound>& outNewClosedParts, TopoDS_Compound& outNewOpenCavity, gp_Dir& outOpenToolDir);
