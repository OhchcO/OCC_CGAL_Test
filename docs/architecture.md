# 项目架构说明

## 1. 当前文件结构

```text
cMake_test/
├── main.cpp
├── cavity_common.h
├── cavity_globals.cpp
├── geom_utils.h
├── geom_utils.cpp
├── convex_hull.h
├── convex_hull.cpp
├── slicing.h
├── slicing.cpp
├── occ_utils.h
├── occ_utils.cpp
├── boolean_ops.h
├── boolean_ops.cpp
├── line_topology.h
├── line_topology.cpp
├── export_utils.h
├── export_utils.cpp
├── closed_cavity.h
├── closed_cavity.cpp
├── open_cavity.h
├── open_cavity.cpp
├── interior_open.h
├── interior_open.cpp
├── visualize_cavities.FCMacro
└── CMakeLists.txt
```

---

## 2. 各模块职责

### `cavity_common.h`
- 存放所有公共类型、结构体、枚举和全局变量声明。
- 核心数据结构：
  - `Point3D`
  - `OneLine`
  - `Loop`
  - `Face2D`
  - `CavityFeature`
  - `HullItem`
  - `HullGroup`
  - `OpenCavityFilterParams`
  - `OpenCavitySplitParams`

### `cavity_globals.cpp`
- 定义全局变量：
  - `inputPath`
  - `savePath`
  - `SCALE`
  - `gBuilder`
  - `gAllFacesCompound`
  - `gAllLinesCompound`
  - `gCompoundsInitialized`

### `geom_utils.h / geom_utils.cpp`
底层 2D 几何工具。
主要函数：
- `IsPointEqual`
- `Cross2D`
- `IsValueBetween`
- `IsPointOnSegment2D`
- `SegmentsIntersect2D`
- `PointSegmentDistance2D`
- `SegmentSegmentDistance2D`
- `LoopDistance2D`
- `CalculateArea`
- `IsPointInLoop`
- `GetFace2DCenter`
- `IsFaceInsideFace`
- `GetLoopBBox`
- `BBoxDistance2D`
- `TranslateFace2DZ`

### `convex_hull.h / convex_hull.cpp`
凸包计算与合并。
主要函数：
- `ConvertFaceToPoints2D_EK`
- `ConvertPointsToHullFace_EK`
- `ComputeConvexHullFace`
- `ComputeMergedHullFace`
- `CollectFaceOuterPoints`
- `BuildMergedHullGroups`

### `slicing.h / slicing.cpp`
Z 轴分析与切片逻辑。
主要函数：
- `GetFaceZRange`
- `GetShapeZRange`
- `GetExtremaZOfFace`
- `GetSplitPointsAlongZ`
- `SliceModelAtZ`

### `occ_utils.h / occ_utils.cpp`
OpenCascade 辅助工具。
主要函数：
- `GetCoplanarFaces`
- `IsConcave`
- `GetConcaveEdges`
- `IsConvex`
- `GetConvexEdgesFromOppositeFaces`
- `ClassifyAndGetEdgesToRemove`
- `SubtractLinesAndMapFaces`
- `BuildCompoundFromOneEdgeVector`
- `InitGlobalCompounds`
- `IsHorizontalFace`
- `GetFaceNormal`
- `GetHorizontalFaceZ`
- `CalculateFaceAreaOCC`
- `CountFacesInCompound`
- `CollectCavityFaceIds`
- `GetFacesByFaceIds`
- `GetCavityCapFaces`
- `GetCavityCompound`
- `AddCompoundFaces`
- `GetCompoundZRange`

### `boolean_ops.h / boolean_ops.cpp`
基于 Clipper2 的 2D 布尔运算。
主要函数：
- `LoopToPath`
- `PathToLoop`
- `FacesToClipperPaths`
- `ExtractFacesFromPolyNode`
- `BooleanFacesSingle`
- `BooleanFaces`
- `RecoverOriginalFaceIdsByGeometry`

### `line_topology.h / line_topology.cpp`
线段处理与拓扑重建。
主要函数：
- `PrintAllLines`
- `PrintEndpointGapDiagnostics`
- `DeduplicateLines`
- `HealOpenEndpointGaps`
- `BuildLoops`
- `BuildTopologyAndExtractFaces`
- `ExtractFacesFromTree`
- `DebugExportLoops`

### `export_utils.h / export_utils.cpp`
导出与调试输出工具。
主要函数：
- `ExportOneLinesToBrep`
- `ExportCavityFaces`
- `ExportFace2DToBrep`
- `ExportSoftEdgesToBrep`

### `closed_cavity.h / closed_cavity.cpp`
封闭型腔处理管线。
主要函数：
- `SplitCavityRecursive`
- `SplitCavity`
- `SeparateDisconnectedCavities`
- `GenerateCavityFeatures`
- `ExtractTrueCavityFaces`
- `ConvertPartsToFeatures`
- `ProcessAndSplitClosedCavityFeatures`

### `open_cavity.h / open_cavity.cpp`
开放型腔处理管线。
主要函数：
- `CalculateFaceCompactness`
- `BuildFilletRadiusMap`
- `IsSmallFilletCavity`
- `CleanAndFilterOpenCavities`
- `CollectOpenCavityFaceIds`
- `GetOpenCavityCapFaces`
- `ExportOpenCavityFaces`
- `GetOpenCavityCompound`
- `GenerateOpenCavityFeatures`
- `ExtractTrueOpenCavityFaces`
- `IsCompoundInsideExact`
- `MergeNestedIslands`
- `SplitCompoundAtZ`
- `RecursiveSplitCavity`
- `MergeCloseSplitZs`
- `CollectOpenCavitySplitZs`
- `GetCompoundZRangeSafe`
- `CollectOpenCavitySplitZsForPart`
- `MergeGeometricallyConnectedOpenParts`
- `SplitOpenCavityByZPlan`
- `SplitOpenCavityConservatively`
- `ProcessAndSplitOpenCavityFeatures`

### `interior_open.h / interior_open.cpp`
中部开放型腔特殊分支。
主要函数：
- `GetFace2DZForInteriorOpen`
- `HasInteriorOpenCavity`
- `FindNextLowerSplitPoint`
- `AddFaceIdsFromFace2D`
- `PrintIdSet`
- `CollectClosedSideFaceIdsAroundOpenBand`
- `ShouldDiscardInteriorOpenRegion`
- `BuildInteriorOpenCavityFeature`
- `ExtractFacesByFeatureSourceIds`
- `GetInteriorOpenBandCapFaces`
- `GetAdjacentHorizontalCapFacesAtZ`
- `GetHorizontalFacesForBand`
- `GetCompoundXYRange`
- `IsFaceInsideXYBox`
- `IsFaceCrossingXYBox`
- `SplitFacesByXYBox`
- `CollectFaceIdsFromFace2D`
- `CollectFaceIdsFromCompound`
- `IsFaceOverlappingFace2D`
- `SplitFacesByClipperRegion`
- `SplitInteriorOpenSourceFacesByBand`
- `ProcessInteriorOpenCavityByBand`

---

## 3. 主调用链

```text
main()
  ├── 读取 STEP
  ├── 构建圆角半径映射 BuildFilletRadiusMap
  ├── 计算切分点 GetSplitPointsAlongZ
  ├── 逐层遍历切分点
  │     └── SliceModelAtZ
  │           ├── GetCoplanarFaces
  │           ├── ClassifyAndGetEdgesToRemove
  │           ├── SubtractLinesAndMapFaces
  │           ├── BuildTopologyAndExtractFaces
  │           │     ├── DeduplicateLines
  │           │     ├── HealOpenEndpointGaps
  │           │     └── BuildLoops
  │           └── ExportFace2DToBrep / ExportOneLinesToBrep
  ├── 构建凸包
  │     ├── ConvertFaceToPoints2D_EK
  │     ├── ComputeConvexHullFace
  │     └── BuildMergedHullGroups
  ├── 开放型腔初步过滤
  │     ├── BooleanFacesSingle
  │     ├── RecoverOriginalFaceIdsByGeometry
  │     └── CleanAndFilterOpenCavities
  ├── HasInteriorOpenCavity
  │     ├── 是 → ProcessInteriorOpenCavityByBand
  │     │           ├── BuildInteriorOpenCavityFeature
  │     │           ├── ExtractFacesByFeatureSourceIds
  │     │           └── GetInteriorOpenBandCapFaces
  │     └── 否
  │           ├── ProcessAndSplitClosedCavityFeatures
  │           │     ├── SeparateDisconnectedCavities
  │           │     ├── SplitCavity
  │           │     ├── GenerateCavityFeatures
  │           │     └── ExtractTrueCavityFaces
  │           └── ProcessAndSplitOpenCavityFeatures
  │                 ├── CollectOpenCavityFaceIds
  │                 ├── GetOpenCavityCompound
  │                 ├── MergeNestedIslands
  │                 ├── SplitOpenCavityConservatively
  │                 └── ExtractTrueOpenCavityFaces
```

---

## 4. 依赖关系

```text
cavity_common.h
    ^
    |-- 所有模块都依赖

geom_utils
    ^
    |-- convex_hull
    |-- occ_utils
    |-- export_utils
    |-- slicing
    |-- line_topology
    |-- closed_cavity
    |-- open_cavity
    |-- interior_open

occ_utils
    ^
    |-- slicing
    |-- export_utils
    |-- closed_cavity
    |-- open_cavity
    |-- interior_open

boolean_ops
    ^
    |-- open_cavity
    |-- interior_open

line_topology
    ^
    |-- slicing

export_utils
    ^
    |-- slicing
    |-- line_topology
    |-- closed_cavity
    |-- open_cavity
    |-- interior_open

slicing
    ^
    |-- closed_cavity
    |-- open_cavity
    |-- interior_open

closed_cavity
open_cavity
interior_open
    ^
    |-- main.cpp
```

---

## 5. 总结

当前架构是一个清晰的分层结构：

1. **公共层**  
   `cavity_common.h` + `cavity_globals.cpp`

2. **工具层**  
   `geom_utils`、`occ_utils`、`export_utils`、`boolean_ops`、`line_topology`

3. **算法层**  
   `convex_hull`、`slicing`

4. **业务管线层**  
   `closed_cavity`、`open_cavity`、`interior_open`

5. **入口层**  
   `main.cpp`

这种结构的优点是：
- 几何工具和业务逻辑分离
- 封闭、开放、中部开放三条管线彼此独立
- 导出和调试工具集中管理
- 主流程只负责流程编排
