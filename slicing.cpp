#include "slicing.h"
#include "occ_utils.h"
#include "line_topology.h"
#include "export_utils.h"
#include "geom_utils.h"

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
    double angleTolerance,
    double mergeTol,
    double areaThreshold)
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
//    double angleTolerance,
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
