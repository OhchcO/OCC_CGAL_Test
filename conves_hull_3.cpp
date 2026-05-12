#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <vector>
#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Point_3                                Point_3;
typedef CGAL::Surface_mesh<Point_3>               Surface_mesh;


#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <BRepTools.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shell.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp_Explorer.hxx>

#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GCPnts_UniformAbscissa.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <STEPControl_Reader.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <TopoDS_Vertex.hxx>

#include <STEPControl_Writer.hxx>
#include <Interface_Static.hxx>

/**
 * @brief 将 OCC Shape 保存为 STEP 文件
 * @param shape 要导出的形状
 * @param filePath 保存路径 (例如 "E:/Result.stp")
 * @return true 导出成功, false 导出失败
 */
bool SaveShapeToStep(const TopoDS_Shape& shape, const std::string& filePath) {
    if (shape.IsNull()) {
        std::cerr << "Error: Shape is null, cannot save to STEP." << std::endl;
        return false;
    }

    STEPControl_Writer writer;

    // 1. 设置 STEP 文件的单位和模式（可选，默认通常是 MM）
    // 例如设置坐标精度或导出模式
    Interface_Static::SetCVal("write.step.unit", "MM");

    // 2. 将 Shape 转换并添加到 Writer 中
    // 第二个参数是导出模式，STEPControl_AsIs 表示按原样导出
    IFSelect_ReturnStatus status = writer.Transfer(shape, STEPControl_AsIs);

    if (status != IFSelect_RetDone) {
        std::cerr << "Error: Failed to transfer shape to STEP writer." << std::endl;
        return false;
    }

    // 3. 写入磁盘
    status = writer.Write(filePath.c_str());

    if (status != IFSelect_RetDone) {
        std::cerr << "Error: Failed to write STEP file to disk." << std::endl;
        return false;
    }

    std::cout << "Successfully saved STEP file to: " << filePath << std::endl;
    return true;
}

/**
 * @brief 从 STEP 模型中高质量提取点云
 * @param shape 输入的 OCC Shape
 * @param deflection 挠度（越小面上的点越密）
 * @param edgeSampleSpacing 边缘采样间距（mm，越小边上的点越密）
 */
std::vector<Point_3> ExtractPointsFromStep(const TopoDS_Shape& shape,
    double deflection = 0.01,
    double edgeSampleSpacing = 0.5) {
    std::vector<Point_3> cgal_points;
    BRep_Builder B;
    TopoDS_Compound pc;
    B.MakeCompound(pc);

    // --- 1. 离散化面 (Facets) ---
    BRepMesh_IncrementalMesh mesher(shape, deflection);

    TopExp_Explorer faceExp(shape, TopAbs_FACE);
    for (; faceExp.More(); faceExp.Next()) {
        TopoDS_Face face = TopoDS::Face(faceExp.Current());
        TopLoc_Location loc;
        Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);
        if (!tri.IsNull()) {
            for (int i = 1; i <= tri->NbNodes(); i++) {
                gp_Pnt p = tri->Node(i).Transformed(loc.Transformation());
                cgal_points.emplace_back(p.X(), p.Y(), p.Z());
            }
        }
    }

    // --- 2. 强制边缘采样 (Edge Sampling) - 关键改进！ ---
    // 解决平面模型只有角点的问题
    TopExp_Explorer edgeExp(shape, TopAbs_EDGE);
    for (; edgeExp.More(); edgeExp.Next()) {
        TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());
        BRepAdaptor_Curve curve(edge);

        // 使用均匀弦长采样
        GCPnts_UniformAbscissa sampler(curve, edgeSampleSpacing);
        if (sampler.IsDone()) {
            for (int i = 1; i <= sampler.NbPoints(); i++) {
                gp_Pnt p = curve.Value(sampler.Parameter(i));
                cgal_points.emplace_back(p.X(), p.Y(), p.Z());
            }
        }
    }

    std::cout << "Total points extracted: " << cgal_points.size() << std::endl;
    return cgal_points;
}

/**
 * @brief 将 CGAL 凸包多面体转换为 OpenCASCADE 干净的实体
 * @param poly 输入的 CGAL Polyhedron_3 对象
 * @return TopoDS_Shape 返回处理后的实体 (Solid)，如果失败则返回 Null Shape
 */
TopoDS_Shape ConvertCgalPolyToCleanSolid(const Polyhedron_3& poly) {
    if (poly.is_empty()) return TopoDS_Shape();

    // --- 1. 基础转换：将 CGAL Facets 转为 OCC Faces ---
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);

    for (auto f = poly.facets_begin(); f != poly.facets_end(); ++f) {
        auto h = f->facet_begin();
        BRepBuilderAPI_MakePolygon polygon;
        do {
            Point_3 p = h->vertex()->point();
            polygon.Add(gp_Pnt(CGAL::to_double(p.x()),
                CGAL::to_double(p.y()),
                CGAL::to_double(p.z())));
        } while (++h != f->facet_begin());

        polygon.Close();
        if (polygon.IsDone()) {
            BRepBuilderAPI_MakeFace faceMaker(polygon.Wire());
            if (faceMaker.IsDone()) {
                builder.Add(compound, faceMaker.Face());
            }
        }
    }

    // --- 2. 缝合 (Sewing)：将散面合成壳 (Shell) ---
    BRepBuilderAPI_Sewing sewer;
    sewer.SetTolerance(1e-6); // 设置缝合精度
    sewer.Add(compound);
    sewer.Perform();
    TopoDS_Shape sewedShape = sewer.SewedShape();

    // --- 3. 实体化 (MakeSolid)：将 Shell 封装成 Solid ---
    TopoDS_Shape finalResult = sewedShape;

    // 自动寻找 Shell（兼容 Compound 嵌套情况）
    TopExp_Explorer shellExp(sewedShape, TopAbs_SHELL);
    if (shellExp.More()) {
        TopoDS_Shell shell = TopoDS::Shell(shellExp.Current());
        BRepBuilderAPI_MakeSolid solidMaker(shell);
        if (solidMaker.IsDone()) {
            finalResult = solidMaker.Solid();
        }
    }

    // --- 4. 净化 (UnifySameDomain)：合并共面三角形，消除切割线 ---
    try {
        // 注意：UnifySameDomain 建议在 Solid 级别操作效果最好
        ShapeUpgrade_UnifySameDomain unifier(finalResult);
        unifier.Build();
        finalResult = unifier.Shape();
    }
    catch (...) {
        std::cerr << "Warning: UnifySameDomain failed, returning sewed shape." << std::endl;
    }

    return finalResult;
}


#if 0
int main() {
    //// 1. 读取点云并计算凸包
    //std::ifstream in((argc > 1) ? argv[1] : CGAL::data_file_path("E:\\soft\\CGAL\\CGAL-5.6.2\\data\\points_3\\cube.xyz"));
    //std::vector<Point_3> points;
    //Point_3 p;
    //while (in >> p) {
    //    points.push_back(p);
    //}
    //Polyhedron_3 poly;


    std::string stepFile = "E:\\soft\\code\\cMake_test\\input\\6_stp.stp";
    std::string savePath = "E:\\soft\\code\\cMake_test\\output\\";

    // 2. [OCC] 读取 STEP 模型
    STEPControl_Reader reader;
    if (reader.ReadFile(stepFile.c_str()) != IFSelect_RetDone) {
        std::cerr << "Error: Cannot read STEP file!" << std::endl;
        return -1;
    }
    reader.TransferRoots();
    TopoDS_Shape shape = reader.OneShape();

    // 3. 调用提取函数
    // 参数说明：shape, 挠度(0.01mm), 边缘采样间距(0.5mm)
    std::cout << "Starting point extraction..." << std::endl;
    std::vector<Point_3> cgal_points = ExtractPointsFromStep(shape, 0.001, 0.5);

    // 4. [可选] 将提取的点云保存为 BREP 验证
    BRep_Builder builder;
    TopoDS_Compound pcCompound;
    builder.MakeCompound(pcCompound);
    for (const auto& p : cgal_points) {
        gp_Pnt occPnt(p.x(), p.y(), p.z());
        // 1. 创建顶点
        TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(occPnt).Vertex();
        // 2. 将顶点加入 Compound
        builder.Add(pcCompound, v);
    }
    BRepTools::Write(pcCompound, (savePath + "2_Point_Cloud.brep").c_str());

    // 1. 计算凸包
    Polyhedron_3 poly;
    CGAL::convex_hull_3(cgal_points.begin(), cgal_points.end(), poly);

    // 2. 一键调用封装函数
    TopoDS_Shape myCleanSolid = ConvertCgalPolyToCleanSolid(poly);

    // 3. 检查并保存结果
    if (!myCleanSolid.IsNull() && myCleanSolid.ShapeType() == TopAbs_SOLID) {
        std::cout << "Success! Created a clean solid without extra edges." << std::endl;
        BRepTools::Write(myCleanSolid, (savePath + "Final_Clean_Solid.brep").c_str());
    }
    else {
        std::cerr << "Error: Could not create a valid solid." << std::endl;
    }

	// 4. 将最终结果保存为 STEP 文件
    SaveShapeToStep(myCleanSolid, (savePath + "Final_Product.stp").c_str());

    system("pause");
    return 0;
}
#endif