// 1. 屏蔽 Windows 的 min/max 干扰
#ifndef NOMINMAX
#define NOMINMAX
#endif

// 2. 先包含标准库，奠定基础环境
#include <iostream>
#include <vector>
#include <mutex>
#include <string>

// 3. 包含 CGAL 头文件（趁环境还干净，先让 CGAL 解析完）
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

// 4. [核心技巧] 定义 OCC 宏保护，防止其 Handle 干扰全局
// 有时候需要暂时取消某些宏定义，但在包含 OCC 前通常不需要。
// 包含 OCC 核心头文件
#include <Standard_Handle.hxx> // 先单独引入 OCC 的句柄系统
#include <STEPControl_Reader.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Face.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <gp_Pnt.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepLib.hxx>
// 定义 CGAL 内核
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;

#if 0
int main() {
    std::string savePath = "E:\\soft\\code\\cMake_test\\output\\"; 
    try {
        // 1. [OCC] 读取 STEP 模型
        std::string fileName = "E:\\soft\\code\\Project1\\input\\01.stp"; // 替换为你的文件名
        STEPControl_Reader reader;
        if (reader.ReadFile(fileName.c_str()) != IFSelect_RetDone) {
            std::cerr << "Error: Cannot read STEP file." << std::endl;
            // 演示目的：如果文件不存在，我们手动创建一个方块点云继续演示
        }
        reader.TransferRoots();
        TopoDS_Shape shape = reader.OneShape();

        // 2. [OCC] 将模型离散化 (Mesh)
        // 参数 0.1 是挠度值，越小网格越密
        BRepMesh_IncrementalMesh mesher(shape, 0.0001);
        BRepTools::Write(shape, (savePath + "1_Mesh_Model.brep").c_str());

        // 3. [OCC -> CGAL] 提取 Mesh 顶点并转换为 CGAL 点
        std::vector<Point_3> cgal_points;
        TopExp_Explorer faceExp(shape, TopAbs_FACE);
        BRep_Builder B;
        TopoDS_Compound pointCloudCompound;
        B.MakeCompound(pointCloudCompound);
        for (; faceExp.More(); faceExp.Next()) {
            TopoDS_Face face = TopoDS::Face(faceExp.Current());
            TopLoc_Location loc;
            Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, loc);

            if (!triangulation.IsNull()) {
                for (int i = 1; i <= triangulation->NbNodes(); i++) {
                    gp_Pnt p = triangulation->Node(i).Transformed(loc.Transformation());
                    cgal_points.emplace_back(p.X(), p.Y(), p.Z());
                    B.Add(pointCloudCompound, BRepBuilderAPI_MakeVertex(p).Vertex());
                }
            }
        }
        BRepTools::Write(pointCloudCompound, (savePath + "2_Point_Cloud.brep").c_str());
        std::cout << "Extracted " << cgal_points.size() << " points from OCC Mesh." << std::endl;


        // 在 cgal_points 提取循环结束后，强制加入 8 个标准角点
        cgal_points.clear(); // 先清空提取到的点
        cgal_points.emplace_back(0, 0, 0);
        cgal_points.emplace_back(10, 0, 0);
        cgal_points.emplace_back(0, 10, 0);
        cgal_points.emplace_back(10, 10, 0);
        cgal_points.emplace_back(0, 0, 10);
        cgal_points.emplace_back(10, 0, 10);
        cgal_points.emplace_back(0, 10, 10);
        cgal_points.emplace_back(10, 10, 10);

        std::cout << "DEBUG: Manually injected 8 corners of a cube." << std::endl;
        // 后面继续跑 hull 计算和导出
        // 4. [CGAL] 计算三维凸包
        Polyhedron hull;
        CGAL::convex_hull_3(cgal_points.begin(), cgal_points.end(), hull);
        std::cout << "CGAL Convex Hull computed. Vertices: " << hull.size_of_vertices() << std::endl;
        TopoDS_Compound hullFacesCompound;
        B.MakeCompound(hullFacesCompound);

        for (auto f = hull.facets_begin(); f != hull.facets_end(); ++f) {
            auto h = f->facet_begin();
            gp_Pnt p1(h->vertex()->point().x(), h->vertex()->point().y(), h->vertex()->point().z());
            gp_Pnt p2((++h)->vertex()->point().x(), (++h)->vertex()->point().y(), (++h)->vertex()->point().z());
            gp_Pnt p3((++h)->vertex()->point().x(), (++h)->vertex()->point().y(), (++h)->vertex()->point().z());

            // 关键：Polygon (线) -> Wire (封闭线) -> Face (面)
            BRepBuilderAPI_MakePolygon poly(p1, p2, p3, Standard_True);
            if (poly.IsDone()) {
                BRepBuilderAPI_MakeFace mf(poly.Wire());
                if (mf.IsDone()) {
                    B.Add(hullFacesCompound, mf.Face()); // 存入散面集合
                }
            }
        }
        BRepTools::Write(hullFacesCompound, (savePath + "3_Hull_Faces.brep").c_str());
        // 5. [CGAL -> OCC] 重建凸包实体 (以最简单的三角形面为例)
        BRepBuilderAPI_Sewing sewer(1e-6);
        sewer.Add(hullFacesCompound); // 直接把刚才的散面集合丢进去
        sewer.Perform();
        TopoDS_Shape sewedShape = sewer.SewedShape();
        BRepTools::Write(sewedShape, (savePath + "4_Sewed_Shell.brep").c_str());

        // 5. 实体化 (Solid)
        if (sewedShape.ShapeType() == TopAbs_SHELL) {
            BRepBuilderAPI_MakeSolid ms(TopoDS::Shell(sewedShape));
            if (ms.IsDone()) {
                BRepTools::Write(ms.Solid(), (savePath + "5_Final_Solid.brep").c_str());
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
	system("pause"); // Windows 下暂停，方便查看输出
    return 0;
}
#endif