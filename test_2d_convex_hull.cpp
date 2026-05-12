#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <vector>
#include <iostream>
#include <iomanip>
// --- C++ 标准库 (解决 cout, endl) ---
#include <iostream>
#include <BRepTools.hxx>
// --- OpenCASCADE 核心拓扑 (解决 TopoDS_Edge, TopoDS_Compound) ---
#include <TopoDS_Edge.hxx>
#include <TopoDS_Compound.hxx>

// --- OpenCASCADE 基础几何 (解决 gp_Pnt) ---
#include <gp_Pnt.hxx>

// --- OpenCASCADE 建模算法 (解决 BRep_Builder, BRepBuilderAPI_MakeEdge) ---
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
using namespace std;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;

// --- 你的结构体定义 ---
struct Point3D { double x, y, z; };
struct OneLine { Point3D start; Point3D end; int faceId; };
typedef std::vector<OneLine> Loop;

struct Face2D {
    Loop outerLoop;
    // 这里暂时不处理 innerLoops
};

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

// --- 你要求的转换函数 ---
Loop ComputeConvexHullFromFace(const Face2D& face) {
    std::vector<Point_2> points_2d;
    if (face.outerLoop.empty()) return Loop();

    double currentZ = face.outerLoop[0].start.z;

    // 1. 收集所有端点
    for (const auto& line : face.outerLoop) {
        points_2d.emplace_back(line.start.x, line.start.y);
        points_2d.emplace_back(line.end.x, line.end.y);
    }

    // 2. 计算凸包
    std::vector<Point_2> result_points;
    CGAL::convex_hull_2(points_2d.begin(), points_2d.end(), std::back_inserter(result_points));

    // 3. 还原为 OneLine 闭合环
    Loop hullLoop;
    for (size_t i = 0; i < result_points.size(); ++i) {
        OneLine line;
        line.start = { CGAL::to_double(result_points[i].x()), CGAL::to_double(result_points[i].y()), currentZ };
        size_t nextIdx = (i + 1) % result_points.size();
        line.end = { CGAL::to_double(result_points[nextIdx].x()), CGAL::to_double(result_points[nextIdx].y()), currentZ };
        line.faceId = -1;
        hullLoop.push_back(line);
    }
    return hullLoop;
}

#if 0
int main() {
    string savePath = "E:\\soft\\code\\cMake_test\\output\\";
    // 1. 构建一个“开口”的 C 字形型腔（Z=5.0 平面）
    // 点：(0,0), (10,0), (10,10), (0,10) -> 此时 (0,10) 到 (0,0) 是断开的
    Face2D myFace;
    double Z = 5.0;
    myFace.outerLoop.push_back({ {0, 0, Z}, {10, 0, Z}, 1 });   // 底边
    myFace.outerLoop.push_back({ {10, 0, Z}, {10, 10, Z}, 1 }); // 右边
    myFace.outerLoop.push_back({ {10, 10, Z}, {0, 10, Z}, 1 });  // 顶边

    std::cout << "--- 原始数据：开口型腔 (3 条边) ---" << std::endl;
    std::cout << "Z 高度: " << Z << std::endl;

    // 2. 运行凸包补全
    Loop closedLoop = ComputeConvexHullFromFace(myFace);
    ExportOneLinesToBrep(closedLoop, savePath + "discretized.brep");
    // 3. 验证结果
    std::cout << "\n--- 补全后数据：凸包闭合环 (" << closedLoop.size() << " 条边) ---" << std::endl;
    std::cout << std::fixed << std::setprecision(1);

    bool hasCheckZ = true;
    for (size_t i = 0; i < closedLoop.size(); ++i) {
        const auto& l = closedLoop[i];
        std::cout << "线段 " << i << ": (" << l.start.x << ", " << l.start.y << ", " << l.start.z << ") -> ("
            << l.end.x << ", " << l.end.y << ", " << l.end.z << ")" << std::endl;

        if (l.start.z != Z || l.end.z != Z) hasCheckZ = false;
    }

    if (hasCheckZ) {
        std::cout << "\n验证成功：所有线段的 Z 值完美保持为 " << Z << std::endl;
    }
    else {
        std::cout << "\n验证失败：Z 值发生偏移！" << std::endl;
    }
    system("pause");
    return 0;
}
#endif