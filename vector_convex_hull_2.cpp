#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>

#include <vector>

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

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef std::vector<Point_2> Points;

using namespace std;

//结构体定义
struct Point3D {
    double x, y, z;
};

struct OneLine {
    Point3D start;
    Point3D end;
    int faceId;
};


/**
 * @brief 将 CGAL 凸包顶点序列转换为有序线段集
 * @param hullPoints CGAL 计算出的有序顶点 (vector<Point_2>)
 * @param z 该层在 3D 空间中的高度
 * @param faceId 默认赋予的面标识符，默认为 -1
 * @return std::vector<OneLine> 首尾相连的线段集合
 */
std::vector<OneLine> ConvertHullToOneLines(const std::vector<Point_2>& hullPoints, double z, int faceId = -1) {
    std::vector<OneLine> lines;

    // 凸包至少需要 2 个点才能构成线段，3 个点才能构成封闭区域
    if (hullPoints.size() < 2) {
        return lines;
    }

    for (size_t i = 0; i < hullPoints.size(); ++i) {
        OneLine line;

        // 1. 设置起点坐标
        line.start.x = CGAL::to_double(hullPoints[i].x());
        line.start.y = CGAL::to_double(hullPoints[i].y());
        line.start.z = z;

        // 2. 设置终点坐标
        // 使用取模运算 (i + 1) % size，当 i 是最后一个点时，nextIdx 会指向 0，从而闭合环路
        size_t nextIdx = (i + 1) % hullPoints.size();
        line.end.x = CGAL::to_double(hullPoints[nextIdx].x());
        line.end.y = CGAL::to_double(hullPoints[nextIdx].y());
        line.end.z = z;

        // 3. 设置标识符
        line.faceId = faceId;

        lines.push_back(line);
    }

    return lines;
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

#if 0
int main()
{
string savePath = "E:\\soft\\code\\cMake_test\\output\\";
  Points points, result;
  points.push_back(Point_2(0,0));
  points.push_back(Point_2(10,0));
  points.push_back(Point_2(10,10));
  points.push_back(Point_2(6,5));
  points.push_back(Point_2(4,1));


  CGAL::convex_hull_2( points.begin(), points.end(), back_inserter(result) );
  cout << result.size() << " points on the convex hull" << endl;

  vector<OneLine> lines = ConvertHullToOneLines(result, 0.0);
  ExportOneLinesToBrep(lines, savePath + "discretized.brep");
  system("pause");
  return 0;
}
#endif