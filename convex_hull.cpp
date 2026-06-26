#include "convex_hull.h"
#include "geom_utils.h"

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
