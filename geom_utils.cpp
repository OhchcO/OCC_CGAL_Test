#include "geom_utils.h"

bool IsPointEqual(const Point3D& p1, const Point3D& p2, double tol) {
    return fabs(p1.x - p2.x) < tol &&
        fabs(p1.y - p2.y) < tol;
}
// ================== 封闭型腔特征提取函数 ==================

// 从所有分层切片中提取出 CAVITY 面的 faceId 集合

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

bool IsValueBetween(double value, double a, double b, double tol) {
    return value >= std::min(a, b) - tol && value <= std::max(a, b) + tol;
}

bool IsPointOnSegment2D(const Point3D& p, const Point3D& a, const Point3D& b, double tol) {
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

