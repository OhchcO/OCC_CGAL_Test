#include "line_topology.h"
#include "export_utils.h"
#include "geom_utils.h"

// --- PrintAllLines ---
void PrintAllLines(const std::vector<OneLine>& lines)
{
    printf("========== ALL LINES (%d) ==========\n", (int)lines.size());
    for (int i = 0; i < lines.size(); i++) {
        const auto& l = lines[i];
        printf("[%4d] (%.6f, %.6f) -> (%.6f, %.6f)\n",
            i,
            l.start.x, l.start.y,
            l.end.x, l.end.y);
    }
    printf("=====================================\n");
}
// 获取面的 Z 范围

// --- 线段去重/修复 ---
void PrintEndpointGapDiagnostics(const std::vector<OneLine>& lines, double snapTol) {
    struct EndpointInfo {
        Point3D p;
        int lineIndex = -1;
        bool isStart = true;
        int degree = 0;
    };

    std::vector<EndpointInfo> endpoints;
    endpoints.reserve(lines.size() * 2);
    for (int i = 0; i < (int)lines.size(); ++i) {
        endpoints.push_back({ lines[i].start, i, true, 0 });
        endpoints.push_back({ lines[i].end, i, false, 0 });
    }

    for (size_t i = 0; i < endpoints.size(); ++i) {
        for (size_t j = i + 1; j < endpoints.size(); ++j) {
            double dx = endpoints[i].p.x - endpoints[j].p.x;
            double dy = endpoints[i].p.y - endpoints[j].p.y;
            double dist = std::sqrt(dx * dx + dy * dy);
            if (dist < snapTol) {
                endpoints[i].degree++;
                endpoints[j].degree++;
            }
        }
    }

    int openEndpointCount = 0;
    double maxNearestGap = 0.0;
    std::cout << "  [BuildLoops诊断] 未用线段端点断口:" << std::endl;
    for (size_t i = 0; i < endpoints.size(); ++i) {
        if (endpoints[i].degree > 0) continue;
        openEndpointCount++;

        double bestDist = std::numeric_limits<double>::max();
        int bestIndex = -1;
        for (size_t j = 0; j < endpoints.size(); ++j) {
            if (i == j) continue;
            double dx = endpoints[i].p.x - endpoints[j].p.x;
            double dy = endpoints[i].p.y - endpoints[j].p.y;
            double dist = std::sqrt(dx * dx + dy * dy);
            if (dist < bestDist) {
                bestDist = dist;
                bestIndex = (int)j;
            }
        }

        maxNearestGap = std::max(maxNearestGap, bestDist);
        std::cout << "    line=" << endpoints[i].lineIndex
            << (endpoints[i].isStart ? ".start" : ".end")
            << " p=(" << endpoints[i].p.x << "," << endpoints[i].p.y << "," << endpoints[i].p.z << ")"
            << " nearest=" << bestDist;
        if (bestIndex >= 0) {
            std::cout << " -> line=" << endpoints[bestIndex].lineIndex
                << (endpoints[bestIndex].isStart ? ".start" : ".end")
                << " p=(" << endpoints[bestIndex].p.x << "," << endpoints[bestIndex].p.y << "," << endpoints[bestIndex].p.z << ")";
        }
        std::cout << std::endl;
    }

    std::cout << "  [BuildLoops诊断] openEndpointCount=" << openEndpointCount
        << " maxNearestGap=" << maxNearestGap
        << " snapTol=" << snapTol << std::endl;
}

int HealOpenEndpointGaps(std::vector<OneLine>& lines, double healTol) {
    struct EndpointRef {
        int lineIndex = -1;
        bool isStart = true;
        Point3D p;
        int degree = 0;
    };

    auto Distance2DLocal = [](const Point3D& a, const Point3D& b) {
        double dx = a.x - b.x;
        double dy = a.y - b.y;
        return std::sqrt(dx * dx + dy * dy);
    };

    std::vector<EndpointRef> endpoints;
    endpoints.reserve(lines.size() * 2);
    for (int i = 0; i < (int)lines.size(); ++i) {
        endpoints.push_back({ i, true, lines[i].start, 0 });
        endpoints.push_back({ i, false, lines[i].end, 0 });
    }

    for (size_t i = 0; i < endpoints.size(); ++i) {
        for (size_t j = i + 1; j < endpoints.size(); ++j) {
            if (Distance2DLocal(endpoints[i].p, endpoints[j].p) < 1e-3) {
                endpoints[i].degree++;
                endpoints[j].degree++;
            }
        }
    }

    int healedCount = 0;
    std::vector<bool> endpointHealed(endpoints.size(), false);
    while (true) {
        int bestA = -1;
        int bestB = -1;
        double bestDist = healTol;

        for (int a = 0; a < (int)endpoints.size(); ++a) {
            if (endpointHealed[a] || endpoints[a].degree > 0) continue;
            for (int b = a + 1; b < (int)endpoints.size(); ++b) {
                if (endpointHealed[b] || endpoints[b].degree > 0) continue;
                if (endpoints[a].lineIndex == endpoints[b].lineIndex) continue;

                double dist = Distance2DLocal(endpoints[a].p, endpoints[b].p);
                if (dist <= bestDist) {
                    bestDist = dist;
                    bestA = a;
                    bestB = b;
                }
            }
        }

        if (bestA < 0 || bestB < 0) break;

        Point3D mergedPoint;
        mergedPoint.x = (endpoints[bestA].p.x + endpoints[bestB].p.x) * 0.5;
        mergedPoint.y = (endpoints[bestA].p.y + endpoints[bestB].p.y) * 0.5;
        mergedPoint.z = (endpoints[bestA].p.z + endpoints[bestB].p.z) * 0.5;

        auto ApplyEndpoint = [&](const EndpointRef& endpoint) {
            if (endpoint.isStart) {
                lines[endpoint.lineIndex].start = mergedPoint;
            }
            else {
                lines[endpoint.lineIndex].end = mergedPoint;
            }
        };

        ApplyEndpoint(endpoints[bestA]);
        ApplyEndpoint(endpoints[bestB]);
        endpointHealed[bestA] = true;
        endpointHealed[bestB] = true;
        healedCount++;
    }

    return healedCount;
}

// 核心成环算法：只找大环，自动忽略碎线
vector<Loop> BuildLoops(const vector<OneLine>& inputLines) {
    vector<Loop> loops;
    vector<OneLine> lines = DeduplicateLines(inputLines);
    int healedCount = HealOpenEndpointGaps(lines, 0.02);
    if (healedCount > 0) {
        cout << "  [BuildLoops] 端点断口愈合数量: " << healedCount
            << " (healTol=0.02mm)" << endl;
    }
    //debug 
    //std::string Name = savePath + "Slice_debug_line.brep";
    //ExportOneLinesToBrep(lines, Name);

    if (lines.empty()) return loops;

    vector<bool> visited(lines.size(), false);
    const double snapTol = 0.001;

    auto Snap = [](const Point3D& p, double tol) -> std::pair<int64_t, int64_t> {
        return { (int64_t)round(p.x / tol), (int64_t)round(p.y / tol) };
    };

    auto MakeSnapKey = [](int64_t x, int64_t y) -> int64_t {
        return x * 1000000007LL + y;
    };

    std::unordered_map<int64_t, std::vector<size_t>> endMap;
    for (size_t i = 0; i < lines.size(); ++i) {
        auto sk = Snap(lines[i].start, snapTol);
        auto ek = Snap(lines[i].end, snapTol);
        int64_t startKey = MakeSnapKey(sk.first, sk.second);
        int64_t endKey = MakeSnapKey(ek.first, ek.second);
        endMap[startKey].push_back(i);
        if (startKey != endKey) {
            endMap[endKey].push_back(i);
        }
    }

    for (size_t i = 0; i < lines.size(); ++i) {
        if (visited[i]) continue;

        Loop currentLoop;
        std::vector<size_t> trialIndices;
        std::vector<bool> trialUsed(lines.size(), false);

        currentLoop.push_back(lines[i]);
        trialIndices.push_back(i);
        trialUsed[i] = true;
        Point3D startPt = lines[i].start;
        Point3D currEnd = lines[i].end;

        while (true) {
            auto key = Snap(currEnd, snapTol);

            size_t bestIdx = SIZE_MAX;
            double bestDist = 1e9;
            bool bestReverse = false;

            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    int64_t mapKey = MakeSnapKey(key.first + dx, key.second + dy);
                    auto it = endMap.find(mapKey);
                    if (it == endMap.end()) continue;

                    for (size_t idx : it->second) {
                        if (visited[idx] || trialUsed[idx]) continue;
                        const auto& l = lines[idx];
                        double dS = fabs(currEnd.x - l.start.x) + fabs(currEnd.y - l.start.y);
                        double dE = fabs(currEnd.x - l.end.x) + fabs(currEnd.y - l.end.y);
                        if (dS < snapTol && dS < bestDist) {
                            bestDist = dS; bestIdx = idx; bestReverse = false;
                        }
                        if (dE < snapTol && dE < bestDist) {
                            bestDist = dE; bestIdx = idx; bestReverse = true;
                        }
                    }
                }
            }

            if (bestIdx == SIZE_MAX) break;

            trialUsed[bestIdx] = true;
            trialIndices.push_back(bestIdx);
            if (bestReverse) {
                OneLine rev = lines[bestIdx];
                swap(rev.start, rev.end);
                currentLoop.push_back(rev);
                currEnd = rev.end;
            } else {
                currentLoop.push_back(lines[bestIdx]);
                currEnd = lines[bestIdx].end;
            }

            if (IsPointEqual(currEnd, startPt, snapTol)) break;
        }

        if (IsPointEqual(currentLoop.back().end, currentLoop.front().start, snapTol)
            && currentLoop.size() >= 3)
        {
            for (size_t idx : trialIndices) {
                visited[idx] = true;
            }
            loops.push_back(currentLoop);
        }
    }

    cout << "\n================ 最终结果 =================" << endl;
    cout << " 成功生成环数量：" << loops.size() << endl;
    for (size_t i = 0; i < loops.size(); i++) {
        cout << "  环 " << i + 1 << "：" << loops[i].size() << " 条线段" << endl;
    }
    int usedCount = 0;
    for (bool v : visited) if (v) usedCount++;
    cout << "  输入线段: " << lines.size() << " 已用: " << usedCount << " 未用: " << lines.size() - usedCount << endl;
    if (usedCount < (int)lines.size()) {
        std::vector<OneLine> unusedLines;
        for (size_t i = 0; i < lines.size(); ++i) {
            if (!visited[i]) {
                unusedLines.push_back(lines[i]);
            }
        }
        std::string unusedName = savePath + "BuildLoops_UnusedLines_" + std::to_string((int)lines.size())
            + "_unused_" + std::to_string((int)unusedLines.size()) + ".brep";
        ExportOneLinesToBrep(unusedLines, unusedName);
        cout << "  [BuildLoops] 未用线段已导出: " << unusedName << endl;
        PrintEndpointGapDiagnostics(unusedLines, snapTol);
    }
    cout << "==========================================\n" << endl;

    return loops;
}



// 计算多边形面积（鞋带公式 Shoelace Formula）

// --- 拓扑构建 ---
void ExtractFacesFromTree(LoopNode* node, int depth, std::vector<Face2D>& faces) {
    // 根节点（depth == 0）是虚拟节点，没有自身的线段，直接跳过其自身的构造
    if (depth > 0) {
        Face2D newFace;
        // 奇数层代表实体，偶数层代表空腔
        newFace.type = (depth % 2 != 0) ? FaceType::SOLID : FaceType::CAVITY;
        newFace.outerLoop = node->loopLines;
        
        // 它的直接子节点就是它的内环
        for (LoopNode* child : node->children) {
            newFace.innerLoops.push_back(child->loopLines);
        }
        faces.push_back(newFace);
    }
    
    // 继续递归遍历子节点
    for (LoopNode* child : node->children) {
        ExtractFacesFromTree(child, depth + 1, faces);
    }
}

// 对每个环可视化导出为独立的 BREP 文件
void DebugExportLoops(const std::vector<Loop>& loops, const std::string& prefix, const std::string& savePath) {
    for (size_t i = 0; i < loops.size(); ++i) {
        BRep_Builder builder;
        TopoDS_Compound comp;
        builder.MakeCompound(comp);
        
        for (const auto& line : loops[i]) {
            gp_Pnt p1(line.start.x, line.start.y, line.start.z);
            gp_Pnt p2(line.end.x, line.end.y, line.end.z);
            if (!p1.IsEqual(p2, 1e-7)) {
                TopoDS_Edge anEdge = BRepBuilderAPI_MakeEdge(p1, p2);
                builder.Add(comp, anEdge);
            }
        }
        
        double area = CalculateArea(loops[i]);
        std::string fileName = savePath + prefix + "_Loop_" + std::to_string(i) + "_Area" + std::to_string((int)area) + ".brep";
        BRepTools::Write(comp, fileName.c_str());
        std::cout << "  环 " << i << " 面积=" << area << " -> " << fileName << std::endl;
    }
}

// 构建嵌套树并提取 Face2D
std::vector<Face2D> BuildTopologyAndExtractFaces(const std::vector<OneLine>& layerLines) {
    // 1. 组装成环
    //PrintAllLines(layerLines);
    std::vector<Loop> loops = BuildLoops(layerLines);
    
    // 调试：导出每个环
    //DebugExportLoops(loops, "DebugLayer", savePath);
    
    // 2. 创建树节点并计算面积
    std::vector<LoopNode*> nodes;
    for (const auto& loop : loops) {
        LoopNode* node = new LoopNode();
        node->loopLines = loop;
        node->area = CalculateArea(loop);
        nodes.push_back(node);
    }
    
    // 3. 按面积从大到小排序
    std::sort(nodes.begin(), nodes.end(), [](LoopNode* a, LoopNode* b) {
        return a->area > b->area;
    });
    
    // 4. 构建嵌套树
    LoopNode root; // 虚拟根节点（第 0 层）

    for (size_t i = 0; i < nodes.size(); ++i) {
        // 取当前环的重心作为测试点（比只取第一个点更稳健）
        Point3D center = {0.0, 0.0, 0.0};
        int ptCount = 0;
        for (const auto& l : nodes[i]->loopLines) {
            center.x += l.start.x;
            center.y += l.start.y;
            center.z += l.start.z;
            ptCount++;
        }
        if (ptCount > 0) {
            center.x /= ptCount;
            center.y /= ptCount;
            center.z /= ptCount;
        }

        bool foundParent = false;
        // 从比它大一点的环开始往上找（即从 i-1 倒序遍历到 0）
        // 找到的第一个包含它的环，就是它的直接父节点（面积最小的包围者）
        for (int j = (int)i - 1; j >= 0; --j) {
            // 用重心测试：重心在多边形内部的概率远高于任意一个顶点
            if (IsPointInLoop(center, nodes[j]->loopLines)) {
                nodes[j]->children.push_back(nodes[i]);
                foundParent = true;
                break;
            }
        }
        // 如果没有环包含它，它就是最外层，挂在虚拟根节点下
        if (!foundParent) {
            root.children.push_back(nodes[i]);
        }
    }
    
    // 5. 提取 Face2D (全空间剖分)
    std::vector<Face2D> finalFaces;
    // 从虚拟根节点开始提取，根节点深度为0
    ExtractFacesFromTree(&root, 0, finalFaces); 
    
    // 清理内存
    for (LoopNode* node : nodes) {
        delete node;
    }
    
    return finalFaces;
}

// 将提取出的 Face2D 转换回 OCCT 的 TopoDS_Face 并保存为 BREP
