#include "boolean_ops.h"


// 将自定义的 Loop 转为 Clipper 的 Path64
Path64 LoopToPath(const Loop& loop) {
    Path64 path;
    for (const auto& line : loop) {
        path.push_back(Point64(line.start.x * SCALE, line.start.y * SCALE));
    }
    return path;
}

// 将 Clipper 的 Path64 还原为 Loop
Loop PathToLoop(const Path64& path, double zHeight) {
    Loop loop;
    if (path.empty()) return loop;

    for (size_t i = 0; i < path.size(); ++i) {
        size_t next_i = (i + 1) % path.size();
        OneLine line;
        line.start.x = static_cast<double>(path[i].x) / SCALE;
        line.start.y = static_cast<double>(path[i].y) / SCALE;
        line.start.z = zHeight;

        line.end.x = static_cast<double>(path[next_i].x) / SCALE;
        line.end.y = static_cast<double>(path[next_i].y) / SCALE;
        line.end.z = zHeight;

        line.faceId = -1; // 布尔运算后的新线段丢失了原始面ID
        loop.push_back(line);
    }
    return loop;
}

// 将 Face2D 转换为 Clipper 的输入路径集合
Paths64 FacesToClipperPaths(const std::vector<Face2D>& faces) {
    Paths64 result;
    for (const auto& face : faces) {
        // 只提取实体面参与运算
        if (face.outerLoop.empty()) continue;

        Path64 outer = LoopToPath(face.outerLoop);
        if (!IsPositive(outer)) std::reverse(outer.begin(), outer.end()); // 外环必须正向
        result.push_back(outer);

        for (const auto& inner : face.innerLoops) {
            Path64 innerPath = LoopToPath(inner);
            if (IsPositive(innerPath)) std::reverse(innerPath.begin(), innerPath.end()); // 内环必须反向
            result.push_back(innerPath);
        }
    }
    return result;
}

// ================== 3. 结果还原：从 PolyTree 重建 Face2D ==================

// 递归遍历 PolyTree，还原出全空间剖分的 Face2D（和之前的手写树逻辑完全一致！）
void ExtractFacesFromPolyNode(const PolyPath64* node, int depth, double zHeight, std::vector<Face2D>& faces) {
    if (depth > 0) {
        // 在布尔运算的结果中，奇数层(depth 1, 3, 5...)代表正向实体区域
        // 偶数层(depth 2, 4, 6...)代表孔洞区域
        // 因为 Face2D 结构已经包含了 innerLoops 来存储孔洞
        // 所以我们只需要在奇数层创建 Face2D，并将其直接子节点作为孔洞加入
        if (depth % 2 != 0) {
            Face2D newFace;
            newFace.type = FaceType::SOLID; // 提取出来的都是正向实体
            newFace.outerLoop = PathToLoop(node->Polygon(), zHeight);

            // 它的直接子节点就是它的孔洞 (depth + 1, 为偶数层)
            for (size_t i = 0; i < node->Count(); ++i) {
                newFace.innerLoops.push_back(PathToLoop(node->Child(i)->Polygon(), zHeight));
            }
            faces.push_back(newFace);
        }
    }

    // 递归处理子节点
    for (size_t i = 0; i < node->Count(); ++i) {
        ExtractFacesFromPolyNode(node->Child(i), depth + 1, zHeight, faces);
    }
}

// ================== 4. 终极封装：面面布尔运算 ==================

/**
 * 执行单个面与单个面的布尔运算（避免多个 clips 路径方向相互抵消）
 * @param subject 目标面 (A)
 * @param clip 裁剪面 (B)
 * @param clipType 运算类型：ClipType::Intersection(交集), Union(并集), Difference(A-B)
 * @param zHeight 输出面的 Z 高度
 */
std::vector<Face2D> BooleanFacesSingle(
    const Face2D& subject,
    const Face2D& clip,
    ClipType clipType,
    double zHeight)
{
    Clipper64 clipper;

    // 1. 添加单个运算主体 A 和单个裁剪面 B
    std::vector<Face2D> subjects = { subject };
    std::vector<Face2D> clips = { clip };
    clipper.AddSubject(FacesToClipperPaths(subjects));
    clipper.AddClip(FacesToClipperPaths(clips));

    // 2. 执行布尔运算，要求输出保留树状结构的 PolyTree
    PolyTree64 solutionTree;
    // 使用 NonZero 填充规则，完美识别由于方向导致的实体与孔洞
    clipper.Execute(clipType, FillRule::NonZero, solutionTree);

    // 3. 从 PolyTree 完美重建所有的 Face2D
    std::vector<Face2D> resultFaces;
    ExtractFacesFromPolyNode(&solutionTree, 0, zHeight, resultFaces);

    return resultFaces;
}

/**
 * 执行面与面的布尔运算（支持多个 subjects 和多个 clips）
 * @param subjects 目标面集合 (A)
 * @param clips 裁剪面集合 (B)
 * @param clipType 运算类型：ClipType::Intersection(交集), Union(并集), Difference(A-B)
 * @param zHeight 输出面的 Z 高度
 */
std::vector<Face2D> BooleanFaces(
    const std::vector<Face2D>& subjects,
    const std::vector<Face2D>& clips,
    ClipType clipType,
    double zHeight)
{
    // 如果 clips 只有一个 Face2D，直接使用原始方法
    if (clips.size() <= 1) {
        Clipper64 clipper;
        clipper.AddSubject(FacesToClipperPaths(subjects));
        clipper.AddClip(FacesToClipperPaths(clips));
        
        PolyTree64 solutionTree;
        clipper.Execute(clipType, FillRule::NonZero, solutionTree);
        
        std::vector<Face2D> resultFaces;
        ExtractFacesFromPolyNode(&solutionTree, 0, zHeight, resultFaces);
        return resultFaces;
    }
    
    // 如果 clips 有多个 Face2D，逐个处理，避免路径方向相互抵消
    std::vector<Face2D> result = subjects;
    for (const auto& clip : clips) {
        std::vector<Face2D> newResult;
        for (const auto& subject : result) {
            std::vector<Face2D> diffResult = BooleanFacesSingle(subject, clip, clipType, zHeight);
            newResult.insert(newResult.end(), diffResult.begin(), diffResult.end());
        }
        result = newResult;
    }
    return result;
}

/**
 * @brief 【特征追溯引擎】空间比对碰撞法：找回被 Clipper2 弄丢的原始面 ID
 * @param pocketResults Clipper2 布尔运算直接吐出来的原始型腔面集合（此时 faceId 均为 -1）
 * @param currentSolidFaces 这一层最原始的工件实体面集合（内部保留了 STEP 模型的真实 faceId）
 * @return 恢复了原始 faceId（硬边界）且将空气边严格标记为 -1（软边界）的全新 Face2D 型腔集合
 */
std::vector<Face2D> RecoverOriginalFaceIdsByGeometry(
    const std::vector<Face2D>& pocketResults,
    const std::vector<Face2D>& currentSolidFaces)
{
    std::vector<Face2D> recoveredCavities = pocketResults; // 拷贝一份准备改写
    const double match_tol = 0.05; // 50微米空间碰撞重合容差

    // 1. 遍历每一个型腔区域
    for (auto& cavityFace : recoveredCavities) {

        // 2. 盘查当前型腔区域的外环线段
        for (auto& cLine : cavityFace.outerLoop) {
            gp_Pnt cStart(cLine.start.x, cLine.start.y, cLine.start.z);
            gp_Pnt cEnd(cLine.end.x, cLine.end.y, cLine.end.z);

            bool isMatched = false;

            // 3. 去最原始的实体面阵营里进行空间高精碰撞
            for (const auto& sFace : currentSolidFaces) {
                for (const auto& sLine : sFace.outerLoop) {
                    gp_Pnt sStart(sLine.start.x, sLine.start.y, sLine.start.z);
                    gp_Pnt sEnd(sLine.end.x, sLine.end.y, sLine.end.z);

                    // 🎯 空间几何对碰成功（支持正向重合或反向首尾重合）
                    if ((cStart.Distance(sStart) < match_tol && cEnd.Distance(sEnd) < match_tol) ||
                        (cStart.Distance(sEnd) < match_tol && cEnd.Distance(sStart) < match_tol))
                    {
                        // 🟢 【黄金继承】：将型腔这根线的 ID，完美恢复成它亲生父母在 STEP 里的原始面 ID！
                        cLine.faceId = sLine.faceId;
                        isMatched = true;
                        break;
                    }
                }
                if (isMatched) break;
            }

            // 4. 🔴 如果遍历了所有的实体边界都碰不上，铁证如山：它就是悬空的【虚拟软边界】
            if (!isMatched) {
                cLine.faceId = -1; // 强制给它盖章为软边界标记
            }
        }

        // 5. 同样的逻辑，顺手盘查可能存在的内孔边界（如果有的话）
        for (auto& innerLoop : cavityFace.innerLoops) {
            for (auto& cLine : innerLoop) {
                gp_Pnt cStart(cLine.start.x, cLine.start.y, cLine.start.z);
                gp_Pnt cEnd(cLine.end.x, cLine.end.y, cLine.end.z);
                bool isMatched = false;

                for (const auto& sFace : currentSolidFaces) {
                    for (const auto& sLine : sFace.outerLoop) {
                        gp_Pnt sStart(sLine.start.x, sLine.start.y, sLine.start.z);
                        gp_Pnt sEnd(sLine.end.x, sLine.end.y, sLine.end.z);

                        if ((cStart.Distance(sStart) < match_tol && cEnd.Distance(sEnd) < match_tol) ||
                            (cStart.Distance(sEnd) < match_tol && cEnd.Distance(sStart) < match_tol))
                        {
                            cLine.faceId = sLine.faceId;
                            isMatched = true;
                            break;
                        }
                    }
                    if (isMatched) break;
                }
                if (!isMatched) {
                    cLine.faceId = -1;
                }
            }
        }
    }

    return recoveredCavities;
}

/**
 * @brief 【调试专用】将型腔集合中所有 faceId == -1 的虚拟软边界单独导出为 BREP 线框
 * @param cavityFaces 已经过 ID 恢复（RecoverOriginalFaceIdsByGeometry）的型腔面集合
 * @param filePath 导出的文件路径（例如 "D:/debug_soft_edges.brep"）
 */
