#include "closed_cavity.h"
#include "slicing.h"
#include "occ_utils.h"
#include "geom_utils.h"
#include "export_utils.h"
#include "open_cavity.h"

// --- 核心递归切分 ---
void SplitCavityRecursive(const std::set<int>& currentFaceIds,
    TopTools_DataMapOfIntegerShape& idToFace,
    int& nextId,
    const std::vector<double>& splitZs,
    int zIndex,
    std::vector<std::set<int>>& result) {

    // 终止条件：没有切分点或没有待处理的面
    if (zIndex >= splitZs.size() || currentFaceIds.empty()) {
        if (!currentFaceIds.empty()) result.push_back(currentFaceIds);
        return;
    }

    double curSplitZ = splitZs[zIndex];
    std::set<int> aboveIds, belowIds;
    std::set<int> facesToSplit;

    // --- 1. 初步筛选：区分哪些面需要被“切一刀” ---
    for (int id : currentFaceIds) {
        TopoDS_Face face = TopoDS::Face(idToFace.Find(id));

        // 水平面分类
        if (IsHorizontalFace(face)) {
            double z = GetHorizontalFaceZ(face);

            // 如果面不在当前切分高度，按常规高度判断
            if (std::abs(z - curSplitZ) > 1e-4) {
                if (z > curSplitZ) aboveIds.insert(id);
                else belowIds.insert(id);
                continue;
            }

            // --- 关键修改：处理刚好在切分点处的水平面 ---
            gp_Dir normal = GetFaceNormal(face);
            double dot = normal.Z(); // 法线与 [0,0,1] 的点积

            if (dot < -0.5) {
                // 1. 法线向下 [0,0,-1]，与进刀方向相同
                // 它是下一层的“盖子”，传给下方
                belowIds.insert(id);
                cout << "  [水平面判定] ID:" << id << " 法线向下 -> 归入下方 (盖子)" << endl;
            }
            else {
                // 2. 法线向上 [0,0,1]，与进刀方向相反
                // 它是这一层的“底面”，归入上方
                aboveIds.insert(id);
                cout << "  [水平面判定] ID:" << id << " 法线向上 -> 归入上方 (底面)" << endl;
            }
            continue;
        }

        // 侧面根据包围盒分类
        double zmin, zmax;
        GetFaceZRange(face, zmin, zmax);

        if (zmin >= curSplitZ - 1e-4) {
            aboveIds.insert(id);
        }
        else if (zmax <= curSplitZ + 1e-4) {
            belowIds.insert(id);
        }
        else {
            facesToSplit.insert(id); // 跨越切分点的侧面
        }
    }

    // --- 2. 构造切分工具：无限大平面 ---
    gp_Pln cuttingPlane(gp_Pnt(0, 0, curSplitZ), gp_Dir(0, 0, 1));
    TopoDS_Face planeFace = BRepBuilderAPI_MakeFace(cuttingPlane);

    // --- 3. 核心切割循环 ---
    for (int id : facesToSplit) {
        TopoDS_Face faceToCut = TopoDS::Face(idToFace.Find(id));

        // A. 求取交线 (Section)
        BRepAlgoAPI_Section section(faceToCut, planeFace);
        section.Build();

        if (!section.IsDone() || !TopExp_Explorer(section.Shape(), TopAbs_EDGE).More()) {
            belowIds.insert(id); // 异常情况：归入下方
            continue;
        }

        // B. 使用 SplitShape 物理分割面
        BRepFeat_SplitShape splitter(faceToCut);
        TopExp_Explorer edgeExp(section.Shape(), TopAbs_EDGE);
        for (; edgeExp.More(); edgeExp.Next()) {
            splitter.Add(TopoDS::Edge(edgeExp.Current()), faceToCut);
        }
        splitter.Build();

        if (splitter.IsDone()) {
            // C. 收集并分类新产生的所有“子面”
            BRep_Builder debugBuilder;
            TopoDS_Compound splitResultComp;
            debugBuilder.MakeCompound(splitResultComp);

            TopExp_Explorer faceExp(splitter.Shape(), TopAbs_FACE);
            for (; faceExp.More(); faceExp.Next()) {
                TopoDS_Face subFace = TopoDS::Face(faceExp.Current());

                double smin, smax;
                GetFaceZRange(subFace, smin, smax);
                double zMid = (smin + smax) / 2.0; // 用中点判定归属最稳健

                int newId = nextId++;
                idToFace.Bind(newId, subFace);
                debugBuilder.Add(splitResultComp, subFace);

                if (zMid > curSplitZ) {
                    aboveIds.insert(newId);
                }
                else {
                    belowIds.insert(newId);
                }
            }

            // D. 可视化导出：这一层当前面被切开后的结果
            std::string debugPath = savePath + "Split_Z" + std::to_string((int)curSplitZ) + "_ID" + std::to_string(id) + ".brep";
            BRepTools::Write(splitResultComp, debugPath.c_str());
        }
        else {
            belowIds.insert(id); // 切割失败保留原始面
        }
    }

    // --- 4. 分类结果可视化 (Above vs Below) ---
    if (!aboveIds.empty() || !belowIds.empty()) {
        BRep_Builder resultBuilder;

        // 导出“上方集合” (当前层的结果)
        TopoDS_Compound compAbove;
        resultBuilder.MakeCompound(compAbove);
        for (int id : aboveIds) {
            resultBuilder.Add(compAbove, idToFace.Find(id));
        }
        //std::string abovePath = savePath + "Layer_Z" + std::to_string((int)curSplitZ) + "_ABOVE_FINAL.brep";
        //BRepTools::Write(compAbove, abovePath.c_str());

        // 导出“下方集合” (留给下一层递归的形状)
        TopoDS_Compound compBelow;
        resultBuilder.MakeCompound(compBelow);
        for (int id : belowIds) {
            resultBuilder.Add(compBelow, idToFace.Find(id));
        }
        //std::string belowPath = savePath + "Layer_Z" + std::to_string((int)curSplitZ) + "_BELOW_REMAINING.brep";
        //BRepTools::Write(compBelow, belowPath.c_str());

        //cout << "  [Visualize] Z=" << curSplitZ << " 层提取完成: " << endl;
        //cout << "    - 上方导出: " << abovePath << " (" << aboveIds.size() << " 个面)" << endl;
        //cout << "    - 下方导出: " << belowPath << " (" << belowIds.size() << " 个面)" << endl;
    }

    // --- 4. 结果收集与递归 ---
    // 将当前层（切分线以上的所有面）存入结果
    if (!aboveIds.empty()) {
        result.push_back(aboveIds);
    }

    // 将切分线以下的所有面交给下一个切分点递归处理
    SplitCavityRecursive(belowIds, idToFace, nextId, splitZs, zIndex + 1, result);
}

// 入口函数
// 将返回值改为 Compound 数组
// ============================================================
// 入口函数：针对单个独立型腔的自适应切分（智能台阶过滤版）
// ============================================================
std::vector<TopoDS_Compound> SplitCavity(const TopoDS_Compound& cavity) {
    TopTools_DataMapOfIntegerShape idToFace;
    std::set<int> allIds;
    int nextId = 1;

    // 1. 初始化 ID 到面的映射
    TopExp_Explorer exp(cavity, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        int id = nextId++;
        idToFace.Bind(id, TopoDS::Face(exp.Current()));
        allIds.insert(id);
    }

    // ==============================================================
    // 🌟 新增：建立 Edge -> Face 的拓扑邻接图，用于打探邻居面的情报
    // ==============================================================
    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(cavity, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);

    // 2. 获取并排序切分点
    std::set<double> zSet;
    for (int id : allIds) {
        TopoDS_Face f = TopoDS::Face(idToFace.Find(id));

        // 我们只考察水平平面是否能成为切分点
        if (IsHorizontalFace(f)) {
            double zFace = GetHorizontalFaceZ(f);

            // 默认假设它是真正的底面（不论是封闭底还是开放悬空底），默认不切！
            bool isIntermediateStep = false;

            // 扫荡这个平面的所有边界边，看看邻居是往上长还是往下走
            TopExp_Explorer edgeExp(f, TopAbs_EDGE);
            for (; edgeExp.More(); edgeExp.Next()) {
                const TopoDS_Edge& edge = TopoDS::Edge(edgeExp.Current());

                // 如果这条边有相邻的面
                if (edgeToFaces.Contains(edge)) {
                    const TopTools_ListOfShape& neighbors = edgeToFaces.FindFromKey(edge);
                    TopTools_ListIteratorOfListOfShape it(neighbors);
                    for (; it.More(); it.Next()) {
                        TopoDS_Face nFace = TopoDS::Face(it.Value());

                        // 排除自己，只看真正的邻居面
                        if (!nFace.IsSame(f)) {
                            // 探测邻居面在 Z 轴上的极限范围
                            double nZmin, nZmax;
                            GetFaceZRange(nFace, nZmin, nZmax);

                            // 🎯 核心判据：如果发现任何一个邻居面比我还低（容差 1e-3 防止浮点抖动）
                            // 铁证如山：我下方还有深坑（比如小腔体侧壁）或者外围有更低的侧壁！
                            // 我绝不是底面，我是一个“中间台阶面 (Step)”，必须在这里切一刀！
                            if (nZmin < zFace - 1e-3) {
                                isIntermediateStep = true;
                                break;
                            }
                        }
                    }
                }
                // 只要发现一个向下的邻居，就已经实锤是台阶了，立刻停止排查这层面的其他边，提升效率
                if (isIntermediateStep) break;
            }

            // 判决时刻：只有被判定为“中间台阶”的平面，才允许登记为切分点
            if (isIntermediateStep) {
                zSet.insert(zFace);
            }
            else {
                std::cout << "  [智能拦截] 侦测到水平面 (Z=" << zFace
                    << ") 所有邻面均朝上或为悬空边界，确认为极限大底面，已免于切分！" << std::endl;
            }
        }
    }

    std::vector<double> splitZs(zSet.begin(), zSet.end());
    // 从大到小降序排列（从型腔最顶端往下排）
    std::sort(splitZs.begin(), splitZs.end(), std::greater<double>());

    // 【已被移除】：if (!splitZs.empty()) splitZs.pop_back(); 
    // 原因：上面的智能拦截网已经把真正的闭合底面和开放底面全部保护起来了，
    // zSet 里剩下的【全都是必须切分的中间台阶】，如果硬弹出一个，反而会漏切！

    // 3. 带着极其精准的切分高度列表，执行递归切分
    std::vector<std::set<int>> idGroups; // 存储各组切分后的面 ID
    SplitCavityRecursive(allIds, idToFace, nextId, splitZs, 0, idGroups);

    // 4. 【关键步骤】在 idToFace 销毁前，将打散的 ID 转换回真正的 3D 几何壳体
    std::vector<TopoDS_Compound> finalShapes;
    BRep_Builder builder;

    for (const auto& group : idGroups) {
        TopoDS_Compound comp;
        builder.MakeCompound(comp);
        for (int id : group) {
            if (idToFace.IsBound(id)) {
                builder.Add(comp, idToFace.Find(id));
            }
        }
        finalShapes.push_back(comp);
    }

    // 返回由物理切分产生的一组实实在在的加工层级壳体
    return finalShapes;
}
// 将一个包含多个独立区域的 Compound 拆分为多个独立的 Compound
std::vector<TopoDS_Compound> SeparateDisconnectedCavities(const TopoDS_Compound& cavityCompound) {
    std::vector<TopoDS_Compound> individualCavities;

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopExp::MapShapesAndAncestors(cavityCompound, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

    // 额外建立 Vertex -> Face 的映射，处理"只共点不共边"的情况
    TopTools_IndexedDataMapOfShapeListOfShape vertexFaceMap;
    TopExp::MapShapesAndAncestors(cavityCompound, TopAbs_VERTEX, TopAbs_FACE, vertexFaceMap);

    TopTools_MapOfShape visitedFaces;
    TopExp_Explorer exp(cavityCompound, TopAbs_FACE);

    for (; exp.More(); exp.Next()) {
        TopoDS_Face startFace = TopoDS::Face(exp.Current());
        if (visitedFaces.Contains(startFace)) continue;

        BRep_Builder builder;
        TopoDS_Compound singleCavity;
        builder.MakeCompound(singleCavity);

        TopTools_ListOfShape queue;
        queue.Append(startFace);
        visitedFaces.Add(startFace);

        while (!queue.IsEmpty()) {
            TopoDS_Face currentFace = TopoDS::Face(queue.First());
            queue.RemoveFirst();
            builder.Add(singleCavity, currentFace);

            // 通过共享边查找邻接面
            TopExp_Explorer edgeExp(currentFace, TopAbs_EDGE);
            for (; edgeExp.More(); edgeExp.Next()) {
                const TopoDS_Shape& edge = edgeExp.Current();
                if (edgeFaceMap.Contains(edge)) {
                    const TopTools_ListOfShape& adjFaces = edgeFaceMap.FindFromKey(edge);
                    for (TopTools_ListIteratorOfListOfShape it(adjFaces); it.More(); it.Next()) {
                        const TopoDS_Shape& adjFace = it.Value();
                        if (!visitedFaces.Contains(adjFace)) {
                            visitedFaces.Add(adjFace);
                            queue.Append(adjFace);
                        }
                    }
                }
            }

            // 通过共享顶点查找邻接面（处理切分后只共点不共边的情况）
            TopExp_Explorer vertExp(currentFace, TopAbs_VERTEX);
            for (; vertExp.More(); vertExp.Next()) {
                const TopoDS_Shape& vertex = vertExp.Current();
                if (vertexFaceMap.Contains(vertex)) {
                    const TopTools_ListOfShape& adjFaces = vertexFaceMap.FindFromKey(vertex);
                    for (TopTools_ListIteratorOfListOfShape it(adjFaces); it.More(); it.Next()) {
                        const TopoDS_Shape& adjFace = it.Value();
                        if (!visitedFaces.Contains(adjFace)) {
                            visitedFaces.Add(adjFace);
                            queue.Append(adjFace);
                        }
                    }
                }
            }
        }
        individualCavities.push_back(singleCavity);
    }
    return individualCavities;
}


/**
 * 功能：将物理分离的型腔实体与所有切片的 2D 环进行匹配，重组为 CavityFeature 列表
 * 修正：bottomZ 自动向匹配下方最近的一个总切分点对齐
 */
std::vector<CavityFeature> GenerateCavityFeatures(
    const std::vector<TopoDS_Compound>& isolatedCavities,
    const std::vector<std::vector<Face2D>>& allLayerFaces,
    const std::vector<double>& sortedSplitPoints, // 传入降序排列的总切分点
    const TopTools_DataMapOfShapeInteger& faceToIdMap)
{
    std::vector<CavityFeature> finalFeatures;

    for (size_t i = 0; i < isolatedCavities.size(); ++i) {
        CavityFeature feat;
        feat.featureId = (int)i;
        feat.type = CavityType::CLOSED;
        feat.topZ = -1e9;   // 初始化极小值
        feat.bottomZ = 1e9; // 初始化极大值，稍后修正为下一个切分点
        double lastLoopZ = 1e9; // 记录探测到的物理环的最低 Z

        // --- 步骤 A: 收集当前型腔实体的面 ID 集合 ---
        TopExp_Explorer exp(isolatedCavities[i], TopAbs_FACE);
        for (; exp.More(); exp.Next()) {
            if (faceToIdMap.IsBound(exp.Current())) {
                feat.sourceFaceIds.insert(faceToIdMap.Find(exp.Current()));
            }
        }

        // --- 步骤 B: 空间 XY 包围盒计算 ---
        Bnd_Box box;
        BRepBndLib::Add(isolatedCavities[i], box);
        double x1, y1, zmin_box, x2, y2, zmax_box;
        box.Get(x1, y1, zmin_box, x2, y2, zmax_box);

        // --- 步骤 C: 全层扫描与匹配 ---
        for (const auto& layer : allLayerFaces) {
            for (const auto& f2d : layer) {
                // 1. 类型过滤：只看空腔类型的切片
                if (f2d.type != FaceType::CAVITY || f2d.outerLoop.empty()) continue;

                // 2. ID 匹配：投票法统计外环所有线段的 faceId，取出现频率最高的
                //    比只看第一个线段更稳健，避免被虚拟软边界(faceId==-1)干扰
                std::map<int, int> idVoteCount;
                for (const auto& l : f2d.outerLoop) {
                    if (l.faceId >= 0) {
                        idVoteCount[l.faceId]++;
                    }
                }
                int dominantFaceId = -1;
                int maxVotes = 0;
                for (const auto& [fid, cnt] : idVoteCount) {
                    if (cnt > maxVotes) {
                        maxVotes = cnt;
                        dominantFaceId = fid;
                    }
                }

                if (dominantFaceId >= 0 && feat.sourceFaceIds.count(dominantFaceId)) {

                    // 3. 空间匹配：检查切片环的点是否在型腔实体的 XY 范围内
                    Point3D p = f2d.outerLoop[0].start;
                    if (p.x >= x1 - 0.5 && p.x <= x2 + 0.5 && p.y >= y1 - 0.5 && p.y <= y2 + 0.5) {

                        double curZ = p.z;

                        // 转换外环
                        CavityLoop out;
                        out.lines = f2d.outerLoop;
                        out.zHeight = curZ;
                        out.isOuter = true;
                        feat.stepLoops.push_back(out);

                        // 转换内环（孔洞/孤岛）
                        for (const auto& inner : f2d.innerLoops) {
                            CavityLoop in;
                            in.lines = inner;
                            in.zHeight = curZ;
                            in.isOuter = false;
                            feat.stepLoops.push_back(in);
                        }

                        // 更新探测到的 Z 极限
                        feat.topZ = std::max(feat.topZ, curZ);
                        lastLoopZ = std::min(lastLoopZ, curZ);
                    }
                }
            }
        }

        // --- 步骤 D: 特征收尾与底面修正 ---
        if (!feat.stepLoops.empty()) {
            // 修正 bottomZ：寻找检测到的最低环 Z 对应的下一个切分点
            feat.bottomZ = lastLoopZ; // 默认值
            for (size_t k = 0; k < sortedSplitPoints.size(); ++k) {
                // 如果当前点与最低环高度重合（带容差）
                if (std::abs(sortedSplitPoints[k] - lastLoopZ) < 1e-4) {
                    // 如果存在更低的一个切分点，则那个点才是型腔的物理底面
                    if (k + 1 < sortedSplitPoints.size()) {
                        feat.bottomZ = sortedSplitPoints[k + 1];
                    }
                    break;
                }
            }

            feat.totalDepth = std::abs(feat.topZ - feat.bottomZ);

            // 按高度从高到低排序，模拟加工顺序
            std::sort(feat.stepLoops.begin(), feat.stepLoops.end(),
                [](const CavityLoop& a, const CavityLoop& b) { return a.zHeight > b.zHeight; });

            finalFeatures.push_back(feat);
        }
    }

    // ==================== 详细信息打印开始 ====================

    for (const auto& feat : finalFeatures) {
        cout << ">>> [型腔 ID: " << feat.featureId << "]" << endl;
        cout << "    - 空间范围: Z_Top = " << feat.topZ << " / Z_Floor = " << feat.bottomZ << endl;
        cout << "    - 加工深度: " << feat.totalDepth << endl;
        cout << "    - 进刀方向: [0, 0, -1]" << endl;

        // 打印原始面 ID 清单
        cout << "    - 关联原始面 ID (" << feat.sourceFaceIds.size() << "个): ";
        for (int oid : feat.sourceFaceIds) cout << oid << " ";
        cout << endl;

        // 打印层级环信息
        cout << "    - 层级切片环清单 (按 Z 降序):" << endl;

        double lastZ = 999999.9;
        for (size_t k = 0; k < feat.stepLoops.size(); ++k) {
            const auto& loop = feat.stepLoops[k];

            if (fabs(loop.zHeight - lastZ) > 1e-4) {
                cout << "      [高度 Z = " << loop.zHeight << "]" << endl;
                lastZ = loop.zHeight;
            }

            string typeStr = loop.isOuter ? "外边界 (Outer)" : "孤岛/内孔 (Inner)";
            cout << "        * " << typeStr << ": 线段数 " << loop.lines.size();

            if (!loop.lines.empty()) {
                const auto& startPt = loop.lines.front().start;
                const auto& endPt = loop.lines.back().end;
                cout << " | 闭合性检校: 起(" << startPt.x << "," << startPt.y
                    << ") -> 终(" << endPt.x << "," << endPt.y << ")";
            }
        }
    }

    return finalFeatures;
}


// 根据识别到的型腔特征，从原始模型中提取真正属于型腔的面，形成新的 Compound
TopoDS_Compound ExtractTrueCavityFaces(
    const CavityFeature& feat,
    const TopoDS_Shape& mainShape,
    const TopTools_DataMapOfShapeInteger& faceToIdMap)
{
    BRep_Builder builder;
    TopoDS_Compound refinedResult;
    builder.MakeCompound(refinedResult);

    // 1. 准备“开口环（最上方环）”用于子面的水平位置校验
    TopoDS_Wire topWire;
    if (!feat.stepLoops.empty()) {
        const auto& tLoop = feat.stepLoops.front();
        BRepBuilderAPI_MakePolygon poly;
        for (const auto& l : tLoop.lines) poly.Add(gp_Pnt(l.start.x, l.start.y, l.start.z));
        const auto& lastL = tLoop.lines.back();
        poly.Add(gp_Pnt(lastL.end.x, lastL.end.y, lastL.end.z));
        if (poly.IsDone()) topWire = poly.Wire();
    }

    if (topWire.IsNull()) return refinedResult;

    TopExp_Explorer exp(mainShape, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Face originalFace = TopoDS::Face(exp.Current());
        if (!faceToIdMap.IsBound(originalFace)) continue;

        int oid = faceToIdMap.Find(originalFace);
        // 只处理属于该型腔特征关联的面
        if (feat.sourceFaceIds.count(oid)) {

            double zmin, zmax;
            GetFaceZRange(originalFace, zmin, zmax);

            // --- 策略 A：处理跨越切分点(topZ)的面 ---
            if (zmax > feat.topZ + 1e-4 && zmin < feat.topZ - 1e-4) {

                gp_Pln splitPln(gp_Pnt(0, 0, feat.topZ), gp_Dir(0, 0, 1));
                TopoDS_Face algoPlane = BRepBuilderAPI_MakeFace(splitPln);
                BRepAlgoAPI_Section section(originalFace, algoPlane);
                section.Build();

                if (section.IsDone()) {
                    BRepFeat_SplitShape splitter(originalFace);
                    TopExp_Explorer eExp(section.Shape(), TopAbs_EDGE);
                    for (; eExp.More(); eExp.Next())
                        splitter.Add(TopoDS::Edge(eExp.Current()), originalFace);
                    splitter.Build();

                    if (splitter.IsDone()) {
                        TopExp_Explorer subExp(splitter.Shape(), TopAbs_FACE);
                        for (; subExp.More(); subExp.Next()) {
                            TopoDS_Face subF = TopoDS::Face(subExp.Current());

                            double sMin, sMax;
                            GetFaceZRange(subF, sMin, sMax);

                            // 校验 1: 子面必须位于 topZ 以下
                            if (sMax <= feat.topZ + 1e-4) {
                                // 校验 2: 子面必须在水平方向上靠近型腔开口环
                                BRepExtrema_DistShapeShape distCalc(subF, topWire);
                                if (distCalc.IsDone() && distCalc.Value() < 1e-2) {
                                    builder.Add(refinedResult, subF);
                                }
                            }
                        }
                    }
                }
            }
            // --- 策略 B：处理完全在 topZ 以下的面 (不在待切分列表中的面) ---
            // 只要面在 topZ 以下且属于 sourceFaceIds，则全部直接加入
            else if (zmax <= feat.topZ + 1e-4) {
                builder.Add(refinedResult, originalFace);
                // 打印 Debug 信息，确认这包含水平面和其他非跨层侧面
                // cout << ">>>> [直接加入] 面 ID: " << oid << " 完全在高度范围内。" << endl;
            }
        }
    }

    return refinedResult;
}



// 用于面与面布尔运算的工业级引擎，完美支持交集、并集、差集，输出全空间剖分的 Face2D 结果

// --- ConvertPartsToFeatures ---
std::vector<CavityFeature> ConvertPartsToFeatures(
    const std::vector<TopoDS_Compound>& parts,
    const CavityFeature& parentFeat,
    const TopTools_DataMapOfShapeInteger& faceToIdMap)
{
    std::vector<CavityFeature> subFeatures;

    for (const auto& part : parts) {
        CavityFeature subFeat;
        subFeat.featureId = -1; // 建议在外部循环中赋值
        subFeat.type = parentFeat.type;
        subFeat.toolDirection = parentFeat.toolDirection;

        // 1. 获取物理边界 & ID 映射
        double subZmin = 1e9, subZmax = -1e9;
        TopExp_Explorer exp(part, TopAbs_FACE);
        for (; exp.More(); exp.Next()) {
            TopoDS_Face f = TopoDS::Face(exp.Current());
            double fzmin, fzmax;
            GetFaceZRange(f, fzmin, fzmax);
            subZmin = std::min(subZmin, fzmin);
            subZmax = std::max(subZmax, fzmax);

            if (faceToIdMap.IsBound(f)) {
                subFeat.sourceFaceIds.insert(faceToIdMap.Find(f));
            }
        }
        subFeat.topZ = subZmax;
        subFeat.bottomZ = subZmin;
        subFeat.totalDepth = std::abs(subFeat.topZ - subFeat.bottomZ);

        // 2. 几何嵌套过滤：从父特征中继承属于本零件的轮廓线段
        Bnd_Box subBox;
        BRepBndLib::Add(part, subBox);
        double bx1, by1, bz1, bx2, by2, bz2;
        subBox.Get(bx1, by1, bz1, bx2, by2, bz2);

        for (const auto& loop : parentFeat.stepLoops) {
            // 高度重合检测
            if (loop.zHeight <= subFeat.topZ + 1e-3 && loop.zHeight >= subFeat.bottomZ - 1e-3) {
                if (!loop.lines.empty()) {
                    Point3D pt = loop.lines.front().start;
                    // XY 包围盒包含检测
                    if (pt.x >= bx1 - 1e-3 && pt.x <= bx2 + 1e-3 &&
                        pt.y >= by1 - 1e-3 && pt.y <= by2 + 1e-3) {
                        subFeat.stepLoops.push_back(loop);
                    }
                }
            }
        }
        subFeatures.push_back(subFeat);
    }
    return subFeatures;
}




// --- 主流程 ---
void ProcessAndSplitClosedCavityFeatures(
    const std::vector<std::vector<Face2D>>& allLayerFaces,
    const TopoDS_Shape& mainShape,
    const std::vector<double>& splitPoints,
    const TopTools_DataMapOfShapeInteger& faceToIdMap,
    const std::string& savePath,
    std::vector<TopoDS_Compound>& outTrueClosedCavities,
    std::vector<CavityFeature>& outFeatures)
{
    std::vector<TopoDS_Compound> trueClosedCavities;

    // 1. 提取所有腔体面 ID 集合
    std::set<int> cavityFaceIds = CollectCavityFaceIds(allLayerFaces);
    if (cavityFaceIds.empty()) {
        return;
    }

    // 2. 得到纯粹的“型腔壳体”
    std::string cavityFacesFileName = savePath + "CavityWallFaces.brep";
    ExportCavityFaces(mainShape, cavityFaceIds, faceToIdMap, cavityFacesFileName);

    // 获取型腔合并面 Compound，供后续所有三维切割和特征工程使用
    TopoDS_Compound cavityCompound = GetCavityCompound(mainShape, cavityFaceIds, faceToIdMap);

    // 3. 利用 BFS（广度优先搜索）拓扑连通算法，把“拼在一起”的壳体按物理区域拆开
    std::vector<TopoDS_Compound> isolatedCavities = SeparateDisconnectedCavities(cavityCompound);
    cout << "    [拓扑分析] 成功将合并壳体拆分为 " << isolatedCavities.size() << " 个空间独立的型腔连通区域。" << endl;

    // 4. 生成高级特征表达
    std::vector<CavityFeature> cavityFeatures = GenerateCavityFeatures(isolatedCavities, allLayerFaces, splitPoints, faceToIdMap);
    
    for (const auto& feat : cavityFeatures) {
        // 智能剪裁：利用 topZ 解析切割，剔除超出开粗顶部的侧壁面（还原真正的工艺加工面）
        TopoDS_Compound trueCavity = ExtractTrueCavityFaces(feat, mainShape, faceToIdMap);
		trueClosedCavities.push_back(trueCavity);
        string filePath = savePath + "Final_True_ClosedCavity_" + to_string(feat.featureId) + ".brep";
        BRepTools::Write(trueCavity, filePath.c_str());
    }

    // 5. 沿z轴分割
    int totalPartCount = 0;
    std::vector<CavityFeature> finalMachinableClosedFeatures;
    for (size_t i = 0; i < trueClosedCavities.size(); ++i) {
        cout << "在处理第 [" << i << "] 个独立型腔区域..." << endl;

        // 获取当前父型腔沿 Z 轴递推切分后的所有层
        std::vector<TopoDS_Compound> parts = RecursiveSplitCavity(trueClosedCavities[i]);

        std::vector<CavityFeature> converted = ConvertPartsToFeatures(parts, cavityFeatures[i], faceToIdMap);

        for (size_t k = 0; k < converted.size(); ++k) {
            converted[k].featureId = totalPartCount;

            std::string fileName = savePath + "Final_CAM_ClosedPart_" + std::to_string(totalPartCount) + ".brep";
            BRepTools::Write(parts[k], fileName.c_str());

            finalMachinableClosedFeatures.push_back(converted[k]);

            totalPartCount++;
        }

    }

    for (const auto& feat : finalMachinableClosedFeatures) {
        std::cout << "\n>>> [独立可加工闭合型腔特征 ID: " << feat.featureId << "]" << std::endl;
        std::cout << "    - 空间范围: Z_Top = " << feat.topZ << " / Z_Floor = " << feat.bottomZ << std::endl;
        std::cout << "    - 深度: " << feat.totalDepth << " mm" << std::endl;
        std::cout << "    - 包含 2D 层级切片数量: " << feat.stepLoops.size() << " 圈" << std::endl;
    }

    outTrueClosedCavities = trueClosedCavities;
    outFeatures = finalMachinableClosedFeatures;
}

