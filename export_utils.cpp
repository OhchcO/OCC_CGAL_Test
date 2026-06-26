#include "export_utils.h"
#include "geom_utils.h"
#include "occ_utils.h"

// --- ExportOneLinesToBrep + ExportCavityFaces ---
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

// 保存整个封闭型腔（侧壁 + 顶/底面）为 BREP
void ExportCavityFaces(const TopoDS_Shape& solid,
                        const std::set<int>& cavityFaceIds,
                        const TopTools_DataMapOfShapeInteger& faceToIdMap,
                        const std::string& fileName) {
    // 1. 获取侧壁面
    TopoDS_Compound wallFaces = GetFacesByFaceIds(solid, cavityFaceIds, faceToIdMap);
    
    // 2. 获取顶/底面/中间平台面
    TopoDS_Compound capFaces = GetCavityCapFaces(solid, cavityFaceIds, faceToIdMap);
    
    // 3. 合并
    BRep_Builder builder;
    TopoDS_Compound result;
    builder.MakeCompound(result);
    
    TopExp_Explorer exp(wallFaces, TopAbs_FACE);
    for (; exp.More(); exp.Next()) builder.Add(result, exp.Current());
    
    TopExp_Explorer expCap(capFaces, TopAbs_FACE);
    for (; expCap.More(); expCap.Next()) builder.Add(result, expCap.Current());
    
    BRepTools::Write(result, fileName.c_str());
    std::cout << "  完整封闭型腔已保存至: " << fileName << std::endl;
}

// 获取整个封闭型腔的合并面（侧壁 + 顶/底面/中间平台），返回 Compound

// --- ExportFace2DToBrep ---
void ExportFace2DToBrep(const std::vector<Face2D>& faces, const std::string& fileName) {
    if (faces.empty()) return;

    BRep_Builder builder;
    TopoDS_Compound comp;
    builder.MakeCompound(comp);

    for (const auto& face2d : faces) {
        try {
            if (face2d.outerLoop.size() < 3) continue;

            // 1. 构建外环实心面
            BRepBuilderAPI_MakePolygon outerPoly;
            for (const auto& line : face2d.outerLoop) {
                outerPoly.Add(gp_Pnt(line.start.x, line.start.y, line.start.z));
            }
            outerPoly.Close(); 
            if (!outerPoly.IsDone()) continue;
            
            TopoDS_Wire outerWire = outerPoly.Wire();
            if (outerWire.IsNull() || !outerWire.Closed()) continue;

            BRepBuilderAPI_MakeFace outerFaceMaker(outerWire);
            if (!outerFaceMaker.IsDone()) continue;
            TopoDS_Shape currentShape = outerFaceMaker.Face();

            // 2. 用布尔减法挖去所有内环（孔洞）
            for (const auto& innerLoop : face2d.innerLoops) {
                if (innerLoop.size() < 3) continue;
                
                BRepBuilderAPI_MakePolygon innerPoly;
                for (const auto& line : innerLoop) {
                    innerPoly.Add(gp_Pnt(line.start.x, line.start.y, line.start.z));
                }
                innerPoly.Close();
                
                if (innerPoly.IsDone()) {
                    TopoDS_Wire innerWire = innerPoly.Wire();
                    if (!innerWire.IsNull() && innerWire.Closed()) {
                        BRepBuilderAPI_MakeFace innerFaceMaker(innerWire);
                        if (innerFaceMaker.IsDone()) {
                            // 核心修复：使用极其稳定的布尔减法，大面减小面
                            BRepAlgoAPI_Cut cutAlgo(currentShape, innerFaceMaker.Face());
                            cutAlgo.Build();
                            if (cutAlgo.IsDone()) {
                                currentShape = cutAlgo.Shape();
                            }
                        }
                    }
                }
            }

            // 3. 将最终生成好的面加入 Compound
            builder.Add(comp, currentShape);
        }
        catch (const Standard_Failure& e) {
            std::cerr << "  [警告] 生成面时发生 OCCT 异常: " << e.GetMessageString() << std::endl;
        }
        catch (...) {
            std::cerr << "  [警告] 生成面时发生未知异常" << std::endl;
        }
    }

    BRepTools::Write(comp, fileName.c_str());
    std::cout << "  重建的面已保存至: " << fileName << std::endl;
}



//获取切分结果：切分面 + 切分线（减去凹边）

// --- ExportSoftEdgesToBrep ---
void ExportSoftEdgesToBrep(const std::vector<Face2D>& cavityFaces, const std::string& filePath)
{
    std::vector<OneLine> softLines;

    // 1. 搜刮所有型腔面的外环
    for (const auto& face : cavityFaces) {
        for (const auto& line : face.outerLoop) {
            if (line.faceId == -1) {
                softLines.push_back(line);
            }
        }

        // 2. 顺手搜刮可能存在的内孔（理论上凸包求差的内孔一般贴着实体，但安全起见也扫一遍）
        for (const auto& innerLoop : face.innerLoops) {
            for (const auto& line : innerLoop) {
                if (line.faceId == -1) {
                    softLines.push_back(line);
                }
            }
        }
    }

    if (softLines.empty()) {
        std::cout << "⚠️ [软边界可视化] Z 轴该层未检测到任何 faceId == -1 的软边界线段！" << std::endl;
        return;
    }

    // 3. 借用你原有的 ExportOneLinesToBrep 刚性写入磁盘
    ExportOneLinesToBrep(softLines, filePath);
    std::cout << "🚀 [软边界可视化成功] 共提取出 " << softLines.size() << " 条软边界线段 -> " << filePath << std::endl;
}
//开放型腔清洗流程
/**
 * @brief 【特征量化】已知 Face2D 的外环计算其几何紧实度 (Isoperimetric Quotient)
 * @param face 输入的 2D 面结构
 * @return 紧实度值，范围 (0, 1]。越接近 1 越胖（趋近于圆/正方形），越接近 0 越骨感（趋近于长毛刺）
 */
