#include "cavrity_mfr.h"

int main(int argc, char* argv[]) {
    // 命令行参数: CgalApp.exe <inputPath> <savePath> <stepfile>
    // 若无参数则使用默认值，保持向后兼容
    if (argc >= 4) {
        inputPath = argv[1];
        savePath  = argv[2];
    } else {
        inputPath = "E:\\soft\\code\\cMake_test\\input\\";
        savePath  = "E:\\soft\\code\\cMake_test\\output\\";
    }
    // 确保路径以反斜杠结尾
    if (!inputPath.empty() && inputPath.back() != '\\' && inputPath.back() != '/')
        inputPath += "\\";
    if (!savePath.empty() && savePath.back() != '\\' && savePath.back() != '/')
        savePath += "\\";

    std::string stepfile;
    if (argc >= 4) {
        stepfile = argv[3];
    } else {
        stepfile = "4_stp_stp.stp";
    }

    // 确保输出目录存在
    _mkdir(savePath.c_str());

    STEPControl_Reader reader;
    std::string inputFileName = inputPath + stepfile;
    if (reader.ReadFile(inputFileName.c_str()) != IFSelect_RetDone) {
        cout << "无法读取 STEP 文件！" << endl;
        return 1;
    }
    reader.TransferRoots();
    TopoDS_Shape mainShape = reader.OneShape();
    
    cout << "模型已加载，开始自适应获取切分点..." << endl;
    
    // 为模型中的每个面分配一个唯一的 ID
    TopTools_DataMapOfShapeInteger faceToIdMap;
    int currentFaceId = 1;
    TopExp_Explorer expFace(mainShape, TopAbs_FACE);
    for (; expFace.More(); expFace.Next()) {
        if (!faceToIdMap.IsBound(expFace.Current())) {
            faceToIdMap.Bind(expFace.Current(), currentFaceId++);
        }
    }
    cout << "已为 " << currentFaceId - 1 << " 个面分配了 ID。" << endl;

    std::map<int, double> filletMap = BuildFilletRadiusMap(mainShape, faceToIdMap);
    cout << "已建立圆角半径映射，共 " << filletMap.size() << " 个面。" << endl;

    OpenCavityFilterParams openCavityParams;

    std::vector<double> splitPoints = GetSplitPointsAlongZ(mainShape);
    
    cout << "检测到 " << splitPoints.size() << " 个切分点:" << endl;
    for (int i = 0; i < splitPoints.size(); i++) {
        cout << "  [" << i << "] Z = " << splitPoints[i] << endl;
    }

    cout << "\n开始逐层切分模型 (从上往下)..." << endl;
    
    
    // 分离存储：封闭型腔和开放型腔的数据使用独立容器，避免混合污染
    std::vector<std::vector<Face2D>> allClosedLayerFaces;  // 封闭型腔专用
    std::vector<std::vector<Face2D>> allOpenLayerFaces;    // 开放型腔专用

    for (int i = (int)splitPoints.size() - 1; i >= 0; i--) {
        double splitZ = splitPoints[i];
        int reversedIndex = splitPoints.size() - 1 - i;
        cout << "\n切分 (倒序 #" << reversedIndex << "/" << splitPoints.size() << ") (Z = " << splitZ << ")" << endl;
        
        if ( 2 )
        {
            // 提取本层的面（包含 SOLID + CAVITY）
            std::vector<Face2D> currentLayerFaces = SliceModelAtZ(mainShape, splitZ, i, splitPoints.size()-1, faceToIdMap);

            // 将封闭型腔相关的切片数据保存到专用容器
            allClosedLayerFaces.push_back(currentLayerFaces);

            // 1. 分离出实体面参与运算
            std::vector<Face2D> currentSolidFaces;
            std::vector<Face2D> currentConvexHullFaces;
            for (const auto& f : currentLayerFaces) {
                if (f.type == FaceType::SOLID) {
                    currentSolidFaces.push_back(f);

                    // 1. 转成Point_2 的格式，准备计算二维凸包
                    vector<Point_2> points_2d = ConvertFaceToPoints2D_EK(f);

                    // 2. 计算 2D 凸包
                    std::vector<Point_2> hull_points;
                    CGAL::convex_hull_2(points_2d.begin(), points_2d.end(), std::back_inserter(hull_points));

                    // 3. 还原为 OneLine 闭合环并存入 outerLoop
                    Face2D HullFace = ConvertPointsToHullFace_EK(hull_points, splitZ);
                    currentConvexHullFaces.push_back(HullFace);
                }
            }

            std::vector<Face2D> currentCavityFaces;
            for (const auto& f : currentLayerFaces) {
                if (f.type == FaceType::CAVITY) {
                    currentCavityFaces.push_back(f);
                }
            }

            std::vector<HullItem> hullItems;
            for (size_t h = 0; h < currentConvexHullFaces.size() && h < currentSolidFaces.size(); ++h) {
                HullItem item;
                item.hull = currentConvexHullFaces[h];
                item.solid = currentSolidFaces[h];
                item.index = (int)h;
                hullItems.push_back(item);
            }

            std::vector<HullGroup> mergedHullGroups = BuildMergedHullGroups(hullItems, openCavityParams);
            std::vector<Face2D> currentMergedHullFaces;

            cout << "  [开放型腔凸包合并] 原始凸包数量: " << currentConvexHullFaces.size()
                << " -> 合并后大凸包数量: " << mergedHullGroups.size()
                << " (阈值=" << openCavityParams.hullMergeDistance << "mm)" << endl;

            std::vector<Face2D> pocketResults;

            for (size_t groupIdx = 0; groupIdx < mergedHullGroups.size(); ++groupIdx) {
                const HullGroup& group = mergedHullGroups[groupIdx];
                Face2D mergedHullFace = ComputeMergedHullFace(hullItems, group.memberIndices, splitZ);
                if (mergedHullFace.outerLoop.empty()) {
                    continue;
                }
                currentMergedHullFaces.push_back(mergedHullFace);

                cout << "    - 大凸包组[" << groupIdx << "] 成员数: " << group.memberIndices.size()
                    << " 触发距离: " << group.triggerDistance << endl;

                std::vector<Face2D> solidsInsideHull;
                for (int memberIndex : group.memberIndices) {
                    if (memberIndex >= 0 && memberIndex < (int)hullItems.size()) {
                        solidsInsideHull.push_back(hullItems[memberIndex].solid);
                    }
                }

                std::vector<Face2D> cavitiesInsideHull;
                for (const auto& cavityFace : currentCavityFaces) {
                    bool shouldSubtract = IsFaceInsideFace(cavityFace, mergedHullFace);
                    if (!shouldSubtract) {
                        std::vector<Face2D> intersection =
                            BooleanFacesSingle(cavityFace, mergedHullFace, ClipType::Intersection, splitZ);
                        shouldSubtract = !intersection.empty();
                    }
                    if (shouldSubtract) {
                        cavitiesInsideHull.push_back(cavityFace);
                    }
                }

                std::vector<Face2D> hullResult = { mergedHullFace };
                for (const auto& solidFace : solidsInsideHull) {
                    std::vector<Face2D> newResult;
                    for (const auto& resultFace : hullResult) {
                        std::vector<Face2D> diffResult = BooleanFacesSingle(resultFace, solidFace, ClipType::Difference, splitZ);
                        newResult.insert(newResult.end(), diffResult.begin(), diffResult.end());
                    }
                    hullResult = newResult;
                }

                for (const auto& cavityFace : cavitiesInsideHull) {
                    std::vector<Face2D> newResult;
                    for (const auto& resultFace : hullResult) {
                        std::vector<Face2D> diffResult = BooleanFacesSingle(resultFace, cavityFace, ClipType::Difference, splitZ);
                        newResult.insert(newResult.end(), diffResult.begin(), diffResult.end());
                    }
                    hullResult = newResult;
                }

                pocketResults.insert(pocketResults.end(), hullResult.begin(), hullResult.end());
            }
            
            std::vector<Face2D> tempOpenCavityFaces = RecoverOriginalFaceIdsByGeometry(pocketResults, currentSolidFaces);
            for (auto& f : tempOpenCavityFaces) {
                f.type = FaceType::OPENCAVITY;
            }

			//到处软边界线段，检查 ID 恢复和软边界标记是否正确（faceId == -1 的线段应该就是软边界）
            std::string softEdgesName = savePath + "Slice_" + std::to_string(i) + "_ONLY_SoftEdges_Z" + std::to_string(splitZ) + ".brep";
            ExportSoftEdgesToBrep(tempOpenCavityFaces, softEdgesName);

            vector<Face2D> currentOpenCavityFaces = CleanAndFilterOpenCavities(tempOpenCavityFaces, filletMap, openCavityParams);

            //debug 
            for(auto& f : currentOpenCavityFaces) {
				vector<Face2D> singleFaceVec = { f };
				std::string singleFaceName = savePath + "Slice_" + "SingleOpenCavityFace_ID" + ".brep";
				ExportFace2DToBrep(singleFaceVec, singleFaceName);
                int cc = 0;
            }

            // 将开放型腔数据保存到专用容器
            allOpenLayerFaces.push_back(currentOpenCavityFaces);



			//DEBUG: 导出当前层的实体面和腔面，检查切分结果
            std::string currentLayerFacesName = savePath + "Slice_" + std::to_string(i) + "currentLayerFaces.brep";
            ExportFace2DToBrep(currentLayerFaces, currentLayerFacesName);
            std::string currentSmallHullFaceName = savePath + "Slice_" + std::to_string(i) + "_IndependentSmallHulls_Z" + std::to_string(splitZ) + ".brep";
            ExportFace2DToBrep(currentConvexHullFaces, currentSmallHullFaceName);
            std::cout << "  [DEBUG] 本层独立小凸包已保存: " << currentSmallHullFaceName
                << " 数量=" << currentConvexHullFaces.size() << std::endl;

            std::string currentMergedHullFaceName = savePath + "Slice_" + std::to_string(i) + "_MergedBigHulls_Z" + std::to_string(splitZ) + ".brep";
            ExportFace2DToBrep(currentMergedHullFaces, currentMergedHullFaceName);
            std::cout << "  [DEBUG] 本层合并大凸包已保存: " << currentMergedHullFaceName
                << " 数量=" << currentMergedHullFaces.size() << std::endl;
            std::string currentSolidFaceName = savePath + "Slice_" + std::to_string(i) + "currentSolidFace.brep";
            ExportFace2DToBrep(currentSolidFaces, currentSolidFaceName);
            std::string currentCavityFaceName = savePath + "Slice_" + std::to_string(i) + "currentCavityFacee.brep";
            ExportFace2DToBrep(currentCavityFaces, currentCavityFaceName);
            std::string tempOpenCavityFaceName = savePath + "Slice_" + std::to_string(i) + "temp_OpenCavityFaces_Z" + std::to_string(splitZ) + ".brep";
            ExportFace2DToBrep(tempOpenCavityFaces, tempOpenCavityFaceName);
            std::string currentOpenCavityFaceName = savePath + "Slice_" + std::to_string(i) + "CGAL_OpenCavityFaces_Z" + std::to_string(splitZ) + ".brep";
            ExportFace2DToBrep(currentOpenCavityFaces, currentOpenCavityFaceName);
            int aaa = 0;
        }

    }
    double modelMinZ = 0.0, modelMaxZ = 0.0;
    GetShapeZRange(mainShape, modelMinZ, modelMaxZ);
    bool isInteriorOpen = HasInteriorOpenCavity(allOpenLayerFaces, allClosedLayerFaces, modelMaxZ, splitPoints, 0.5);
    std::vector<TopoDS_Compound> trueClosedCavities;
    std::vector<CavityFeature> openCavityFeatures;
    std::vector<CavityFeature> closedCavityFeatures;
    std::vector<TopoDS_Compound> interiorNewClosedParts;
    TopoDS_Compound interiorNewOpenCavity;
    gp_Dir interiorOpenToolDir;

    if (isInteriorOpen) {
        // ================= 中部开放型腔独立旁路 =================
        ProcessInteriorOpenCavityByBand(allClosedLayerFaces, allOpenLayerFaces, mainShape, splitPoints, faceToIdMap, savePath, interiorNewClosedParts, interiorNewOpenCavity, interiorOpenToolDir);
    } else {
        // ================= 提取封闭型腔特征面 =================
        ProcessAndSplitClosedCavityFeatures(allClosedLayerFaces, mainShape, splitPoints, faceToIdMap, savePath, trueClosedCavities, closedCavityFeatures);

        // ================= 提取开放型腔特征面 =================
        openCavityFeatures = ProcessAndSplitOpenCavityFeatures(allOpenLayerFaces, allClosedLayerFaces, trueClosedCavities, mainShape, splitPoints, faceToIdMap, savePath);
    }

    cout << "\n";
    cout << "======================================================================\n";
    cout << "                      特征提取结果汇总\n";
    cout << "======================================================================\n\n";

    if (isInteriorOpen) {
        cout << "【分支类型】中部开放型腔分支 (InteriorOpen)\n\n";

        cout << "--- 重组封闭型腔 (Z轴分割后) ---\n";
        cout << "  数量: " << interiorNewClosedParts.size() << " 个\n";
        for (size_t i = 0; i < interiorNewClosedParts.size(); ++i) {
            Bnd_Box box;
            BRepBndLib::Add(interiorNewClosedParts[i], box);
            double cxmin, cymin, czmin, cxmax, cymax, czmax;
            box.Get(cxmin, cymin, czmin, cxmax, cymax, czmax);
            cout << "  Part_" << i
                << "  类型=OTHER"
                << "  进刀方向=(0,0,-1)"
                << "  Z范围=[" << czmax << "," << czmin << "]"
                << "  深度=" << (czmax - czmin) << "mm\n";
        }

        cout << "\n--- 不可加工区域 (侧向加工) ---\n";
        Bnd_Box openBox;
        BRepBndLib::Add(interiorNewOpenCavity, openBox);
        double oxmin, oymin, ozmin, oxmax, oymax, ozmax;
        openBox.Get(oxmin, oymin, ozmin, oxmax, oymax, ozmax);
        cout << "  类型=OTHER"
            << "  进刀方向=(" << interiorOpenToolDir.X() << "," << interiorOpenToolDir.Y() << "," << interiorOpenToolDir.Z() << ")"
            << "  Z范围=[" << ozmax << "," << ozmin << "]"
            << "  深度=" << (ozmax - ozmin) << "mm\n";
    } else {
        cout << "【分支类型】正常分支\n\n";

        cout << "--- 封闭型腔 (Z轴分割后) ---\n";
        cout << "  数量: " << closedCavityFeatures.size() << " 个\n";
        for (const auto& feat : closedCavityFeatures) {
            cout << "  ID=" << feat.featureId
                << "  类型=CLOSED  进刀方向=(0,0,-1)"
                << "  Z范围=[" << feat.topZ << "," << feat.bottomZ << "]"
                << "  深度=" << feat.totalDepth << "mm"
                << "  切片层数=" << feat.stepLoops.size() << "\n";
        }

        cout << "\n--- 开放型腔 (加工可行性赛选后) ---\n";
        cout << "  保留: " << openCavityFeatures.size() << " 个\n";
        for (const auto& feat : openCavityFeatures) {
            cout << "  ID=" << feat.featureId
                << "  类型=" << (feat.type == CavityType::OPEN ? "OPEN" : feat.type == CavityType::CLOSED ? "CLOSED" : "OTHER")
                << "  进刀方向=(" << feat.toolDirection.X() << "," << feat.toolDirection.Y() << "," << feat.toolDirection.Z() << ")"
                << "  Z范围=[" << feat.topZ << "," << feat.bottomZ << "]"
                << "  深度=" << feat.totalDepth << "mm";
            if (feat.toolDirection.IsEqual(gp_Dir(0, 0, -1), 1e-6)) {
                cout << "  Z轴分割=" << feat.stepLoops.size() << "层";
            } else {
                cout << "  跳过Z轴分割(侧向)";
            }
            cout << "\n";
        }
    }

    cout << "\n======================================================================\n";
    cout << "生成的 BREP 文件保存在 " << savePath << " 目录。\n";
    cout << "======================================================================\n";

    system("pause");
    return 0;
}