#include "regression_test.h"
#include "geom_utils.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <set>

static std::string TrimCopy(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    return s.substr(start, end - start + 1);
}

static std::string EscapeJsonString(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s) {
        switch (c) {
            case '\"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:   out += c; break;
        }
    }
    return out;
}

static std::string FeatureTypeToString(CavityType t) {
    if (t == CavityType::OPEN) return "OPEN";
    if (t == CavityType::CLOSED) return "CLOSED";
    return "OTHER";
}

static std::string MakeFeatureSignature(const RegressionFeatureRecord& r) {
    std::ostringstream oss;
    oss << r.type << "|" << r.topZ << "|" << r.bottomZ << "|";
    std::vector<int> ids = r.originalFaceIds;
    std::sort(ids.begin(), ids.end());
    for (size_t i = 0; i < ids.size(); ++i) {
        if (i) oss << "-";
        oss << ids[i];
    }
    return oss.str();
}

static std::set<int> CollectAllFaceIdsFromLayers(const std::vector<std::vector<Face2D>>& allLayerFaces) {
    std::set<int> ids;
    for (const auto& layer : allLayerFaces) {
        for (const auto& f : layer) {
            for (const auto& seg : f.outerLoop) ids.insert(seg.faceId);
            for (const auto& loop : f.innerLoops) {
                for (const auto& seg : loop) ids.insert(seg.faceId);
            }
        }
    }
    return ids;
}

static std::map<int, int> CountPartNumberPerOriginalFaceInFeature(const CavityFeature& feat, const std::set<int>& allFaceIds) {
    std::map<int, int> partCountMap;
    std::set<int> featureFaceIds = feat.sourceFaceIds;
    for (const auto& loop : feat.stepLoops) {
        for (const auto& seg : loop.lines) {
            if (seg.faceId > 0 && featureFaceIds.count(seg.faceId)) {
                partCountMap[seg.faceId] += 1;
            }
        }
    }
    for (int fid : featureFaceIds) {
        if (partCountMap.find(fid) == partCountMap.end()) {
            partCountMap[fid] = 1;
        }
    }
    return partCountMap;
}

static std::vector<std::string> ReadFileLines(const std::string& path) {
    std::ifstream ifs(path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line)) lines.push_back(line);
    return lines;
}

static bool StartsWith(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

static std::string ExtractStringAfter(const std::string& s, const std::string& key) {
    size_t p = s.find(key);
    if (p == std::string::npos) return "";
    p += key.size();
    size_t q = s.find_first_of(",}", p);
    if (q == std::string::npos) q = s.size();
    return TrimCopy(s.substr(p, q - p));
}

static double ExtractDoubleAfter(const std::string& s, const std::string& key, double fallback = 0.0) {
    std::string val = ExtractStringAfter(s, key);
    if (val.empty()) return fallback;
    try { return std::stod(val); } catch (...) { return fallback; }
}

static int ExtractIntAfter(const std::string& s, const std::string& key, int fallback = 0) {
    std::string val = ExtractStringAfter(s, key);
    if (val.empty()) return fallback;
    try { return std::stoi(val); } catch (...) { return fallback; }
}

static std::vector<int> ExtractIntListAfter(const std::string& s, const std::string& key) {
    std::vector<int> out;
    size_t p = s.find(key);
    if (p == std::string::npos) return out;
    p = s.find('[', p);
    size_t q = s.find(']', p);
    if (p == std::string::npos || q == std::string::npos) return out;
    std::string inner = s.substr(p + 1, q - p - 1);
    std::istringstream iss(inner);
    std::string token;
    while (std::getline(iss, token, ',')) {
        token = TrimCopy(token);
        if (token.empty()) continue;
        try { out.push_back(std::stoi(token)); } catch (...) {}
    }
    return out;
}

static bool ContainsString(const std::string& s, const std::string& target) {
    return s.find(target) != std::string::npos;
}

void SaveModelRegressionResult(
    const std::string& outputPath,
    const std::string& modelName,
    const std::string& runId,
    const std::vector<CavityFeature>& features,
    const std::vector<std::vector<Face2D>>& allLayerFaces) {

    RegressionModelResult result;
    result.model = modelName;
    result.runId = runId;

    std::set<int> allFaceIds = CollectAllFaceIdsFromLayers(allLayerFaces);
    std::map<int, int> globalReferenceCount;

    for (const auto& feat : features) {
        RegressionFeatureRecord rec;
        rec.featureId = std::to_string(feat.featureId);
        rec.type = FeatureTypeToString(feat.type);
        rec.topZ = feat.topZ;
        rec.bottomZ = feat.bottomZ;
        rec.originalFaceIds.assign(feat.sourceFaceIds.begin(), feat.sourceFaceIds.end());
        std::sort(rec.originalFaceIds.begin(), rec.originalFaceIds.end());

        std::map<int, int> partCountMap = CountPartNumberPerOriginalFaceInFeature(feat, allFaceIds);
        for (int fid : rec.originalFaceIds) {
            globalReferenceCount[fid] += 1;
            RegressionOriginalFaceRecord r;
            r.originalFaceId = fid;
            r.partCount = partCountMap[fid];
            r.isSplit = r.partCount > 1;
            rec.originalFaces.push_back(r);
            if (r.isSplit) rec.splitOriginalFaces.push_back(fid);
        }
        rec.hasFaceSplit = !rec.splitOriginalFaces.empty();
        rec.signature = MakeFeatureSignature(rec);
        result.features.push_back(rec);
    }

    result.featureCount = (int)result.features.size();
    result.splitOriginalFaceCount = 0;
    for (const auto& f : result.features) result.splitOriginalFaceCount += (int)f.splitOriginalFaces.size();
    result.sharedOriginalFaceCount = 0;
    for (const auto& kv : globalReferenceCount) if (kv.second > 1) result.sharedOriginalFaceCount += 1;

    {
        std::ostringstream typeSig;
        for (size_t i = 0; i < result.features.size(); ++i) {
            if (i) typeSig << "|";
            typeSig << result.features[i].type;
        }
        result.featureSignature = typeSig.str();
    }

    std::ostringstream oss;
    oss << "{\n";
    oss << "  \"model\": \"" << EscapeJsonString(result.model) << "\",\n";
    oss << "  \"runId\": \"" << EscapeJsonString(result.runId) << "\",\n";

    oss << "  \"features\": [\n";
    for (size_t fi = 0; fi < result.features.size(); ++fi) {
        const auto& f = result.features[fi];
        oss << "    {\n";
        oss << "      \"featureId\": \"" << EscapeJsonString(f.featureId) << "\",\n";
        oss << "      \"type\": \"" << f.type << "\",\n";
        oss << "      \"topZ\": " << f.topZ << ",\n";
        oss << "      \"bottomZ\": " << f.bottomZ << ",\n";

        oss << "      \"originalFaceIds\": [";
        for (size_t i = 0; i < f.originalFaceIds.size(); ++i) {
            if (i) oss << ", ";
            oss << f.originalFaceIds[i];
        }
        oss << "],\n";

        oss << "      \"hasFaceSplit\": " << (f.hasFaceSplit ? "true" : "false") << ",\n";

        oss << "      \"splitOriginalFaces\": [";
        for (size_t i = 0; i < f.splitOriginalFaces.size(); ++i) {
            if (i) oss << ", ";
            oss << f.splitOriginalFaces[i];
        }
        oss << "],\n";

        oss << "      \"originalFaces\": [\n";
        for (size_t i = 0; i < f.originalFaces.size(); ++i) {
            const auto& o = f.originalFaces[i];
            oss << "        { \"originalFaceId\": " << o.originalFaceId
                << ", \"isSplit\": " << (o.isSplit ? "true" : "false")
                << ", \"partCount\": " << o.partCount << " }";
            if (i + 1 < f.originalFaces.size()) oss << ",";
            oss << "\n";
        }
        oss << "      ],\n";

        oss << "      \"signature\": \"" << EscapeJsonString(f.signature) << "\"\n";
        oss << "    }";
        if (fi + 1 < result.features.size()) oss << ",";
        oss << "\n";
    }
    oss << "  ],\n";

    oss << "  \"summary\": {\n";
    oss << "    \"featureCount\": " << result.featureCount << ",\n";
    oss << "    \"splitOriginalFaceCount\": " << result.splitOriginalFaceCount << ",\n";
    oss << "    \"sharedOriginalFaceCount\": " << result.sharedOriginalFaceCount << ",\n";
    oss << "    \"featureSignature\": \"" << EscapeJsonString(result.featureSignature) << "\"\n";
    oss << "  }\n";
    oss << "}\n";

    std::ofstream ofs(outputPath);
    ofs << oss.str();
}

RegressionDiffReport CompareModelRegressionResults(
    const std::string& baselineJsonPath,
    const std::string& currentJsonPath,
    double zTol) {

    RegressionDiffReport report;
    report.overallResult = RegressionCompareResult::PASS;

    auto baselineLines = ReadFileLines(baselineJsonPath);
    auto currentLines = ReadFileLines(currentJsonPath);

    auto JoinLines = [](const std::vector<std::string>& v) {
        std::ostringstream oss;
        for (const auto& s : v) oss << s << "\n";
        return oss.str();
    };

    std::string baselineText = JoinLines(baselineLines);
    std::string currentText = JoinLines(currentLines);

    auto UpgradeResult = [&](RegressionDiffReport& r, RegressionCompareResult next) {
        if (static_cast<int>(next) > static_cast<int>(r.overallResult)) r.overallResult = next;
    };

    int baselineCount = ExtractIntAfter(baselineText, "\"featureCount\":");
    int currentCount = ExtractIntAfter(currentText, "\"featureCount\":");
    if (baselineCount != currentCount) {
        UpgradeResult(report, RegressionCompareResult::REGRESSION);
        report.reasons.push_back("featureCount changed: " + std::to_string(baselineCount) + " -> " + std::to_string(currentCount));
    }

    std::string baselineSig = ExtractStringAfter(baselineText, "\"featureSignature\":");
    std::string currentSig = ExtractStringAfter(currentText, "\"featureSignature\":");
    baselineSig = TrimCopy(baselineSig);
    currentSig = TrimCopy(currentSig);
    baselineSig.erase(std::remove(baselineSig.begin(), baselineSig.end(), '\"'), baselineSig.end());
    currentSig.erase(std::remove(currentSig.begin(), currentSig.end(), '\"'), currentSig.end());

    if (baselineSig != currentSig) {
        UpgradeResult(report, RegressionCompareResult::REGRESSION);
        report.reasons.push_back("featureSignature changed: " + baselineSig + " -> " + currentSig);
    }

    auto CollectFeatureSignatures = [](const std::string& text) {
        std::vector<std::string> sigs;
        std::istringstream iss(text);
        std::string line;
        while (std::getline(iss, line)) {
            line = TrimCopy(line);
            if (ContainsString(line, "\"signature\":")) {
                size_t p = line.find(':');
                if (p != std::string::npos) {
                    std::string val = TrimCopy(line.substr(p + 1));
                    val.erase(std::remove(val.begin(), val.end(), '\"'), val.end());
                    val.erase(std::remove(val.begin(), val.end(), ','), val.end());
                    sigs.push_back(val);
                }
            }
        }
        return sigs;
    };

    std::vector<std::string> baselineSigs = CollectFeatureSignatures(baselineText);
    std::vector<std::string> currentSigs = CollectFeatureSignatures(currentText);

    size_t commonCount = std::min(baselineSigs.size(), currentSigs.size());
    for (size_t i = 0; i < commonCount; ++i) {
        if (baselineSigs[i] != currentSigs[i]) {
            UpgradeResult(report, RegressionCompareResult::REGRESSION);
            report.reasons.push_back("feature[" + std::to_string(i) + "] signature changed: " + baselineSigs[i] + " -> " + currentSigs[i]);
        }
    }
    if (currentSigs.size() > baselineSigs.size()) {
        UpgradeResult(report, RegressionCompareResult::CHANGE);
        report.reasons.push_back("new features added: " + std::to_string(currentSigs.size() - baselineSigs.size()));
    } else if (currentSigs.size() < baselineSigs.size()) {
        UpgradeResult(report, RegressionCompareResult::REGRESSION);
        report.reasons.push_back("features removed: " + std::to_string(baselineSigs.size() - currentSigs.size()));
    }

    auto CountMatches = [](const std::string& text, const std::string& key) {
        int count = 0;
        size_t pos = 0;
        while ((pos = text.find(key, pos)) != std::string::npos) {
            ++count;
            pos += key.size();
        }
        return count;
    };

    int baselineSplitCount = CountMatches(baselineText, "\"isSplit\": true");
    int currentSplitCount = CountMatches(currentText, "\"isSplit\": true");
    if (baselineSplitCount != currentSplitCount) {
        UpgradeResult(report, RegressionCompareResult::CHANGE);
        report.reasons.push_back("isSplit count changed: " + std::to_string(baselineSplitCount) + " -> " + std::to_string(currentSplitCount));
    }

    return report;
}
