#include "cavity_common.h"

// ======================== Global Variable Definitions ========================

std::string inputPath;
std::string savePath;

const double SCALE = 1000000.0;

BRep_Builder gBuilder;
TopoDS_Compound gAllFacesCompound;
TopoDS_Compound gAllLinesCompound;
bool gCompoundsInitialized = false;
