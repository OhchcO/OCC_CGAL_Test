cmake_minimum_required(VERSION 3.10)
project(OCC_CGAL_Test)

# 1. 强制要求 C++17 标准 (与你 VS 里的设置保持一致)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 2. 设置各个库的根路径
set(BOOST_ROOT "E:/soft/boost_1_74_0")
set(CGAL_DIR "E:/soft/CGAL/CGAL-5.6.2")
# 指向包含 inc 和 win64 的那个总目录
set(OCC_BASE_DIR "C:/OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0")

# 3. 寻找依赖项
find_package(Boost REQUIRED)
find_package(CGAL REQUIRED)

# 4. 包含头文件路径
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CGAL_INCLUDE_DIRS})
include_directories("${OCC_BASE_DIR}/inc")
include_directories("E:/soft/CGAL/CGAL-5.6.2/auxiliary/gmp/include")

# 5. 添加可执行程序
add_executable(CgalApp
	main.cpp
	cavity_common.h
	cavity_globals.cpp
	geom_utils.h
	geom_utils.cpp
	convex_hull.h
	convex_hull.cpp
	slicing.h
	slicing.cpp
	occ_utils.h
	occ_utils.cpp
	boolean_ops.h
	boolean_ops.cpp
	line_topology.h
	line_topology.cpp
	export_utils.h
	export_utils.cpp
	closed_cavity.h
	closed_cavity.cpp
	open_cavity.h
	open_cavity.cpp
	interior_open.h
	interior_open.cpp
	regression_test.h
	regression_test.cpp
	clipper.engine.cpp
	clipper.offset.cpp
	clipper.rectclip.cpp
)

# 6. 链接库文件 (这里直接给绝对路径，不再依赖脆弱的 find_package)
target_link_libraries(CgalApp PRIVATE 
    CGAL::CGAL 
    ${Boost_LIBRARIES}
    "E:/soft/CGAL/CGAL-5.6.2/auxiliary/gmp/lib/libgmp-10.lib"
    "E:/soft/CGAL/CGAL-5.6.2/auxiliary/gmp/lib/libmpfr-4.lib"
    # 直接链接 OCC 的 vc14 库文件，确保 v145 翻译官能直接找到它们
    "${OCC_BASE_DIR}/win64/vc14/lib/TKernel.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKMath.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKBRep.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKGeomAlgo.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKTopAlgo.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKPrim.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKG3d.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKSTEP.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKSTEPBase.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKSTEPAttr.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKXSBase.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKMesh.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKTopAlgo.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKShHealing.lib"
    "${OCC_BASE_DIR}/win64/vc14/lib/TKBO.lib"      
    "${OCC_BASE_DIR}/win64/vc14/lib/TKPrim.lib"    
    "${OCC_BASE_DIR}/win64/vc14/lib/TKTopAlgo.lib" 
    "${OCC_BASE_DIR}/win64/vc14/lib/TKFeat.lib"   
    "${OCC_BASE_DIR}/win64/vc14/lib/TKBool.lib" 
    "${OCC_BASE_DIR}/win64/vc14/lib/TKGeomBase.lib"
)