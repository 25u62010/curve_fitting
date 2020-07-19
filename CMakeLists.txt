cmake_minimum_required( VERSION 2.8 )
project( g2o_curve_fitting )

set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-std=c++11 -O3" )

# 添加cmake模块以使用ceres库
list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules )

# 寻找G2O
find_package( G2O REQUIRED )
include_directories( 
    ${G2O_INCLUDE_DIRS}
    "/usr/include/eigen3"
    ${PROJECT_SOURCE_DIR}/include
)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
# OpenCV
find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_DIRS} )
add_library(${PROJECT_NAME} SHARED
src/curve_fitting.cpp
)
add_executable( curve_fitting example.cpp )
# 与G2O和OpenCV链接
target_link_libraries( curve_fitting 
	${PROJECT_NAME}
	${PROJECT_SOURCE_DIR}/lib
    ${OpenCV_LIBS}
    g2o_core g2o_stuff
)
