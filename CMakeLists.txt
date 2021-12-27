cmake_minimum_required(VERSION 3.5)
project(basic_utils LANGUAGES CXX)
add_definitions(-std=c++14 -w)

SET(CMAKE_BUILD_TYPE "Relese")  
#SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -g2 -ggdb")  
# SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall") 

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen3)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Filter)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LQR)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LQR/riccati)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(linearizer_test examples/linearizerTest.cpp)
add_executable(lqr_test examples/lqrTest.cpp)
add_executable(kalman_test examples/kalmanTest.cpp)
add_executable(ori_test examples/oriTest.cpp)
add_executable(pseudoinv_test examples/pseudoinverseTest.cpp)