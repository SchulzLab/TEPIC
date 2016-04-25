cmake_minimum_required(VERSION 3.0.0)
project(TEPIC_HIC)

include_directories(../include/)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES hicModule.cpp readMatrix.cpp readMatrix.hpp)
add_executable(TEPIC_HIC ${SOURCE_FILES})