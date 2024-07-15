cmake_minimum_required(VERSION 3.00)
project(hls)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)



add_subdirectory(least_square)

add_executable(${PROJECT_NAME} test.cpp matrix.cpp)
add_executable(ba ba.cpp)
