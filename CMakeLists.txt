cmake_minimum_required(VERSION 3.24)
project(CompMath)

set(CMAKE_CXX_STANDARD 23)

add_subdirectory(tests)
add_subdirectory(src)
add_subdirectory(lib)

find_package (Eigen3 3.4 REQUIRED NO_MODULE)