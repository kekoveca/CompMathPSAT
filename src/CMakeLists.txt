cmake_minimum_required(VERSION 3.24)
project(CompMath)

include(CMakePrintHelpers)
set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS_RELEASE -Wall -Wextra -O3)
set(CMAKE_CXX_FLAGS_DEBUG -g)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DMY_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -DMY_RELEASE")

file(GLOB_RECURSE source *.hpp *.h *.cpp *.cc)
add_library(CompMath INTERFACE ${source})