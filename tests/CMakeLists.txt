include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_subdirectory(googletest)

file(GLOB files "*.cpp" "*.cc")

foreach(filename ${files})
    get_filename_component (TName ${filename} NAME_WE)
    add_executable("${TName}" ${filename})
    add_test(NAME "${TName}" COMMAND ${filename})
    target_link_libraries(${TName} CompMath gtest gtest_main Eigen3::Eigen)
endforeach()