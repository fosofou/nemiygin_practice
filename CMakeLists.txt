cmake_minimum_required(VERSION 3.0.0)
project(laba1 VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()



add_executable(laba1 main.cpp matrix.cpp openmp_settings.cpp)

target_include_directories(laba1 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(laba1 OpenMP::OpenMP_CXX)
endif()
target_compile_options(laba1 PRIVATE -fopenmp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
