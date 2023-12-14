cmake_minimum_required(VERSION 3.0.0)
project(laba1 VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

add_executable(laba1 main.cpp matrix.cpp)

target_include_directories(laba1 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
