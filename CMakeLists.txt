cmake_minimum_required(VERSION 3.20)

project(
    SliQEC_v2
    LANGUAGES CXX
    DESCRIPTION "SliQEC_v2 - A BDD-based Quantum Circuit Equivalence Checker")

file(GLOB SRC "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(SliQEC_v2 ${SRC})

target_include_directories(SliQEC_v2
    PUBLIC
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/cudd/cudd"
        "${PROJECT_SOURCE_DIR}/cudd/util"
    )

target_compile_options(SliQEC_v2 PUBLIC -Wall -O3 -DNDEBUG)

target_link_libraries(SliQEC_v2 ${CMAKE_SOURCE_DIR}/cudd/cudd/.libs/libcudd.a -lboost_program_options)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
