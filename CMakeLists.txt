cmake_minimum_required(VERSION 3.20)

project(
    SliQEC_v2
    LANGUAGES CXX
    DESCRIPTION "SliQEC_v2 - A BDD-based Quantum Circuit Equivalence Checker")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

add_subdirectory(src)
