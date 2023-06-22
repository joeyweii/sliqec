add_library(
    ${PROJECT_NAME} STATIC
    ${PROJECT_SOURCE_DIR}/include/bddSystem.h
    ${PROJECT_SOURCE_DIR}/include/qCkt.h
    gateOpe.cpp
    propCheck.cpp
    misc.cpp
    qasmParser.cpp)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/cudd/cudd"
        "${PROJECT_SOURCE_DIR}/cudd/util"
    )

target_compile_options(${PROJECT_NAME} PUBLIC -Wall -O3 -DNDEBUG)

target_link_libraries(${PROJECT_NAME} 
    ${PROJECT_SOURCE_DIR}/cudd/cudd/.libs/libcudd.a)
