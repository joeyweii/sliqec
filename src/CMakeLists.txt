add_executable(${PROJECT_NAME}
    main.cpp
    gateOpe.cpp
    memMeasure.cpp
    misc.cpp
    propCheck.cpp
    qasmParser.cpp
    qcCheck.cpp)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/cudd/cudd"
        "${PROJECT_SOURCE_DIR}/cudd/util"
    )

target_compile_options(${PROJECT_NAME} PUBLIC -Wall -O3 -DNDEBUG)

target_link_libraries(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/cudd/cudd/.libs/libcudd.a -lboost_program_options)
