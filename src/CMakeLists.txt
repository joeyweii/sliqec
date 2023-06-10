add_subdirectory("${PROJECT_SOURCE_DIR}/lib/sliqsys" "lib/sliqsys" EXCLUDE_FROM_ALL)

add_executable(${PROJECT_NAME}
    main.cpp
    memMeasure.cpp
    qcCheck.cpp)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        "${PROJECT_SOURCE_DIR}/include"
    )

target_compile_options(${PROJECT_NAME} PUBLIC -Wall -O3 -DNDEBUG)

target_link_libraries(${PROJECT_NAME} PUBLIC sliqsys)
