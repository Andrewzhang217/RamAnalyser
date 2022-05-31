cmake_minimum_required(VERSION 3.21)
project(RamAnalyser)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)

include(FetchContent)

if (NOT biosoup_FOUND)
    if (ram_install)
        message(FATAL_ERROR "Missing package biosoup!")
    endif ()

    FetchContent_Declare(
            biosoup
            GIT_REPOSITORY https://github.com/rvaser/biosoup
            GIT_TAG 0.10.0)

    FetchContent_GetProperties(biosoup)
    if (NOT biosoup_POPULATED)
        FetchContent_Populate(biosoup)
        add_subdirectory(
                ${biosoup_SOURCE_DIR}
                ${biosoup_BINARY_DIR}
                EXCLUDE_FROM_ALL)
    endif ()
endif ()

if (NOT thread_pool_FOUND)
    if (ram_install)
        message(FATAL_ERROR "Missing package thread_pool!")
    endif ()

    FetchContent_Declare(
            thread_pool
            GIT_REPOSITORY https://github.com/rvaser/thread_pool
            GIT_TAG 4.0.0)

    FetchContent_GetProperties(thread_pool)
    if (NOT thread_pool_POPULATED)
        FetchContent_Populate(thread_pool)
        add_subdirectory(
                ${thread_pool_SOURCE_DIR}
                ${thread_pool_BINARY_DIR}
                EXCLUDE_FROM_ALL)
    endif ()
endif ()

add_executable(RamAnalyser src/main.cpp src/minimizer_engine.cpp src/minimizer_engine.h)

target_link_libraries(RamAnalyser
        biosoup::biosoup
        thread_pool::thread_pool)
