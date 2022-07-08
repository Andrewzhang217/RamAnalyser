cmake_minimum_required(VERSION 3.21)
project(RamAnalyser)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)

include(FetchContent)

find_package(biosoup 0.10.0 QUIET)
find_package(thread_pool 4.0.0 QUIET)
find_package(bioparser 3.0.13 QUIET)
find_package(ram 2.2.0 QUIET)

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

if (NOT bioparser_FOUND)
    FetchContent_Declare(
            bioparser
            GIT_REPOSITORY https://github.com/rvaser/bioparser
            GIT_TAG 3.0.13)

    FetchContent_GetProperties(bioparser)
    if (NOT bioparser_POPULATED)
        FetchContent_Populate(bioparser)
        add_subdirectory(
                ${bioparser_SOURCE_DIR}
                ${bioparser_BINARY_DIR}
                EXCLUDE_FROM_ALL)
    endif ()
endif ()

if (NOT ram_FOUND)
    FetchContent_Declare(
            ram
            GIT_REPOSITORY https://github.com/lbcb-sci/ram.git
            GIT_TAG 2.2.0)

    FetchContent_GetProperties(ram)
    if (NOT ram_POPULATED)
        FetchContent_Populate(ram)
        add_subdirectory(
                ${ram_SOURCE_DIR}
                ${ram_BINARY_DIR}
                EXCLUDE_FROM_ALL)
    endif ()
endif ()

file(GLOB_RECURSE srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/**.cpp")
file(GLOB_RECURSE headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/**.h")
add_library(RamAnalyser ${srcs} ${headers} src/alias.h)

target_include_directories(RamAnalyser
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>")
target_link_libraries(RamAnalyser
        biosoup::biosoup
        thread_pool::thread_pool
        bioparser::bioparser
        ram::ram)
target_compile_options(RamAnalyser PUBLIC
        $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>
        $<$<CONFIG:Debug,RelWithDebInfo>:-fno-omit-frame-pointer>)
target_link_options(RamAnalyser PUBLIC
        $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>)
include(${CMAKE_CURRENT_LIST_DIR}/test/RamAnalyserTest.cmake)

## this makes the headers accessible for other projects which uses spa lib
#target_include_directories(spa PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
