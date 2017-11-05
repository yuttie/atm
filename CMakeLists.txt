cmake_minimum_required(VERSION 3.9)

project(atm
    VERSION 1.2.0
    DESCRIPTION "An implementation of the purity measure and more."
    LANGUAGES CXX)

set(default_build_type "Release")
if(CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}'.")
else()
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
endif()

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(PROJECT_BINARY_DIR ${CMAKE_BINARY_DIR})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS         "-Wall")
set(CMAKE_EXE_LINKER_FLAGS  "-static")

set(Boost_USE_STATIC_LIBS    ON)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost
    REQUIRED
    COMPONENTS timer chrono system)
include_directories(${Boost_INCLUDE_DIRS})

configure_file(
    "${PROJECT_SOURCE_DIR}/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h"
    )
include_directories("${PROJECT_BINARY_DIR}")

add_executable(atm ${PROJECT_SOURCE_DIR}/main.cpp)
include_directories(${PROJECT_SOURCE_DIR})
target_link_libraries(atm ${Boost_LIBRARIES})

install(TARGETS atm DESTINATION bin)
