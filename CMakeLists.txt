cmake_minimum_required(VERSION 3.1)
project(mSofaPlugin)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

file(GLOB_RECURSE HEADER_FILES
    "src/*.h"
    "src/*.inl"
)

file(GLOB_RECURSE SOURCE_FILES
    "initPlugin.cpp"
    "src/*.cpp"
)

file(GLOB_RECURSE SCENES_FILES
    "scenes/*.scn"
    "scenes/*.xml"
)

find_package(SofaFramework REQUIRED)
find_package(Sofa.Simulation.Core REQUIRED)
find_package(SofaBaseLinearSolver REQUIRED)
find_package(Metis REQUIRED)
include_directories(${METIS_INCLUDE_DIRS})

set_source_files_properties(${IGNORED_FILES} PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties(${DEPRECATED_FILES} PROPERTIES HEADER_FILE_ONLY TRUE)

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${HEADER_FILES} ${DEPRECATED_FILES} ${IGNORED_FILES} ${SCENES_FILES} ${README_FILES} )
set(${PROJECT_NAME}_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DPLUGIN_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\"")

target_link_libraries(${PROJECT_NAME}
    Sofa.Core
    Sofa.Simulation.Core
    SofaBaseLinearSolver
    metis
)

## Install rules for the library and headers; CMake package configurations files
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${${PROJECT_NAME}_INCLUDE_DIRS}>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>")
sofa_create_package_with_targets(PACKAGE_NAME ${PROJECT_NAME} PACKAGE_VERSION 0.1 TARGETS ${PROJECT_NAME} INCLUDE_INSTALL_DIR ${${PROJECT_NAME}_INCLUDE_DIRS})
