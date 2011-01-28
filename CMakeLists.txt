cmake_minimum_required(VERSION 2.8)

project(MilianWolffBachelor)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(OptimizeForArchitecture)

find_package(GSL REQUIRED)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(hs_hossam hamiltonian+selection.cpp)
OptimizeForArchitecture(hs_hossam)
target_link_libraries(hs_hossam ${GSL_LIBRARIES})

MACRO(BUILD_HS PROTONS)
# add_definitions(-DPROTONS=${PROTONS})
set(HS_NAME hs-${PROTONS})
add_executable(${HS_NAME} main.cpp)
set_target_properties(${HS_NAME} PROPERTIES COMPILE_FLAGS -DPROTONS=${PROTONS})
OptimizeForArchitecture(${HS_NAME})
target_link_libraries(${HS_NAME})
# remove_definitions(-DPROTONS=${PROTONS})
ENDMACRO(BUILD_HS)

BUILD_HS(2)
BUILD_HS(5)
BUILD_HS(7)
