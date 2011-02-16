cmake_minimum_required(VERSION 2.8)

project(MilianWolffBachelor)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(OptimizeForArchitecture)

find_package(GSL REQUIRED)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Qt4 REQUIRED)
include_directories(${QT_INCLUDE_DIR})

find_package(OpenMP)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_executable(hs_hossam hamiltonian+selection.cpp)
OptimizeForArchitecture(hs_hossam)
target_link_libraries(hs_hossam ${GSL_LIBRARIES})

add_executable(hs main.cpp)
OptimizeForArchitecture(hs)
target_link_libraries(hs ${QT_QTCORE_LIBRARY})

add_executable(convolute convolute.cpp)
OptimizeForArchitecture(convolute)
target_link_libraries(convolute ${QT_QTCORE_LIBRARY})
