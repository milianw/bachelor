cmake_minimum_required(VERSION 2.8)

project(MilianWolffBachelor)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(OptimizeForArchitecture)

find_package(GSL REQUIRED)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Qt4 REQUIRED)
include_directories(${QT_INCLUDE_DIR})

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.42.0 COMPONENTS mpi serialization program_options regex system filesystem REQUIRED )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${Boost_INCLUDE_DIRS} /usr/include/openmpi)

add_subdirectory(spinlib)
add_subdirectory(mpi)
add_subdirectory(non-mpi)
add_subdirectory(convolute)
add_subdirectory(dbg)

add_executable(test test.cpp)
OptimizeForArchitecture(test)
target_link_libraries(test hslib ${Boost_LIBRARIES})