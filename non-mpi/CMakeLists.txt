find_package(OpenMP)

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

  add_executable(hs main.cpp)
  OptimizeForArchitecture(hs)
  target_link_libraries(hs ${QT_QTCORE_LIBRARY} hslib)
endif(OPENMP_FOUND)