add_executable(convolute main.cpp)
OptimizeForArchitecture(convolute)
target_link_libraries(convolute ${QT_QTCORE_LIBRARY})