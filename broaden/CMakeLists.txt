add_executable(broaden broaden.c)
OptimizeForArchitecture(broaden)
target_link_libraries(broaden ${FFTW_LIBRARIES})