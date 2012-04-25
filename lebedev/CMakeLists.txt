add_executable(lebedev lebedev.c)
OptimizeForArchitecture(lebedev)
#target_link_libraries(lebedev ${FFTW_LIBRARIES})