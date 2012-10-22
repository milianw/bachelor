add_executable(powder powder.c)
OptimizeForArchitecture(powder)
target_link_libraries(powder ${FFTW_LIBRARIES})