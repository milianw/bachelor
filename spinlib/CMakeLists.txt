add_library(hslib spinhamiltonian.cpp resonancefield.cpp orcaparser.cpp experiment.cpp helpers)
target_link_libraries(hslib  ${Boost_LIBRARIES})
OptimizeForArchitecture(hslib)