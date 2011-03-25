add_library(orca orca.cpp)
OptimizeForArchitecture(orca hslib)

add_library(memcheck memcheck.cpp)
OptimizeForArchitecture(memcheck hslib)

add_library(resfield resfield.cpp)
OptimizeForArchitecture(resfield hslib)

add_library(hamiltonian orca.cpp)
OptimizeForArchitecture(hamiltonian hslib)