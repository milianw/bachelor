#!/bin/bash
for f in $1; do
  FREQUENCY_THRESHOLD=1E-6 mpirun -np 3 ./mpi/hs-mpi -s $f -m 9.5 -i 2048 -o /local_scratch/milianw/$(basename $f);
done
