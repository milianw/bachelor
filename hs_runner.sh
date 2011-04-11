#!/bin/bash
if [[ "$2" != "" ]]; then
  mw=$2
else
  mw=9.5
fi
../run_and_plot_intensity.sh ./non-mpi/hs -s $1 -i 0-0:auto:$mw
