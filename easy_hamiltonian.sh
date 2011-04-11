#!/bin/bash
./dbg/orca2easy $1 > /tmp/test.m || exit
t=$(basename $1);
f="/tmp/$t.easy.dat";
if [[ "$2" != "" ]]; then
  B=$2
else
  B=0.3
fi
echo "hamiltonian2(Sys, $B * 1000 * B_direction, 1);" >> /tmp/test.m
echo "exit;" >> /tmp/test.m
matlab -nojvm -nosplash -r "run('/tmp/test.m');"
