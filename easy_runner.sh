#!/bin/bash
./dbg/orca2easy $1 > /tmp/test.m || exit
t=$(basename $1);
f="/tmp/$t.easy.dat";
echo "saveSpec(B, Spec, '$f');" >> /tmp/test.m
echo "exit;" >> /tmp/test.m
matlab -nojvm -nosplash -r "run('/tmp/test.m');" || exit
echo "set title '$t (easyspin)'; plot '$f' w lines notitle;" | gnuplot -persist
