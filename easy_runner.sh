#!/bin/bash
./dbg/orca2easy $1 > /tmp/test.m || exit
t=$(basename $1);
f="/tmp/$t.easy.dat";
if [[ "$2" != "" ]]; then
  mw=$2
else
  mw=9.5
fi
echo "Exp.mwFreq = $mw;" >> /tmp/test.m
echo "[B, Spec] = pepper(Sys, Exp, Opt);" >> /tmp/test.m
echo "saveSpec(B, Spec, '$f');" >> /tmp/test.m
echo "exit;" >> /tmp/test.m
matlab -nojvm -nosplash -r "run('/tmp/test.m');" || exit

echo "plotting files from: $f"

./convolute/convolute $f > ${f}.conv

t="Spectrum for $f, mwFreq = $mw GHz (EasySpin)"
echo "$GNUPLOT_CMD set title '$t';
      set xlabel 'static B-field in mT';
      set ylabel 'intensity in a.u.';
      set grid;
      plot '$f' w impulses notitle" | gnuplot -persist
