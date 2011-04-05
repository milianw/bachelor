#!/bin/bash
./dbg/orca2easy $1 > /tmp/test.m || exit
t=$(basename $1);
f="/tmp/$t.easy.dat";
echo "saveSpec(B, Spec, '$f');" >> /tmp/test.m
echo "exit;" >> /tmp/test.m
matlab -nojvm -nosplash -r "run('/tmp/test.m');" || exit

echo "plotting files from: $f"

./convolute/convolute $f > ${f}.conv

t="Spectrum for $f, mwFreq = 9.5 GHz (EasySpin)"
pltraw=", '$f' w impulses notitle";
if [[ "$NOPLOT_RAW" == "1" ]]; then
  pltraw=""
fi
echo "$GNUPLOT_CMD set title '$t';
      set xlabel 'static B-field in Tesla';
      set ylabel 'intensity in a.u.';
      set grid;
      plot '${f}.conv' w lines notitle$pltraw" | gnuplot -persist
