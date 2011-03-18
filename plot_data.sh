#!/bin/bash
echo "plotting files from: $1"

if [ -d $1 ]; then
  cat $1/* | sort -g > /tmp/.plt
else
  cat $1 | sort -g > /tmp/.plt
fi

./convolute/convolute $1 > /tmp/.plt_conv
a=($(basename $1 | tr ":" "\n"))
t="Spectrum for ${a[0]} 1H, ${a[1]} 14H, mwFreq = ${a[4]} GHz, B-Range: ${a[2]} (${a[3]} steps) ${a[5]}"
echo "$GNUPLOT_CMD set title '$t';
      set xlabel 'static B-field in Tesla';
      set ylabel 'intensity in a.u.';
      set grid;
      plot '/tmp/.plt_conv' w lines notitle, '/tmp/.plt' w impulses notitle" | gnuplot -persist
