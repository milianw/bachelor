#!/bin/bash
echo "plotting files from: $1"
# cat $1/* | sort -g > /tmp/.plt
./convolute $1 > /tmp/.plt
a=($(basename $1 | tr ":" "\n"))
t="Spectrum for ${a[0]} protons, mwFreq = ${a[3]} GHz, B-Range: ${a[1]} (${a[2]} steps) ${a[4]}"
echo "$GNUPLOT_CMD set title '$t'; set xlabel 'static B-field in Tesla'; set ylabel 'intensity in a.u.'; set grid; plot '/tmp/.plt' w lines notitle" | gnuplot -persist
