#!/bin/bash
echo "plotting files from: $1"

plt_base=/tmp/.plt_$(basename $1)_

if [ -d $1 ]; then
  cat $1/* | sort -g > ${plt_base}raw
else
  cat $1 | sort -g > ${plt_base}raw
fi

./convolute/convolute $1 > ${plt_base}conv
if [[ "$(basename $1)" == "intensity.data" ]]; then
  # mpi
  a=($(basename $(dirname $1) | tr ":" "\n"))
else
  # non-mpi
  a=($(basename $1 | tr ":" "\n"))
fi

if [[ "${a[1]}" == "-1" ]]; then
  steps="auto"
else
  steps=${a[1]}
fi

if [[ "$TITLE_ADD" != "" ]]; then
  TITLE_ADD="\\n$TITLE_ADD"
fi

t="Spectrum for ${a[2]} 1H, ${a[3]} 14H, mwFreq = ${a[4]} GHz\\nB-Range: ${a[0]} ($steps steps) ${a[5]}$TITLE_ADD"

pltraw="'${plt_base}conv' using 1:(\$3 == 1 ? \$2 : 1/0) w impulses notitle,";

if [[ "$PLOT_NO_RAW" != "" ]]; then
  pltraw=""
fi

echo "$GNUPLOT_CMD set title \"$t\";
      set xlabel 'static B-field in Tesla';
      set ylabel 'intensity in a.u.';
      set grid;
      plot $pltraw \\
      '${plt_base}conv' using 1:(\$3 == 0 ? \$2 : 1/0) w lines notitle lt 7" | gnuplot -persist
