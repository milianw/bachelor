#!/bin/bash
echo "plotting files from: $1"

if [[ "$(basename "$1")" == "intensity.data" ]]; then
  # mpi
  basestring=$(basename "$(dirname "$1")")
else
  # non-mpi
  basestring=$(basename "$1")
fi

plt_base=/tmp/.plt_${basestring}_

if [ -d "$1" ]; then
  cat "$1"/* | sort -g > "${plt_base}raw"
else
  cat "$1" | sort -g > "${plt_base}raw"
fi

./convolute/convolute "$1" > "${plt_base}conv"

date=${basestring/=*/}
exp=${basestring/*=/}

a=($(echo $exp | tr ":" "\n"))
if [[ "${a[1]}" == "-1" ]]; then
  steps="auto"
else
  steps=${a[1]}
fi

if [[ "$TITLE_ADD" != "" ]]; then
  TITLE_ADD="\\n$TITLE_ADD"
fi

t="Spectrum for ${a[2]}x 1H, ${a[3]}x 14H, mwFreq = ${a[4]} GHz\\nB-Range: ${a[0]} ($steps steps) ${a[5]} $date $TITLE_ADD"

echo "$GNUPLOT_CMD set title \"$t\";
      set xlabel 'static B-field in Tesla';
      set ylabel 'intensity in a.u.';
      set grid;
      plot '${plt_base}conv' using 1:(\$3 == 0 ? \$2 : 1/0) w lines notitle, \\
      '${plt_base}conv' using 1:(\$3 == 1 ? \$2 : 1/0) w impulses notitle" | gnuplot -persist
