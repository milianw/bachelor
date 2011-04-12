#!/bin/bash
if [[ "$TERM" == "" ]]; then
  TERM=pdf
fi

for f in $1; do
  n=$(basename $(dirname $(dirname $f)));
  GNUPLOT_CMD="set term $TERM; set output '$n.$TERM';" CONVOLUTE_WIDTH=5E-5 TITLE_ADD=$n ../plot_data.sh $f;
 done
