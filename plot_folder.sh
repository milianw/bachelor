#!/bin/bash
echo "plotting files from: $1"
cat $1/* | sort -g > /tmp/.plt
t=$(basename $1)
echo "$GNUPLOT_CMD set grid; plot '/tmp/.plt' w lines title '$t'" | gnuplot -persist
