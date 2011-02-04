#!/bin/bash
echo "plotting files from: $1"
cat $1/* | sort -g > /tmp/.plt
echo "$GNUPLOT_CMD set grid; plot '/tmp/.plt' w lines" | gnuplot -persist
