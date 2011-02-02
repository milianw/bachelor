#!/bin/bash
echo "running: $@"
$@ | tee /tmp/.plt
echo
base=$(tail -n 1 /tmp/.plt)
echo "plotting files from: $base"
cat $base/* | sort -n > /tmp/.plt
echo "plot '/tmp/.plt' w lines" | gnuplot -persist
