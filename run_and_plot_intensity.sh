#!/bin/bash
echo "running: $@"
#$@ | tee /tmp/.plt
echo
base=$(tail -n 1 /tmp/.plt)
$(dirname $0)/plot_folder.sh "$base"
