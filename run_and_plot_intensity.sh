#!/bin/bash
echo "running: $@"
$@ | tee /tmp/.plt
echo
base=$(tail -n 1 /tmp/.plt)
if [[ -d "$base" || -f "$base" ]]; then
  $(dirname $0)/plot_data.sh "$base"
else
  echo aborted
  exit 1
fi
