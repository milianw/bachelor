#!/bin/bash
echo "running: $@"
$@ | tee /tmp/.plt
echo
base=$(tail -n 1 /tmp/.plt)
if [[ -d "$base" ]]; then
  $(dirname $0)/plot_folder.sh "$base"
else
  echo aborted
  exit 1
fi
