#!/bin/bash
echo $@
$@ > /tmp/.plt
echo "plot '/tmp/.plt' w lines" | gnuplot -persist
