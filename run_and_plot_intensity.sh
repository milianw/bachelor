#!/bin/bash
$@ > /tmp/.plt && (echo "plot '/tmp/.plt' w lines" | gnuplot -persist)
