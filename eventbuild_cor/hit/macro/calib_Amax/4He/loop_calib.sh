#!/bin/bash

for run in `seq 2400 2413`
do
    for i in `seq 0 15`
    do
	for j in `seq 0 15`
	do
            root -l -q 'calib.C('$run','$i','$j')'
	done
    done
done

