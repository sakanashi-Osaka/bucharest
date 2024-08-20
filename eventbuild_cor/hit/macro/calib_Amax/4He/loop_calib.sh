#!/bin/bash

#for run in `seq 2236 2249`
#for run in `seq 2250 2269`
#for run in `seq 2270 2289`
for run in `seq 2290 2304`
do
    for i in `seq 0 15`
    do
	for j in `seq 0 15`
	do
            root -l -q 'calib.C('$run','$i','$j')'
	done
    done
done

