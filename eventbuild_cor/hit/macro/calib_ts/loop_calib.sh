#!/bin/bash
for i in `seq 0 15`
do
    for j in `seq 0 15`
    do
        root -l -q 'calib.C(2305,'$i','$j')'
    done
done

