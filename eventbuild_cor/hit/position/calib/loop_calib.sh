#!/bin/bash
for i in `seq 2305 2413` 
do
    root -l -q 'calib.C('$i')'
done
