#!/bin/bash
for i in `seq 2307 2413` 
do
    root -l -q 'calib.C('$i')'
done
