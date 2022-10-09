#!/bin/bash
for i in `seq 2380 2411` #2300-2411 finished
do
    root -l -q 'timing.C('$i')'
done
