#!/bin/bash
#for i in `seq 2305 2339` #2300-2411 finished
#for i in `seq 2340 2379` #2300-2411 finished
for i in `seq 2380 2411` #2300-2411 finished
do
    root -l -q 'timing.C('$i')'
done
