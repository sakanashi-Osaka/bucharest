#!/bin/bash
for i in `seq 2380 2380` #test2300-2379 finished, hit2290-2379 finished
do
    root -l -q 'hit.C+('$i')'
done
