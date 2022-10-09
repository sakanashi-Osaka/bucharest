#!/bin/bash
for i in `seq 2370 2411` #2300-2411 finished
do
    root -l -q 'hit.C+('$i')'
done
