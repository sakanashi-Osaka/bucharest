#!/bin/bash
for i in `seq 2306 2413` #2305-2413 finished
do
    root -l -q 'position.C('$i')'
done
