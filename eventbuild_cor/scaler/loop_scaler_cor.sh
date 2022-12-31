#!/bin/bash
for i in `seq 2305 2413` # finished
do
    root -l -q 'scaler.C+('$i')'
done
