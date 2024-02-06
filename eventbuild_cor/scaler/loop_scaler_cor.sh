#!/bin/bash
#for i in `seq 2236 2269` # finished
for i in `seq 2270 2304` # finished
#for i in `seq 2305 2413` # finished
do
    root -l -q 'scaler.C+('$i')'
done
