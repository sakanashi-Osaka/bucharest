#!/bin/bash
for i in `seq 2272 2272` # finished

#for i in `seq 2287 2295` #2305-2413 finished
#for i in `seq 2296 2304` #2305-2413 finished

#for i in `seq 2305 2339` #2305-2413 finished
#for i in `seq 2340 2379` #2305-2413 finished
#for i in `seq 2380 2413` #2305-2413 finished
do
    root -l -q 'sorter_cor.C+('$i')'
done
