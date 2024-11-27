#!/bin/bash
#for i in `seq 2236 2269`
for i in `seq 2270 2304`


#for i in `seq 2305 2339` #2300-2411 finished
#for i in `seq 2340 2379` #2300-2411 finished
#for i in `seq 2380 2411` #2300-2411 finished
#for i in `seq 2412 2413` #2300-2411 finished
do
    root -l -q 'ene_gamma.C('$i')'
done
