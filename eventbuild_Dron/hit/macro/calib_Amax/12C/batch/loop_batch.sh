#!/bin/bash
#for i in `seq 2236 2259`
#for i in `seq 2260 2279`
#for i in `seq 2280 2304`

for i in `seq 2236 2413`
#for i in `seq 2305 2339`
#for i in `seq 2340 2369`
#for i in `seq 2370 2399`
#for i in `seq 2400 2413`
do
    qsub "batch"$i".sh"
done
