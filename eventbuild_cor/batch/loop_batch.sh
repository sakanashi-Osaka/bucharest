#!/bin/bash
for i in `seq 2350 2413` #run2305-2413 finished
#for i in `seq 2305 2413` #run2305-2413 finished
#for i in `seq 2287 2304`
#for i in `seq 2236 2286`
do
    qsub "batch"$i".sh"
done
