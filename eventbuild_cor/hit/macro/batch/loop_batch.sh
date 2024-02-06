#!/bin/bash
#for i in `seq 2236 2304`
#for i in `seq 2305 2354`
#for i in `seq 2355 2413`
#for i in `seq 2305 2413`
for i in `seq 2305 2339`
do
    qsub "batch"$i".sh"
done
