#!/bin/bash
#for i in `seq 2350 2413`
#for i in `seq 2305 2413` 
#for i in `seq 2287 2304`

for i in `seq 2296 2304`
#for i in `seq 2286 2295`
#for i in `seq 2276 2285`
#for i in `seq 2266 2275`
#for i in `seq 2256 2265`
#for i in `seq 2246 2255`
#for i in `seq 2236 2245`
do
    qsub "batch"$i".sh"
done
