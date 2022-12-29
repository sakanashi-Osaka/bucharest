#!/bin/bash
#for i in `seq 2305 2339` #hit finished
#for i in `seq 2340 2379` #hit finished
for i in `seq 2380 2413` #hit finished
do
    root -l -q 'hit.C+('$i')'
done
