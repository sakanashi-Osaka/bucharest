#!/bin/bash
for i in `seq 2280 2304` #2305-2413 finished
do
    root -l -q 'position.C+('$i')'
done
