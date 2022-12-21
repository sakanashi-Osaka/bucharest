#!/bin/bash
for i in `seq 2353 2353` #2305-2413 finished
do
    root -l -q 'sorter_cor.C+('$i')'
done
