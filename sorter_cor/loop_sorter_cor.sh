#!/bin/bash
for i in `seq 2400 2413` #2300-2413 finished
do
    root -l -q 'sorter_cor.C+('$i')'
done
