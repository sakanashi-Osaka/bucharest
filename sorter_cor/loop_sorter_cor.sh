#!/bin/bash
for i in `seq 2380 2380` #2300-2411 finished
do
    root -l -q 'sorter_cor.C+('$i')'
done
