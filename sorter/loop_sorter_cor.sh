#!/bin/bash
for i in `seq 2400 2425` #2305-2425 finished
do
    root -l -q 'sorter_cor.C+('$i')'
done
