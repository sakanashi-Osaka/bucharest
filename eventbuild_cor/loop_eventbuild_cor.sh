#!/bin/bash
for i in `seq 2380 2380` #2290-2379 finished
do
    root -l -q 'eventbuild_cor.C+('$i')'
done
