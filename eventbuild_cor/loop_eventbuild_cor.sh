#!/bin/bash
for i in `seq 2290 2300` #2300-2411 finished
do
    root -l -q 'eventbuild_cor.C+('$i')'
done
