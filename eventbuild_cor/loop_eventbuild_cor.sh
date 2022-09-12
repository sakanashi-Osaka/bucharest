#!/bin/bash
for i in `seq 2385 2399` #2300-2399 finished
do
    root -l -q 'eventbuild_cor.C+('$i')'
done
