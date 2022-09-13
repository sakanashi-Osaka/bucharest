#!/bin/bash
for i in `seq 2331 2331` #2300-2413 finished
do
    root -l -q 'eventbuild_cor.C+('$i')'
done
