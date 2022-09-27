#!/bin/bash
for i in `seq 2269 2269` #2300-2330 finished
do
    root -l -q 'eventbuild_cor.C+('$i')'
done
