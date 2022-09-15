#!/bin/bash
for i in `seq 2300 2330` #2300-2330 finished
do
    root -l -q 'eventbuild_cor.C+('$i')'
done
