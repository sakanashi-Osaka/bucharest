#!/bin/bash
for i in `seq 2300 2389` #2300-2389 finished
do
    root -l -q 'hit.C+('$i')'
done
