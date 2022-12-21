#!/bin/bash
for i in `seq 2360 2389` #hit finished
do
    root -l -q 'hit.C+('$i')'
done
