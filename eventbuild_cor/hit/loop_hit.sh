#!/bin/bash
for i in `seq 2300 2349` #2300-2349 finished
do
    root -l -q 'hit.C+('$i')'
done
