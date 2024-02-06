#!/bin/bash
for i in `seq 2250 2304` #2200- 
do
    root -l -q 'tmp_time.C+('$i')'
done
