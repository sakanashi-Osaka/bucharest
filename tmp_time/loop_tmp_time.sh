#!/bin/bash
for i in `seq 2420 2425` # 2300-2425 finished
do
    root -l -q 'tmp_time.C+('$i')'
done
