#!/bin/bash
for i in `seq 2280 2304` 
do
    root -l -q 'calib.C('$i')'
done
