#!/bin/bash
for i in `seq 2300 2399`
do
    qsub "batch"$i".sh"
done
