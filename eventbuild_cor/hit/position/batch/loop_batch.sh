#!/bin/bash
for i in `seq 2290 2304`
do
    qsub "batch"$i".sh"
done
