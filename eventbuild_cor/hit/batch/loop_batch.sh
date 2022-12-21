#!/bin/bash
for i in `seq 2350 2413`
do
    qsub "batch"$i".sh"
done
