#!/bin/bash
for i in `seq 2301 2309`
do
    root -l -q 'test.cpp+('$i')'
done
