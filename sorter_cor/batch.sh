#! /bin/bash
#PBS -q AL -l select=1:ncpus=1:mem=80gb
cd /home/sakra/exp/Bucharest2022/sorter_cor
./test 2300 #2290-2389 finished
