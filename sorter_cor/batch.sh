#! /bin/bash
#PBS -q AL -l select=1:ncpus=1:mem=80gb
cd /home/sakra/data/Bucharest2022/sorter_cor
./test 2407 #2300-2407 finished
