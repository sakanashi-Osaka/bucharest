#! /bin/bash
#PBS -q AL -l select=1:ncpus=1:mem=80gb
cd /home/sakra/exp/Bucharest2022/sorter_cor
./test 2269 #2269,2300-2413 finished
