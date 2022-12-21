#! /bin/bash
#PBS -q AL -l select=1:ncpus=1:mem=128gb
cd /home/sakra/exp/Bucharest2022/sorter_cor
#2305-2369 2nd sort finished  ##2305-2413 1st sort finished
./test 2410
./test 2411
./test 2412
./test 2413
