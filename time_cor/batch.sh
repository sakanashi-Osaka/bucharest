#! /bin/bash
#PBS -q AL -l select=1:ncpus=1:mem=8gb
cd /home/sakra/exp/Bucharest2022/time_cor
./test 2260
