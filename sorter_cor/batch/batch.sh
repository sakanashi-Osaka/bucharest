#! /bin/bash
#PBS -q AL -l select=1:ncpus=1:mem=128gb
#PBS -M sakanashi@ne.phys.sci.osaka-u.ac.jp
cd /home/sakra/exp/Bucharest2022/sorter_cor
#2305-2369 2nd sort finished  ##2305-2413 1st sort finished
./test 2306
