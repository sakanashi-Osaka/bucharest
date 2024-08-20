#! /bin/bash
#PBS -q AM -l select=1:ncpus=1:mem=128gb
#PBS -M sakanashi@ne.phys.sci.osaka-u.ac.jp
cd /home/sakra/exp/Bucharest2022/sorter_cor
## 1st sort finished
##run2236-2239,2270- finishe
./test 2411
