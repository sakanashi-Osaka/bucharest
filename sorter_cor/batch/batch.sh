#! /bin/bash
#PBS -q AM -l select=1:ncpus=1:mem=128gb
#PBS -M sakanashi@ne.phys.sci.osaka-u.ac.jp
cd /home/sakra/exp/Bucharest2022/sorter_cor
## 1st sort finished
##run2236-2239,2270- finished
./test 2340
./test 2341
./test 2342
./test 2343
./test 2344
./test 2345
./test 2346
./test 2347
./test 2348
./test 2349
