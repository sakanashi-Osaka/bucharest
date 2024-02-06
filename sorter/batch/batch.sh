#! /bin/bash
#PBS -q AM -l select=1:ncpus=1:mem=128gb
#PBS -M sakanashi@ne.phys.sci.osaka-u.ac.jp
cd /home/sakra/exp/Bucharest2022/sorter
#2305--2425 sort finished
#2280-2304 sort finished
#2230-2379 sort finished
./test 2239
./test 2238
./test 2237
./test 2236
./test 2235
./test 2234
./test 2233
./test 2232
./test 2231
./test 2230
