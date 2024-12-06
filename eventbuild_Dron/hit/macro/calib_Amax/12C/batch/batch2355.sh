#! /bin/bash
#PBS -q AM -l select=1:ncpus=1:mem=2gb
#PBS -M sakanashi@ne.phys.sci.osaka-u.ac.jp
cd /home/sakra/exp/Bucharest2022/eventbuild_Dron/hit/macro/calib_Amax/12C
./calib 2355