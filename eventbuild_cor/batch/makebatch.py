#for i in range(2305,2414):
#for i in range(2287,2305):
for i in range(2236,2287):
    f=open("batch"+str(i)+".sh","a")
    f.write("#! /bin/bash\n")
    f.write("#PBS -q AL -l select=1:ncpus=1:mem=2gb\n")
    f.write("#PBS -M sakanashi@ne.phys.sci.osaka-u.ac.jp\n")
    f.write("cd /home/sakra/exp/Bucharest2022/eventbuild_cor\n")
    f.write("./eventbuild_cor "+str(i))