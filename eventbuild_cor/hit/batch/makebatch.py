for i in range(2287,2305):
    f=open("batch"+str(i)+".sh","a")
    f.write("#! /bin/bash\n")
    f.write("#PBS -q AL -l select=1:ncpus=1:mem=4gb\n")
    f.write("#PBS -M sakanashi@ne.phys.sci.osaka-u.ac.jp\n")
    f.write("cd /home/sakra/exp/Bucharest2022/eventbuild_cor/hit\n")
    f.write("./hit "+str(i))
