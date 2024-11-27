for i in range(2236,2414):
    f=open("batch"+str(i)+".sh","a")
    f.write("#! /bin/bash\n")
    f.write("#PBS -q AM -l select=1:ncpus=1:mem=128gb\n")
    f.write("#PBS -M sakanashi@ne.phys.sci.osaka-u.ac.jp\n")
    f.write("cd /home/sakra/exp/Bucharest2022/sorter_cor\n")
    f.write("./test "+str(i))
