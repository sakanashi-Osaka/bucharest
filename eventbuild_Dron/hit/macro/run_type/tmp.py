import numpy as np

runinfo = []
with open("run.txt", "r") as f:
    for line in f:
        line = line.strip()
        line = line.split(" ")

        for num in range(2):
            line[num] = int(line[num])
        runinfo.append(line)
arr = np.array(runinfo) 

current = []
for run in runinfo:
    with open("/home/sakra/exp/Bucharest2022/eventbuild_cor/scaler/log/log"+str(run[0])+".txt", "r") as f1:
        for index,line1 in enumerate(f1):
            line1 = line1.strip()
            if index == 61:
                current.append(int(line1))
                break

arr_c = np.array(current)
arr_c = arr_c.reshape(len(current),1)
arr = np.hstack((arr, arr_c))  # arr[run_number][run_type][current]
            

count = [0,0,0,0]
count_c = [0,0,0,0]
for run in arr:
    count[run[1]] +=1
    count_c[run[1]] += run[2]
    
print(count)
print(count_c)
print(count_c[0]/count[0],count_c[1]/count[1],count_c[2]/count[2])
