import os
import shutil

def analizeFile(sourceFile,resultFile):#analize the source file and save the information in the result file
    f=open(sourceFile,"r")#open file and read all its lines
    lines=f.readlines()
    data={}
    num={}
    for line in lines:#for each line get the informations
        type,omp,mpi,creationTime,communicationTime,executionTime,totalTime=line.strip()[:-1].split(";")
        if type not in data:
            data[type]={}
            num[type]={}
        if omp not in data[type]:
            data[type][omp]={}
            num[type][omp]={}
        if mpi not in data[type][omp]:
            data[type][omp][mpi]=[float(creationTime),float(communicationTime),float(executionTime),float(totalTime)]
            num[type][omp][mpi]=1
        else:
            data[type][omp][mpi][0]+=float(creationTime)
            data[type][omp][mpi][1]+=float(communicationTime)
            data[type][omp][mpi][2]+=float(executionTime)
            data[type][omp][mpi][3]+=float(totalTime)
            num[type][omp][mpi]+=1
    f.close()#close the source file and open the result file
    f=open(resultFile,"w")
    
    f.write("Modality;OMP;MPI;tree_creation_time;communication_time;rb_search_time;program_execution_time;speedup;efficiency\n")
    for typeKey in data.keys():#for each element in the dictionary calculate the mean and save it in the result file
        ompList=list(data[typeKey].keys())
        ompList.sort(key=int)
        for ompKey in ompList:
            mpiList=list(data[typeKey][ompKey].keys())
            mpiList.sort()
            for mpiKey in mpiList:
                it=num[typeKey][ompKey][mpiKey]
                data[typeKey][ompKey][mpiKey][0]/=it
                data[typeKey][ompKey][mpiKey][1]/=it
                data[typeKey][ompKey][mpiKey][2]/=it
                data[typeKey][ompKey][mpiKey][3] /= it
                t_serial = data["Sequential"]["0"]["0"][3]
                t_parallel_single = data["OMP+MPI"]["1"]["1"][3]
                speedup = t_parallel_single / data[typeKey][ompKey][mpiKey][3]
                if typeKey != "Sequential":
                    efficiency = t_parallel_single/(int(ompKey)*int(mpiKey)*data[typeKey][ompKey][mpiKey][3])*100
                else:
                    t_parallel_single = data["OMP+MPI"]["1"]["1"][3]/it
                    efficiency = (t_parallel_single/t_serial)*100
                formatted_line = (
                    f"{typeKey};{ompKey};{mpiKey};"
                    f"{data[typeKey][ompKey][mpiKey][0]:.8f};"
                    f"{data[typeKey][ompKey][mpiKey][1]:.8f};"
                    f"{data[typeKey][ompKey][mpiKey][2]:.8f};"
                    f"{data[typeKey][ompKey][mpiKey][3]:.8f};"
                    f"{speedup:.8f};{efficiency:.4f}\n"
                )
                f.write(formatted_line)
    f.close()#close the result file

for path, currentDirectory, files in os.walk("Informations/OMP_MPI"):#for each file in the directory Informations, analize them and save the results in the relative path in Results directory
    files = sorted(files)
    for file in files:
        resultPath=path.replace("Informations","Results")
        if not(os.path.exists(resultPath)):
            os.makedirs(resultPath)
        analizeFile(os.path.join(path, file),os.path.join(resultPath,file))            