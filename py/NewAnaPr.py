import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

data=[]

def Data_zero():
    for eachCol in data:
        dfirst = eachCol[0]
        for eachEle in eachCol:
            eachCol-=dfirst

def Data_proc():
    for eachClo in data:
        for eachEle in eachClo:
            eachEle=-eachEle

def Plt3(data,indx,indy,indz,color,ax):
    x0,y0,z0 = [],[],[]
    x1,y1,z1 = [],[],[]
    # x.append(0)
    # y.append(0)
    # z.append(0)
    for each in data:
        if each[color]=='b':
            x0.append(each[indx])
            y0.append(each[indy])
            z0.append(each[indz])
        if each[color]=='r':
            x1.append(each[indx])
            y1.append(each[indy])
            z1.append(each[indz])
    ax.scatter(x0,y0,z0,c='b')
    ax.scatter(x1,y1,z1,c='r')





def ReadData(path, colNum,begin):
    file  = open(path)
    txyz, tlla = [],[]
    for ind in range(colNum):
        data.append([])
    k=0
    for each in file.readlines():
        k+=1
        if(k<begin):continue
        dataN = each.split(',')
        if(len(dataN)<colNum):
            continue

        for ind in range(colNum):
            data[ind].append(float(dataN[ind]))

def AnaSv(file,colNum,begin):
    ReadData(file,colNum,begin)
    plt.figure("anaPr")
    plt.scatter(data[0], data[1])
    plt.scatter(data[0], data[2])

    plt.figure("anadiff")
    plt.scatter(data[0], data[3])
    plt.scatter(data[0], data[4])

    # plt.figure("anaAdd")
    # plt.plot(data[0],data[6])
    # plt.plot(data[0],data[7])
    # plt.figure("anastd")
    # plt.scatter(data[0],data[8])
    # plt.scatter(data[0],data[9])

    plt.show()

def Anaxyz(file,colNum,begin):
    ReadData(file,colNum,begin)
    # plt.figure("anaXyz")
    # plt.scatter(data[0], data[1])
    # plt.scatter(data[0], data[2])

    plt.plot(data[2], data[1])
    data.clear()

def Anacol(file,width, colX,colY, begin):
    ReadData(file, width, begin)

    # plt.plot(data[colX], data[colY])
    plt.scatter(data[colX], data[colY])
    data.clear()



if __name__=="__main__":
    # ss = "../log/1215_07_56xyz"
    ss = "../log/1219_07_41xyz"
    width = 4;

    # plt.figure(ss)
    # # Anaxyz(ss+"RAC.txt",width,10)
    # Anaxyz(ss+"RTK.txt",4,10)
    # # Anaxyz(ss+"SIG.txt",4,10000)
    #
    # Anaxyz(ss+"UBX.txt",width,10)
    # Anaxyz(ss+"KAL2.txt",width,10)
    # Anaxyz(ss+"NVT.txt",width,10)


    # Anaxyz(ss+"KAL.txt",width,10)
    # AnaSv("../log/SV/3_30.txt",5,1)
    #


    file_dir="../log/SV/"
    width = 7;
    plt.figure(ss)

    # for root, dirs, files in os.walk(file_dir):
    #     print(root) #当前目录路径
    #     # print(dirs) #当前路径下所有子目录
    #     print(files) #当前路径下所有非目录子文件
    #
    #     for file in files:
    #         ss=root+file
    #         Anacol(ss,width,0,1,10);
    #
    ss = "/home/alan/projects/GPS/log/logDebug.txt"
    # Anacol(ss,width,0,1,1);

    #
    ss = "../log/SV/3_06.txt"
    Anacol(ss,width,0,1,16);
    Anacol(ss,width,0,2,16);
    Anacol(ss,width,0,3,16);
    Anacol(ss,width,0,4,16);
    Anacol(ss,width,0,5,16);
    # Anacol(ss,width,0,6,16);


    plt.show()






