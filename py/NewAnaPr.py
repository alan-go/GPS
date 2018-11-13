import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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


if __name__=="__main__":
    # ReadData("../log/logtu.txt",2,121)
    # ReadData("../log/xyzOf1030_08_44RAC.txt",4,1)
    # ReadData("../log/xyzOf1030_08_44SIG.txt",4,1)
    # ReadData("../log/xyzOf1030_08_44RTK.txt",4,1)
    # ReadData("../log/xyzOf1030_08_44KAL.txt",4,1)
    # ReadData("../log/SV/0_22.txt",13,1)
    ReadData("../log/SV/3_30.txt",13,1)
    plt.figure("ana")
    # plt.scatter(data[0],data[4])
    # plt.scatter(data[0],data[5])
    plt.plot(data[0],data[6])
    # plt.scatter(data[0],data[7])
    plt.plot(data[0],data[8])
    plt.figure("ana2")
    plt.scatter(data[0],data[7])



    # plt.scatter(data[0],data[9])
    # plt.scatter(data[0],data[10])
    # plt.scatter(data[0],data[11])

    # plt.plot(data[0],data[1],"r")
    # plt.plot(data[0],data[2],"g")
    # plt.plot(data[0],data[4],"b")
    # plt.scatter(data[0],data[3])
    # plt.scatter(data[0],data[6])
    # plt.plot(data[2],data[1])
    # plt.plot(data[0],data[2])
    # plt.plot(data[0],data[3])
    # plt.plot(data[0],data[4])

    # plt.plot(data[0],data[5])
    # plt.plot(data[0],data[6])
    # plt.plot(data[0],data[10])


    plt.show()

