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

    plt.figure("file")
    plt.plot(data[2], data[1])


    data.clear()


if __name__=="__main__":
    # AnaSv("../log/SV/3_30.txt",5,1)
    #
    # Anaxyz("../log/xyzSIG.txt",4,1)
    # Anaxyz("../log/xyzSIGNew.txt",4,1)
    Anaxyz("../log/xyzRTK.txt",4,1)
    Anaxyz("../log/xyzKAL.txt",4,1)
    # Anaxyz("../log/xyzBX.txt",4,1)
    Anaxyz("../log/xyzRAC.txt",4,1)
    plt.show()






