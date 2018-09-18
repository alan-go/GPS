import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def Plt2(data,indx,indy,figNme):
    x,y = [],[]
    for each in data:
        x.append(each[indx])
        y.append(each[indy])
    plt.figure(figNme)
    # plt.scatter(x,y)
    plt.plot(x,y)

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

def Test():
    a=0


def AnaData(path):
    file  = open(path)
    txyz, tlla = [],[]
    svDatas, svNames = [],[]
    for each in file.readlines():
        if 'xyz' in each:
            txyz0 = each.split(',')
            time = float(txyz0[1])
            x = float(txyz0[3])
            y = float(txyz0[4])
            z = float(txyz0[5])
            txyz.append([time,x,y,z])
        if 'LLA' in each:
            tlla0 = each.split(',')
            tlla.append([float(tlla0[1]),float(tlla0[3]),float(tlla0[4]),float(tlla0[5])])


        if 'svs' in each:
            color = 'g'
            if 'Sp3'in each:
                color = 'b'
            if 'Brt'in each:
                color = 'r'
            sv = each.split(',')
            data = [float(sv[1]),float(sv[4]),float(sv[6]), float(sv[8]),float(sv[9]),float(sv[10]),float(sv[11]),color]   # time,pr,prres,norm,x,y,z

            if not sv[2] in svNames:
                svNames.append(sv[2])
                svDatas.append([])
                print("add sv",sv[2],"time =", sv[1])

            ind = svNames.index(sv[2])
            svDatas[ind].append(data)


    Plt2(txyz,0,1,'pos')
    Plt2(tlla,0,3,'height')

    #####################
    # fig = plt.figure('svs')
    # ax = plt.subplot(111, projection='3d')
    for i,eachsv in enumerate(svDatas):
        fig = plt.figure(svNames[i])
        ax = fig.add_subplot(111, projection='3d')
        Plt3(eachsv,4,5,6,7,ax)
        # Plt2(eachsv,0,1,svNames[i])
    #####################

    plt.show()

if __name__=="__main__":
    AnaData("../log/log.txt")
    Test()