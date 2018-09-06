import numpy as np
import matplotlib.pyplot as plt
import math

def RMSE(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())
def Analysis(path):
    file = open(path)
    listx,listy,listz = [],[],[]
    for each in file.readlines():
        if 'xyz' in each:
            position = each.split(',')
            listx.append(float(position[3]))
            listy.append(float(position[4]))
            listz.append(float(position[5]))
    targetsx = [-2171471.444215721 for i in range(len(listx))]
    targetsy = [4386134.629738129 for i in range(len(listy))]
    targetsz = [4076258.730295825 for i in range(len(listz))]
    RMSE_X = RMSE(np.array(targetsx),np.array(listx))
    RMSE_Y = RMSE(np.array(targetsy),np.array(listy))
    RMSE_Z = RMSE(np.array(targetsz),np.array(listz))
    print("RMSE_X: ",RMSE_X,"RMSE_Y: ",RMSE_Y,"RMSE_Z: ",RMSE_Z,)

    plt.plot(listx)
    plt.show()

if __name__=="__main__":
    path = "../log/log.txt"
    Analysis(path)
    # Analysis("../log/logbg.txt")


