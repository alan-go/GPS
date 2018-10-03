import numpy as np
import matplotlib.pyplot as plt
import math

def ana(data):
    sum2 = []
    mean = np.average(data)
    for each in data:
        sum2.append((each-mean)*(each-mean))
    sigma2 = np.average(sum2)
    sigma = np.sqrt(sigma2)
    print("sigma = ",sigma2,sigma)



if __name__=="__main__":
    path = "../log/logDebug.txt"
    print("read",path)

    file = open(path)
    data = []

    for each in file.readlines():
        temp = each.split(',')
        data.append(float(temp[2]))
    print(len(data))
    ana(data)
    x = range(len(data))
    plt.plot(x,data)
    plt.show()



