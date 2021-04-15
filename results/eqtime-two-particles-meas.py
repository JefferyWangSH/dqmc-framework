import numpy as np
from matplotlib import pyplot as plt


def readDate(filename):

    px = []
    MomentumDist = []

    dirtData = {}
    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split() 
                dirtData[float(strList[7])] = float(strList[5])
    file.close()
    
    px = list(dirtData.keys())
    px.sort(reverse = False)
    for key in dirtData.keys():
        MomentumDist.append(dirtData[key])
        
    return px, MomentumDist


def plotFigure(px, obs, label):
    plt.figure()
    plt.grid(linestyle='-.')

    plt.plot(px, obs, 's', ms=4, label=label, ls=":")

    plt.xlabel('$( {p}_{x}, {p}_{y} )$\n')
    plt.ylabel(label)

    plt.legend()
    plt.show()


if __name__ == "__main__":

    px, MomentumDist = readDate("eq-u-4.txt")

    plotFigure(px, MomentumDist, '${n}_{k}$')
