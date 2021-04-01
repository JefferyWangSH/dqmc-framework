import numpy as np
from matplotlib import pyplot as plt


def readDate(filename):

    u = []
    DoubleOccu = []
    KineticEnergy = []
    StructFactor = []

    dirtData = {}
    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split()
                data = [float(strList[i+1]) for i in range(len(strList)-1)]
                dirtData[float(strList[0])] = data
    file.close()
    
    u = list(dirtData.keys())
    u.sort(reverse = False)
    for key in dirtData.keys():
        DoubleOccu.append(dirtData[key][0])
        KineticEnergy.append(dirtData[key][1])
        StructFactor.append(dirtData[key][2])

    return u, DoubleOccu, KineticEnergy, StructFactor


def plotFigure(u, obs, label):
    plt.figure()
    plt.grid(linestyle='-.')

    plt.plot(u, obs, 's', ms=4, label=label, ls=":")

    plt.xlabel('${U/t}$\n')
    plt.ylabel(label)

    plt.legend()
    plt.show()


if __name__ == "__main__":

    u, doubleoccu, kineticenergy, structfactor = readDate("meas-eqtime.txt")

    plotFigure(u, doubleoccu, 'Double Occu')
    plotFigure(u, kineticenergy, 'Kinetic Energy')
    plotFigure(u, structfactor, "Struct Factor")
