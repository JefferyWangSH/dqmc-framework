import numpy as np
from matplotlib import pyplot as plt


def readDate(filename):

    T = []
    DoubleOccu = []
    errDoubleOccu = []
    KineticEnergy = []
    errKineticEnergy = []
    StructFactor = []
    errStructFactor = []

    dirtData = {}
    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split()
                data = [float(strList[i+2]) for i in range(5)]
                err = [float(strList[i+7]) for i in range(5)]
                dirtData[1/float(strList[1])] = [data, err]
    file.close()
    
    T = [abs(key) for key in dirtData.keys()]
    T.sort(reverse = False)
    for key in dirtData.keys():
        DoubleOccu.append(dirtData[key][0][0])
        errDoubleOccu.append(dirtData[key][1][0])
        KineticEnergy.append(dirtData[key][0][1])
        errKineticEnergy.append(dirtData[key][1][1])
        StructFactor.append(dirtData[key][0][2])
        errStructFactor.append(dirtData[key][1][2])

    return T, [DoubleOccu, errDoubleOccu], [KineticEnergy, errKineticEnergy], [StructFactor, errStructFactor]


def plotFigure(T, data, label):
    plt.figure()
    plt.grid(linestyle='-.')

    obs, err = data

    # plt.plot(u, obs, 's', ms=4, label=label, ls=":")
    plt.errorbar(T, obs, err, label=label, ms=4, fmt='o', ecolor='r', color='b', elinewidth=1.5, capsize=4)

    plt.xlabel('${T}$\n')
    plt.ylabel(label)

    plt.legend()
    plt.show()


if __name__ == "__main__":

    T, doubleoccu, kineticenergy, structfactor = readDate("eq-u-4.txt")

    plotFigure(T, doubleoccu, 'Double Occu')
    # plotFigure(T, kineticenergy, 'Kinetic Energy')
    # plotFigure(T, structfactor, "Struct Factor")
