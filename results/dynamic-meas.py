import numpy as np
from matplotlib import pyplot as plt


def readDate(filename):

    tau = []
    Matsubara = []
    errMatsubara = []

    dirtData = {}
    n = 0
    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0 | n != 0:
                strList = line.split()
                dirtData[float(strList[0])] = [float(strList[1]), float(strList[2])]
            n = n + 1
    file.close()
    
    tau = list(dirtData.keys())
    tau.sort(reverse = False)
    for key in dirtData.keys():
        Matsubara.append(dirtData[key][0])
        errMatsubara.append(dirtData[key][1])

    return tau, [Matsubara, errMatsubara]


def plotFigure(u, data):
    plt.figure()
    plt.grid(linestyle='-.')

    obs, err = data

    # plt.plot(u, obs, 's', ms=4, label=label, ls=":")
    plt.errorbar(u, obs, err, label='${G(\\tau)}$', ms=4, fmt='o', ecolor='r', color='b', elinewidth=1.5, capsize=4)

    plt.xlabel('${\\tau / \\Delta\\tau}$\n')
    plt.ylabel('${G(\\tau)}$')
    plt.title('${k = (\\pi/2, \\pi/2)}$')

    plt.legend()
    plt.show()


if __name__ == "__main__":

    tau, matsubara = readDate("dynamic-u4-attractive.txt")

    plotFigure(tau, matsubara)
