import numpy as np
from matplotlib import pyplot as plt


def readDate(filename):

    u = []
    DoubleOccu = []
    errDoubleOccu = []
    KineticEnergy = []
    errKineticEnergy = []
    StructFactor = []
    errStructFactor = []
    localCorr = []
    errlocalCorr = []

    dirtData = {}
    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split()
                data = [float(strList[i+2]) for i in range(5)]
                err = [float(strList[i+7]) for i in range(5)]
                dirtData[float(strList[0])] = [data, err]
    file.close()
    
    u = [abs(key) for key in dirtData.keys()]
    u.sort(reverse = False)
    for key in dirtData.keys():
        DoubleOccu.append(dirtData[key][0][0])
        errDoubleOccu.append(dirtData[key][1][0])
        KineticEnergy.append(dirtData[key][0][1])
        errKineticEnergy.append(dirtData[key][1][1])
        StructFactor.append(dirtData[key][0][2])
        errStructFactor.append(dirtData[key][1][2])
        localCorr.append(dirtData[key][0][4])
        errlocalCorr.append(dirtData[key][1][4])

    return u, [DoubleOccu, errDoubleOccu], [KineticEnergy, errKineticEnergy], [StructFactor, errStructFactor], [localCorr, errlocalCorr]


def plotFigure(u, data, label):
    plt.figure()
    plt.grid(linestyle='-.')

    obs, err = data

    # plt.plot(u, obs, 's', ms=4, label=label, ls=":")
    plt.errorbar(u, obs, err, label=label, ms=4, fmt='o:', elinewidth=1.5, capsize=4)

    plt.xlabel('${U/t}$\n', fontsize = 13)
    plt.ylabel(label, fontsize = 13)

    plt.legend(fontsize = 12)
    plt.show()


if __name__ == "__main__":

    u_b4, doubleoccu_b4, kineticenergy_b4, structfactor_b4, localCorr_b4 = readDate("eqtime_L4_beta4_repulsive.txt")
    u_b6, doubleoccu_b6, kineticenergy_b6, structfactor_b6, localCorr_b6 = readDate("eqtime_L4_beta6_repulsive.txt")
    u_b8, doubleoccu_b8, kineticenergy_b8, structfactor_b8, localCorr_b8 = readDate("eqtime_L4_beta8_repulsive.txt")
    u_b12, doubleoccu_b12, kineticenergy_b12, structfactor_b12, localCorr_b12 = readDate("eqtime_L4_beta12_repulsive.txt")

    plt.figure()
    plt.grid(linestyle='-.')

    mean_effective_t_b4 = [k/kineticenergy_b4[0][0] for k in kineticenergy_b4[0]]
    err_effective_t_b4 = [k/kineticenergy_b4[0][0] for k in kineticenergy_b4[1]]
    mean_effective_t_b6 = [k/kineticenergy_b6[0][0] for k in kineticenergy_b6[0]]
    err_effective_t_b6 = [k/kineticenergy_b6[0][0] for k in kineticenergy_b6[1]]
    mean_effective_t_b8 = [k/kineticenergy_b8[0][0] for k in kineticenergy_b8[0]]
    err_effective_t_b8 = [k/kineticenergy_b8[0][0] for k in kineticenergy_b8[1]]
    mean_effective_t_b12 = [k/kineticenergy_b12[0][0] for k in kineticenergy_b12[0]]
    err_effective_t_b12 = [k/kineticenergy_b12[0][0] for k in kineticenergy_b12[1]]

    effective_t_b4 = [mean_effective_t_b4, err_effective_t_b4]
    effective_t_b6 = [mean_effective_t_b6, err_effective_t_b6]
    effective_t_b8 = [mean_effective_t_b8, err_effective_t_b8]
    effective_t_b12 = [mean_effective_t_b12, err_effective_t_b12]

    plt.errorbar(u_b4, structfactor_b4[0], structfactor_b4[1], label="${\\beta = 4.0}$", ms=4, fmt='o:', elinewidth=1.5, capsize=4)
    plt.errorbar(u_b6, structfactor_b6[0], structfactor_b6[1], label="${\\beta = 6.0}$", ms=4, fmt='s:', elinewidth=1.5, capsize=4)
    plt.errorbar(u_b8, structfactor_b8[0], structfactor_b8[1], label="${\\beta = 8.0}$", ms=4, fmt='^:', elinewidth=1.5, capsize=4)
    # plt.errorbar(u_b12, structfactor_b12[0], structfactor_b12[1], label="${\\beta = 12.0}$", ms=4, fmt='D:', elinewidth=1.5, capsize=4)

    plt.xlabel('${U/t}$\n', fontsize = 16)
    plt.ylabel("${S^{zz}(\\pi,\\pi)}$ / ${N}$", fontsize = 16)


    plt.legend(fontsize = 13)

    plt.savefig("structure_factor_repulsive.pdf")
    plt.show()
