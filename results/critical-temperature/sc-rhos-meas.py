import numpy as np
from matplotlib import pyplot as plt

def analyseDataFromBins(filename):

    nn = 0
    tmp_rho_s = 0.0
    tmp_rho_s_2 = 0.0
    mean_rho_s = 0.0
    err_rho_s = 0.0

    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split()
                tmp_rho_s += float(strList[1])
                tmp_rho_s_2 += float(strList[1])**2
                nn += 1
    file.close()

    mean_rho_s = tmp_rho_s / nn
    err_rho_s = (tmp_rho_s_2 /nn - mean_rho_s**2)**0.5 / (nn -1)**0.5
    
    return mean_rho_s, err_rho_s


def analyseDataFromMeans(infile, outfile):

    data = {}

    with open(infile, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split()

                if ( float(strList[0]) not in data.keys()):
                    # data structure: { beta: [num, mean, error] }
                    data[float(strList[0])] = [1.0, float(strList[2]), float(strList[2])**2]  
                else:
                    data[float(strList[0])][0] += 1.0
                    data[float(strList[0])][1] += float(strList[2])
                    data[float(strList[0])][2] += float(strList[2])**2
    file.close()


    for beta in data.keys():
        data[beta][1] =  data[beta][1] / data[beta][0]
        if ( abs(data[beta][0] - 1) < 10**-3 ):
            data[beta][2] = 0
        else:
            data[beta][2] =  ( data[beta][2] / data[beta][0] - data[beta][1]**2 )**0.5 / (data[beta][0] - 1)**0.5

    # output
    with open(outfile, mode='w') as file:
        list_beta = list(data.keys())
        list_beta.sort(reverse=False)
        for beta in list_beta:
            file.write('{0:>15} {1:>15.5f} {2:>15.8f} {3:>15.8f} \n'.format(beta, 1/beta, data[beta][1], data[beta][2]))

    return data


def readDate(filename):

    temp = []
    rhos = []
    errRhos = []

    dirtData = {}
    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split()
                data = float(strList[2])
                err = float(strList[3])
                dirtData[float(strList[1])] = [data, err]
    file.close()
    
    temp = [key for key in dirtData.keys()]
    temp.sort(reverse = False)
    for key in temp:
        rhos.append(dirtData[key][0])
        errRhos.append(dirtData[key][1])
    
    return temp, [rhos, errRhos]


def plotFigure(x, data, label, title):
    plt.figure()
    plt.grid(linestyle='-.')

    obs, err = data

    plt.errorbar(x, obs, err, label=label, ms=4, fmt='o', ecolor='r', color='b', elinewidth=1.5, capsize=4)

    # bench line to determine the critical temperature of BKT transition 
    # benchmarkTemp = np.arange(0, 0.5, 50)
    benchmarkLine = [ 2 * x / np.pi for x in x]
    plt.plot(x, benchmarkLine, ms=2, label="${2T/\\pi}$", ls="dashed")

    plt.xlabel('${T}$\n')
    plt.ylabel(label)
    plt.title(title)

    plt.legend()
    plt.show()


if __name__ == "__main__":

    data = analyseDataFromMeans(infile="../rhos_L_6/sc_rhos_L_6_u_-4.0.txt", outfile="../rhos_L_6/rhosForPlot_L_6_u_-4.0.txt")
    data = analyseDataFromMeans(infile="../rhos_L_8/sc_rhos_L_8_u_-4.0.txt", outfile="../rhos_L_8/rhosForPlot_L_8_u_-4.0.txt")

    # data_mean, data_err = analyseDataFromBins("bins_L_8_u_-4.0_b_6.0.txt")
    # print(data_mean, data_err)

    # temp, rhos = readData("...")
    # plotFigure(x=temp, data=rhos, label="${\\rho_{s}}$", 
    #             title="${ \\langle n \\rangle = 0.5 \\quad L = 6 \\quad u = - 4.0}$")

    temp_L_6, rhos_L_6 = readDate(filename="../rhos_L_6/rhosForPlot_L_6_u_-4.0.txt")
    temp_L_8, rhos_L_8 = readDate(filename="../rhos_L_8/rhosForPlot_L_8_u_-4.0.txt")

    plt.figure()

    # error-bar plot
    # color: lightcoral, orangered, blue, royalbl
    plt.errorbar(temp_L_6, rhos_L_6[0], rhos_L_6[1], label="${ L = 6 }$", ms=5, fmt='^:', ecolor='lightcoral', color='orangered', elinewidth=1.5, capsize=4)
    plt.errorbar(temp_L_8, rhos_L_8[0], rhos_L_8[1], label="${ L = 8 }$", ms=4, fmt='o:', ecolor='blue', color='royalblue', elinewidth=1.5, capsize=4)
    
    # benchmark line 
    benchmarkTemp = [ 1/16, 1/4]
    benchmarkRhos = [ 2 * temp / np.pi for temp in benchmarkTemp]
    plt.plot(benchmarkTemp, benchmarkRhos, ms=2, alpha=0.85, label="${ 2T / \\pi }$", ls="-")

    plt.title("${ \\langle n \\rangle = 0.5 \\quad u = - 4.0}$")
    plt.xlabel("${ T }$ \n")
    plt.ylabel("${ \\rho_{s} }$")

    plt.legend()
    plt.grid(linestyle='-.')
    plt.show()
