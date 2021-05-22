import numpy as np
from matplotlib import pyplot as plt


def readDate(filename):

    l = 0
    data_dirt = {}
    
    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split()
                tuple_key = (int(strList[0]), int(strList[1]))
                data_dirt[tuple_key] = float(strList[4])
    file.close()
    
    for key in data_dirt.keys():
        l = max(max(key[0], key[1]), l)
    
    arr = np.zeros((l+1, l+1))
    
    for i in range(l+1):
        for j in range(l+1):
            if (i, j) in data_dirt.keys():
                arr[i][j] = data_dirt[(i, j)]
            
            elif (l - i, l - j) in data_dirt.keys():
                # time reverse symmetry G(k, beta/2) = G(-k, beta/2)
                arr[i][j] = data_dirt[(l - i, l - j)]
            else:
                # filled with 0.0
                data_dirt[(i, j)] = 0.0
            

    return arr


def plotFigure(arr, figname, interpolation, title):

    plt.imshow(arr, interpolation=interpolation, origin='lower', vmin=0, vmax=0.5, cmap='viridis')
    # interpolation: none, nearest, bilinear, bicubic, spline16, spline36, hanning, hamming, hermite,
    #                kaiser, quadric, catrom, gaussian, bessel, mitchell, sinc, lanczos.

    plt.colorbar()
    plt.xticks(())
    plt.yticks(())
    plt.title(title, y=-0.22, fontsize=20)

    plt.savefig(figname)
    plt.show()


if __name__ == "__main__":

    arr = readDate('fermi_surface_L8_beta4_u-1.0.txt')

    plotFigure(arr, figname="fermi_surface_L8_b4_u-1.0.pdf", interpolation="none", title="${\\mathrm{U} = -1.0}$\n")
    plotFigure(arr, figname="fermi_surface_L8_b4_u-1.0_ip.pdf", interpolation="hamming", title="${\\mathrm{U} = -1.0}$\n")