import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np


def read_file(path):
    dataSet = nc.Dataset(path, "r")
    data = dataSet.variables["data"][:, :]
    return data


def main():
    FOLDER_PATH = "Tests\\Data\\"
    fig = plt.figure(figsize=(4, 4), dpi=80)
    for i in range(1, 6):
        rho_data = np.array(read_file(FOLDER_PATH + "rho_" + str(i) + ".nc"))
        plt.cla()
        plt.imshow(rho_data.T)
        plt.clim(0.8, 2.2)
        ax = plt.gca()
        ax.invert_yaxis()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_aspect('equal')
        plt.pause(0.001)
    plt.savefig("fv.png", dpi=240)
    plt.show()


main()
