import os
import numpy as np
import matplotlib.pyplot as plt


def plot_cell_size(simulation_name):
    dir_name = simulation_name+"/results_from_time_0"
    os.chdir(dir_name)
    if "ib_cell_size.dat" in os.listdir():
        data = np.loadtxt("ib_cell_size.dat")
        time = data[:, 0]
        cell_size = data[:, 5::5].transpose()
        plt.figure(figsize=(7, 4))
        plt.xlabel("Time")
        plt.ylabel("Cell size")
        for i in range(cell_size.shape[0]):
            plt.plot(time, cell_size[i])
        plt.savefig("cell_size.png")

    else:
        print("ib_cell_size.dat not found the specified directory!")
        quit()


if __name__ == '__main__':
    os.chdir("../testoutput/")
    sims = os.listdir()
    sims.sort()
    print("Pick simulation:")
    for index in range(len(sims)):
        print(sims[index]+" ({})".format(index))
    sim_index = -1
    while 0 > sim_index or sim_index > len(sims):
        sim_index = int(input("Simulation index: "))
    plot_cell_size(sims[sim_index])
