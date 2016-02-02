import multiprocessing
import subprocess
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

# Globally accessible directory paths and names
path_to_exec   = os.environ.get('CHASTE_BUILD_DIR') + '/projects/ImmersedBoundary/apps/'
path_to_output = os.environ.get('CHASTE_TEST_OUTPUT') + '/'
exec_name = 'Exe_VaryAdhesionOfSingleCell'

# Range for two simulation parameters
num_global_consts = 5
num_local_consts  = 11
num_kicks_per_sim = 11

def main():
    run_simulations()
    combine_output()
    plot_results()


# Create a list of commands and pass them to separate processes
def run_simulations():

    # if os.path.isdir(path_to_output + exec_name):
    #     print("Py: Deleting previous output")
    #     shutil.rmtree(path_to_output + exec_name)

    print("Py: Starting simulations")

    # Make a list of calls to a Chaste executable
    command_list = []
    for global_const in range(num_global_consts):
        for local_const in range(num_local_consts):
            for kick in range(num_kicks_per_sim):
                command = path_to_exec + exec_name \
                          + ' --G ' + str(global_const) \
                          + ' --L ' + str(local_const) \
                          + ' --K ' + str(kick)
                command_list.append(command)

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool

    pool.map(execute_command, command_list)



# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)


# Gather the output from all simulations and put it in the same file
def combine_output():

    print("Py: Combining output")

    combined_results = open(path_to_output + exec_name + '/combined_results.dat', 'w')
    combined_results.write("global_const,local_const,kick,change_in_centroid\n")

    for global_const in range(num_global_consts):
        for local_const in range(num_local_consts):
            for kick in range(num_kicks_per_sim):
                file_name = path_to_output + exec_name + '/' \
                                           + str(global_const) + '_' \
                                           + str(local_const) + '_' \
                                           + str(kick) + '/results.dat'

                if os.path.isfile(file_name):
                    local_results = open(file_name, 'r')
                    combined_results.write(local_results.readline())
                    local_results.close()
                    combined_results.write("\n")

    combined_results.close()



def plot_results():

    print("Py: Plotting results coming soon")

    from numpy import genfromtxt
    my_data = genfromtxt(path_to_output + exec_name + '/combined_results.dat', delimiter=',', skip_header=1)

    # x_vals = my_data[:,0]
    # y_vals = my_data[:,1]
    #
    # plt.scatter(x_vals, y_vals)
    # plt.savefig(path_to_output + exec_name + '/figure.pdf')

if __name__ == "__main__":
    main()