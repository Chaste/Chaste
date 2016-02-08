import multiprocessing
import subprocess
import os
import shutil

# Globally accessible directory paths and names
path_to_exec   = '/scratch/Cooper/chaste_build/projects/ImmersedBoundary/apps/'#os.environ.get('CHASTE_BUILD_DIR') + '/projects/ImmersedBoundary/apps/'
path_to_output = '/scratch/Cooper/chaste_test_output/'#os.environ.get('CHASTE_TEST_OUTPUT') + '/'
sim_name = 'CornerTurning'
exec_name = 'Exe_' + sim_name

# Range for two simulation parameters
num_global_consts = 4
num_local_consts  = 15
num_kicks_per_sim = 6

def main():
    #run_simulations()
    #combine_output()
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
    combined_results.write("global_const,local_const,kick,asymmetry\n")

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

    print("Py: Plotting results")

    import numpy as np
    import matplotlib.pyplot as plt

    # Set LaTeX font rendering
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    my_data = np.genfromtxt(path_to_output + exec_name + '/combined_results.dat', delimiter=',', skip_header=1)

    x_vals = my_data[:,1] # local const
    y_vals = my_data[:,3] # summary statistic

    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)

    x_range, y_range = x_max - x_min, y_max - y_min
    margin = 0.05;

    x_lims = [x_min - margin * x_range, x_max + margin * x_range]
    y_lims = [y_min - margin * y_range, y_max + margin * y_range]

    num_sims_per_plot = num_local_consts * num_kicks_per_sim

    for plot in range(1):

        plt.clf()

        plt.xlabel(r'Local spring constant multiple', fontsize=12)
        plt.ylabel(r'Proportion of right-to-left asymmetry', fontsize=12)

        title = 'Global cell-cell spring constant scaled by ' + str(0.1 * (1+plot))

        plt.title(title, fontsize=12)

        # highlight_x_vals = my_data[plot * num_sims_per_plot : (plot+1) * num_sims_per_plot,1]
        # highlight_y_vals = my_data[plot * num_sims_per_plot : (plot+1) * num_sims_per_plot,3]

        # plt.scatter(x_vals, y_vals, color = '0.75')
        plt.scatter(x_vals, y_vals, color = '#F39200')

        plt.xlim(x_lims)
        plt.ylim(y_lims)

        plt.savefig(path_to_output + exec_name + '/Fig_' + sim_name + str(plot) + '.pdf')

if __name__ == "__main__":
    main()
