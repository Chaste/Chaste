import multiprocessing
import os
import subprocess

import numpy as np
import matplotlib.pyplot as plt

# import immersed_boundary as ib

# Globally accessible directory paths and names
path_to_exec   = os.environ.get('CHASTE_BUILD_DIR') + '/projects/ImmersedBoundary/apps/Exe_ConvergenceWithMesh'
path_to_output = os.environ.get('CHASTE_TEST_OUTPUT') + '/convergence/node_spacing/'

# Range for simulation parameters
num_sims = 4

def main():
    run_simulations()
    combine_output()
    plot_results()


# Create a list of commands and pass them to separate processes
def run_simulations():

    # Make a list of calls to a Chaste executable
    command_list = []
    for sim_id in range(num_sims):

        num_nodes = int(round(2**(5 + 0.5 * sim_id)))

        command = 'nice -n 19 ' + path_to_exec \
                  + ' --ID ' + str(sim_id) \
                  + ' --NN ' + str(num_nodes)
        command_list.append(command)

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Py: Starting simulations with " + str(count) + " processes")

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

    combined_results = open(path_to_output + '/combined_results.dat', 'w')
    combined_results.write("simulation_id,dl_at_end,esf_at_end\n")

    for sim_id in range(num_sims):
        file_name = path_to_output + '/sim/' + str(sim_id) + '/results.dat'

        if os.path.isfile(file_name):
            local_results = open(file_name, 'r')
            combined_results.write(local_results.readline())
            local_results.close()
            combined_results.write("\n")

    combined_results.close()


def plot_results():

    print("Py: Plotting results")

    ################################################################
    # Global options to ensure consistency between different plots #
    ################################################################

    bg_gray = '0.75'
    bg_line_width = 0.5
    font_size = 8

    mycol_gr = '#006532' #Green
    mycol_or = '#F39200' #Orange
    mycol_bl = '#2C2E83' #Blue
    mycol_lb = '#0072bd' #Light blue

    # Useful to convert picas to inches or cm
    # SIAM max fig width = 31pi = 13.12cm = 5.17in
    pica_to_cm = 127.0/300.0
    pica_to_in = 1.0/6.0
    golden_rat = 1.618

    fig_width = 14 * pica_to_in # num picas, but need it in inches

    # Set LaTeX font rendering
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', serif='computer modern roman')
    plt.rc('figure', figsize=[fig_width, fig_width / golden_rat])
    plt.rc('axes', linewidth=bg_line_width, edgecolor=bg_gray, axisbelow=True)
    plt.rc('xtick', labelsize=font_size)
    plt.rc('ytick', labelsize=font_size)
    plt.rc('xtick.major', size=0, pad=4)
    plt.rc('ytick.major', size=0, pad=4)


    # Data for esf vs time plot
    esf_data = np.genfromtxt(path_to_output + '/combined_results.dat',
                             delimiter=',',
                             skip_header=1)

    node_spacing = -1.0 * np.log2(esf_data[:, 1])
    esf_at_end = esf_data[:, 2]

    # Calculate limits to give a comfy margin
    x_min, x_max = min(node_spacing), max(node_spacing)
    y_min, y_max = min(esf_at_end), max(esf_at_end)

    x_range, y_range = x_max - x_min, y_max - y_min
    margin = 0.05

    x_lims = [x_min - margin * x_range, x_max + margin * x_range]
    y_lims = [y_min - margin * y_range, y_max + margin * y_range]

    # Plot the two continuous lines
    plt.plot(node_spacing, esf_at_end, marker='o',
             linestyle='-',
             linewidth=0.75,
             color=mycol_or,
             markerfacecolor='#FFFFFF',
             markeredgecolor=mycol_or,
             markersize=4.0,
             markeredgewidth=0.75)

    # Set the axis limits
    plt.xlim(x_lims)
    plt.ylim(y_lims)

    # Label the axes
    plt.xlabel(r'$-\log_2(\Delta l)$', fontsize=font_size, labelpad=1.0*font_size)
    plt.ylabel(r'ESF', fontsize=font_size, labelpad=1.0*font_size)

    # Customise axes
    ax = plt.gca()
    ax.grid(b=True, which='major', color=bg_gray, linestyle='dotted', dash_capstyle='round')

    # Export figure
    plt.savefig(path_to_output + 'ConvergenceWithNodeSpacing.pdf', bbox_inches='tight', pad_inches=0.0)

if __name__ == "__main__":
    main()
