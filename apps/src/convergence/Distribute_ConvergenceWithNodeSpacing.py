import multiprocessing
import os
import subprocess
import time

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck


# Globally accessible directory paths, names, and variables
chaste_build_dir = os.environ.get('CHASTE_BUILD_DIR')
executable = os.path.join(chaste_build_dir, 'projects/ImmersedBoundary/apps', 'Exe_ConvergenceWithNodeSpacing')

if not(os.path.isfile(executable)):
    raise Exception('Py: Could not find executable: ' + executable)

chaste_test_dir = os.environ.get('CHASTE_TEST_OUTPUT')
path_to_output = os.path.join(chaste_test_dir, 'convergence', 'num_nodes')

pdf_name = os.path.join(path_to_output, 'ConvergenceWithNodeSpacing.pdf')
results_header_string = 'simulation_id,dl_at_start,esf_at_end\n'

num_sims = 16


def main():
    run_simulations()
    combine_output()
    plot_results()
    compress_output()


# Create a list of commands and pass them to separate processes
def run_simulations():

    # 2^7 (128) up to 2^12 (4096) at third-integer powers, rounded to integer value
    param_values = (2**(np.linspace(7, 12, num=num_sims))).astype(int)

    # Make a list of calls to a Chaste executable
    command_list = []

    for idx, param_val in enumerate(param_values):

        command = 'nice -n 19 ' + executable \
                  + ' --ID ' + str(idx) \
                  + ' --NN ' + str(param_val)
        command_list.append(command)

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Py: Starting simulations with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool

    pool.map(execute_command, reversed(command_list))


# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)


# Gather the output from all simulations and put it in the same file
def combine_output():

    print("Py: Combining output")

    if not(os.path.isdir(path_to_output)):
        raise Exception('Py: Could not find output directory: ' + path_to_output)

    combined_results = open(os.path.join(path_to_output, 'combined_results.dat'), 'w')
    combined_results.write(results_header_string)

    for sim_id in range(num_sims):
        file_name = os.path.join(path_to_output, 'sim', str(sim_id), 'results.dat')

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

    # SIAM max fig width = 31pi = 13.12cm = 5.17in

    # latex line width in inches
    latex_line_width = 5.126

    # Potentially useful pica conversions
    pica_to_cm = 127.0/300.0
    pica_to_in = 1.0/6.0

    golden_rat = 1.618

    fig_width = 0.45 * latex_line_width

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
    esf_data = np.genfromtxt(os.path.join(path_to_output, 'combined_results.dat'),
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
             markerfacecolor=mycol_lb,
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
    # ax.get_yaxis().get_major_formatter()(False)

    y_formatter = tck.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)

    ax.grid(b=True, which='major', color=bg_gray, linestyle='dotted', dash_capstyle='round')

    # Export figure
    plt.savefig(pdf_name, bbox_inches='tight', pad_inches=0.0)


# Compress output and suffix with date run
def compress_output():

    # Check that output directory exists
    if not (os.path.isdir(path_to_output)):
        raise Exception('Py: Could not find output directory: ' + path_to_output)

    # Change cwd to one above output path
    os.chdir(os.path.join(path_to_output, '..'))

    simulation_name = os.path.basename(os.path.normpath(path_to_output))

    today = time.strftime('%Y-%m-%d')

    # Compress entire output folder and append with the the simulations were started
    os.system('tar -zcf ' + simulation_name + '_' + today + '.tar.gz ' + simulation_name)

if __name__ == "__main__":
    main()
