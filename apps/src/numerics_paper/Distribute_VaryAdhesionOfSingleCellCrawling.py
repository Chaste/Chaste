import multiprocessing
import os
import subprocess

import numpy as np
import matplotlib.pyplot as plt

# import immersed_boundary as ib

# Globally accessible directory paths and names
path_to_exec   = os.environ.get('CHASTE_BUILD_DIR') + '/projects/ImmersedBoundary/apps/'
path_to_output = os.environ.get('CHASTE_TEST_OUTPUT') + '/numerics_paper/'
sim_name = 'VaryAdhesionOfSingleCellCrawling'
exec_name = 'Exe_' + sim_name

# Fixed params
fix_crl = 0.25
fix_csc = 1e7
fix_trl = 0.01
fix_di = 0.02
fix_ts = 20000

# Range for two simulation parameters
num_tsc = 20
num_ad = 20

def main():
    run_simulations()
    # make_movies_parallel()
    # combine_output()
    # plot_results()


# Create a list of commands and pass them to separate processes
def run_simulations():

    # Make a list of calls to a Chaste executable
    command_list = []
    for tsc in range(num_tsc):
        fix_tsc = 1e6 * (1 + 0.2 * (1 + tsc))

        for ad in range(num_ad):
            fix_ad = 1.0 + 0.05 * (-8 + tsc)

            command = 'nice -n 19 ' + path_to_exec + exec_name \
                      + ' --ID ' + str(tsc) + '_' + str(ad) \
                      + ' --CRL ' + str(fix_crl) \
                      + ' --CSC ' + str(fix_csc) \
                      + ' --TRL ' + str(fix_trl) \
                      + ' --DI ' + str(fix_di) \
                      + ' --TS ' + str(fix_ts) \
                      + ' --TSC ' + str(fix_tsc) \
                      + ' --AD ' + str(fix_ad)
            command_list.append(command)

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Py: Starting simulations with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool

    pool.map(execute_command, command_list)



# Make an mp4 movie from each pvd file
def make_movies_parallel():

    print("Py: Combining chaste output to movies")

    command_list = []
    for global_const in range(num_global_consts):
        for local_const in range(num_local_consts):
            for kick in range(num_kicks_per_sim):
                dir_name = path_to_output + exec_name + '/sim/' \
                           + str(global_const) + '_' \
                           + str(local_const) + '_' \
                           + str(kick) + '/'

                string_of_dir_name = '"' + dir_name + '"'

                command_list.append("""python -c 'import immersed_boundary as ib; ib.pvd_to_mp4(""" + string_of_dir_name + """)'""")

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Py: Creating movies from simulation output with " + str(count) + " processes")

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
    combined_results.write("global_const,local_const,kick,change_in_centroid,right_symmetry\n")

    for global_const in range(num_global_consts):
        for local_const in range(num_local_consts):
            for kick in range(num_kicks_per_sim):
                file_name = path_to_output + exec_name + '/sim/' \
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

    bg_gray = '0.75'
    bg_line_width = 0.5

    # Set LaTeX font rendering
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', serif='computer modern roman')
    plt.rc('figure', figsize=[10,6.18]) # inches
    plt.rc('axes', linewidth=bg_line_width, edgecolor=bg_gray, axisbelow=True)
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rc('xtick.major', size=0, pad=4)
    plt.rc('ytick.major', size=0, pad=4)


    my_data = np.genfromtxt(path_to_output + exec_name + '/combined_results.dat', delimiter=',', skip_header=1)

    col_to_plot = 4 # the column to plot from the data file

    x_vals = my_data[:,1] # local const
    y_vals = my_data[:,col_to_plot] # summary statistic

    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)

    x_range, y_range = x_max - x_min, y_max - y_min
    margin = 0.05;

    x_lims = [x_min - margin * x_range, x_max + margin * x_range]
    y_lims = [y_min - margin * y_range, y_max + margin * y_range]

    num_sims_per_plot = num_local_consts * num_kicks_per_sim

    for plot in range(num_global_consts):

        plt.clf()

        plt.xlabel(r'Local spring constant multiple', fontsize=12, labelpad=20)
        plt.ylabel(r'Right asymmetry', fontsize=12, labelpad=20)

        title = 'Global cell-cell spring constant scaled by ' + str(0.5 + 0.25 * plot)

        plt.title(title, fontsize=14, y = 1.05)

        highlight_x_vals = my_data[plot * num_sims_per_plot : (plot+1) * num_sims_per_plot,1]
        highlight_y_vals = my_data[plot * num_sims_per_plot : (plot+1) * num_sims_per_plot,col_to_plot]

        plt.scatter(x_vals, y_vals, color = '0.5', marker = 'o', facecolors='none', edgecolors='0.5')
        plt.scatter(highlight_x_vals, highlight_y_vals, color = '#F39200')

        plt.xlim(x_lims)
        plt.ylim(y_lims)


        ax = plt.gca()
        ax.grid(b=True, which='major', color=bg_gray, linestyle='dotted', dash_capstyle='round')

        plt.savefig(path_to_output + exec_name + '/Fig_' + exec_name + str(plot) + '_pdf.pdf', bbox_inches='tight', pad_inches=0.4)
        plt.savefig(path_to_output + exec_name + '/Fig_' + exec_name + str(plot) + '_eps.eps', bbox_inches='tight', pad_inches=0.5)

if __name__ == "__main__":
    main()