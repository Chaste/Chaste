# This script generates and submits jobs for a parameter sweep of our new force
import multiprocessing
import subprocess

def main():
    run_simulations()
    #combine_output()
    plot_results()

def run_simulations():
    # generate a list of bash commands
    command_list = []

    # ... and populate it with calls to a Chaste executable
    for simulation_id in range(3):
        for spring_const in range(3):
            command = "../../../../cmake_build/projects/ImmersedBoundary/apps/Exe_VaryCellCellAdhesion"\
                  + " --ID " + str(simulation_id)\
                  +  " --K " + str(spring_const)
            command_list.append(command)

    # use `count' no of processes
    count = multiprocessing.cpu_count()

    # generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # ... and pass the list of bash commands to the pool
    print("Starting...")
    pool.map(execute_command, command_list)
    print("Done.")

# This is a helper function that allows to run a bash command in a separate process
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

def combine_output():

    fd = open('./testoutput/Exe_VaryCellCellAdhesion/combined_results.csv','a')

    for simulation_id in range(3):
        for spring_const in range(3):
            file_name = "./testoutput/Exe_VaryCellCellAdhesion/" + str(spring_const) + "_" + str(simulation_id) + "/results.dat"

            csvfile = open(file_name, 'r')
            fd.write(csvfile.readline())
            csvfile.close()
            fd.write("\n")

    fd.close()

def plot_results():
    print("not yet plotting")

if __name__ == "__main__":
    main()