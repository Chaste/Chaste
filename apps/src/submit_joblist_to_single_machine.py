# This script generates and submits jobs for a parameter sweep of our new force
import multiprocessing
import subprocess

# generate a list of bash commands
command_list = []

# ... and populate it with calls to a Chaste executable
for simulation_id in range(100):
    command = "./ExampleCellBasedExecutable --ID " + str(simulation_id)
    command_list.append(command)
    
# This is a helper function that allows to run a bash command in a separate process
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

# use `count' no of processes 
count = 3 # use multiprocessing.cpu_count() for the number of cores on your machine

# generate a pool of workers
pool = multiprocessing.Pool(processes=count)

# ... and pass the list of bash commands to the pool
pool.map(execute_command, command_list)


