# This script generates and submits jobs for a parameter sweep of our new force
import multiprocessing
import subprocess

# generate a list of bash commands
command_list = []

# ... and populate it with calls to a Chaste executable
for simulation_id in range(10):
    command = "../../../../chaste_build/projects/ImmersedBoundary/apps/Exe_VaryReynoldsNumber"\
              + " --ID " + str(simulation_id)\
              +  " --R " + str(1e-4 * (1 + simulation_id))
    command_list.append(command)
    
# This is a helper function that allows to run a bash command in a separate process
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

# use `count' no of processes 
count = multiprocessing.cpu_count()

# generate a pool of workers
pool = multiprocessing.Pool(processes=count)

# ... and pass the list of bash commands to the pool
pool.map(execute_command, command_list)


