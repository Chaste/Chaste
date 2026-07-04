#!/bin/bash

# Chaste template script
# Use this template to create a script for running tests on supercomputers
#  * First edit this file to suit your system and job scheduler
#  * Then compile on the head node

# EDIT HERE - commands for the scheduler
# Here is an example to run in the current directory (works on GridEngine)
#$ -cwd
# GridEngine join output and error files
#$ -j oe
# PBS join output and error files
#PBS -j oe

# EDIT HERE - launch command
# Define the command and any other machine files
# Here is an example for some supercomputer from 7 years ago:
# export MPI_LAUNCH_COMMAND="mpirun -np $NSLOTS -machinefile $TMPDIR/machines"
export MPI_LAUNCH_COMMAND=""

# EDIT HERE - any other miscellaneous commands that need to be run before the tests
# For example in PBS:
# cd $PBS_O_WORKDIR

