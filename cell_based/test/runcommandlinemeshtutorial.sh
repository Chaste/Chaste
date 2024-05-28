#!/bin/bash
# This bash script is utilised alongside the TestCommandLineVertexMeshTutorial.hpp file.
# First, you should build Chaste outside of the main chaste folder in a directory named "build".
# If your build folder has a different name such as "chaste_build" e.t.c, then replace the "build"
# on line 21 with your folders name. Second, you should also investigate the accompanying 
# TestCommandLineVertexMeshTutorial.hpp to see how the simulations are executed

# This bash script will be utilised to executed the TestCommandLineVertexMeshTutorial multiple times.
# Creating 9 total simulations with a tissue height and width from 2->10.

# To execute this bash script, after first building Chaste, open a new terminal and navigate to Chaste/cell_based/test/ .
# Now enter into the terminal "bash runcommandlinemeshtutorial.sh" and press enter. You should see multiple tests being executed
# in succession.

# N is the max number of cells we wish to have output. 
N=10

# This for loop will run the TestCommandLineVertexMeshTutorial with an initial size of 2 to a maximum size of N.
# Due to how the test is set up a different output will be created for each mesh size.
for ((i = 2; i <= N; i += 1)); do
    ~/build/cell_based/test/TestCommandLineVertexMeshTutorial -opt $i &
done
echo All done
