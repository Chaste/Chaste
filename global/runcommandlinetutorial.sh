#!/bin/bash

# This bash script accompanies that TestCommandLineArgumentsTutorial.
# First, you should build Chaste outside of the main chaste folder in a directory named "build".
# If youre build folder has a different name such as "chaste_build" e.t.c, then replace the "build"
# on line 19 with your folders name. Second, you should also investigate the accompanying 
# TestCommandLineArgumentsTutorial.hpp to see how the simulations are executed

# This bash script will be utilised to executed the TestCommandLineArgumentsTutorial multiple times.
# Outputting the sum of the command lines arguments at each iteration to the terminal.

# To execute this bash script, after first building Chaste, open a new terminal and navigate to Chaste/global/ .
# Now enter into the terminal "bash runcommandlinetutorial.sh" and press enter. You should see multiple tests being executed
# in succession.

# Here we will declare some values we wish to later pass to a for loop.  
N=2
L=3
M=4

# Here we set up a simple for loop over variables i,j and k based on the values of N,L and M.
for ((i = 0; i <= N; i += 1)); do
  for ((j = 1; j <= L; j += 1)); do
    for ((k = 2; k <= M; k += 1)); do
    # Each loop runs an instance of the TestCommandLineArgumentsTutorial with opt1,opt2 and opt3 taking on the 
    # values of i,j and k resepctivley.      
     ~/build/global/test/TestCommandLineArgumentsTutorial -opt1 $i -opt2 $j -opt3 $k &
    done
  done
done     
