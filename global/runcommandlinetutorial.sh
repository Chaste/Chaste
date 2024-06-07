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
# As bash script cannot directly handle double arithmetic we will First
# declare these as larger variables to be handled by another programme later.
N=2000
L=3000
M=4000

# Here we set up a simple for loop over variables i,j and k based on the values of N,L and M.
for ((i = 1; i <= N; i += 1000)); do
  for ((j = 1001; j <= L; j += 1000)); do
    for ((k = 2001; k <= M; k += 1000)); do
    # As bash cannot directly handle double floating point arithmetic we will utilise
    # awk to divide our varaibles and convert them to doubles.
    idouble=$(awk "BEGIN {printf \"%.2f\",$i/120}")
    jdouble=$(awk "BEGIN {printf \"%.2f\",$j/150}")
    kdouble=$(awk "BEGIN {printf \"%.2f\",$k/180}")

    # Each loop runs an instance of the TestCommandLineArgumentsTutorial with a vector
    # containing the double corrected version of our varaibles.
     ~/build/global/test/TestCommandLineArgumentsTutorial --my-vector-of-arguments $idouble $jdouble $kdouble &
    done
  done
done
