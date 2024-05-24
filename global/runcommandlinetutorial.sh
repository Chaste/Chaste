#!/bin/bash

# This bash script accompanies that TestCommandLineArgumentsTutorial
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