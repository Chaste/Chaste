#!/bin/bash
#
# This script compiles the output from several different runs, defined by filling
# in the SETUP section below. It prints out one line of visulizer output out of
# every 25 and sticks files together to make movies more simply.
#

# BEGIN SETUP
results=/tmp/pmxaw/testoutput/Noddy_WNT_No_Area_No_Length
first_experiment=500
last_experiment=600
# END OF SETUP

echo "Cleaning up previous results files"

echo "Compiling results files"
for (( i=$first_experiment ; i<=$last_experiment ; i=i+5 ))
do
	echo "results $i"
	more $results/results_from_time_"$i"/results.vizvoronoi >> $results/All_results
done
