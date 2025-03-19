/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTCOMMANDLINEARGUMENTSTUTORIAL_HPP_
#define TESTCOMMANDLINEARGUMENTSTUTORIAL_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <iomanip>
/* Most Chaste code uses PETSc to solve linear algebra problems.  This involves starting PETSc at the beginning of a test-suite
 * and closing it at the end.  (If you never run code in parallel then it is safe to replace PetscSetupAndFinalize.hpp with FakePetscSetup.hpp)
 */
#include "PetscSetupAndFinalize.hpp"
#include "CommandLineArguments.hpp"

/**
 * @file
 *
 * This is a tutorial for providing command line arguements to a test from a bash script.
 * The example bash script also includes a for loop as an example of how such input could
 * be used to run multiple tests with varying parameter inputs.
 */

 // NOTE: These tests will run if executed from terminal directly with their default varaible values. However, these
 // tutorial tests are intended to be executed through the use of bash/shell scripts (.sh files).
class TestCommandLineArgumentsTutorial : public CxxTest::TestSuite
{
public:

    /* To script the running of this test with lots of different arguments, copy and paste the following code into a bash script (`.sh`) file.
       This bash script file will need to be saved in your chaste build folder.
    ```sh
      #!/bin/bash

      # This bash script accompanies that TestCommandLineArgumentsTutorial. This script will execute the test TestCommandLineArgumentsTutorial
      # multiple times with several command line arguments as test variables.
      # Here we will declare some values we wish to later pass to a for loop.

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
          ./global/test/TestCommandLineArgumentsTutorial -opt1 $i -opt2 $j -opt3 $k &
          done
        done
      done
    ```
    */
    void TestCommandLineDefaultTutorial()
    {
        // First, we set up some variables. These variables should later be assigned our command line argument values.
        // However we will give them default values to ensure that this test can always run.
        unsigned outp1 = 1;
        unsigned outp2 = 2;
        unsigned outp3 = 3;

        // Here an if statement is utilised to check that we have passed some command line arguments to our test.
        if (CommandLineArguments::Instance()->OptionExists("-opt1") && CommandLineArguments::Instance()->OptionExists("-opt2") && CommandLineArguments::Instance()->OptionExists("-opt3"))
        {
           // If command line arguments have been passed to our test, we utilise CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption()
           // to take in our command line arguements.
           outp1 = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-opt1");
           outp2 = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-opt2");
           outp3 = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-opt3");
        }

        unsigned sum = outp1 + outp2 + outp3;

        std::cout << "When we add "<< outp1 << " ,"<<  outp2 << " and "<< outp3 <<" we get " << sum <<"\n";
    }

    /*
      In addition to using unsigned integers as command line arguments we can also pass in both doubles and a vector of arguments.
      This could be very useful for several reasons, such as having a large list of parameters one may wish to pass into a test.

      To script the running of this test with lots of different arguments, copy and paste the following code into a bash script (`.sh`) file.
      ```sh
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
           ./global/test/TestCommandLineArgumentsTutorial --my-vector-of-arguments $idouble $jdouble $kdouble &
        done
       done
      done
      ```
    */
    void TestCommandLineDoubleTutorial()
    {
        // Here an if statement is utilised to check that we have passed some command line arguments to our test.
        // To ensure the test can still run in the absence of command line arguments we will initialise the vector
        // with some default values if no command line option exists.
        if (CommandLineArguments::Instance()->OptionExists("--my-vector-of-arguments"))
        {
           // First, we need to setup our vector based on the vector passed in as our command line arguement.
           std::vector<double> vector_of_doubles = CommandLineArguments::Instance()->GetDoublesCorrespondingToOption("--my-vector-of-arguments");
           double sum = 0.0;
           // Now summing over our vector.
           for(double i = 0; i < vector_of_doubles.size(); i++){
            sum += vector_of_doubles[i];
           }
           // We will finally cout our result for the summed vector components.
           // Here we have chosen to set our precision to 4 to ensure the correct number of significant figures.
           std::cout << std::setprecision(5) << "When we add "<< vector_of_doubles[0] << ", " << vector_of_doubles[1] << " and " << vector_of_doubles[2] << " together from our vector we get " << sum  <<std::endl;

        }
    }
};

#endif /*TESTCOMMANDLINEARGUMENTSTUTORIAL_HPP_*/
