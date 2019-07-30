/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef TESTDEBUG_HPP_
#define TESTDEBUG_HPP_

#include <cxxtest/TestSuite.h>
#include "Debug.hpp"

#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestDebug : public CxxTest::TestSuite
{
private:
    void AnotherFunctionOnTheStack()
    {
        STACK;
    }

public:

    // Can't really test these other than that they compile and visually looking at output
    void TestDebugMacros()
    {
        TRACE("Some trace");

        // Note that these macros do nothing in NDEBUG -- we should ensure that the variables get used anyway
        double use_vars = 0.0;

        unsigned my_var = 3141;
        PRINT_VARIABLE(my_var);
        use_vars += (double) my_var;

        double another_var = 2.81;
        PRINT_2_VARIABLES(my_var, another_var);
        use_vars += another_var;

        double cancer_curing_constant = 0.053450242435;
        PRINT_3_VARIABLES(my_var, another_var, cancer_curing_constant);
        use_vars += cancer_curing_constant;

        double heart_disease_ending_constant = -3e-141;
        PRINT_4_VARIABLES(my_var, another_var, cancer_curing_constant, heart_disease_ending_constant);
        use_vars += heart_disease_ending_constant;

        std::cout << "\n\n";

        for (unsigned i=0; i<10; i++)
        {
            HOW_MANY_TIMES_HERE("inside for loop");

            for (unsigned j=0; j<2; j++)
            {
                HOW_MANY_TIMES_HERE("nested loop");
            }
        }

        for (unsigned j=0; j<10 /*change to 11 and it should quit*/; j++)
        {
            QUIT_AFTER_N_VISITS(11);
        }

        std::cout << "\n\n\n";

        for (unsigned j=0; j<3; j++)
        {
            TRACE_FROM_NTH_VISIT("hello",2);
        }

        std::vector<double> vec(4);
        vec[0] = 0.0;
        vec[1] = 1.0;
        vec[2] = 2.7;
        vec[3] = 3.1;
        PRINT_VECTOR(vec);

        MARK; // Something like: "DEBUG: ./global/test/TestDebug.hpp at line 110"
        MARK_IN_ORDER; // Something like: "DEBUG: proc 0 ./global/test/TestDebug.hpp at line 111" etc
        MARK_ON_PROCESS(1); // In serial there is no output.  In parallel, something like: "DEBUG: proc 1: ./global/test/TestDebug.hpp at line 112"

        // Check interaction with process isolation
        PetscTools::IsolateProcesses();
        TRACE("Processes are isolated.");
        PetscTools::IsolateProcesses(false);

        AnotherFunctionOnTheStack();

        MARK_MEMORY;   // Should print that PRINT_MEMORY will now be relative to here
        PRINT_MEMORY;  // Should be 0 Mb
        UNMARK_MEMORY; // Should print that PRINT_MEMORY will now be the absolute footprint
        PRINT_MEMORY;  // Should be > 0 Mb

        MARK_MEMORY;
        std::vector<double> poinless_vector;
        poinless_vector.resize(100000);
        PRINT_MEMORY;                        // Should be x Mb
        poinless_vector.resize(200000);
        PRINT_MEMORY;                        // Should be at least 2x Mb
    }
};

#endif /*TESTDEBUG_HPP_*/
