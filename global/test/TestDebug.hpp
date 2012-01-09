/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef TESTDEBUG_HPP_
#define TESTDEBUG_HPP_

#include <cxxtest/TestSuite.h>
#include "Debug.hpp"

// Interestingly, this test won't fork if you attempt to run in parallel unless this is included
#include "PetscSetupAndFinalize.hpp"

class TestDebug : public CxxTest::TestSuite
{
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

        MARK; // Something like: "DEBUG: ./global/test/TestDebug.hpp at line 95"
    }
};

#endif /*TESTDEBUG_HPP_*/
