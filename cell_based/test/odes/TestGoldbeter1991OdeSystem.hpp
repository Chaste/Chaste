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

#ifndef TESTGOLDBETER1991ODESYSTEM_HPP_
#define TESTGOLDBETER1991ODESYSTEM_HPP_


#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <ctime>
#include <vector>
#include <iostream>

#include "OutputFileHandler.hpp"
#include "Goldbeter1991OdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"

//#include "RungeKutta4IvpOdeSolver.hpp"
//#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
//#include "BackwardEulerIvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/*
 * Basic tests only - no archiving test
 *
 */
class TestGoldbeter1991OdeSystem : public CxxTest::TestSuite
{
public:

    void TestGoldbeter1991Equation()
    {
        Goldbeter1991OdeSystem ode_system;

        double time = 0.0;
        std::vector<double> initial_conditions;
        initial_conditions.push_back(0.01);
        initial_conditions.push_back(0.01);
        initial_conditions.push_back(0.01);

        std::vector<double> derivs(initial_conditions.size());
        ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct
        TS_ASSERT_DELTA(derivs[0], 0.0240, 1e-4);
        TS_ASSERT_DELTA(derivs[1], -0.9414, 1e-4);
        TS_ASSERT_DELTA(derivs[2], -0.3233, 1e-4);
    }


    void TestGoldbeter1991OSolver()
    {
        Goldbeter1991OdeSystem ode_system;

        RungeKutta4IvpOdeSolver ode_solver;

        OdeSolution solutions;

        std::vector<double> initial_conditions = ode_system.GetInitialConditions();
        double start_time = 0.0;
        double end_time = 100.0;
        double h_value = 0.01; // 1.0 // maximum tolerance

        //Test the hard coded ics
        TS_ASSERT_DELTA(initial_conditions[0], 0.01, 1e-6);
        TS_ASSERT_DELTA(initial_conditions[1], 0.01, 1e-6);
        TS_ASSERT_DELTA(initial_conditions[2], 0.01, 1e-6);

        double cpu_start_time = (double) std::clock();
        solutions = ode_solver.Solve(&ode_system, initial_conditions, start_time, end_time, h_value, h_value);
        double cpu_end_time = (double) std::clock();
        double cpu_elapsed_time = (cpu_end_time - cpu_start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "1. Solver Elapsed time = " << cpu_elapsed_time << "\n";

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 7.8 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], end_time, 1e-2);
        //std::cout <<  "End time = " << solutions.rGetTimes()[end] << "\n";
        // Decent results - checked with numpy # [ 0.54706214  0.29369527  0.00678837]
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 0.5470, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 0.2936, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 0.0067, 1e-4);
    }

    void TestArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "gb1991_ode.arch";

        {

            std::vector<double> state_variables;
            state_variables.push_back(3.0);
            state_variables.push_back(4.0);
            state_variables.push_back(5.0);

            Goldbeter1991OdeSystem ode_system(state_variables);

            ode_system.SetDefaultInitialCondition(2, 3.25);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 3u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.01, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[1], 0.01, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[2], 3.2500, 1e-4);

            double var1 = ode_system.GetStateVariable(0);
            double var2 = ode_system.GetStateVariable(1);
            double var3 = ode_system.GetStateVariable(2);

            TS_ASSERT_DELTA(var1, 3.0, 1e-3);
            TS_ASSERT_DELTA(var2, 4.0, 1e-3);
            TS_ASSERT_DELTA(var3, 5.0, 1e-3);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem* const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        {
            AbstractOdeSystem* p_ode_system;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_ode_system;

            // Check that archiving worked correctly
            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 3u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.01, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[1], 0.01, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[2], 0.01, 1e-4);

            double var1 = p_ode_system->GetStateVariable(0);
            double var2 = p_ode_system->GetStateVariable(1);
            double var3 = p_ode_system->GetStateVariable(2);

            TS_ASSERT_DELTA(var1, 3.0, 1e-3);
            TS_ASSERT_DELTA(var2, 4.0, 1e-3);
            TS_ASSERT_DELTA(var3, 5.0, 1e-3);

            // Tidy up
            delete p_ode_system;
        }
    }
};

#endif /* TESTGOLDBETER1991ODESYSTEM_HPP_ */
