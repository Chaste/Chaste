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


#ifndef _TESTMOCKEULERIVPODESOLVER_HPP_
#define _TESTMOCKEULERIVPODESOLVER_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "MockEulerIvpOdeSolver.hpp"
#include "Ode1.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestMockEulerIvpOdeSolver: public CxxTest::TestSuite
{
public:

    void TestMockEulerSolver()
    {
        Ode1 ode_system;

        // Initialising the instance of our solver class
        MockEulerIvpOdeSolver euler_solver;

        // Initialising the instance of our solution class
        OdeSolution solutions;

        // Solving the ODE problem and writing to solution
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.001, 2.0);

        int last = solutions.GetNumberOfTimeSteps();

        // Test to see if this worked
        double testvalue = solutions.rGetSolutions()[last][0];

        TS_ASSERT_DELTA(testvalue, 2.0, 0.01);

        TS_ASSERT_EQUALS(euler_solver.GetCallCount(), 1U);

        ode_system.SetDefaultInitialCondition(0, 0.0);

        state_variables = ode_system.GetInitialConditions();
        solutions = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.001, 2.0);

        last = solutions.GetNumberOfTimeSteps();

        // Test to see if this worked
        testvalue = solutions.rGetSolutions()[last][0];

        TS_ASSERT_DELTA(testvalue, 2.0, 0.01);

        TS_ASSERT_EQUALS(euler_solver.GetCallCount(), 2u);
    }

    void TestArchivingMockEulerSolver()
    {
        EXIT_IF_PARALLEL;

        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "mock_euler_solver.arch";

        // Create and archive simulation time
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set up a solver
            AbstractIvpOdeSolver* const p_mock_euler_ivp_ode_solver = new MockEulerIvpOdeSolver;


            Ode1 ode_system;
            p_mock_euler_ivp_ode_solver->SolveAndUpdateStateVariable(&ode_system, 0, 1, 0.01);
            TS_ASSERT_DELTA(ode_system.rGetStateVariables()[0], 1.0, 1e-2);

            // Should always archive a pointer
            output_arch << p_mock_euler_ivp_ode_solver;

            // Change stimulus a bit
            delete p_mock_euler_ivp_ode_solver;
        }

        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Create a pointer
            AbstractIvpOdeSolver* p_mock_euler;
            input_arch >> p_mock_euler;

            Ode1 ode_system;
            p_mock_euler->SolveAndUpdateStateVariable(&ode_system, 0, 1, 0.01);
            TS_ASSERT_DELTA(ode_system.rGetStateVariables()[0], 1.0, 1e-2);

            TS_ASSERT_EQUALS(static_cast<MockEulerIvpOdeSolver&> (*p_mock_euler).GetCallCount(), 2u);

            delete p_mock_euler;
        }
    }
};

#endif //_TESTMOCKEULERIVPODESOLVER_HPP_
