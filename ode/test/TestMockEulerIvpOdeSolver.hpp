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


#ifndef _TESTMOCKEULERIVPODESOLVER_HPP_
#define _TESTMOCKEULERIVPODESOLVER_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <cxxtest/TestSuite.h>

#include "MockEulerIvpOdeSolver.hpp"
#include "Ode1.hpp"
#include "OutputFileHandler.hpp"

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

    void TestArchivingMockEulerSolver() throw(Exception)
    {
        OutputFileHandler handler("archive",false);
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
