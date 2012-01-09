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


#ifndef TESTSOLVINGSTIFFODESYSTEMS_HPP_
#define TESTSOLVINGSTIFFODESYSTEMS_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include <iostream>

#include "BackwardEulerIvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "FieldNoyesReactionSystem.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestSolvingStiffOdeSystems : public CxxTest::TestSuite
{
public:

    /**
     *  Solve the Field-Noyes system.
     *  This is a stiff ODE so Runge-Kutta won't work - program hangs if
     *  end time > 0.01 and variables go out of bounds.
     *  Can be solved ok with backward Euler.
     */
    void TestFieldNoyesReactionSystem()
    {
        FieldNoyesReactionSystem ode_system;
        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution ode_solution;
        std::vector<double> ode_state_variables = ode_system.GetInitialConditions();

        // note small end time
        ode_solution = rk4_solver.Solve(&ode_system, ode_state_variables, 0.0, 0.01, 0.001, 0.001);

        int last = ode_solution.GetNumberOfTimeSteps();
        for (int i=0; i<=last; i++)
        {
            double x = ode_solution.rGetSolutions()[i][0];
            TS_ASSERT_LESS_THAN(x, 9.7e3);
            TS_ASSERT_LESS_THAN(0.99, x);
        }

        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system.GetNumberOfStateVariables());
        OdeSolution ode_solution2;
        std::vector<double> ode_state_variables2 = ode_system.GetInitialConditions();
        ode_solution2 = backward_euler_solver.Solve(&ode_system, ode_state_variables2, 0.0, 0.1, 0.001, 0.001);

        last = ode_solution2.GetNumberOfTimeSteps();
        for (int i=0; i<=last; i++)
        {
            double x = ode_solution2.rGetSolutions()[i][0];
            TS_ASSERT_LESS_THAN(x, 9.7e3);
            TS_ASSERT_LESS_THAN(0.99, x);
        }
    }
};

#endif /*TESTSOLVINGSTIFFODESYSTEMS_HPP_*/
