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
