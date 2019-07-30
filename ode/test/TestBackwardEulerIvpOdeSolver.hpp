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


#ifndef TESTBACKWARDEULERIVPODESOLVER_HPP_
#define TESTBACKWARDEULERIVPODESOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "CheckpointArchiveTypes.hpp"

#include "OdeThirdOrder.hpp"
#include "OdeThirdOrderWithEvents.hpp"
#include "Ode4.hpp"
#include "Ode5.hpp"
#include "Ode5Jacobian.hpp"
#include "VanDerPolOde.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "OutputFileHandler.hpp"
#include "ArchiveLocationInfo.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestBackwardEulerIvpOdeSolver: public CxxTest::TestSuite
{
public:
    void TestBackwardEulerSystemOf3Equations()
    {
        OdeThirdOrder ode_system;

        double h_value = 0.01;

        // Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system.GetNumberOfStateVariables());

        // Cover the SetEpsilonForNumericalJacobian() method
        backward_euler_solver.SetEpsilonForNumericalJacobian(1e-6);

        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = backward_euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        unsigned last = solutions.GetNumberOfTimeSteps();

        double numerical_solution[3];
        numerical_solution[0] = solutions.rGetSolutions()[last][0];
        numerical_solution[1] = solutions.rGetSolutions()[last][1];
        numerical_solution[2] = solutions.rGetSolutions()[last][2];

        // The tests
        double analytical_solution[3];

        analytical_solution[0] = -sin(2.0);
        analytical_solution[1] = sin(2.0)+cos(2.0);
        analytical_solution[2] = 2*sin(2.0);

        double global_error_euler = 0.5*2*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(numerical_solution[0], analytical_solution[0], global_error_euler);
        TS_ASSERT_DELTA(numerical_solution[1], analytical_solution[1], global_error_euler);
        TS_ASSERT_DELTA(numerical_solution[2], analytical_solution[2], global_error_euler);
    }

    void TestBackwardEulerNonlinearEquation()
    {
        Ode4 ode_system;

        double h_value = 0.01;

        // Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system.GetNumberOfStateVariables());
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();

        solutions = backward_euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        unsigned last = solutions.GetNumberOfTimeSteps();

        double numerical_solution;
        numerical_solution = solutions.rGetSolutions()[last][0];

        // The tests
        double analytical_solution = 1.0/(1.0+exp(-12.5));

        TS_ASSERT_DELTA(numerical_solution, analytical_solution, 1.0e-4);
    }

    void TestBackwardEulerAnotherNonlinearEquation()
    {
        Ode5 ode_system;

        double h_value = 0.01;
        double end_time = 1.0;

        // Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system.GetNumberOfStateVariables());
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();

        solutions = backward_euler_solver.Solve(&ode_system, state_variables, 0.0, end_time, h_value, h_value);
        unsigned last = solutions.GetNumberOfTimeSteps();

        double numerical_solution;
        numerical_solution = solutions.rGetSolutions()[last][0];

        // The tests
        double analytical_solution = 1.0/(1.0+4.0*exp(-100.0*end_time));

        TS_ASSERT_DELTA(numerical_solution, analytical_solution, 1.0e-3);
    }

    void TestBackwardEulerSystemOf3EquationsWithEvents()
    {
        OdeThirdOrderWithEvents ode_system_with_events;

        double h_value = 0.01;

        // Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system_with_events.GetNumberOfStateVariables());
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system_with_events.GetInitialConditions();
        solutions = backward_euler_solver.Solve(&ode_system_with_events, state_variables, 0.0, 2.0, h_value, h_value);
        unsigned last = solutions.GetNumberOfTimeSteps();

        // Final time should be pi/6 (?)
        TS_ASSERT_DELTA( solutions.rGetTimes()[last], 0.5236, 0.01);

        // Penultimate y0 should be greater than -0.5
        TS_ASSERT_LESS_THAN(-0.5,solutions.rGetSolutions()[last-1][0]);

        // Final y0 should be less than -0.5
        TS_ASSERT_LESS_THAN( solutions.rGetSolutions()[last][0], -0.5);

        // Solver should correctly state the stopping event occurred
        TS_ASSERT_EQUALS(backward_euler_solver.StoppingEventOccurred(), true);
    }

    void TestBackwardEulerAnotherNonlinearEquationAnalytic()
    {
        Ode5Jacobian ode_system;

        double h_value = 0.01;
        double end_time = 1.0;

        // Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system.GetNumberOfStateVariables());
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();

        solutions = backward_euler_solver.Solve(&ode_system, state_variables, 0.0, end_time, h_value, h_value);
        unsigned last = solutions.GetNumberOfTimeSteps();

        double numerical_solution;
        numerical_solution = solutions.rGetSolutions()[last][0];

        // The tests
        double analytical_solution = 1.0/(1.0+4.0*exp(-100.0*end_time));

        TS_ASSERT_DELTA(numerical_solution, analytical_solution, 1.0e-3);
    }

    void TestBackwardEulerVanDerPolOde()
    {
        VanDerPolOde ode_system;

        double h_value = 0.01;
        double end_time = 100.0;

        // Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system.GetNumberOfStateVariables());
        backward_euler_solver.ForceUseOfNumericalJacobian(); // coverage
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();

        solutions = backward_euler_solver.Solve(&ode_system, state_variables, 0.0, end_time, h_value, 5*h_value);
        unsigned last = solutions.GetNumberOfTimeSteps();

//        OutputFileHandler handler("");
//        out_stream rabbit_file=handler.OpenOutputFile("foxrabbit.dat");
//
//        for (unsigned i=0; i<last; i++)
//        {
//            (*rabbit_file) << solutions.rGetSolutions()[i][0] << "\t" << solutions.rGetSolutions()[i][1] << "\n" << std::flush;
//        }
//        rabbit_file->close();

        // assert that we are within a [-2,2] in x and [-2,2] in y (on limit cycle)
        TS_ASSERT_DELTA(solutions.rGetSolutions()[last][0], 0, 2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[last][1], 0, 2);
    }

    void TestArchivingSolver()
    {
        OutputFileHandler handler("archive", false);
        ArchiveLocationInfo::SetArchiveDirectory(handler.FindFile(""));
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("backward_euler_solver.arch");

        VanDerPolOde ode_system;

        double h_value = 0.01;
        double end_time = 100.0;

        // Create and archive simulation time
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set up a solver
            AbstractIvpOdeSolver* const p_backward_euler_solver = new BackwardEulerIvpOdeSolver(ode_system.GetNumberOfStateVariables());

            // Should always archive a pointer
            output_arch << p_backward_euler_solver;

            // Change stimulus a bit
            delete p_backward_euler_solver;
        }

        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Create a pointer
            AbstractIvpOdeSolver* p_backward_euler;
            input_arch >> p_backward_euler;
            OdeSolution solutions;

            std::vector<double> state_variables = ode_system.GetInitialConditions();

            solutions = p_backward_euler->Solve(&ode_system, state_variables, 0.0, end_time, h_value, 5*h_value);
            unsigned last = solutions.GetNumberOfTimeSteps();

            // assert that we are within a [-2,2] in x and [-2,2] in y (on limit cycle)
            TS_ASSERT_DELTA(solutions.rGetSolutions()[last][0], 0, 2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[last][1], 0, 2);

            delete p_backward_euler;
        }
    }
};

#endif /*TESTBACKWARDEULERIVPODESOLVER_HPP_*/
