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


#ifndef TESTBACKWARDEULERIVPODESOLVER_HPP_
#define TESTBACKWARDEULERIVPODESOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "OdeThirdOrder.hpp"
#include "OdeThirdOrderWithEvents.hpp"
#include "Ode4.hpp"
#include "Ode5.hpp"
#include "Ode5Jacobian.hpp"
#include "VanDerPolOde.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "OutputFileHandler.hpp"

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

        analytical_solution[0] = -sin(2);
        analytical_solution[1] = sin(2)+cos(2);
        analytical_solution[2] = 2*sin(2);

        double global_error_euler = 0.5*2*(exp(2)-1)*h_value;
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

    void TestBackwardEulerAnotherNonlinearEquationAnalytic() throw(Exception)
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

    void TestArchivingSolver() throw(Exception)
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "backward_euler_solver.arch";

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
