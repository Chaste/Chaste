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


#ifndef _TESTABSTRACTIVPODESOLVER_HPP_
#define _TESTABSTRACTIVPODESOLVER_HPP_

#include <cxxtest/TestSuite.h>

#include <iostream>
#include <sstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "OdeSolution.hpp"

#include "Ode1.hpp"
#include "Ode2.hpp"
#include "Ode4.hpp"
#include "OdeFirstOrder.hpp"
#include "OdeSecondOrder.hpp"
#include "OdeSecondOrderWithEvents.hpp"
#include "OdeThirdOrder.hpp"
#include "ParameterisedOde.hpp"

#include "FileComparison.hpp"
#include "NumericFileComparison.hpp"

#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "OutputFileHandler.hpp"
#include "ArchiveLocationInfo.hpp"


class TestAbstractIvpOdeSolver: public CxxTest::TestSuite
{
private:

    void MyTestGenericSolver(AbstractIvpOdeSolver& rSolver, double startTime, double endTime, double dt, double samplingTime)
    {
        // Initialise the instances of our ODE system and solution classes
        Ode1 ode_system;

        // Solving the ODE problem. Note that dt and the sampling time are different
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        const OdeSolution solutions = rSolver.Solve(&ode_system, state_variables, startTime, endTime, dt, samplingTime);

        int num_timesteps = solutions.GetNumberOfTimeSteps();

        // The number of timesteps should be (just about) equal to
        // end_time/sampling_time = 2/0.01 = 200
        TS_ASSERT_DELTA(num_timesteps, (endTime-startTime)/samplingTime, 1);

        // Also check the size of the data is correct
        TS_ASSERT_EQUALS(solutions.rGetSolutions().size(), (unsigned) (num_timesteps+1));
        TS_ASSERT_EQUALS(solutions.rGetTimes().size(), (unsigned) (num_timesteps+1));

        int last = num_timesteps;

        // Test the solution is correct
        double testvalue = solutions.rGetSolutions()[last][0];

        // Exact solution of Ode1 is y=t-t0
        TS_ASSERT_DELTA(testvalue, endTime-startTime, 0.01);

        // For coverage of OdeSolution
        std::vector<double> var0 = solutions.GetVariableAtIndex(0);
        TS_ASSERT_DELTA(var0[last], testvalue, 1e-12);

        // Test second version of Solve
        ode_system.SetStateVariables(ode_system.GetInitialConditions());
        state_variables = ode_system.rGetStateVariables();
        rSolver.Solve(&ode_system, state_variables, startTime, endTime, dt);
        TS_ASSERT_DELTA(state_variables[0], endTime-startTime, 0.01);

        // No stopping event was specified in the ODE, so check the
        // solver correctly states it didn't stop due to a
        // stopping event.
        TS_ASSERT_EQUALS(rSolver.StoppingEventOccurred(), false);
    }

    // Test a given solver on an ODE which has a stopping event defined
    void MyTestSolverOnOdesWithEvents(AbstractIvpOdeSolver& rSolver)
    {
        // ODE which has solution y0 = cos(t), and stopping event y0<0,
        // ie should stop when t = pi/2;
        OdeSecondOrderWithEvents ode_with_events;

        OdeSolution solutions;
        std::vector<double> state_variables = ode_with_events.GetInitialConditions();
        solutions = rSolver.Solve(&ode_with_events, state_variables, 0.0, 2.0, 0.001, 0.001);

        int num_timesteps = solutions.GetNumberOfTimeSteps();

        // Final time should be around pi/2
        TS_ASSERT_DELTA( solutions.rGetTimes()[num_timesteps], M_PI_2, 0.01);

        // Penultimate y0 should be greater than zero
        TS_ASSERT_LESS_THAN( 0, solutions.rGetSolutions()[num_timesteps-1][0]);

        // Final y0 should be less than zero
        TS_ASSERT_LESS_THAN( solutions.rGetSolutions()[num_timesteps][0], 0);

        // Solver should correctly state the stopping event occurred
        TS_ASSERT_EQUALS(rSolver.StoppingEventOccurred(), true);

        // This is to cover the exception when a stopping event occurs before the first timestep.
        TS_ASSERT_THROWS_THIS(rSolver.Solve(&ode_with_events, state_variables, 2.0, 3.0, 0.001),
                "(Solve without sampling) Stopping event is true for initial condition");

        ///////////////////////////////////////////////
        // Repeat with sampling time larger than dt
        ///////////////////////////////////////////////

        state_variables = ode_with_events.GetInitialConditions();
        solutions = rSolver.Solve(&ode_with_events, state_variables, 0.0, 2.0, 0.001, 0.01);

        num_timesteps = solutions.GetNumberOfTimeSteps();

        // Final time should be around pi/2
        TS_ASSERT_DELTA( solutions.rGetTimes()[num_timesteps], M_PI_2, 0.01);

        // Penultimate y0 should be greater than zero
        TS_ASSERT_LESS_THAN( 0, solutions.rGetSolutions()[num_timesteps-1][0]);

        // Final y0 should be less than zero
        TS_ASSERT_LESS_THAN( solutions.rGetSolutions()[num_timesteps][0], 0);

        // Solver should correctly state the stopping event occurred
        TS_ASSERT_EQUALS(rSolver.StoppingEventOccurred(), true);

        // Cover the check event isn't initially true exception
        std::vector<double> bad_init_cond;
        bad_init_cond.push_back(-1);  //y0 < 0 so stopping event true
        bad_init_cond.push_back(0.0);
        TS_ASSERT_THROWS_THIS(rSolver.Solve(&ode_with_events, bad_init_cond, 0.0, 2.0, 0.001, 0.01),
                "(Solve with sampling) Stopping event is true for initial condition");
    }

public:

    void TestCoverageOfWriteToFile()
    {
        Ode2 ode_system;
        OdeSolution solutions;
        EulerIvpOdeSolver solver;

        // Solve
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = solver.Solve(&ode_system, state_variables, 0.0, 0.1, 0.1, 0.1);

        // Write
        solutions.WriteToFile("OdeSolution", "Ode2_8", "time");

        // Write at lower precision with derived quantities (this ODE doesn't actually have any derived quantities)
        solutions.CalculateDerivedQuantitiesAndParameters(&ode_system);

        solutions.WriteToFile("OdeSolution", "Ode2_4", "time", 1, false, 4, true);
        PetscTools::Barrier("TestCoverageOfWriteToFile");
        NumericFileComparison comparer(OutputFileHandler::GetChasteTestOutputDirectory() + "OdeSolution/Ode2_4.dat",
                                       "ode/test/data/Ode2_4.dat");

        TS_ASSERT(comparer.CompareFiles(1e-6));
        //The info file should now contain the ODE solver (EulerIvpOdeSolver) identifier
        FileFinder generated_file("OdeSolution/Ode2_4.info", RelativeTo::ChasteTestOutput);
        FileFinder reference_file("ode/test/data/Ode2_4.info", RelativeTo::ChasteSourceRoot);
        FileComparison file_comparer(generated_file, reference_file);
        TS_ASSERT(file_comparer.CompareFiles());
    }

    void TestEulerSolver()
    {
        EulerIvpOdeSolver euler_solver;

        MyTestGenericSolver(euler_solver,  0.0, 2.0, 0.001, 0.001);
        MyTestGenericSolver(euler_solver,  1.0, 2.0, 0.001, 0.01);
        MyTestGenericSolver(euler_solver, -1.0, 2.0, 0.001, 2);
        MyTestGenericSolver(euler_solver,  0.0, 0.4, 0.01,  0.34);

        MyTestSolverOnOdesWithEvents(euler_solver);

        // Test SolveAndUpdateStateVariable()
        Ode1 ode_system;
        euler_solver.SolveAndUpdateStateVariable(&ode_system, 0, 1, 0.01);
        TS_ASSERT_DELTA(ode_system.rGetStateVariables()[0], 1.0, 1e-2);

        // Cover an exception. This throws because SolveAndUpdateStateVar
        // called but the state is not set up in this ODE system.
        OdeSecondOrder ode2;
        TS_ASSERT_THROWS_THIS(euler_solver.SolveAndUpdateStateVariable(&ode2, 0, 1, 0.01),
                "SolveAndUpdateStateVariable() called but the state variable vector in the ODE system is not set up");
    }

    void TestRungeKutta2Solver()
    {
        RungeKutta2IvpOdeSolver rk2_solver;

        MyTestGenericSolver(rk2_solver,  0.0, 2.0, 0.001, 0.001);
        MyTestGenericSolver(rk2_solver,  1.0, 2.0, 0.001, 0.01);
        MyTestGenericSolver(rk2_solver, -1.0, 2.0, 0.001, 2);
        MyTestGenericSolver(rk2_solver,  0.0, 0.4, 0.01,  0.34);

        MyTestSolverOnOdesWithEvents(rk2_solver);

        // Test SolveAndUpdateStateVariable()
        Ode1 ode_system;
        rk2_solver.SolveAndUpdateStateVariable(&ode_system, 0, 1, 0.01);
        TS_ASSERT_DELTA(ode_system.rGetStateVariables()[0], 1.0, 1e-2);

    }

    void TestRungeKutta4Solver()
    {
        RungeKutta4IvpOdeSolver rk4_solver;

        MyTestGenericSolver(rk4_solver,  0.0, 2.0, 0.001, 0.001);
        MyTestGenericSolver(rk4_solver,  1.0, 2.0, 0.001, 0.01);
        MyTestGenericSolver(rk4_solver, -1.0, 2.0, 0.001, 2);
        MyTestGenericSolver(rk4_solver,  0.0, 0.4, 0.01,  0.34);

        MyTestSolverOnOdesWithEvents(rk4_solver);

        // Test SolveAndUpdateStateVariable()
        Ode1 ode_system;
        rk4_solver.SolveAndUpdateStateVariable(&ode_system, 0, 1, 0.01);
        TS_ASSERT_DELTA(ode_system.rGetStateVariables()[0], 1.0, 1e-2);
    }

    void TestWithParameters()
    {
        ParameterisedOde ode; // dy/dt = a, y(0) = 0.
        EulerIvpOdeSolver euler_solver;

        TS_ASSERT_EQUALS(ode.GetParameter(0), 0);

        // Test with a = 0 => y = 0.
        euler_solver.SolveAndUpdateStateVariable(&ode, 0, 1, 0.01);
        TS_ASSERT_DELTA(ode.rGetStateVariables()[0], 0.0, 1e-6);

        // Test with a = 5 => y = 5t.
        ode.SetStateVariables(ode.GetInitialConditions());
        ode.SetParameter(0, 5.0);
        euler_solver.SolveAndUpdateStateVariable(&ode, 0, 1, 0.01);
        TS_ASSERT_DELTA(ode.rGetStateVariables()[0], 5.0, 1e-2);

        // Test with a = 5 => y = 5t, for calculating derived quantities
        std::vector<double> inits = ode.GetInitialConditions();
        OdeSolution solution = euler_solver.Solve(&ode, inits, 0, 1, 0.01, 0.1);

        // Exception coverage
        TS_ASSERT_THROWS_THIS(solution.WriteToFile("OdeSolution", "ParameterisedOde", "seconds", 1, false, 4, true),
                              "You must first call ""CalculateDerivedQuantitiesAndParameters()"" in order to write derived quantities.");

        // Check solution and derived quantity for all times
        for (unsigned i=0; i<solution.rGetSolutions().size(); i++)
        {
            TS_ASSERT_DELTA(solution.rGetSolutions()[i][0], 5.0*solution.rGetTimes()[i], 1e-2);
             // (derived quantity = 2a+y)
            TS_ASSERT_DELTA(solution.rGetDerivedQuantities(&ode)[i][0], 2.0*5.0 + solution.rGetSolutions()[i][0], 1e-2);
        }

        solution.CalculateDerivedQuantitiesAndParameters(&ode);

        // Check that the new methods work correctly...
        std::vector<double> ys = solution.GetAnyVariable("y");
        std::vector<double> as = solution.GetAnyVariable("a");
        std::vector<double> two_a_plus_y = solution.GetAnyVariable("2a_plus_y");

        TS_ASSERT_THROWS_THIS(solution.GetAnyVariable("sausages"),
                              "No state variable, parameter, or derived quantity named \'sausages\'.");
        TS_ASSERT_THROWS_THIS(solution.GetVariableAtIndex(3),
                              "Invalid index passed to GetVariableAtIndex().");
        for (unsigned i=0; i<solution.rGetTimes().size(); ++i)
        {
            TS_ASSERT_DELTA(ys[i], 5*solution.rGetTimes()[i], 1e-9);
            TS_ASSERT_DELTA(as[i], 5, 1e-9);
            TS_ASSERT_DELTA(2*as[i]+ys[i], two_a_plus_y[i], 1e-9);
        }

        // Check the derived quantity is written to the file properly too.
        solution.WriteToFile("OdeSolution", "ParameterisedOde", "seconds", 1, false, 4, true);
        PetscTools::Barrier("TestWithParameters");
        NumericFileComparison comparer(OutputFileHandler::GetChasteTestOutputDirectory() + "OdeSolution/ParameterisedOde.dat",
                                       "ode/test/data/ParameterisedOde.dat");
        TS_ASSERT(comparer.CompareFiles(1e-6));
    }

    void TestLastTimeStep()
    {
        Ode1 ode_system;

        // Initialise the instance of our solver class
        EulerIvpOdeSolver euler_solver;

        // Initialise the instance of our solution class
        OdeSolution solutions;

        // Solving the ODE problem. Note that dt and the sampling time are different
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.000037, 0.000037);

        int last = solutions.GetNumberOfTimeSteps();

        // Test to see if this worked
        double testvalue = solutions.rGetSolutions()[last-1][0];

        TS_ASSERT_DELTA(testvalue, 2.0, 0.001);
    }

    void TestGlobalError()
    {
        OdeFirstOrder ode_system;

        double h_value = 0.01;

        // Euler solver solution worked out
        EulerIvpOdeSolver euler_solver;
        OdeSolution solutions_euler;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions_euler = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last = solutions_euler.GetNumberOfTimeSteps();
        double testvalue_euler = solutions_euler.rGetSolutions()[last][0];

        //R unge Kutta 2 solver solution worked out
        RungeKutta2IvpOdeSolver rk2_solver;
        OdeSolution solutions_rk2;

        state_variables = ode_system.GetInitialConditions();
        solutions_rk2 = rk2_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last2 = solutions_rk2.GetNumberOfTimeSteps();
        double testvalue_rk2 = solutions_rk2.rGetSolutions()[last2][0];

        // Runge Kutta 4 solver solution worked out
        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution solutions_rk4;

        state_variables = ode_system.GetInitialConditions();
        solutions_rk4 = rk4_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last3 = solutions_rk4.GetNumberOfTimeSteps();
        double testvalue_rk4 = solutions_rk4.rGetSolutions()[last3][0];

        // The tests
        double exact_solution = exp(2.0);

        double global_error_euler;
        global_error_euler = 0.5*exp(2.0)*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_euler, exact_solution, global_error_euler);

        double global_error_rk2;
        global_error_rk2 = (1.0/6.0)*h_value*exp(2.0)*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk2, exact_solution, global_error_rk2);

        double global_error_rk4;
        global_error_rk4 = (1.0/24.0)*pow((double)h_value,(double)3)*exp(2.0)*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk4, exact_solution, global_error_rk4);

        //The info files should now contain the ODE solver identifiers
        solutions_rk2.WriteToFile("OdeSolution", "RungeKutta2", "time", 1, false, 4, true);
        solutions_rk4.WriteToFile("OdeSolution", "RungeKutta4", "time", 1, false, 4, true);
        solutions_euler.WriteToFile("OdeSolution", "Euler", "time", 1, false, 4, true);
        PetscTools::Barrier("TestGlobalError"); // Ensure files exist

        {
            FileFinder generated_file("OdeSolution/RungeKutta2.info", RelativeTo::ChasteTestOutput);
            FileFinder reference_file("ode/test/data/RungeKutta2.info", RelativeTo::ChasteSourceRoot);
            FileComparison file_comparer(generated_file, reference_file);
            TS_ASSERT(file_comparer.CompareFiles());
        }
        {
            FileFinder generated_file("OdeSolution/RungeKutta4.info", RelativeTo::ChasteTestOutput);
            FileFinder reference_file("ode/test/data/RungeKutta4.info", RelativeTo::ChasteSourceRoot);
            FileComparison file_comparer(generated_file, reference_file);
            TS_ASSERT(file_comparer.CompareFiles());
        }
        {
            FileFinder generated_file("OdeSolution/Euler.info", RelativeTo::ChasteTestOutput);
            FileFinder reference_file("ode/test/data/Euler.info", RelativeTo::ChasteSourceRoot);
            FileComparison file_comparer(generated_file, reference_file);
            TS_ASSERT(file_comparer.CompareFiles());
        }
    }

    void TestGlobalErrorSystemOf2Equations()
    {
        OdeSecondOrder ode_system;

        double h_value = 0.01;

        // Euler solver solution worked out
        EulerIvpOdeSolver euler_solver;
        OdeSolution solutions_euler;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions_euler = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last = solutions_euler.GetNumberOfTimeSteps();

        double testvalue_euler[2];
        testvalue_euler[0] = solutions_euler.rGetSolutions()[last][0];
        testvalue_euler[1] = solutions_euler.rGetSolutions()[last][1];

        // Runge Kutta 2 solver solution worked out
        RungeKutta2IvpOdeSolver rk2_solver;
        OdeSolution solutions_rk2;

        state_variables = ode_system.GetInitialConditions();
        solutions_rk2 = rk2_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last2 = solutions_rk2.GetNumberOfTimeSteps();

        double testvalue_rk2[2];
        testvalue_rk2[0] = solutions_rk2.rGetSolutions()[last2][0];
        testvalue_rk2[1] = solutions_rk2.rGetSolutions()[last2][1];

        // Runge Kutta 4 solver solution worked out
        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution solutions_rk4;

        state_variables = ode_system.GetInitialConditions();
        solutions_rk4 = rk4_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last3 = solutions_rk4.GetNumberOfTimeSteps();
        double testvalue_rk4[2];
        testvalue_rk4[0] = solutions_rk4.rGetSolutions()[last3][0];
        testvalue_rk4[1] = solutions_rk4.rGetSolutions()[last3][1];

        // The tests
        double exact_solution[2];

        exact_solution[0] = sin(2.0);
        exact_solution[1] = cos(2.0);

        double global_error_euler;
        global_error_euler = 0.5*1*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_euler[0], exact_solution[0], global_error_euler);
        TS_ASSERT_DELTA(testvalue_euler[1], exact_solution[1], global_error_euler);

        double global_error_rk2;
        global_error_rk2 = (1.0/6.0)*h_value*1*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk2[0], exact_solution[0],global_error_rk2);
        TS_ASSERT_DELTA(testvalue_rk2[1], exact_solution[1], global_error_rk2);

        double global_error_rk4;
        global_error_rk4 = (1.0/24.0)*pow((double)h_value,(double)3)*1*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk4[0], exact_solution[0], global_error_rk4);
        TS_ASSERT_DELTA(testvalue_rk4[1], exact_solution[1], global_error_rk4);
    }

    void TestGlobalErrorSystemOf3Equations()
    {
        OdeThirdOrder ode_system;

        double h_value = 0.01;

        // Euler solver solution worked out
        EulerIvpOdeSolver euler_solver;
        OdeSolution solutions_euler;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions_euler = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last = solutions_euler.GetNumberOfTimeSteps();

        double testvalue_euler[3];
        testvalue_euler[0] = solutions_euler.rGetSolutions()[last][0];
        testvalue_euler[1] = solutions_euler.rGetSolutions()[last][1];
        testvalue_euler[2] = solutions_euler.rGetSolutions()[last][2];

        // Runge Kutta 2 solver solution worked out
        RungeKutta2IvpOdeSolver rk2_solver;
        OdeSolution solutions_rk2;

        state_variables = ode_system.GetInitialConditions();
        solutions_rk2 = rk2_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last2 = solutions_rk2.GetNumberOfTimeSteps();

        double testvalue_rk2[3];
        testvalue_rk2[0] = solutions_rk2.rGetSolutions()[last2][0];
        testvalue_rk2[1] = solutions_rk2.rGetSolutions()[last2][1];
        testvalue_rk2[2] = solutions_rk2.rGetSolutions()[last2][2];

        // Runge Kutta 4 solver solution worked out
        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution solutions_rk4;

        state_variables = ode_system.GetInitialConditions();
        solutions_rk4 = rk4_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last3 = solutions_rk4.GetNumberOfTimeSteps();

        double testvalue_rk4[3];
        testvalue_rk4[0] = solutions_rk4.rGetSolutions()[last3][0];
        testvalue_rk4[1] = solutions_rk4.rGetSolutions()[last3][1];
        testvalue_rk4[2] = solutions_rk4.rGetSolutions()[last3][2];

        // The tests
        double exact_solution[3];

        exact_solution[0] = -sin(2.0);
        exact_solution[1] = sin(2.0)+cos(2.0);
        exact_solution[2] = 2*sin(2.0);

        double global_error_euler;
        global_error_euler = 0.5*2*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_euler[0], exact_solution[0], global_error_euler);
        TS_ASSERT_DELTA(testvalue_euler[1], exact_solution[1], global_error_euler);
        TS_ASSERT_DELTA(testvalue_euler[2], exact_solution[2], global_error_euler);

        double global_error_rk2;
        global_error_rk2 = (1.0/6.0)*h_value*2*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk2[0], exact_solution[0], global_error_rk2);
        TS_ASSERT_DELTA(testvalue_rk2[1], exact_solution[1], global_error_rk2);
        TS_ASSERT_DELTA(testvalue_rk2[2], exact_solution[2], global_error_rk2);

        double global_error_rk4;
        global_error_rk4 = (1.0/24.0)*pow((double)h_value,(double)3)*2*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk4[0], exact_solution[0], global_error_rk4);
        TS_ASSERT_DELTA(testvalue_rk4[1], exact_solution[1], global_error_rk4);
        TS_ASSERT_DELTA(testvalue_rk4[2], exact_solution[2], global_error_rk4);
     }

    void TestGlobalError2()
    {
        Ode4 ode_system;

        double h_value = 0.001;

        // Euler solver solution worked out
        EulerIvpOdeSolver euler_solver;
        OdeSolution solutions_euler;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions_euler = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last = solutions_euler.GetNumberOfTimeSteps();
        double testvalue_euler = solutions_euler.rGetSolutions()[last][0];

        // Runge Kutta 2 solver solution worked out
        RungeKutta2IvpOdeSolver rk2_solver;
        OdeSolution solutions_rk2;

        state_variables = ode_system.GetInitialConditions();
        solutions_rk2 = rk2_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last2 = solutions_rk2.GetNumberOfTimeSteps();
        double testvalue_rk2 = solutions_rk2.rGetSolutions()[last2][0];

        // Runge Kutta 4 solver solution worked out
        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution solutions_rk4;

        state_variables = ode_system.GetInitialConditions();
        solutions_rk4 = rk4_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last3 = solutions_rk4.GetNumberOfTimeSteps();
        double testvalue_rk4 = solutions_rk4.rGetSolutions()[last3][0];

        // The tests
        double alpha = 100;
        double exact_solution = 1/(1+exp(-alpha*2));

        double global_error_euler;
        global_error_euler = 0.5*1/(1+exp(-alpha*2))*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_euler, exact_solution, global_error_euler);

        double global_error_rk2;
        global_error_rk2 = (1.0/6.0)*h_value*1/(1+exp(-alpha*2))*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk2, exact_solution, global_error_rk2);

        double global_error_rk4;
        global_error_rk4 = (1.0/24.0)*pow((double)h_value,(double)3)*1/(1+exp(-alpha*2))*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk4, exact_solution, global_error_rk4);
    }

    void TestArchivingSolvers()
    {
        OutputFileHandler handler("archive",false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("ode_solver.arch"); // PLEASE DO NOT COPY-PASTE THIS FILENAME
        // Archive
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set up a solver
            AbstractIvpOdeSolver* const p_euler_ivp_ode_solver = new EulerIvpOdeSolver;
            AbstractIvpOdeSolver* const p_runge_kutta_2_ode_solver = new RungeKutta2IvpOdeSolver;
            AbstractIvpOdeSolver* const p_runge_kutta_4_ode_solver = new RungeKutta4IvpOdeSolver;

            // Should always archive a pointer
            output_arch << p_euler_ivp_ode_solver;
            output_arch << p_runge_kutta_2_ode_solver;
            output_arch << p_runge_kutta_4_ode_solver;

            // Free memory
            delete p_euler_ivp_ode_solver;
            delete p_runge_kutta_2_ode_solver;
            delete p_runge_kutta_4_ode_solver;
        }

        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Create a pointer
            AbstractIvpOdeSolver* p_euler;
            AbstractIvpOdeSolver* p_rk2;
            AbstractIvpOdeSolver* p_rk4;
            input_arch >> p_euler;
            input_arch >> p_rk2;
            input_arch >> p_rk4;

            // Check the solver now has the properties of the one archived above.
            Ode1 ode_system;
            p_euler->SolveAndUpdateStateVariable(&ode_system, 0, 1, 0.01);
            TS_ASSERT_DELTA(ode_system.rGetStateVariables()[0], 1.0, 1e-2);

            Ode1 ode_system_2;
            p_rk2->SolveAndUpdateStateVariable(&ode_system_2, 0, 1, 0.01);
            TS_ASSERT_DELTA(ode_system_2.rGetStateVariables()[0], 1.0, 1e-2);

            Ode1 ode_system_3;
            p_rk4->SolveAndUpdateStateVariable(&ode_system_3, 0, 1, 0.01);
            TS_ASSERT_DELTA(ode_system_3.rGetStateVariables()[0], 1.0, 1e-2);

            delete p_euler;
            delete p_rk2;
            delete p_rk4;
        }
     }
};

#endif //_TESTABSTRACTIVPODESOLVER_HPP_
