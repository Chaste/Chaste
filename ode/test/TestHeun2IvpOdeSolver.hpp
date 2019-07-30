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

#ifndef _TESTHEUN2IVPODESOLVER_HPP_
#define _TESTHEUN2IVPODESOLVER_HPP_

#include <cxxtest/TestSuite.h>

#include <iostream>
#include <sstream>

#include "CheckpointArchiveTypes.hpp"

#include "HeunIvpOdeSolver.hpp"
#include "Ode1.hpp"
#include "Ode2.hpp"
#include "Ode4.hpp"
#include "OdeFirstOrder.hpp"
#include "OdeSecondOrder.hpp"
#include "OdeSecondOrderWithEvents.hpp"
#include "OdeThirdOrder.hpp"
#include "ParameterisedOde.hpp"
#include "OutputFileHandler.hpp"
#include "ArchiveLocationInfo.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestHeun2IvpOdeSolver : public CxxTest::TestSuite
{
private:
      void MyTestGenericSolver(AbstractIvpOdeSolver& rSolver, double startTime,
                               double endTime, double dt, double samplingTime)
      {
        // Initialise the instances of our ODE system and solution classes
        Ode1 ode_system;
        OdeSolution solutions;

        // Solving the ODE problem. Note that dt and the sampling time are different
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = rSolver.Solve(&ode_system, state_variables, startTime, endTime,
            dt, samplingTime);

        unsigned num_timesteps = solutions.GetNumberOfTimeSteps();

        // The number of time steps should be (just about) equal to
        // end_time/sampling_time = 2/0.01 = 200
        TS_ASSERT_DELTA(num_timesteps, (endTime-startTime)/samplingTime, 1);

        // Also check the size of the data is correct
        TS_ASSERT_EQUALS(solutions.rGetSolutions().size(), num_timesteps+1u);

        unsigned last = num_timesteps;

        // Test to solution is correct
        double testvalue = solutions.rGetSolutions()[last][0];

        // Exact solution of Ode1 is y=t-t0
        TS_ASSERT_DELTA(testvalue, endTime-startTime, 0.01);

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
        std::vector<double> state_variables =
            ode_with_events.GetInitialConditions();
        solutions = rSolver.Solve(&ode_with_events, state_variables, 0.0, 2.0,
            0.001, 0.001);

        unsigned num_timesteps = solutions.GetNumberOfTimeSteps();

        // Final time should be around pi/2
        TS_ASSERT_DELTA( solutions.rGetTimes()[num_timesteps], M_PI_2, 0.01);

        // Penultimate y0 should be greater than zero
        TS_ASSERT_LESS_THAN( 0, solutions.rGetSolutions()[num_timesteps-1][0]);

        // Final y0 should be less than zero
        TS_ASSERT_LESS_THAN( solutions.rGetSolutions()[num_timesteps][0], 0);

        // Solver should correctly state the stopping event occurred
        TS_ASSERT_EQUALS(rSolver.StoppingEventOccurred(), true);

        // This is to cover the exception when a stopping event occurs before the first timestep.
        TS_ASSERT_THROWS_ANYTHING(rSolver.Solve(&ode_with_events, state_variables, 2.0, 3.0, 0.001));

        ///////////////////////////////////////////////
        // Repeat with sampling time larger than dt
        ///////////////////////////////////////////////

        state_variables = ode_with_events.GetInitialConditions();
        solutions = rSolver.Solve(&ode_with_events, state_variables, 0.0, 2.0,
            0.001, 0.01);

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
        bad_init_cond.push_back(-1); //y0 < 0 so stopping event true
        bad_init_cond.push_back(0.0);
        TS_ASSERT_THROWS_ANYTHING(rSolver.Solve(&ode_with_events, bad_init_cond, 0.0, 2.0, 0.001, 0.01));
    }

public:
    void TestHeunSolver()
    {
        HeunIvpOdeSolver heun_solver;

        MyTestGenericSolver(heun_solver, 0.0, 2.0, 0.001, 0.001);
        MyTestGenericSolver(heun_solver, 1.0, 2.0, 0.001, 0.01);
        MyTestGenericSolver(heun_solver, -1.0, 2.0, 0.001, 2);
        MyTestGenericSolver(heun_solver, 0.0, 0.4, 0.01, 0.34);

        MyTestSolverOnOdesWithEvents(heun_solver);

        // Test SolveAndUpdateStateVariable()
        Ode1 ode_system;
        heun_solver.SolveAndUpdateStateVariable(&ode_system, 0, 1, 0.01);
        TS_ASSERT_DELTA(ode_system.rGetStateVariables()[0], 1.0, 1e-2);
    }

    void TestGlobalError()
    {
        OdeFirstOrder ode_system;

        double h_value = 0.01;

        //Heun solver solution worked out
        HeunIvpOdeSolver heun_solver;
        OdeSolution solutions_heun;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions_heun = heun_solver.Solve(&ode_system, state_variables, 0.0, 2.0,
            h_value, h_value);
        unsigned last2 = solutions_heun.GetNumberOfTimeSteps();
        double testvalue_heun = solutions_heun.rGetSolutions()[last2][0];

        // The tests
        double exact_solution = exp(2.0);

        double global_error_heun;
        global_error_heun = (1.0 / 6.0) * h_value * exp(2.0) * (exp(2.0) - 1)
            * h_value;
        TS_ASSERT_DELTA(testvalue_heun, exact_solution, global_error_heun);
    }

    void TestGlobalErrorSystemOf2Equations()
    {
        OdeSecondOrder ode_system;

        double h_value = 0.01;

        //Heun solver solution worked out
        HeunIvpOdeSolver heun_solver;
        OdeSolution solutions_heun;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions_heun = heun_solver.Solve(&ode_system, state_variables, 0.0, 2.0,
            h_value, h_value);
        unsigned last2 = solutions_heun.GetNumberOfTimeSteps();

        double testvalue_heun[2];
        testvalue_heun[0] = solutions_heun.rGetSolutions()[last2][0];
        testvalue_heun[1] = solutions_heun.rGetSolutions()[last2][1];

        // The tests
        double exact_solution[2];

        exact_solution[0] = sin(2.0);
        exact_solution[1] = cos(2.0);

        double global_error_heun;
        global_error_heun = (1.0 / 6.0) * h_value * 1 * (exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_heun[0], exact_solution[0],global_error_heun);
        TS_ASSERT_DELTA(testvalue_heun[1], exact_solution[1], global_error_heun);
    }

    void TestGlobalErrorSystemOf3Equations()
    {
        OdeThirdOrder ode_system;

        double h_value = 0.01;

        //Heun solver solution worked out
        HeunIvpOdeSolver heun_solver;
        OdeSolution solutions_heun;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions_heun = heun_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        unsigned last2 = solutions_heun.GetNumberOfTimeSteps();

        double testvalue_heun[3];
        testvalue_heun[0] = solutions_heun.rGetSolutions()[last2][0];
        testvalue_heun[1] = solutions_heun.rGetSolutions()[last2][1];
        testvalue_heun[2] = solutions_heun.rGetSolutions()[last2][2];

        // The tests
        double exact_solution[3];

        exact_solution[0] = -sin(2.0);
        exact_solution[1] = sin(2.0)+cos(2.0);
        exact_solution[2] = 2*sin(2.0);

        double global_error_heun;
        global_error_heun = (1.0/6.0)*h_value*2*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_heun[0], exact_solution[0], global_error_heun);
        TS_ASSERT_DELTA(testvalue_heun[1], exact_solution[1], global_error_heun);
        TS_ASSERT_DELTA(testvalue_heun[2], exact_solution[2], global_error_heun);
    }

    void TestGlobalError2()
    {
        Ode4 ode_system;

        double h_value = 0.001;

        //Heun solver solution worked out
        HeunIvpOdeSolver heun_solver;
        OdeSolution solutions_heun;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions_heun = heun_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        unsigned last2 = solutions_heun.GetNumberOfTimeSteps();
        double testvalue_heun = solutions_heun.rGetSolutions()[last2][0];

        // The tests
        double alpha = 100;
        double exact_solution = 1/(1+exp(-alpha*2));

        double global_error_heun;
        global_error_heun = (1.0/6.0)*h_value*1/(1+exp(-alpha*2))*(exp(2.0)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_heun, exact_solution, global_error_heun);
    }

    // Test the order of the method by comparing two solutions using dt/2 for the second.
    void TestOrderOnSimpleSystem()
    {
        OdeFirstOrder ode_system1, ode_system2;

        double h_value1 = 0.01;
        double h_value2 = 0.005;
        HeunIvpOdeSolver heun_solver1, heun_solver2;
        OdeSolution solutions_heun1,solutions_heun2;

        std::vector<double> state_variables1 = ode_system1.GetInitialConditions();

        // Solution with ordinary time step
        solutions_heun1 = heun_solver1.Solve(&ode_system1, state_variables1, 0.0, 2.0, h_value1, h_value1);
        unsigned last1 = solutions_heun1.GetNumberOfTimeSteps();
        double testvalue_heun1 = solutions_heun1.rGetSolutions()[last1][0];

        std::vector<double> state_variables2 = ode_system2.GetInitialConditions();

        // Solution with halved time step
        solutions_heun2 = heun_solver2.Solve(&ode_system2, state_variables2, 0.0, 2.0, h_value2, h_value1);
        unsigned last2 = solutions_heun2.GetNumberOfTimeSteps();
        double testvalue_heun2 = solutions_heun2.rGetSolutions()[last2][0];

        double heun_order = 2;

        // Test that the error is going down by a factor of 2^order
        double exact_solution = exp(2.0);
        double error1 = testvalue_heun1 - exact_solution;
        double error2 = testvalue_heun2 - exact_solution;
        TS_ASSERT_DELTA(error1/error2, pow(2,heun_order), 1e-1);
    }

    // Test the order of the method by comparing two solutions using dt/2 for the second.
    void TestOrderOnMediumSystem()
    {
        OdeSecondOrder ode_system1, ode_system2;

        double h_value1 = 0.01;
        double h_value2 = 0.005;
        HeunIvpOdeSolver heun_solver1, heun_solver2;
        OdeSolution solutions_heun1,solutions_heun2;

        std::vector<double> state_variables1 = ode_system1.GetInitialConditions();

        // Solution with ordinary time step
        solutions_heun1 = heun_solver1.Solve(&ode_system1, state_variables1, 0.0, 2.0, h_value1, h_value1);
        unsigned last1 = solutions_heun1.GetNumberOfTimeSteps();
        double testvalue_heun1 = solutions_heun1.rGetSolutions()[last1][0];

        std::vector<double> state_variables2 = ode_system2.GetInitialConditions();
        // Solution with halved time step
        solutions_heun2 = heun_solver2.Solve(&ode_system2, state_variables2, 0.0, 2.0, h_value2, h_value1);
        unsigned last2 = solutions_heun2.GetNumberOfTimeSteps();
        double testvalue_heun2 = solutions_heun2.rGetSolutions()[last2][0];

        unsigned heun_order = 2;

        // Test that the error is going down by a factor of 2^order
        double exact_solution = sin(2.0);
        double error1 = testvalue_heun1 - exact_solution;
        double error2 = testvalue_heun2 - exact_solution;
        TS_ASSERT_DELTA(error1/error2, pow(2.0,(double)heun_order), 1e-1);
    }

    // Test the order of the method by comparing two solutions using dt/2 for the second.
    void TestOrderOnHardSystem()
    {
        OdeThirdOrder ode_system1, ode_system2;

        double h_value1 = 0.01;
        double h_value2 = 0.005;
        HeunIvpOdeSolver heun_solver1, heun_solver2;
        OdeSolution solutions_heun1,solutions_heun2;

        std::vector<double> state_variables1 = ode_system1.GetInitialConditions();

        // Solution with ordinary time step
        solutions_heun1 = heun_solver1.Solve(&ode_system1, state_variables1, 0.0, 2.0, h_value1, h_value1);
        unsigned last1 = solutions_heun1.GetNumberOfTimeSteps();
        double testvalue_heun1 = solutions_heun1.rGetSolutions()[last1][0];

        std::vector<double> state_variables2 = ode_system2.GetInitialConditions();

        // Solution with halved time step
        solutions_heun2 = heun_solver2.Solve(&ode_system2, state_variables2, 0.0, 2.0, h_value2, h_value1);
        unsigned last2 = solutions_heun2.GetNumberOfTimeSteps();
        double testvalue_heun2 = solutions_heun2.rGetSolutions()[last2][0];

        unsigned heun_order = 2;

        // Test that the error is going down by a factor of 2^order
        double exact_solution = -sin(2.0);
        double error1 = testvalue_heun1 - exact_solution;
        double error2 = testvalue_heun2 - exact_solution;
        TS_ASSERT_DELTA(error1/error2, pow(2.0,(double)heun_order), 1e-1);
    }

    void TestArchivingSolvers()
    {
        OutputFileHandler handler("archive",false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("heun_ode_solver.arch");

        // Archive
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set up a solver
            AbstractIvpOdeSolver* const p_heun_ode_solver = new HeunIvpOdeSolver;

            // Should always archive a pointer
            output_arch << p_heun_ode_solver;

            // Free memory
            delete p_heun_ode_solver;
        }

        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Create a pointer
            AbstractIvpOdeSolver *p_heun;
            input_arch >> p_heun;

            // Check the solver now has the properties of the one archived above.
            Ode1 ode_system;
            p_heun->SolveAndUpdateStateVariable(&ode_system, 0, 1, 0.01);
            TS_ASSERT_DELTA(ode_system.rGetStateVariables()[0], 1.0, 1e-2);

            Ode1 ode_system_2;
            p_heun->SolveAndUpdateStateVariable(&ode_system_2, 0, 1, 0.01);
            TS_ASSERT_DELTA(ode_system_2.rGetStateVariables()[0], 1.0, 1e-2);

            Ode1 ode_system_3;
            p_heun->SolveAndUpdateStateVariable(&ode_system_3, 0, 1, 0.01);
            TS_ASSERT_DELTA(ode_system_3.rGetStateVariables()[0], 1.0, 1e-2);

            delete p_heun;
        }
    }
};

#endif //_TESTHEUN2IVPODESOLVER_HPP_
