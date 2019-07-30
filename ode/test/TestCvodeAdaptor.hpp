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

#ifndef _TESTCVODEADAPTOR_HPP_
#define _TESTCVODEADAPTOR_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include <iostream>
#include <cmath>

#include "CvodeAdaptor.hpp"
#include "OutputFileHandler.hpp"
#include "OdeSolution.hpp"
#include "PetscTools.hpp"

#include "Ode1.hpp"
#include "OdeFirstOrder.hpp"
#include "OdeSecondOrder.hpp"
#include "OdeSecondOrderWithEvents.hpp"

#include "PetscSetupAndFinalize.hpp"

class OdeWithRootFunction : public OdeSecondOrderWithEvents
{
public:
    OdeWithRootFunction()
        : OdeSecondOrderWithEvents()
    {
    }

    double CalculateRootFunction(double time, const std::vector<double>& rY)
    {
        return rY[0];
    }
};


class ExceptionalOdeWithRootFunction : public OdeSecondOrderWithEvents
{
private:
    bool mNice;
public:
    ExceptionalOdeWithRootFunction()
        : OdeSecondOrderWithEvents()
    {
        mNice = false;
    }

    void BeNice()
    {
        mNice = true;
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] =  rY[1];
        rDY[1] = -rY[0];
        if (!mNice)
        {
            EXCEPTION("I'm feeling nasty!");
        }
    }

    double CalculateRootFunction(double time, const std::vector<double>& rY)
    {
        EXCEPTION("I'm feeling nasty!");
        return rY[0];
    }
};

class TestCvodeAdaptor: public CxxTest::TestSuite
{
private:

#ifdef CHASTE_CVODE
    void HelperTestOde1(double startTime, double endTime, double samplingTime)
    {
        // Initialise
        Ode1 ode_system;
        OdeSolution solutions;
        CvodeAdaptor solver;

        // Solving the ODE problem.
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = solver.Solve(&ode_system, state_variables, startTime, endTime, endTime-startTime, samplingTime);

        int num_timesteps = solutions.GetNumberOfTimeSteps();

        // The number of timesteps should be (just about) equal to sim_time/sampling_time
        TS_ASSERT_DELTA(num_timesteps, (endTime-startTime)/samplingTime, 1);

        // also check the size of the data is correct
        TS_ASSERT_EQUALS(solutions.rGetSolutions().size(),
                         (unsigned) (num_timesteps+1));

        int last = num_timesteps;

        // Test the solution is correct
        double testvalue = solutions.rGetSolutions()[last][0];

        // Exact solution of Ode1 is y=t-t0
        TS_ASSERT_DELTA(testvalue, endTime-startTime, 0.01);

        // Test second version of Solve
        ode_system.SetStateVariables(ode_system.GetInitialConditions());
        state_variables = ode_system.rGetStateVariables();
        solver.Solve(&ode_system, state_variables, startTime, endTime, endTime-startTime);
        TS_ASSERT_DELTA(state_variables[0], endTime-startTime, 0.01);

        // No stopping event was specified in the ODE, so check the
        // solver correctly states it didn't stop due to a
        // stopping event.
        TS_ASSERT_EQUALS(solver.StoppingEventOccurred(), false);
    }
#endif // CHASTE_CVODE

public:
    void TestBasics()
    {
#ifdef CHASTE_CVODE
        CvodeAdaptor solver;
        solver.SetMaxSteps(1000);
        TS_ASSERT_EQUALS(solver.GetMaxSteps(), 1000);

        TS_ASSERT_DELTA(solver.GetRelativeTolerance(), 1e-4, 1e-12);
        TS_ASSERT_DELTA(solver.GetAbsoluteTolerance(), 1e-6, 1e-12);

        solver.SetTolerances(1e-5, 1e-5);
        TS_ASSERT_DELTA(solver.GetRelativeTolerance(), 1e-5, 1e-12);
        TS_ASSERT_DELTA(solver.GetAbsoluteTolerance(), 1e-5, 1e-12);
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestOnOde1()
    {
#ifdef CHASTE_CVODE
        HelperTestOde1(0.0, 2.0, 0.001);
        HelperTestOde1(1.0, 2.0, 0.01);
        HelperTestOde1(-1.0, 2.0, 2);
        HelperTestOde1(0.0, 0.4, 0.34);
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestGlobalError()
    {
#ifdef CHASTE_CVODE
        OdeFirstOrder ode_system;

        double h_value = 0.01;

        CvodeAdaptor solver;
        solver.SetMaxSteps(1000);
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, 0.1);
        int last = solutions.GetNumberOfTimeSteps();
        double testvalue = solutions.rGetSolutions()[last][0];

        // The tests
        double exact_solution = exp(2.0);

        /// \todo: #890 Work out what the global error should be bounded by
        double global_error = 1e-3;

        TS_ASSERT_DELTA(testvalue, exact_solution, global_error);
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestMultipleCalls()
    {
#ifdef CHASTE_CVODE
        OdeFirstOrder ode_system;

        double h_value = 0.01;

        CvodeAdaptor solver;
        solver.SetMaxSteps(1000);
        OdeSolution solutions;
        std::vector<double> state_variables = ode_system.GetInitialConditions();

        double exact_solution = exp(2.0);
        double global_error = 1e-3;

        {   // Solve in two steps with reset
            solver.SetForceReset(true);
            solver.Solve(&ode_system, state_variables, 0.0, 1.0, h_value, 0.1);
            solutions = solver.Solve(&ode_system, state_variables, 1.0, 2.0, h_value, 0.1);
            int last = solutions.GetNumberOfTimeSteps();
            double testvalue = solutions.rGetSolutions()[last][0];

            TS_ASSERT_DELTA(testvalue, exact_solution, global_error);
        }

        state_variables = ode_system.GetInitialConditions();
        solver.ResetSolver();

        {   // Solve in two steps without reset
            solver.Solve(&ode_system, state_variables, 0.0, 1.0, h_value, 0.1);
            solutions = solver.Solve(&ode_system, state_variables, 1.0, 2.0, h_value, 0.1);
            int last = solutions.GetNumberOfTimeSteps();
            double testvalue = solutions.rGetSolutions()[last][0];

            TS_ASSERT_DELTA(testvalue, exact_solution, global_error);
        }

        state_variables = ode_system.GetInitialConditions();
        solver.ResetSolver();

        // Altering the state variables should also force a reset
        {
            solver.Solve(&ode_system, state_variables, 0.0, 1.0, h_value, 0.1);
            state_variables[0] += 0.00001;
            solutions = solver.Solve(&ode_system, state_variables, 1.0, 2.0, h_value, 0.1);
            int last = solutions.GetNumberOfTimeSteps();
            double testvalue = solutions.rGetSolutions()[last][0];

            TS_ASSERT_DELTA(testvalue, exact_solution, global_error);
        }

        state_variables = ode_system.GetInitialConditions();
        solver.SetForceReset(true);

        {   // Solve in two steps with reset again
            solver.Solve(&ode_system, state_variables, 0.0, 1.0, h_value, 0.1);
            solutions = solver.Solve(&ode_system, state_variables, 1.0, 2.0, h_value, 0.1);
            int last = solutions.GetNumberOfTimeSteps();
            double testvalue = solutions.rGetSolutions()[last][0];

            TS_ASSERT_DELTA(testvalue, exact_solution, global_error);
        }

#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestGlobalErrorSystemOf2Equations()
    {
#ifdef CHASTE_CVODE
        OdeSecondOrder ode_system;

        double h_value = 0.01;

        CvodeAdaptor solver;
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, 0.1);
        int last = solutions.GetNumberOfTimeSteps();

        double testvalue[2];
        testvalue[0] = solutions.rGetSolutions()[last][0];
        testvalue[1] = solutions.rGetSolutions()[last][1];

        double exact_solution[2];
        exact_solution[0] = sin(2.0);
        exact_solution[1] = cos(2.0);
        double global_error = 1e-3;

        TS_ASSERT_DELTA(testvalue[0], exact_solution[0], global_error);
        TS_ASSERT_DELTA(testvalue[1], exact_solution[1], global_error);
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestWithStoppingEvent()
    {
#ifdef CHASTE_CVODE
        OdeSecondOrderWithEvents ode_system;
        CvodeAdaptor solver;
        OdeSolution solutions;

        solver.CheckForStoppingEvents();
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        const double sampling_time = 0.01;
        solutions = solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.1, sampling_time);
        int num_timesteps = solutions.GetNumberOfTimeSteps();

        // Final time should be about pi/2
        double tend = solutions.rGetTimes()[num_timesteps];
        double tstop = solver.GetStoppingTime();
        TS_ASSERT_DELTA( tend, M_PI_2, 0.01);
        TS_ASSERT_EQUALS(tend, tstop);

        // Solver should correctly state the stopping event occurred
        TS_ASSERT_EQUALS(solver.StoppingEventOccurred(), true);
        // ODE system should agree that this is a stopping event
        TS_ASSERT(ode_system.CalculateStoppingEvent(tend, solutions.rGetSolutions()[num_timesteps]));
        // i.e. penultimate y0 should be greater than zero
        TS_ASSERT_LESS_THAN( 0, solutions.rGetSolutions()[num_timesteps-1][0]);
        // and final y0 should be less than zero
        TS_ASSERT_LESS_THAN( solutions.rGetSolutions()[num_timesteps][0], 0);
#if CHASTE_SUNDIALS_VERSION >= 20400
        TS_ASSERT_DELTA(solver.GetLastStepSize(), sampling_time, 1e-6);
#else
        // I'm not convinced that this older version of CVODE is doing the correct thing here...
        TS_ASSERT_DELTA(solver.GetLastStepSize(), 0.0999, 1e-4);
#endif
        std::cout << "1st with Exception\n"<< std::flush;
        // If we try to continue, the stopping event is still true, which is an error
        TS_ASSERT_THROWS_THIS(solutions = solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.1, 0.01),
                "(Solve with sampling) Stopping event is true for initial condition");

        // Alternative Solve method
        state_variables = ode_system.GetInitialConditions();

        std::cout << "2nd Solve\n"<< std::flush;
        solver.Solve(&ode_system, state_variables, 0.0, 2.0, sampling_time);
        TS_ASSERT_EQUALS(solver.StoppingEventOccurred(), true);
        TS_ASSERT_DELTA(solver.GetStoppingTime(), M_PI_2, 0.01);
        std::cout << "2nd with Exception\n"<< std::flush;
        // If we try to continue, the stopping event is still true, which is an error
        TS_ASSERT_THROWS_THIS(solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.1),
                "(Solve) Stopping event is true for initial condition");
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestWithRootFunction()
    {
#ifdef CHASTE_CVODE
        OdeWithRootFunction ode_system;
        CvodeAdaptor solver;
        OdeSolution solutions;

        solver.CheckForStoppingEvents();
        std::vector<double> state_variables = ode_system.GetInitialConditions();

        const double sampling_time = 0.01;
        solutions = solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.1, sampling_time);
        int num_timesteps = solutions.GetNumberOfTimeSteps();

        // Final time should be about pi/2
        double tend = solutions.rGetTimes()[num_timesteps];
        double tstop = solver.GetStoppingTime();
        TS_ASSERT_DELTA(tend, M_PI_2, 1e-4);
        TS_ASSERT_EQUALS(tend, tstop);

        // Solver should correctly state the stopping event occurred
        TS_ASSERT_EQUALS(solver.StoppingEventOccurred(), true);
        // ODE system should agree that this is a root
        TS_ASSERT_DELTA(ode_system.CalculateRootFunction(tend, solutions.rGetSolutions()[num_timesteps]),0,1e-12);
        // i.e. penultimate y0 should be greater than zero
        TS_ASSERT_LESS_THAN( 0, solutions.rGetSolutions()[num_timesteps-1][0]);
        // and final y0 should be less than zero
        TS_ASSERT_LESS_THAN( solutions.rGetSolutions()[num_timesteps][0], 0);

#if CHASTE_SUNDIALS_VERSION >= 20400
        TS_ASSERT_DELTA(solver.GetLastStepSize(), sampling_time, 1e-6);
#else
        // I'm not convinced that this older version of CVODE is doing the correct thing here...
        TS_ASSERT_DELTA(solver.GetLastStepSize(), 0.0999, 1e-4);
#endif

        // Alternative Solve method
        state_variables = ode_system.GetInitialConditions();
        solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.1);
        TS_ASSERT_EQUALS(solver.StoppingEventOccurred(), true);
        TS_ASSERT_DELTA(solver.GetStoppingTime(), M_PI_2, 1e-4)
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestExceptions()
    {
#ifdef CHASTE_CVODE
        ExceptionalOdeWithRootFunction ode_system;
        CvodeAdaptor solver;
        OdeSolution solutions;

        // Exception in EvaluateYDerivatives
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        TS_ASSERT_THROWS_THIS(solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.1),
                "CVODE failed to solve system: CV_RHSFUNC_FAIL");

        // Not enough steps
        ode_system.BeNice();
        solver.SetMaxSteps(1);
        state_variables = ode_system.GetInitialConditions();
        TS_ASSERT_THROWS_THIS(solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.1),
                "CVODE failed to solve system: CV_TOO_MUCH_WORK");

        //Try again with time sampling
        TS_ASSERT_THROWS_THIS(solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.1, 0.1),
                "CVODE failed to solve system: CV_TOO_MUCH_WORK");

        // Exception in root function
        solver.CheckForStoppingEvents();
        state_variables = ode_system.GetInitialConditions();
        TS_ASSERT_THROWS_THIS(solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.1),
                "CVODE failed to solve system: CV_RTFUNC_FAIL");
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestArchivingCvodeAdaptorSolver()
    {
        EXIT_IF_PARALLEL;

#ifdef CHASTE_CVODE
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "cvode_adaptor_solver.arch";

        // Create and archive simulation time
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set up a solver
            AbstractIvpOdeSolver* const p_cvode_adaptor_solver = new CvodeAdaptor;


            Ode1 ode_system;
            p_cvode_adaptor_solver->SolveAndUpdateStateVariable(&ode_system, 0, 1, 0.01);
            TS_ASSERT_DELTA(ode_system.rGetStateVariables()[0], 1.0, 1e-2);

            // Should always archive a pointer
            output_arch << p_cvode_adaptor_solver;

            // Change stimulus a bit
            delete p_cvode_adaptor_solver;
        }

        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Create a pointer
            AbstractIvpOdeSolver* p_cvode_solver;
            input_arch >> p_cvode_solver;

            delete p_cvode_solver;
        }
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }
};

#endif //_TESTCVODEADAPTOR_HPP_
