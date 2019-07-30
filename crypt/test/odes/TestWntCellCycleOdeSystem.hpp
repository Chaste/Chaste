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

#ifndef TESTWNTCELLCYCLEODESYSTEM_HPP_
#define TESTWNTCELLCYCLEODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <vector>
#include <iostream>

#include "OutputFileHandler.hpp"
#include "Timer.hpp"

#include "WntCellCycleOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestWntCellCycleOdeSystem : public CxxTest::TestSuite
{
public:

    void TestWntCellCycleEquations()
    {
        double wnt_level = 0.0;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        WntCellCycleOdeSystem wnt_cell_cycle_system(wnt_level, p_state);

        double time = 0.0;
        std::vector<double> initial_conditions = wnt_cell_cycle_system.GetInitialConditions();

        std::vector<double> derivs(initial_conditions.size());
        wnt_cell_cycle_system.EvaluateYDerivatives(time, initial_conditions, derivs);
        TS_ASSERT_DELTA(initial_conditions[6]+initial_conditions[7], 0.0074, 1e-4);

        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from Matlab code)
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2],-9.370533804903016e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4],1.549680000000000e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[5],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8],0.0, 1e-5);

        /**
         * And the same for a high Wnt level
         */
        wnt_level = 1.0;
        boost::shared_ptr<AbstractCellMutationState> p_wt_state(new WildTypeCellMutationState);
        WntCellCycleOdeSystem wnt_cell_cycle_system2(wnt_level, p_wt_state);
        initial_conditions = wnt_cell_cycle_system2.GetInitialConditions();

        TS_ASSERT_DELTA(initial_conditions[6]+initial_conditions[7], 0.6002, 1e-4);

        wnt_cell_cycle_system2.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from Matlab code)
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2],5.845305754146019e+00, 1e-5);
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4],1.549680000000000e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[5],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8],0.0, 1e-5);

        /**
         * A test for the case mutation = 1
         * (An APC +/- mutation)
         */
        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);
        wnt_level = 1.0;
        WntCellCycleOdeSystem wnt_cell_cycle_system3(wnt_level, p_apc1);
        initial_conditions = wnt_cell_cycle_system3.GetInitialConditions();

        TS_ASSERT_DELTA(initial_conditions[6]+initial_conditions[7], 0.750207, 1e-6);

        wnt_cell_cycle_system3.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from Matlab code)
        TS_ASSERT_DELTA(derivs[0], -1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1], -5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2], 7.3260e+00, 1e-4);
        TS_ASSERT_DELTA(derivs[3], -7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4], 1.549680000000000e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[5], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8], 0.0, 1e-5);

        /**
        * A test for the case mutation = 2
        * (A beta-cat delta45 mutation)
        */
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
        wnt_level = 1.0;
        WntCellCycleOdeSystem wnt_cell_cycle_system4(wnt_level, p_bcat1);
        initial_conditions = wnt_cell_cycle_system4.GetInitialConditions();

        TS_ASSERT_DELTA(initial_conditions[6]+initial_conditions[7], 0.8001, 1e-4);

        wnt_cell_cycle_system4.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from Matlab code)
        TS_ASSERT_DELTA(derivs[0], -1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1], -5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2], 7.8190e+00, 1e-4);
        TS_ASSERT_DELTA(derivs[3], -7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4], 1.549680000000000e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[5], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8], 0.0, 1e-5);

        /**
        * A test for the case mutation = 3
        * (An APC -/- mutation)
        */
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        wnt_level = 1.0;
        WntCellCycleOdeSystem wnt_cell_cycle_system5(wnt_level, p_apc2);
        initial_conditions = wnt_cell_cycle_system5.GetInitialConditions();

        TS_ASSERT_DELTA(initial_conditions[6]+initial_conditions[7],1.0,1e-6);

        wnt_cell_cycle_system5.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from Matlab code)
        TS_ASSERT_DELTA(derivs[0], -1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1], -5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2], 9.7928e+00, 1e-4);
        TS_ASSERT_DELTA(derivs[3], -7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4], 1.549680000000000e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[5], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8], 0.0, 1e-5);
    }

    void TestWntCellCycleSolver()
    {
        double wnt_level = 1.0;
        boost::shared_ptr<AbstractCellMutationState> p_wt_state(new WildTypeCellMutationState);
        WntCellCycleOdeSystem wnt_system(wnt_level, p_wt_state);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.

        double h_value_rk4 = 1e-4;

        RungeKutta4IvpOdeSolver rk4_solver;
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        BackwardEulerIvpOdeSolver back_solver(9);

        OdeSolution solutions_rk4;
        OdeSolution solutions_rkf;

        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        Timer::Reset();
        solutions_rk4 = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value_rk4, h_value_rk4);
        Timer::Print("1. Runge-Kutta");

        double h_value_rkf = 0.1;

        initial_conditions = wnt_system.GetInitialConditions();
        Timer::Reset();
        solutions_rkf = rkf_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value_rkf, 1e-4);
        Timer::Print("2. Runge-Kutta-Fehlberg");

//        WntCellCycleOdeSystem wnt_system_2(wnt_level);
//        initial_conditions = wnt_system.GetInitialConditions();
//
//        h_value = 0.001;
//
//        Timer::Reset();
//        solutions = back_solver.Solve(&wnt_system_2, initial_conditions, 0.0, 100.0, h_value, h_value);
//        Timer::Print("1. Backward Euler");

        // Testing RK4 solution
        // Test solutions are OK for a small time increase...
        int end = solutions_rk4.rGetSolutions().size() - 1;

        // Tests the simulation is ending at the right time...(going into S phase at 5.971 hours)
        TS_ASSERT_DELTA(solutions_rk4.rGetTimes()[end], 5.971, 1e-2);

        // Proper values from Matlab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][0], 2.880603485931000e-01, 1e-3);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][1], 1.000220438771564, 1.02e-2);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][2], 2.453870380958196, 1e-3);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][3], 1.446185835615586, 1e-3);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][4], 1.383272155041549e-01, 1e-3);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][5], 4.975124378109454e-03, 1e-3);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][6]+solutions_rk4.rGetSolutions()[end][7], 6.002649406788524e-01, 1e-3);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][8], 1.00, 1e-3);

        // Testing RKF solution
        // Test solutions are OK for a small time increase...
        end = solutions_rkf.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 5.971 hours)
        TS_ASSERT_DELTA(solutions_rkf.rGetTimes()[end], 5.971, 1e-2);
        // Proper values from Matlab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][0], 2.880603485931000e-01, 1e-3);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][1], 1.000220438771564e+00, 1.02e-2);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][2], 2.453870380958196e+00, 1e-3);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][3], 1.446185835615586e+00, 1e-3);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][4], 1.383272155041549e-01, 1e-3);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][5], 4.975124378109454e-03, 1e-3);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][6]+solutions_rkf.rGetSolutions()[end][7], 6.002649406788524e-01, 1e-3);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][8], 1.00, 1e-3);
    }

    void TestWntCellCycleSolverWithAPCSingleHit()
    {
        double wnt_level = 1.0;
        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);
        WntCellCycleOdeSystem wnt_system(wnt_level, p_apc1);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.

        double h_value = 0.0001;

        RungeKutta4IvpOdeSolver rk4_solver;
        BackwardEulerIvpOdeSolver back_solver(9);

        OdeSolution solutions;
        //OdeSolution solutions2;

        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 3.94 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 4.804, 1e-2);
        // Proper values from Matlab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],2.493601889546602e-01, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],1.0, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],2.922278616458120e+00, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],1.913484075138688e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],1.583557222931072e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],4.975124378109454e-03, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6],0.375, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7],0.375, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8],1.00, 1e-3);
    }

    void TestWntCellCycleSolverWithBetaCateninHit()
    {
        double wnt_level = 0.0;
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
        WntCellCycleOdeSystem wnt_system(wnt_level, p_bcat1);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.

        double h_value = 0.001;

        RungeKutta4IvpOdeSolver rk4_solver;
        BackwardEulerIvpOdeSolver back_solver(9);

        OdeSolution solutions;
        //OdeSolution solutions2;

        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        //Timer::Reset();
        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        //Timer::Print("1. Runge-Kutta");

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 7.8 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 7.835, 1e-2);
        // Proper values from Matlab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 3.242663439545868e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 0.999, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 2.153277726022381e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3], 1.145736207245425e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4], 1.234257806221668e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5], 1.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6], 3.707764665147012e-03, 1e-5);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7], 0.5, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8], 0.00, 1e-3);
    }

    void TestWntCellCycleSolverWithAPCDoubleHit()
    {
        double wnt_level = 0.0;
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        WntCellCycleOdeSystem wnt_system(wnt_level, p_apc2);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.

        double h_value = 0.001;

        RungeKutta4IvpOdeSolver rk4_solver;
        BackwardEulerIvpOdeSolver back_solver(9);

        OdeSolution solutions;
        //OdeSolution solutions2;

        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;

        // Tests the simulation is ending at the right time...(going into S phase at 3.94 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 3.9435, 1e-2);

        // Proper values from Matlab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 2.058373151310055e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 0.999, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 3.699024514542648e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3], 2.687523235896298e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4], 1.834144555072084e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5], 0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6], 0.5, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7], 0.5, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8], 0.00, 1e-3);
    }

    void TestArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "wnt_ode.arch";

        {
            double wnt_level = 0.567;
            boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
            WntCellCycleOdeSystem ode_system(wnt_level, p_state);

            TS_ASSERT_DELTA(ode_system.GetWntLevel(), 0.567, 1e-6);
            TS_ASSERT_EQUALS(ode_system.GetMutationState()->IsType<WildTypeCellMutationState>(), true);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 9u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.7357, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[1], 0.1713, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[2], 0.0690, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[3], 0.0033, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[4], 0.0000, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[5], 0.0087, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[6], 0.2304, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[7], 0.2304, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[8], 0.5670, 1e-4);

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
            TS_ASSERT_DELTA(static_cast<WntCellCycleOdeSystem*>(p_ode_system)->GetWntLevel(), 0.567, 1e-6);
            TS_ASSERT_EQUALS(static_cast<WntCellCycleOdeSystem*>(p_ode_system)->GetMutationState()->IsType<WildTypeCellMutationState>(), true);

            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 9u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.7357, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[1], 0.1713, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[2], 0.0690, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[3], 0.0033, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[4], 0.0000, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[5], 0.0087, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[6], 0.2304, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[7], 0.2304, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[8], 0.5670, 1e-4);

            // Tidy up
            delete p_ode_system;
        }
    }
};

#endif /*TESTWNTCELLCYCLEODESYSTEM_HPP_*/
