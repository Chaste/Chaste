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

#ifndef TESTALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_
#define TESTALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <stdio.h>
#include <vector>
#include <iostream>

#include "OutputFileHandler.hpp"

#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "Timer.hpp"

#include "FakePetscSetup.hpp"

/**
 * This class contains tests for Alarcon2004OxygenBasedCellCycleOdeSystem,
 * a system of ODEs that are used by the cell-cycle model
 * Alarcon2004OxygenBasedCellCycleModel to determine when a cell is ready
 * to divide.
 */
class TestAlarcon2004OxygenBasedCellCycleOdeSystem : public CxxTest::TestSuite
{
public:

    /**
     * Test derivative calculations (correct values calculated using Matlab).
     */
    void TestAlarcon2004Equations()
    {
        // Set up
        double time = 0.0;
        double oxygen_concentration = 1.0;

        Alarcon2004OxygenBasedCellCycleOdeSystem normal_system(oxygen_concentration, false);
        Alarcon2004OxygenBasedCellCycleOdeSystem cancer_system(oxygen_concentration, true);

        std::vector<double> initial_conditions = normal_system.GetInitialConditions();

        std::vector<double> normal_derivs(initial_conditions.size());
        std::vector<double> cancer_derivs(initial_conditions.size());
        normal_system.EvaluateYDerivatives(time, initial_conditions, normal_derivs);
        cancer_system.EvaluateYDerivatives(time, initial_conditions, cancer_derivs);

        /**
         * Test derivatives are correct initially
         * (correct values calculated using Matlab code)
         */

        // Normal cell
        TS_ASSERT_DELTA(normal_derivs[0], 455.630699088, 1e-5);
        TS_ASSERT_DELTA(normal_derivs[1], 1.83000000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs[2], 3.00000000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs[3], 1.50000000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs[4], -5.4060000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs[5], 0.00000000000, 1e-5);

        // Cancer cell
        TS_ASSERT_DELTA(cancer_derivs[0], 455.630699088, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs[1], 1.83600000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs[2], 0.42000000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs[3], 1.50000000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs[4], -5.4060000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs[5], 0.00000000000, 1e-5);

        /**
         * Again test derivatives are correct initially, but for
         * different initial conditions (corresponding to a low
         * oxygen concentration). The usual initial condition for
         * z is zero, so we need to change it to see any difference.
         */
        oxygen_concentration = 0.1;

        Alarcon2004OxygenBasedCellCycleOdeSystem normal_system2(oxygen_concentration, false);
        Alarcon2004OxygenBasedCellCycleOdeSystem cancer_system2(oxygen_concentration, true);

        std::vector<double> normal_derivs2(initial_conditions.size());
        normal_system2.SetDefaultInitialCondition(2, 0.1);

        std::vector<double> cancer_derivs2(initial_conditions.size());
        cancer_system2.SetDefaultInitialCondition(2, 0.1);

        std::vector<double> initial_conditions2 = normal_system2.GetInitialConditions();

        normal_system2.EvaluateYDerivatives(time, initial_conditions2, normal_derivs2);
        cancer_system2.EvaluateYDerivatives(time, initial_conditions2, cancer_derivs2);

        // Normal cell
        TS_ASSERT_DELTA(normal_derivs2[0], 455.630699088, 1e-5);
        TS_ASSERT_DELTA(normal_derivs2[1], 1.81500000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs2[2], 2.94545454545, 1e-5);
        TS_ASSERT_DELTA(normal_derivs2[3], 1.50000000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs2[4], -5.4060000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs2[5], 0.00000000000, 1e-5);

        // Cancer cell
        TS_ASSERT_DELTA(cancer_derivs2[0], 455.630699088, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs2[1], 1.82100000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs2[2], 0.36545454545, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs2[3], 1.50000000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs2[4], -5.4060000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs2[5], 0.00000000000, 1e-5);
    }

    /**
     * Test two ODE solvers with this ODE system (correct values calculated using the Matlab solver ode15s).
     *
     */
    void TestAlarcon2004Solver()
    {
        // Set up
        double oxygen_concentration = 1.0;
        Alarcon2004OxygenBasedCellCycleOdeSystem alarcon_system(oxygen_concentration, false);

        // Create ODE solvers
        RungeKutta4IvpOdeSolver rk4_solver;
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        BackwardEulerIvpOdeSolver back_solver(6);

        // Set up for solver
        OdeSolution solutions;
        std::vector<double> initial_conditions = alarcon_system.GetInitialConditions();
        double h_value = 1e-4; // maximum tolerance

        // Solve the ODE system using a Runge Kutta fourth order solver
        Timer::Reset();
        solutions = rk4_solver.Solve(&alarcon_system, initial_conditions, 0.0, 10.0, h_value, h_value);
        Timer::Print("1. Runge-Kutta");

        // Reset maximum tolerance for Runge Kutta Fehlber solver
        h_value = 1e-1;

        // Solve the ODE system using a Runge Kutta Fehlber solver
        initial_conditions = alarcon_system.GetInitialConditions();
        Timer::Reset();
        solutions = rkf_solver.Solve(&alarcon_system, initial_conditions, 0.0, 10.0, h_value, 1e-4);
        Timer::Print("1. Runge-Kutta-Fehlberg");

        // Test that solutions are accurate for a small time increase
        int end = solutions.rGetSolutions().size() - 1;

        // Test that the solver stops at the right time
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 9.286356375, 1e-2);

        // Test solution - note the high tolerances
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 0.004000000000000, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 0.379221366479055, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 0.190488726735972, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3], 9.962110289977730, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4], 0.096476600742599, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5], 1.000000000000000, 1e-3);
    }

    void TestArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "alarcon_ode.arch";

        {
            double oxygen_concentration = 0.7;
            bool is_labelled = true;

            Alarcon2004OxygenBasedCellCycleOdeSystem ode_system(oxygen_concentration, is_labelled);

            TS_ASSERT_DELTA(ode_system.GetOxygenConcentration(), 0.70, 1e-6);
            TS_ASSERT_EQUALS(ode_system.IsLabelled(), true);


            ode_system.SetDefaultInitialCondition(2, 3.25);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 6u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.90, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[1], 0.01, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[2], 3.25, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[3], 5.00, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[4], 1.00, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[5], 0.70, 1e-6);

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
            TS_ASSERT_DELTA(static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(p_ode_system)->GetOxygenConcentration(), 0.70, 1e-6);
            TS_ASSERT_EQUALS(static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(p_ode_system)->IsLabelled(), true);

            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 6u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.90, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[1], 0.01, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[2], 0.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[3], 5.00, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[4], 1.00, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[5], 0.70, 1e-6);

            // Tidy up
            delete p_ode_system;
        }
    }
};

#endif /*TESTALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_*/
