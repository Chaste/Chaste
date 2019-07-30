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

#ifndef TESTTYSONNOVAK2001ODESYSTEM_HPP_
#define TESTTYSONNOVAK2001ODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <vector>
#include <iostream>

#include "TysonNovak2001OdeSystem.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "ColumnDataWriter.hpp"
#include "Timer.hpp"

#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestTysonNovak2001OdeSystem : public CxxTest::TestSuite
{
public:

    void TestTysonNovakEquation()
    {
        TysonNovak2001OdeSystem tyson_novak_system;

        double time = 0.0;
        std::vector<double> initial_conditions;
        initial_conditions.push_back(0.6);
        initial_conditions.push_back(0.1);
        initial_conditions.push_back(1.5);
        initial_conditions.push_back(0.6);
        initial_conditions.push_back(0.6);
        initial_conditions.push_back(0.85);

        std::vector<double> derivs(initial_conditions.size());
        tyson_novak_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct
        // Divided by 60 to change to hours
        TS_ASSERT_DELTA(derivs[0], -4.400000000000000e-02*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[1], -6.047872340425530e+00*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[2], 3.361442884485838e-02*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[3], 4.016602000735009e-02*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[4], 8.400000000000001e-03*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[5], 7.777500000000001e-03*60.0, 1e-5);
    }

    void TestTysonNovakSolver()
    {
        TysonNovak2001OdeSystem tyson_novak_system;

        // Solve system using backward Euler solver

        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits

        double dt = 0.1/60.0;

        //Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(6);

        std::vector<double> state_variables = tyson_novak_system.GetInitialConditions();

        Timer::Reset();
        OdeSolution solutions = backward_euler_solver.Solve(&tyson_novak_system, state_variables, 0.0, 75.8350/60.0, dt, dt);
        Timer::Print("1. Tyson Novak Backward Euler");

        // If you run it up to about 75min the ODE will stop, anything less and it will not and this test will fail
        TS_ASSERT_EQUALS(backward_euler_solver.StoppingEventOccurred(), true);

        unsigned end = solutions.rGetSolutions().size() - 1;

        // The following code provides nice output for gnuplot
        // use the command
        // plot "tyson_novak.dat" u 1:2
        // or
        // plot "tyson_novak.dat" u 1:3 etc. for the various proteins...

//        OutputFileHandler handler("");
//        out_stream file=handler.OpenOutputFile("tyson_novak.dat");
//        for (unsigned i=0; i<=end; i++)
//        {
//            (*file) << solutions.rGetTimes()[i]<< "\t" << solutions.rGetSolutions()[i][0] << "\t" << solutions.rGetSolutions()[i][1] << "\t" << solutions.rGetSolutions()[i][2] << "\t" << solutions.rGetSolutions()[i][3] << "\t" << solutions.rGetSolutions()[i][4] << "\t" << solutions.rGetSolutions()[i][5] << "\n" << std::flush;
//        }
//        file->close();

        ColumnDataWriter writer("TysonNovak", "TysonNovak");
        if (PetscTools::AmMaster()) // if master process
        {
            int step_per_row = 1;
            int time_var_id = writer.DefineUnlimitedDimension("Time", "s");

            std::vector<int> var_ids;
            for (unsigned i=0; i<tyson_novak_system.rGetStateVariableNames().size(); i++)
            {
                var_ids.push_back(writer.DefineVariable(tyson_novak_system.rGetStateVariableNames()[i],
                                                        tyson_novak_system.rGetStateVariableUnits()[i]));
            }
            writer.EndDefineMode();

            for (unsigned i = 0; i < solutions.rGetSolutions().size(); i+=step_per_row)
            {
                writer.PutVariable(time_var_id, solutions.rGetTimes()[i]);
                for (unsigned j=0; j<var_ids.size(); j++)
                {
                    writer.PutVariable(var_ids[j], solutions.rGetSolutions()[i][j]);
                }
                writer.AdvanceAlongUnlimitedDimension();
            }
            writer.Close();
        }
        PetscTools::Barrier();

        // Proper values calculated using the Matlab stiff ODE solver ode15s. Note that
        // large tolerances are required for the tests to pass with both chaste solvers
        // and CVODE.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],1.54216806705641, 1e-1);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],1.40562614481544, 1e-1);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],0.95328206604519, 2e-2);
    }

    void TestArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "tn_ode.arch";

        {
            TysonNovak2001OdeSystem ode_system;

            ode_system.SetDefaultInitialCondition(2, 3.25);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 6u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.0999, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[1], 0.9890, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[2], 3.2500, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[3], 1.4211, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[4], 0.6728, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[5], 0.4854, 1e-4);

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
            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 6u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.0999, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[1], 0.9890, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[2], 3.2500, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[3], 1.4211, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[4], 0.6728, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[5], 0.4854, 1e-4);

            // Tidy up
            delete p_ode_system;
        }
    }
};

#endif /*TESTTYSONNOVAK2001ODESYSTEM_HPP_*/
