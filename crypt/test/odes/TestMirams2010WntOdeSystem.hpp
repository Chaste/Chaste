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

#ifndef TESTMIRAMS2010WNTODESYSTEM_HPP_
#define TESTMIRAMS2010WNTODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include <fstream>

#include "OutputFileHandler.hpp"
#include "Mirams2010WntOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "Timer.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestMirams2010WntOdeSystem : public AbstractCellBasedTestSuite
{
public:

    void TestMirams2010WntOdeSystemSetup()
    {
#ifdef CHASTE_CVODE
        double wnt_level = 0.5;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        Mirams2010WntOdeSystem wnt_system(wnt_level, p_state);
        // Solve system using CVODE solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.
        double h_value = 0.1;

        CvodeAdaptor cvode_solver;

        OdeSolution solutions;
        //OdeSolution solutions2;

        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();
        std::cout << "Timings for 100 hours\n";
        Timer::Reset();
        solutions = cvode_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        Timer::Print("1. Cvode");

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 7.8 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 100, 1e-2);

        // Decent results
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 67.5011, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 67.5011, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], wnt_level, 1e-4);
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestGarysWntOdeSystemApc2Hit()
    {
#ifdef CHASTE_CVODE
        double wnt_level = 0.5;
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        Mirams2010WntOdeSystem wnt_system(wnt_level, p_apc2);

        // Solve system using CVODE solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits
        double h_value = 0.01;
        CvodeAdaptor cvode_solver;
        OdeSolution solutions;
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        Timer::Reset();
        solutions = cvode_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        Timer::Print("1. Cvode");

        // Test solutions are OK for a small time increase
        int end = solutions.rGetSolutions().size() - 1;

        // Test the simulation is ending at the right time (going into S phase at 7.8 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 100, 1e-2);

        // Check results are correct
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 433.114, 2e-3); // Tolerances relaxed for
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 433.114, 2e-3); // different CVODE versions.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], wnt_level, 1e-4);
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestGarysWntOdeSystemBetaCatenin1Hit()
    {
#ifdef CHASTE_CVODE
        double wnt_level = 0.5;
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
        Mirams2010WntOdeSystem wnt_system(wnt_level, p_bcat1);

        // Solve system using CVODE solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits
        double h_value = 0.1;
        CvodeAdaptor cvode_solver;
        OdeSolution solutions;
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        Timer::Reset();
        solutions = cvode_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        Timer::Print("1. Cvode");

        // Test solutions are OK for a small time increase
        int end = solutions.rGetSolutions().size() - 1;

        // Tests the simulation is ending at the right time (going into S phase at 7.8 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 100, 1e-2);

        // Check results are correct
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 67.5011, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 824.0259, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], wnt_level, 1e-4);
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestArchiving()
    {
        EXIT_IF_PARALLEL;
#ifdef CHASTE_CVODE
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "mirams_ode.arch";

        {
            double wnt_level = 0.5;
            boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
            Mirams2010WntOdeSystem ode_system(wnt_level, p_state);

            TS_ASSERT_DELTA(ode_system.GetWntLevel(), 0.50, 1e-6);
            TS_ASSERT_EQUALS(ode_system.GetMutationState()->IsType<WildTypeCellMutationState>(), true);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 3u);
            TS_ASSERT_DELTA(initial_conditions[0], 64.1863, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[1], 64.1863, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[2], 0.5, 1e-6);

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
            TS_ASSERT_DELTA(static_cast<Mirams2010WntOdeSystem*>(p_ode_system)->GetWntLevel(), 0.50, 1e-6);
            TS_ASSERT_EQUALS(static_cast<Mirams2010WntOdeSystem*>(p_ode_system)->GetMutationState()->IsType<WildTypeCellMutationState>(), true);

            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 3u);
            TS_ASSERT_DELTA(initial_conditions[0], 64.1863, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[1], 64.1863, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[2], 0.5, 1e-6);

            // Tidy up
            delete p_ode_system;
        }
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }
};

#endif /* TESTMIRAMS2010WNTODESYSTEM_HPP_ */
