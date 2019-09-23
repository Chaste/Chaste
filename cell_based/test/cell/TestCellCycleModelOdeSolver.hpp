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

#ifndef TESTCELLCYCLEMODELODESOLVER_HPP_
#define TESTCELLCYCLEMODELODESOLVER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <cmath>

#include "CellCycleModelOdeSolver.hpp"

#include "TysonNovakCellCycleModel.hpp"
#include "Alarcon2004OxygenBasedCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "OdeSystemInformation.hpp"
#include "AbstractCellBasedTestSuite.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * Simple ODE system for use in the test suite. Defines the
 * IVP dy/dt = 1, y(0) = 0.
 */
class SimpleOde : public AbstractOdeSystem
{
public:
    SimpleOde() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde>::Instance();
        ResetToInitialConditions();
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] = 1.0;
    }
};

template<>
void OdeSystemInformation<SimpleOde>::Initialise()
{
    this->mVariableNames.push_back("Variable_1");
    this->mVariableUnits.push_back("Units_1");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}

/**
 * Simple ODE system with stopping event for use in the test suite.
 * Defines the IVP dx/dt = y, dy/dt = -x, x(0) = 1, y(0) = 0.
 * Solutions to this system form circles about the origin. The
 * stopping event is x(0) < 0, which should first occur at time
 * t = pi/2.
 */
class OdeSecondOrderWithEvents : public AbstractOdeSystem
{
public:
    OdeSecondOrderWithEvents() : AbstractOdeSystem(2)
    {
        mpSystemInfo = OdeSystemInformation<OdeSecondOrderWithEvents>::Instance();
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] =  rY[1];
        rDY[1] = -rY[0];
    }

    bool CalculateStoppingEvent(double time, const std::vector<double>& rY)
    {
        return (rY[0] < 0);
    }
};

template<>
void OdeSystemInformation<OdeSecondOrderWithEvents>::Initialise()
{
    this->mVariableNames.push_back("Variable_1");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("Variable_2");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}

class TestCellCycleModelOdeSolver : public AbstractCellBasedTestSuite
{
public:

    void TestMethods()
    {
        typedef CellCycleModelOdeSolver<TysonNovakCellCycleModel,RungeKutta4IvpOdeSolver> RkSolver;

        // Check we can create an instance
        boost::shared_ptr<RkSolver> p_solver = RkSolver::Instance();
        TS_ASSERT(p_solver.get() != NULL);

        // Check singleton-ness
        boost::shared_ptr<RkSolver> p_solver2 = RkSolver::Instance();
        TS_ASSERT_EQUALS(p_solver, p_solver2);

        // Coverage
        TS_ASSERT_THROWS_NOTHING(p_solver->Reset());

        p_solver->Initialise();

        // Check the solver can be called for a simple ODE system
        SimpleOde ode;
        double last_time = 0.0;
        double current_time = 2;
        double dt = 1e-5;

        ode.SetStateVariables(ode.GetInitialConditions());
        p_solver->SolveAndUpdateStateVariable(&ode, last_time, current_time, dt);

        // No stopping event is specified, so check the solver did not stop
        TS_ASSERT_EQUALS(p_solver->StoppingEventOccurred(), false);

        // Check the solver can be called for another ODE system, this time with a stopping event
        OdeSecondOrderWithEvents ode_with_events;

        ode_with_events.SetStateVariables(ode_with_events.GetInitialConditions());
        p_solver->SolveAndUpdateStateVariable(&ode_with_events, last_time, current_time, dt);

        // Check the solver stopped at the correct time
        TS_ASSERT_EQUALS(p_solver->StoppingEventOccurred(), true);
        TS_ASSERT_DELTA(p_solver->GetStoppingTime(), M_PI_2, 1e-4);
    }

    void TestWithBackwardEulerIvpOdeSolver()
    {
        typedef CellCycleModelOdeSolver<TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver> EulerSolver;

        // Check we can create an instance
        boost::shared_ptr<EulerSolver> p_solver = EulerSolver::Instance();
        TS_ASSERT(p_solver.get() != NULL);

        // Check singleton-ness
        boost::shared_ptr<EulerSolver> p_solver2 = EulerSolver::Instance();
        TS_ASSERT_EQUALS(p_solver, p_solver2);

        TS_ASSERT_THROWS_THIS(p_solver->Initialise(), "SetSizeOfOdeSystem() must be called before calling Initialise()");
        p_solver->SetSizeOfOdeSystem(1);
        p_solver->Initialise();

        // Check the solver can be called for a simple ODE system
        SimpleOde ode;
        double last_time = 0.0;
        double current_time = 2;
        double dt = 1e-4;

        ode.SetStateVariables(ode.GetInitialConditions());
        p_solver->SolveAndUpdateStateVariable(&ode, last_time, current_time, dt);

        // No stopping event is specified, so check the solver did not stop
        TS_ASSERT_EQUALS(p_solver->StoppingEventOccurred(), false);

        p_solver->Reset();
        p_solver->SetSizeOfOdeSystem(2);
        p_solver->Initialise();

        TS_ASSERT_EQUALS(p_solver->GetSizeOfOdeSystem(), 2u);

        // Check the solver can be called for another ODE system, this time with a stopping event
        OdeSecondOrderWithEvents ode_with_events;

        ode_with_events.SetStateVariables(ode_with_events.GetInitialConditions());
        p_solver->SolveAndUpdateStateVariable(&ode_with_events, last_time, current_time, dt);

        // Check the solver stopped at the correct time
        TS_ASSERT_EQUALS(p_solver->StoppingEventOccurred(), true);

        TS_ASSERT_DELTA(p_solver->GetStoppingTime(), M_PI_2, 1e-2);
    }

    void TestWithCvodeAdaptor()
    {
#ifdef CHASTE_CVODE
        typedef CellCycleModelOdeSolver<TysonNovakCellCycleModel, CvodeAdaptor> CvodeSolver;

        // Check we can create an instance
        boost::shared_ptr<CvodeSolver> p_solver = CvodeSolver::Instance();
        TS_ASSERT(p_solver.get() != NULL);

        // Check singleton-ness
        boost::shared_ptr<CvodeSolver> p_solver2 = CvodeSolver::Instance();
        TS_ASSERT_EQUALS(p_solver, p_solver2);

        p_solver->Initialise();

        TS_ASSERT_THROWS_NOTHING(p_solver->CheckForStoppingEvents());
        TS_ASSERT_THROWS_NOTHING(p_solver->SetMaxSteps(1000));
        TS_ASSERT_THROWS_NOTHING(p_solver->SetTolerances(1e-5, 1e-5));
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ode_solver.arch";

        typedef CellCycleModelOdeSolver<Alarcon2004OxygenBasedCellCycleModel, BackwardEulerIvpOdeSolver> EulerSolver;

        // Create an output archive
        {
            boost::shared_ptr<EulerSolver> p_solver = EulerSolver::Instance();
            p_solver->SetSizeOfOdeSystem(4);
            p_solver->Initialise();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_solver;

            TS_ASSERT_EQUALS(p_solver->GetSizeOfOdeSystem(), 4u);
            p_solver->Reset();
            TS_ASSERT_EQUALS(p_solver->GetSizeOfOdeSystem(), UNSIGNED_UNSET);
        }

        {
            boost::shared_ptr<EulerSolver> p_solver = EulerSolver::Instance();

            TS_ASSERT_EQUALS(p_solver->GetSizeOfOdeSystem(), UNSIGNED_UNSET);

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_solver;

            TS_ASSERT_EQUALS(p_solver->GetSizeOfOdeSystem(), 4u);
        }
    }
};

#endif /*TESTCELLCYCLEMODELODESOLVER_HPP_*/
