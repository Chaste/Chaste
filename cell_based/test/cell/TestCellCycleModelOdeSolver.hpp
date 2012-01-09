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

#ifndef TESTCELLCYCLEMODELODESOLVER_HPP_
#define TESTCELLCYCLEMODELODESOLVER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "CellCycleModelOdeSolver.hpp"

#include "TysonNovakCellCycleModel.hpp"
#include "Alarcon2004OxygenBasedCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "OdeSystemInformation.hpp"
#include "AbstractCellBasedTestSuite.hpp"

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

    void TestMethods() throw(Exception)
    {
        typedef CellCycleModelOdeSolver<TysonNovakCellCycleModel,RungeKutta4IvpOdeSolver> RkSolver;

        // Check we can create an instance
        boost::shared_ptr<RkSolver> p_solver = RkSolver::Instance();
        TS_ASSERT(p_solver.get() != NULL);

        // Check singleton-ness
        boost::shared_ptr<RkSolver> p_solver2 = RkSolver::Instance();
        TS_ASSERT_EQUALS(p_solver, p_solver2);

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

    void TestWithBackwardEulerIvpOdeSolver() throw(Exception)
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

    void TestWithCvodeAdaptor() throw(Exception)
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
#endif // CHASTE_CVODE
    }

    void TestArchiving() throw(Exception)
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
