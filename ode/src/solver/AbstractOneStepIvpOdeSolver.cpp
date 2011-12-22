/*

Copyright (C) University of Oxford, 2005-2011

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

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "TimeStepper.hpp"
#include "Exception.hpp"
#include <cmath>

OdeSolution AbstractOneStepIvpOdeSolver::Solve(AbstractOdeSystem* pOdeSystem,
                                               std::vector<double>& rYValues,
                                               double startTime,
                                               double endTime,
                                               double timeStep,
                                               double timeSampling)
{
    assert(rYValues.size()==pOdeSystem->GetNumberOfStateVariables());
    assert(endTime > startTime);
    assert(timeStep > 0.0);
    assert(timeSampling >= timeStep);

    mStoppingEventOccurred = false;
    if ( pOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true )
    {
        EXCEPTION("(Solve with sampling) Stopping event is true for initial condition");
    }
    TimeStepper stepper(startTime, endTime, timeSampling);

    // setup solutions if output is required
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(stepper.EstimateTimeSteps());
    solutions.rGetSolutions().push_back(rYValues);
    solutions.rGetTimes().push_back(startTime);
    solutions.SetOdeSystemInformation(pOdeSystem->GetSystemInformation());
    solutions.SetSolverName( GetIdentifier() );

    mWorkingMemory.resize(rYValues.size());

    // Solve the ODE system
    while ( !stepper.IsTimeAtEnd() && !mStoppingEventOccurred )
    {
        InternalSolve(pOdeSystem, rYValues, mWorkingMemory, stepper.GetTime(), stepper.GetNextTime(), timeStep);
        stepper.AdvanceOneTimeStep();
        // write current solution into solutions
        solutions.rGetSolutions().push_back(rYValues);
        // Push back new time into the time solution vector
        if ( mStoppingEventOccurred )
        {
            solutions.rGetTimes().push_back(mStoppingTime);
        }
        else
        {
            solutions.rGetTimes().push_back(stepper.GetTime());
        }
    }

    // stepper.EstimateTimeSteps may have been an overestimate...
    solutions.SetNumberOfTimeSteps(stepper.GetTotalTimeStepsTaken());
    return solutions;
}

void AbstractOneStepIvpOdeSolver::Solve(AbstractOdeSystem* pOdeSystem,
                                        std::vector<double>& rYValues,
                                        double startTime,
                                        double endTime,
                                        double timeStep)
{
    assert(rYValues.size()==pOdeSystem->GetNumberOfStateVariables());
    assert(endTime > startTime);
    assert(timeStep > 0.0);

    mStoppingEventOccurred = false;
    if ( pOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true )
    {
        EXCEPTION("(Solve without sampling) Stopping event is true for initial condition");
    }

    // Perhaps resize working memory
    mWorkingMemory.resize(rYValues.size());
    // And solve...
    InternalSolve(pOdeSystem, rYValues, mWorkingMemory, startTime, endTime, timeStep);
}

void AbstractOneStepIvpOdeSolver::InternalSolve(AbstractOdeSystem* pOdeSystem,
                                                std::vector<double>& rYValues,
                                                std::vector<double>& rWorkingMemory,
                                                double startTime,
                                                double endTime,
                                                double timeStep)
{
    TimeStepper stepper(startTime, endTime, timeStep);
    // Solve the ODE system

    // Which of our vectors holds the current solution?
    // If this is true, it's in rYValues, otherwise it's in rWorkingMemory.
    bool curr_is_curr = false;

    // should never get here if this bool has been set to true;
    assert(!mStoppingEventOccurred);
    while ( !stepper.IsTimeAtEnd() && !mStoppingEventOccurred )
    {
        curr_is_curr = not curr_is_curr;
        // Function that calls the appropriate one-step solver
        CalculateNextYValue(pOdeSystem,
                            stepper.GetNextTimeStep(),
                            stepper.GetTime(),
                            curr_is_curr ? rYValues : rWorkingMemory,
                            curr_is_curr ? rWorkingMemory : rYValues);
        stepper.AdvanceOneTimeStep();
        if ( pOdeSystem->CalculateStoppingEvent(stepper.GetTime(),
                                                curr_is_curr ? rWorkingMemory : rYValues) == true )
        {
            mStoppingTime = stepper.GetTime();
            mStoppingEventOccurred = true;
        }
    }
    // Final answer must be in rYValues
    if (curr_is_curr)
    {
        rYValues.assign(rWorkingMemory.begin(), rWorkingMemory.end());
    }
}
