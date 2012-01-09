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

#include "AbstractCellCycleModelOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "Exception.hpp"

AbstractCellCycleModelOdeSolver::AbstractCellCycleModelOdeSolver()
    : mSizeOfOdeSystem(UNSIGNED_UNSET)
{
}

AbstractCellCycleModelOdeSolver::~AbstractCellCycleModelOdeSolver()
{
}

void AbstractCellCycleModelOdeSolver::Reset()
{
}

void AbstractCellCycleModelOdeSolver::SolveAndUpdateStateVariable(AbstractOdeSystem* pAbstractOdeSystem,
                                                                  double startTime,
                                                                  double endTime,
                                                                  double timeStep)
{
    assert(IsSetUp());
    mpOdeSolver->SolveAndUpdateStateVariable(pAbstractOdeSystem, startTime, endTime, timeStep);
}

bool AbstractCellCycleModelOdeSolver::StoppingEventOccurred()
{
    assert(IsSetUp());
    return mpOdeSolver->StoppingEventOccurred();
}

double AbstractCellCycleModelOdeSolver::GetStoppingTime()
{
    assert(IsSetUp());
    return mpOdeSolver->GetStoppingTime();
}

void AbstractCellCycleModelOdeSolver::SetSizeOfOdeSystem(unsigned sizeOfOdeSystem)
{
    mSizeOfOdeSystem = sizeOfOdeSystem;
}

unsigned AbstractCellCycleModelOdeSolver::GetSizeOfOdeSystem()
{
    return mSizeOfOdeSystem;
}

void AbstractCellCycleModelOdeSolver::CheckForStoppingEvents()
{
#ifdef CHASTE_CVODE
    assert(IsSetUp());
    if (boost::dynamic_pointer_cast<CvodeAdaptor>(mpOdeSolver))
    {
        (boost::static_pointer_cast<CvodeAdaptor>(mpOdeSolver))->CheckForStoppingEvents();
    }
#endif //CHASTE_CVODE
}

void AbstractCellCycleModelOdeSolver::SetMaxSteps(long int numSteps)
{
#ifdef CHASTE_CVODE
    assert(IsSetUp());
    if (boost::dynamic_pointer_cast<CvodeAdaptor>(mpOdeSolver))
    {
        (boost::static_pointer_cast<CvodeAdaptor>(mpOdeSolver))->SetMaxSteps(numSteps);
    }
#endif //CHASTE_CVODE
}

void AbstractCellCycleModelOdeSolver::SetTolerances(double relTol, double absTol)
{
#ifdef CHASTE_CVODE
    assert(IsSetUp());
    if (boost::dynamic_pointer_cast<CvodeAdaptor>(mpOdeSolver))
    {
        (boost::static_pointer_cast<CvodeAdaptor>(mpOdeSolver))->SetTolerances(relTol, absTol);
    }
#endif //CHASTE_CVODE
}

bool AbstractCellCycleModelOdeSolver::IsAdaptive()
{
    bool adaptive = false;
#ifdef CHASTE_CVODE
    assert(IsSetUp());
    if (boost::dynamic_pointer_cast<CvodeAdaptor>(mpOdeSolver))
    {
        adaptive = true;
    }
#endif //CHASTE_CVODE
    return adaptive;
}
