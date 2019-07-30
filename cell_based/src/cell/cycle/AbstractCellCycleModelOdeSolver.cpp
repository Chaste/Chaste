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

#include "AbstractCellCycleModelOdeSolver.hpp"
#include "CvodeAdaptor.hpp"

AbstractCellCycleModelOdeSolver::AbstractCellCycleModelOdeSolver()
    : mSizeOfOdeSystem(UNSIGNED_UNSET)
{
}

AbstractCellCycleModelOdeSolver::~AbstractCellCycleModelOdeSolver()
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
