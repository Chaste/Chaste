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

#include "OdeLinearSystemSolver.hpp"

OdeLinearSystemSolver::OdeLinearSystemSolver(unsigned systemSize, double timeStep)
    : mLinearSystem(systemSize)
{
    assert(timeStep > 0.0);
    mTimeStep = timeStep;

    // Initialise vectors to zero
    mCurrentSolution = PetscTools::CreateAndSetVec(systemSize, 0.0);
    mForceVector = PetscTools::CreateAndSetVec(systemSize, 0.0);
}

OdeLinearSystemSolver::~OdeLinearSystemSolver()
{
    PetscTools::Destroy(mCurrentSolution);
    PetscTools::Destroy(mForceVector);
}

double OdeLinearSystemSolver::GetTimeStep()
{
    return mTimeStep;
}

Mat& OdeLinearSystemSolver::rGetLhsMatrix()
{
    return mLinearSystem.rGetLhsMatrix();
}

Vec& OdeLinearSystemSolver::rGetForceVector()
{
    return mForceVector;
}

void OdeLinearSystemSolver::SetInitialConditionVector(Vec initialConditionsVector)
{
    VecCopy(initialConditionsVector, mCurrentSolution);
}

Vec OdeLinearSystemSolver::SolveOneTimeStep()
{
    // Compute the product of the LHS matrix and the current solution vector,
    // setting the answer to be the RHS vector
    MatMult(mLinearSystem.rGetLhsMatrix(), mCurrentSolution, mLinearSystem.rGetRhsVector());

    // Add timestep multipled by force vector
    PetscVecTools::AddScaledVector(mLinearSystem.rGetRhsVector(), mForceVector, mTimeStep);

    // avoid memory leaks
    PetscTools::Destroy(mCurrentSolution);

    // Having constructed the RHS vector, solve the resulting linear system...
    mCurrentSolution = mLinearSystem.Solve();

    return mCurrentSolution;
}
