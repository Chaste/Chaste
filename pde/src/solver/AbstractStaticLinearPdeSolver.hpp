
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

#ifndef ABSTRACTSTATICLINEARPDESOLVER_HPP_
#define ABSTRACTSTATICLINEARPDESOLVER_HPP_

#include "AbstractLinearPdeSolver.hpp"

/**
 * Abstract class for static linear PDE solves.
 * This class defines the Solve() method. The concrete class should implement
 * the SetupLinearSystem() method (defined in AbstractLinearPdeSolver), based
 * on the PDE being solved and the numerical method.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractStaticLinearPdeSolver : public AbstractLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
public:

    /**
     * Constructor.
     *
     * @param pMesh the mesh
     */
    AbstractStaticLinearPdeSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
        : AbstractLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(pMesh)
    {
    }

    /**
     * Static solve method.
     *
     * @param initialGuess optional initial guess for passing into the linear solve method
     * @return the solution vector
     */
    Vec Solve(Vec initialGuess=nullptr);
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractStaticLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::Solve(Vec initialGuess)
{
    // Set up the linear system
    this->InitialiseForSolve(initialGuess);

    // This method can be overloaded if necessary
    this->PrepareForSetupLinearSystem(nullptr);

    // This method should be implemented by the concrete class
    this->SetupLinearSystem(nullptr, true);

    this->FinaliseLinearSystem(nullptr);

    // Solve the linear system
    Vec solution = this->mpLinearSystem->Solve(initialGuess);

    this->FollowingSolveLinearSystem(solution);

    return solution;
}

#endif /*ABSTRACTSTATICLINEARPDESOLVER_HPP_*/
