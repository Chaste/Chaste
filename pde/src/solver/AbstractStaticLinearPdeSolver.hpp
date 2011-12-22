
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
     */
    Vec Solve(Vec initialGuess=NULL);
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractStaticLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::Solve(Vec initialGuess)
{
    // Set up the linear system
    this->InitialiseForSolve(initialGuess);

    // This method can be overloaded if neccessary
    this->PrepareForSetupLinearSystem(NULL);

    // This method should be implemented by the concrete class
    this->SetupLinearSystem(NULL, true);

    this->FinaliseLinearSystem(NULL);

    // Solve the linear system
    Vec solution = this->mpLinearSystem->Solve(initialGuess);

    this->FollowingSolveLinearSystem(solution);

    return solution;
}

#endif /*ABSTRACTSTATICLINEARPDESOLVER_HPP_*/
