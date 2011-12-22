
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

#ifndef ABSTRACTASSEMBLERSOLVERHYBRID_HPP_
#define ABSTRACTASSEMBLERSOLVERHYBRID_HPP_

#include "AbstractFeVolumeIntegralAssembler.hpp"
#include "AbstractLinearPdeSolver.hpp"
#include "NaturalNeumannSurfaceTermAssembler.hpp"

/**
 * A class which inherits from AbstractFeVolumeIntegralAssembler and implements a method SetupGivenLinearSystem(), which sets up
 * the given linear system using the assembler part of this class, which can be called by SetUpLinearSystem() on a
 * concrete solver.
 *
 * It assumes natural Neumann boundary conditions are needed and uses a NaturalNeumannSurfaceTermAssembler for this
 * part of the vector.
 *
 * See SimpleLinearEllipticSolver for an example of a concrete class
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, InterpolationLevel INTERPOLATION_LEVEL>
class AbstractAssemblerSolverHybrid
   : public AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, true, true, INTERPOLATION_LEVEL>
{
protected:

    /** An assembler for Neumann surface integrals, which are assumed to arise from natural Neumann boundary
     *  conditions, ie such that this surface integral is (for a 1-unknown problem) integral(g phi_i dS),
     *  where g is the Neumann boundary condition function
     */
    NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> mNaturalNeumannSurfaceTermAssembler;

    /** Boundary conditions container */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* mpBoundaryConditions;


public:

    /**
     * Constructor.
     *
     * @param pMesh pointer to the mesh
     * @param pBoundaryConditions pointer to the boundary conditions.
     * @param numQuadPoints number of quadrature points in each dimension to use per element (defaults to 2)
     */
    AbstractAssemblerSolverHybrid(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                  BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* pBoundaryConditions,
                                  unsigned numQuadPoints=2)
        : AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, true, true, INTERPOLATION_LEVEL>(pMesh,numQuadPoints),
          mNaturalNeumannSurfaceTermAssembler(pMesh,pBoundaryConditions),
          mpBoundaryConditions(pBoundaryConditions)
    {
        assert(pMesh);
        assert(pBoundaryConditions);
    }

    /**
     * Destructor.
     */
    virtual ~AbstractAssemblerSolverHybrid()
    {
    }

    /**
     * Implementation of AbstractLinearPdeSolver::SetupLinearSystem, using the assembler that this class
     * also inherits from. Concrete classes inheriting from both this class and
     * AbstractLinearPdeSolver can then have a one-line implementation of
     * AbstractLinearPdeSolver::SetupLinearSystem which calls this method.
     *
     * @param currentSolution The current solution which can be used in setting up
     *  the linear system if needed (NULL if there isn't a current solution)
     * @param computeMatrix Whether to compute the LHS matrix of the linear system
     *  (mainly for dynamic solves)
     * @param pLinearSystem  The linear system to set up.
     */
    void SetupGivenLinearSystem(Vec currentSolution, bool computeMatrix, LinearSystem* pLinearSystem);
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, INTERPOLATION_LEVEL>::SetupGivenLinearSystem(Vec currentSolution,
                                                                                                                     bool computeMatrix,
                                                                                                                     LinearSystem* pLinearSystem)
{
    assert(pLinearSystem->rGetLhsMatrix() != NULL);
    assert(pLinearSystem->rGetRhsVector() != NULL);

    // Assemble the matrix and vector calling methods on AbstractFeVolumeIntegralAssembler
    this->SetMatrixToAssemble(pLinearSystem->rGetLhsMatrix());
    this->SetVectorToAssemble(pLinearSystem->rGetRhsVector(), true);

    if (currentSolution != NULL)
    {
        this->SetCurrentSolution(currentSolution);
    }

    if (computeMatrix)
    {
        this->Assemble();
    }
    else
    {
        this->AssembleVector();
    }

    // Add the Neumann boundary conditions. The boundary conditions put into the BoundaryConditionsContainer
    // are assumed to be natural Neumann BCs.
    mNaturalNeumannSurfaceTermAssembler.SetVectorToAssemble(pLinearSystem->rGetRhsVector(), false);
    mNaturalNeumannSurfaceTermAssembler.Assemble();

    pLinearSystem->FinaliseRhsVector();
    pLinearSystem->SwitchWriteModeLhsMatrix();

    // add Dirichlet BCs
    mpBoundaryConditions->ApplyDirichletToLinearProblem(*pLinearSystem, true);

    pLinearSystem->FinaliseRhsVector();
    pLinearSystem->FinaliseLhsMatrix();
}

#endif /*ABSTRACTASSEMBLERSOLVERHYBRID_HPP_*/
