
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
     */
    AbstractAssemblerSolverHybrid(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                  BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* pBoundaryConditions)
        : AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, true, true, INTERPOLATION_LEVEL>(pMesh),
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
    assert(pLinearSystem->rGetLhsMatrix() != nullptr);
    assert(pLinearSystem->rGetRhsVector() != nullptr);

    // Assemble the matrix and vector calling methods on AbstractFeVolumeIntegralAssembler
    this->SetMatrixToAssemble(pLinearSystem->rGetLhsMatrix());
    this->SetVectorToAssemble(pLinearSystem->rGetRhsVector(), true);

    if (currentSolution != nullptr)
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

//// #2033 - see Test2dHeatEquationWithPeriodicBcs in TestSimpleLinearEllipticSolver.hpp
    //mpBoundaryConditions->ApplyPeriodicBcsToLinearProblem(*pLinearSystem, true);

    pLinearSystem->FinaliseRhsVector();
    pLinearSystem->FinaliseLhsMatrix();
}

#endif /*ABSTRACTASSEMBLERSOLVERHYBRID_HPP_*/
