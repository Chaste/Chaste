/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef LINEARELLIPTICPDESYSTEMSOLVER_HPP_
#define LINEARELLIPTICPDESYSTEMSOLVER_HPP_

#include "AbstractAssemblerSolverHybrid.hpp"
#include "AbstractStaticLinearPdeSolver.hpp"
#include "AbstractLinearEllipticPdeSystem.hpp"

///\todo #2930
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class LinearEllipticPdeSystemSolver
    : public AbstractAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NORMAL>,
      public AbstractStaticLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
protected:

    /** The PDE system to be solved. */
    AbstractLinearEllipticPdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpEllipticPdeSystem;

    /**
     * ///\todo #2930
     * @return the term to be added to the element stiffness matrix - see AbstractFeVolumeIntegralAssembler
     *
     * grad_phi[row] . ( pde_diffusion_term * grad_phi[col])
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double, PROBLEM_DIM>& rU,
        c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM, SPACE_DIM>* pElement);

    /**
     * ///\todo #2930
     * @return the term arising from boundary conditions to be added to the element
     * stiffness vector - see AbstractFeVolumeIntegralAssembler
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double, PROBLEM_DIM>& rU,
        c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM, SPACE_DIM>* pElement);


    // Note: does not have to provide a ComputeVectorSurfaceTerm for surface integrals,
    // the parent AbstractAssemblerSolverHybrid assumes natural Neumann BCs and uses a
    // NaturalNeumannSurfaceTermAssembler for assembling this part of the vector.

    /**
     * ///\todo #2930
     * Delegate to AbstractAssemblerSolverHybrid::SetupGivenLinearSystem.
     *
     * @param currentSolution The current solution which can be used in setting up
     *   the linear system if needed (NULL if there isn't a current solution)
     * @param computeMatrix Whether to compute the LHS matrix of the linear system
     *   (mainly for dynamic solves).
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {
        this->SetupGivenLinearSystem(currentSolution, computeMatrix, this->mpLinearSystem);
    }

public:

    /**
     * Constructor.
     *
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE system
     * @param pBoundaryConditions pointer to the boundary conditions
     */
    LinearEllipticPdeSystemSolver(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                  AbstractLinearEllipticPdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pPdeSystem,
                                  BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions);

    /**
     * Overloaded InitaliseForSolve() which just calls the base class but also
     * sets the matrix as symmetric and sets Conjugate Gradients as the solver
     *
     * @param initialSolution initialSolution (used in base class version of this method)
     */
    void InitialiseForSolve(Vec initialSolution = nullptr);
};

/*
 * As this class is templated over PROBLEM_DIM, we put the implementation
 * in the header file to avoid explicit instantiation.
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
LinearEllipticPdeSystemSolver<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::LinearEllipticPdeSystemSolver(
    AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
    AbstractLinearEllipticPdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pPdeSystem,
    BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions)
    : AbstractAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NORMAL>(pMesh, pBoundaryConditions),
      AbstractStaticLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(pMesh)
{
    mpEllipticPdeSystem = pPdeSystem;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)>LinearEllipticPdeSystemSolver<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ComputeMatrixTerm(
    c_vector<double, ELEMENT_DIM+1>& rPhi,
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
    ChastePoint<SPACE_DIM>& rX,
    c_vector<double, PROBLEM_DIM>& rU,
    c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
    Element<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> matrix_term = zero_matrix<double>(PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1));

    // Loop over PDEs and populate matrix_term
    for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
        c_matrix<double, 1*(ELEMENT_DIM+1), 1*(ELEMENT_DIM+1)> this_stiffness_matrix = zero_matrix<double>(1*(ELEMENT_DIM+1), 1*(ELEMENT_DIM+1));
        c_matrix<double, SPACE_DIM, SPACE_DIM> this_pde_diffusion_term = mpEllipticPdeSystem->ComputeDiffusionTerm(rX, pde_index);

        // This if statement just saves computing phi*phi^T if it is to be multiplied by zero
        double this_source_term = mpEllipticPdeSystem->ComputeLinearInUCoeffInSourceTerm(rX, pde_index, pElement);
        if (this_source_term != 0)
        {
            this_stiffness_matrix = prod(trans(rGradPhi), c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>(prod(this_pde_diffusion_term, rGradPhi))) - this_source_term*outer_prod(rPhi,rPhi);
        }
        else
        {
            this_stiffness_matrix = prod(trans(rGradPhi), c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>(prod(this_pde_diffusion_term, rGradPhi)));
        }

        for (unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                matrix_term(i*PROBLEM_DIM + pde_index, j*PROBLEM_DIM + pde_index) = this_stiffness_matrix(i,j);
            }
        }
    }
    return matrix_term;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> LinearEllipticPdeSystemSolver<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ComputeVectorTerm(
    c_vector<double, ELEMENT_DIM+1>& rPhi,
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
    ChastePoint<SPACE_DIM>& rX,
    c_vector<double, PROBLEM_DIM>& rU,
    c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
    Element<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> vector_term = zero_vector<double>(PROBLEM_DIM*(ELEMENT_DIM+1));

    // Loop over PDEs and populate vector_term
    for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
        double this_source_term = mpEllipticPdeSystem->ComputeConstantInUSourceTerm(rX, pde_index, pElement);
        c_vector<double, ELEMENT_DIM+1> this_vector_term = this_source_term * rPhi;

        for (unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            vector_term(i*PROBLEM_DIM + pde_index) = this_vector_term(i);
        }
    }

    return vector_term;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void LinearEllipticPdeSystemSolver<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::InitialiseForSolve(Vec initialSolution)
{
    AbstractLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::InitialiseForSolve(initialSolution);
    assert(this->mpLinearSystem);
    this->mpLinearSystem->SetMatrixIsSymmetric(true);
    this->mpLinearSystem->SetKspType("cg");
}

#endif /*LINEARELLIPTICPDESYSTEMSOLVER_HPP_*/
