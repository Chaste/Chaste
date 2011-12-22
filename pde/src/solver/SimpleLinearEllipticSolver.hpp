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

#ifndef SIMPLELINEARELLIPTICSOLVER_HPP_
#define SIMPLELINEARELLIPTICSOLVER_HPP_

#include "AbstractAssemblerSolverHybrid.hpp"
#include "AbstractStaticLinearPdeSolver.hpp"
#include "AbstractLinearEllipticPde.hpp"

/**
 * SimpleLinearEllipticSolver.
 *
 * Solver for solving AbstractLinearEllipticPdes.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SimpleLinearEllipticSolver
    : public AbstractAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, 1, NORMAL>,
      public AbstractStaticLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, 1>
{
protected:

    /** The PDE to be solved. */
    AbstractLinearEllipticPde<ELEMENT_DIM,SPACE_DIM>* mpEllipticPde;

    /**
     * The term to be added to the element stiffness matrix - see AbstractFeVolumeIntegralAssembler
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
    virtual c_matrix<double, 1*(ELEMENT_DIM+1), 1*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);

    /**
     * The term arising from boundary conditions to be added to the element
     * stiffness vector - see AbstractFeVolumeIntegralAssembler
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);


    // Note: does not have to provide a ComputeVectorSurfaceTerm for surface integrals,
    // the parent AbstractAssemblerSolverHybrid assumes natural Neumann BCs and uses a
    // NaturalNeumannSurfaceTermAssembler for assembling this part of the vector.



    /**
     * Delegate to AbstractAssemblerSolverHybrid::SetupGivenLinearSystem.
     *
     * @param currentSolution The current solution which can be used in setting up
     *   the linear system if needed (NULL if there isn't a current solution)
     * @param computeMatrix Whether to compute the LHS matrix of the linear system
     *   (mainly for dynamic solves).
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {
        SetupGivenLinearSystem(currentSolution, computeMatrix, this->mpLinearSystem);
    }

public:

    /**
     * Constructor.
     *
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE
     * @param pBoundaryConditions pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    SimpleLinearEllipticSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                               AbstractLinearEllipticPde<ELEMENT_DIM,SPACE_DIM>* pPde,
                               BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                               unsigned numQuadPoints = 2);

    /**
     * Overloaded InitaliseForSolve() which just calls the base class but also
     * sets the matrix as symmetric and sets Conjugate Gradients as the solver
     *
     * @param initialSolution initialSolution (used in base class version of this method)
     */
    void InitialiseForSolve(Vec initialSolution = NULL);
};

#endif /*SIMPLELINEARELLIPTICSOLVER_HPP_*/
