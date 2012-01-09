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

#ifndef CELLBASEDPDESOLVER_HPP_
#define CELLBASEDPDESOLVER_HPP_

#include "SimpleLinearEllipticSolver.hpp"
#include "TetrahedralMesh.hpp"
#include "GaussianQuadratureRule.hpp"

/**
 * A purpose-made elliptic solver that interpolates the source terms from node onto
 * Gauss points, as for a cell-based simulation with PDEs the source will only be
 * known at the cells (nodes), not the Gauss points.
 */
template<unsigned DIM>
class CellBasedPdeSolver
    : public SimpleLinearEllipticSolver<DIM, DIM>
{
private:

    /** The constant in u part of the source term, interpolated onto the current point. */
    double mConstantInUSourceTerm;

    /** The linear in u part of the source term, interpolated onto the current point. */
    double mLinearInUCoeffInSourceTerm;

protected:

    /**
     * The SimpleLinearEllipticSolver version of this method is
     * overloaded using the interpolated source term.
     *
     * @param rPhi
     * @param rGradPhi
     * @param rX
     * @param rU
     * @param rGradU
     * @param pElement
     *
     * @return vector term.
     */
    virtual c_vector<double, 1*(DIM+1)> ComputeVectorTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double, 1>& rU,
        c_matrix<double, 1, DIM>& rGradU /* not used */,
        Element<DIM, DIM>* pElement);

    /**
     * The SimpleLinearEllipticSolver version of this method is
     * overloaded using the interpolated source term.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     *
     * @return The matrix term for the stiffness matrix
     */
    virtual c_matrix<double, 1*(DIM+1), 1*(DIM+1)> ComputeMatrixTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double, 1>& rU,
        c_matrix<double, 1, DIM>& rGradU,
        Element<DIM, DIM>* pElement);

    /**
     * Overridden ResetInterpolatedQuantities() method.
     */
    void ResetInterpolatedQuantities();

    /**
     * Overridden IncrementInterpolatedQuantities() method.
     *
     * @param phiI
     * @param pNode
     */
    void IncrementInterpolatedQuantities(double phiI, const Node<DIM>* pNode);

    /**
     * Create the linear system object if it hasn't been already.
     * Can use an initial solution as PETSc template, or base it on the mesh size.
     *
     * @param initialSolution an initial guess
     */
    void InitialiseForSolve(Vec initialSolution);

public:

    /**
     * Constructor.
     *
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE
     * @param pBoundaryConditions pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    CellBasedPdeSolver(TetrahedralMesh<DIM,DIM>* pMesh,
                                       AbstractLinearEllipticPde<DIM,DIM>* pPde,
                                       BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions,
                                       unsigned numQuadPoints=2);

    /**
     * Destructor.
     */
    ~CellBasedPdeSolver();
};

#endif /*CELLBASEDPDESOLVER_HPP_*/
