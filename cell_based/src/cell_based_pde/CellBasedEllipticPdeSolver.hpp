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

#ifndef CELLBASEDELLIPTICPDESOLVER_HPP_
#define CELLBASEDELLIPTICPDESOLVER_HPP_

#include "SimpleLinearEllipticSolver.hpp"
#include "TetrahedralMesh.hpp"

/**
 * A purpose-made elliptic solver that interpolates the source terms from node onto
 * Gauss points, as for a cell-based simulation with PDEs the source will only be
 * known at the cells (nodes), not the Gauss points.
 */
template<unsigned DIM>
class CellBasedEllipticPdeSolver : public SimpleLinearEllipticSolver<DIM, DIM>
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
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
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
     */
    CellBasedEllipticPdeSolver(TetrahedralMesh<DIM,DIM>* pMesh,
                       AbstractLinearEllipticPde<DIM,DIM>* pPde,
                       BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions);

    /**
     * Destructor.
     */
    virtual ~CellBasedEllipticPdeSolver();
};

#endif /*CELLBASEDELLIPTICPDESOLVER_HPP_*/
