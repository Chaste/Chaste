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

#ifndef CELLBASEDPARABOLICPDESYSTEMSOLVER_HPP_
#define CELLBASEDPARABOLICPDESYSTEMSOLVER_HPP_

#include "LinearParabolicPdeSystemSolver.hpp"
#include "TetrahedralMesh.hpp"

/**
 * \todo #2930 document class
 * A purpose-made parabolic solver that interpolates the source terms from node onto
 * Gauss points, as for a cell-based simulation with PDEs the source will only be
 * known at the cells (nodes), not the Gauss points.
 */
template<unsigned DIM, unsigned PROBLEM_DIM=1>
class CellBasedParabolicPdeSystemSolver : public LinearParabolicPdeSystemSolver<DIM, DIM, PROBLEM_DIM>
{
private:

    /** Vector storing each PDE's source term, interpolated onto the current point. */
    std::vector<double> mInterpolatedSourceTerm;

protected:

    /**
     * \todo #2930 fix comments
     * The SimpleLinearParabolicSolver version of this method is
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
    virtual c_vector<double, PROBLEM_DIM*(DIM+1)> ComputeVectorTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double, PROBLEM_DIM>& rU,
        c_matrix<double, PROBLEM_DIM, DIM>& rGradU /* not used */,
        Element<DIM, DIM>* pElement);

    /**
     * Overridden ResetInterpolatedQuantities() method.
     */
    void ResetInterpolatedQuantities();

    /**
     * \todo #2930 fix comments
     * Overridden IncrementInterpolatedQuantities() method.
     *
     * @param phiI
     * @param pNode
     */
    void IncrementInterpolatedQuantities(double phiI, const Node<DIM>* pNode);

public:

    /**
     * Constructor.
     *
     * @param pMesh pointer to the mesh
     * @param pPdeSystem pointer to the PDE system
     * @param pBoundaryConditions pointer to the boundary conditions
     */
    CellBasedParabolicPdeSystemSolver(TetrahedralMesh<DIM, DIM>* pMesh,
                                      AbstractLinearParabolicPdeSystem<DIM, DIM, PROBLEM_DIM>* pPdeSystem,
                                      BoundaryConditionsContainer<DIM, DIM, PROBLEM_DIM>* pBoundaryConditions);

    /**
     * Destructor.
     */
    virtual ~CellBasedParabolicPdeSystemSolver();
};

/*
 * As this class is templated over PROBLEM_DIM, we put the implementation
 * in the header file to avoid explicit instantiation.
 */

template<unsigned DIM, unsigned PROBLEM_DIM>
CellBasedParabolicPdeSystemSolver<DIM, PROBLEM_DIM>::CellBasedParabolicPdeSystemSolver(
    TetrahedralMesh<DIM,DIM>* pMesh,
    AbstractLinearParabolicPdeSystem<DIM, DIM, PROBLEM_DIM>* pPdeSystem,
    BoundaryConditionsContainer<DIM, DIM, PROBLEM_DIM>* pBoundaryConditions)
    : LinearParabolicPdeSystemSolver<DIM, DIM, PROBLEM_DIM>(pMesh, pPdeSystem, pBoundaryConditions)
{
    mInterpolatedSourceTerm.clear();
    mInterpolatedSourceTerm.resize(PROBLEM_DIM);
    std::fill(mInterpolatedSourceTerm.begin(), mInterpolatedSourceTerm.end(), 0);
}

template<unsigned DIM, unsigned PROBLEM_DIM>
CellBasedParabolicPdeSystemSolver<DIM, PROBLEM_DIM>::~CellBasedParabolicPdeSystemSolver()
{
}

template<unsigned DIM, unsigned PROBLEM_DIM>
c_vector<double, PROBLEM_DIM*(DIM+1)> CellBasedParabolicPdeSystemSolver<DIM, PROBLEM_DIM>::ComputeVectorTerm(
    c_vector<double, DIM+1>& rPhi,
    c_matrix<double, DIM, DIM+1>& rGradPhi,
    ChastePoint<DIM>& rX,
    c_vector<double, PROBLEM_DIM>& rU,
    c_matrix<double, PROBLEM_DIM, DIM>& rGradU /* not used */,
    Element<DIM, DIM>* pElement)
{
    assert(mInterpolatedSourceTerm.size() == rU.size());

    double timestep_inverse = PdeSimulationTime::GetPdeTimeStepInverse();
    c_vector<double, PROBLEM_DIM*(DIM+1)> vector_term = zero_vector<double>(PROBLEM_DIM*(DIM+1));

    // Loop over PDEs and populate vector_term
    for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
        double this_dudt_coefficient = this->mpParabolicPdeSystem->ComputeDuDtCoefficientFunction(rX, pde_index);
        c_vector<double, DIM+1> this_vector_term = (mInterpolatedSourceTerm[pde_index] + timestep_inverse*this_dudt_coefficient*rU(pde_index))* rPhi;

        for (unsigned i=0; i<DIM+1; i++)
        {
            vector_term(i*PROBLEM_DIM + pde_index) = this_vector_term(i);
        }
    }

    return vector_term;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void CellBasedParabolicPdeSystemSolver<DIM, PROBLEM_DIM>::ResetInterpolatedQuantities()
{
    std::fill(mInterpolatedSourceTerm.begin(), mInterpolatedSourceTerm.end(), 0);
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void CellBasedParabolicPdeSystemSolver<DIM, PROBLEM_DIM>::IncrementInterpolatedQuantities(double phiI, const Node<DIM>* pNode)
{
    c_vector<double,PROBLEM_DIM> u_at_node;
    for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
        u_at_node[pde_index] = this->GetCurrentSolutionOrGuessValue(pNode->GetIndex(), pde_index);
    }
    for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
        mInterpolatedSourceTerm[pde_index] += phiI * this->mpParabolicPdeSystem->ComputeSourceTermAtNode(*pNode, u_at_node, pde_index);
    }
}

#endif /*CELLBASEDPARABOLICPDESYSTEMSOLVER_HPP_*/