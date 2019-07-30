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

#include "CellBasedEllipticPdeSolver.hpp"

template<unsigned DIM>
CellBasedEllipticPdeSolver<DIM>::CellBasedEllipticPdeSolver(TetrahedralMesh<DIM,DIM>* pMesh,
                              AbstractLinearEllipticPde<DIM,DIM>* pPde,
                              BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions)
    : SimpleLinearEllipticSolver<DIM, DIM>(pMesh, pPde, pBoundaryConditions)
{
}

template<unsigned DIM>
CellBasedEllipticPdeSolver<DIM>::~CellBasedEllipticPdeSolver()
{
}

template<unsigned DIM>
c_vector<double, 1*(DIM+1)> CellBasedEllipticPdeSolver<DIM>::ComputeVectorTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double, 1>& rU,
        c_matrix<double, 1, DIM>& rGradU /* not used */,
        Element<DIM, DIM>* pElement)
{
    return mConstantInUSourceTerm * rPhi;
}

template<unsigned DIM>
c_matrix<double, 1*(DIM+1), 1*(DIM+1)> CellBasedEllipticPdeSolver<DIM>::ComputeMatrixTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double, 1>& rU,
        c_matrix<double, 1, DIM>& rGradU,
        Element<DIM, DIM>* pElement)
{
    c_matrix<double, DIM, DIM> pde_diffusion_term = this->mpEllipticPde->ComputeDiffusionTerm(rX);

    // This if statement just saves computing phi*phi^T if it is to be multiplied by zero
    if (mLinearInUCoeffInSourceTerm != 0)
    {
        return   prod( trans(rGradPhi), c_matrix<double, DIM, DIM+1>(prod(pde_diffusion_term, rGradPhi)) )
               - mLinearInUCoeffInSourceTerm * outer_prod(rPhi,rPhi);
    }
    else
    {
        return   prod( trans(rGradPhi), c_matrix<double, DIM, DIM+1>(prod(pde_diffusion_term, rGradPhi)) );
    }
}

template<unsigned DIM>
void CellBasedEllipticPdeSolver<DIM>::ResetInterpolatedQuantities()
{
    mConstantInUSourceTerm = 0;
    mLinearInUCoeffInSourceTerm = 0;
}

template<unsigned DIM>
void CellBasedEllipticPdeSolver<DIM>::IncrementInterpolatedQuantities(double phiI, const Node<DIM>* pNode)
{
    mConstantInUSourceTerm += phiI * this->mpEllipticPde->ComputeConstantInUSourceTermAtNode(*pNode);
    mLinearInUCoeffInSourceTerm += phiI * this->mpEllipticPde->ComputeLinearInUCoeffInSourceTermAtNode(*pNode);
}

template<unsigned DIM>
void CellBasedEllipticPdeSolver<DIM>::InitialiseForSolve(Vec initialSolution)
{
    // Linear system created here
    SimpleLinearEllipticSolver<DIM,DIM>::InitialiseForSolve(initialSolution);

    this->mpLinearSystem->SetMatrixIsSymmetric(true);
}

// Explicit instantiation
template class CellBasedEllipticPdeSolver<1>;
template class CellBasedEllipticPdeSolver<2>;
template class CellBasedEllipticPdeSolver<3>;
