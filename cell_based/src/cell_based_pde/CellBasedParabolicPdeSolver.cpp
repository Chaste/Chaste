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

#include "CellBasedParabolicPdeSolver.hpp"

template<unsigned DIM>
CellBasedParabolicPdeSolver<DIM>::CellBasedParabolicPdeSolver(TetrahedralMesh<DIM,DIM>* pMesh,
                              AbstractLinearParabolicPde<DIM,DIM>* pPde,
                              BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions)
     : SimpleLinearParabolicSolver<DIM, DIM>(pMesh, pPde, pBoundaryConditions)
{
}

template<unsigned DIM>
CellBasedParabolicPdeSolver<DIM>::~CellBasedParabolicPdeSolver()
{
}

template<unsigned DIM>
c_vector<double, 1*(DIM+1)> CellBasedParabolicPdeSolver<DIM>::ComputeVectorTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double, 1>& rU,
        c_matrix<double, 1, DIM>& rGradU /* not used */,
        Element<DIM, DIM>* pElement)
{
  return (mInterpolatedSourceTerm
        + PdeSimulationTime::GetPdeTimeStepInverse() * this->mpParabolicPde->ComputeDuDtCoefficientFunction(rX) * rU(0)) * rPhi;
}

template<unsigned DIM>
c_matrix<double, 1*(DIM+1), 1*(DIM+1)> CellBasedParabolicPdeSolver<DIM>::ComputeMatrixTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double, 1>& rU,
        c_matrix<double, 1, DIM>& rGradU,
        Element<DIM, DIM>* pElement)
{

     c_matrix<double, DIM, DIM> pde_diffusion_term = this->mpParabolicPde->ComputeDiffusionTerm(rX, pElement);

     return    prod( trans(rGradPhi), c_matrix<double, DIM, DIM+1>(prod(pde_diffusion_term, rGradPhi)) )
                + PdeSimulationTime::GetPdeTimeStepInverse() * this->mpParabolicPde->ComputeDuDtCoefficientFunction(rX) * outer_prod(rPhi, rPhi);
}

template<unsigned DIM>
void CellBasedParabolicPdeSolver<DIM>::ResetInterpolatedQuantities()
{
    mInterpolatedSourceTerm = 0;
}

template<unsigned DIM>
void CellBasedParabolicPdeSolver<DIM>::IncrementInterpolatedQuantities(double phiI, const Node<DIM>* pNode)
{
    unsigned index_of_unknown = 0;
    double u_at_node = this->GetCurrentSolutionOrGuessValue(pNode->GetIndex(), index_of_unknown);

    mInterpolatedSourceTerm += phiI*this->mpParabolicPde->ComputeSourceTermAtNode(*pNode,u_at_node);
}

// Explicit instantiation
template class CellBasedParabolicPdeSolver<1>;
template class CellBasedParabolicPdeSolver<2>;
template class CellBasedParabolicPdeSolver<3>;
