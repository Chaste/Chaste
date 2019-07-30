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

#include "SimpleLinearEllipticSolver.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double, 1*(ELEMENT_DIM+1), 1*(ELEMENT_DIM+1)>SimpleLinearEllipticSolver<ELEMENT_DIM,SPACE_DIM>:: ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    c_matrix<double, SPACE_DIM, SPACE_DIM> pde_diffusion_term = mpEllipticPde->ComputeDiffusionTerm(rX);

    // This if statement just saves computing phi*phi^T if it is to be multiplied by zero
    if (mpEllipticPde->ComputeLinearInUCoeffInSourceTerm(rX,pElement)!=0)
    {
        return   prod( trans(rGradPhi), c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) )
               - mpEllipticPde->ComputeLinearInUCoeffInSourceTerm(rX,pElement)*outer_prod(rPhi,rPhi);
    }
    else
    {
        return   prod( trans(rGradPhi), c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) );
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,1*(ELEMENT_DIM+1)> SimpleLinearEllipticSolver<ELEMENT_DIM,SPACE_DIM>::ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    return mpEllipticPde->ComputeConstantInUSourceTerm(rX, pElement) * rPhi;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SimpleLinearEllipticSolver<ELEMENT_DIM,SPACE_DIM>::SimpleLinearEllipticSolver(
                                  AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                  AbstractLinearEllipticPde<ELEMENT_DIM,SPACE_DIM>* pPde,
                                  BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions)
        : AbstractAssemblerSolverHybrid<ELEMENT_DIM,SPACE_DIM,1,NORMAL>(pMesh,pBoundaryConditions),
          AbstractStaticLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>(pMesh)
{
    mpEllipticPde = pPde;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SimpleLinearEllipticSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    AbstractLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>::InitialiseForSolve(initialSolution);
    assert(this->mpLinearSystem);
    this->mpLinearSystem->SetMatrixIsSymmetric(true);
    this->mpLinearSystem->SetKspType("cg");
}

// Explicit instantiation
template class SimpleLinearEllipticSolver<1,1>;
template class SimpleLinearEllipticSolver<1,2>;
template class SimpleLinearEllipticSolver<1,3>;
template class SimpleLinearEllipticSolver<2,2>;
template class SimpleLinearEllipticSolver<3,3>;
