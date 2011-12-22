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
                                  BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                                  unsigned numQuadPoints)
        : AbstractAssemblerSolverHybrid<ELEMENT_DIM,SPACE_DIM,1,NORMAL>(pMesh,pBoundaryConditions,numQuadPoints),
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

//////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////

template class SimpleLinearEllipticSolver<1,1>;
template class SimpleLinearEllipticSolver<1,2>;
template class SimpleLinearEllipticSolver<1,3>;
template class SimpleLinearEllipticSolver<2,2>;
template class SimpleLinearEllipticSolver<3,3>;
