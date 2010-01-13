/*

Copyright (C) University of Oxford, 2005-2010

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

#include "TissueSimulationWithNutrientsAssembler.hpp"
#include "TetrahedralMesh.hpp"
#include "SimpleLinearEllipticAssembler.hpp"
#include "GaussianQuadratureRule.hpp"


template<unsigned DIM>
TissueSimulationWithNutrientsAssembler<DIM>::TissueSimulationWithNutrientsAssembler(TetrahedralMesh<DIM,DIM>* pMesh,
                              AbstractLinearEllipticPde<DIM,DIM>* pPde,
                              BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions,
                              unsigned numQuadPoints) :
        BaseClassType(pMesh, pPde, pBoundaryConditions, numQuadPoints)
{
}

template<unsigned DIM>
TissueSimulationWithNutrientsAssembler<DIM>::~TissueSimulationWithNutrientsAssembler()
{
}

template<unsigned DIM>
c_vector<double, 1*(DIM+1)> TissueSimulationWithNutrientsAssembler<DIM>::ComputeVectorTerm(
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
c_matrix<double, 1*(DIM+1), 1*(DIM+1)> TissueSimulationWithNutrientsAssembler<DIM>::ComputeMatrixTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double, 1>& rU,
        c_matrix<double, 1, DIM>& rGradU,
        Element<DIM, DIM>* pElement)
{
    c_matrix<double, DIM, DIM> pde_diffusion_term = this->mpEllipticPde->ComputeDiffusionTerm(rX);

    // if statement just saves computing phi*phi^T if it is to be multiplied by zero
    if (mLinearInUCoeffInSourceTerm!=0)
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
void TissueSimulationWithNutrientsAssembler<DIM>::ResetInterpolatedQuantities()
{
    mConstantInUSourceTerm = 0;
    mLinearInUCoeffInSourceTerm = 0;
}

template<unsigned DIM>
void TissueSimulationWithNutrientsAssembler<DIM>::IncrementInterpolatedQuantities(double phiI, const Node<DIM>* pNode)
{
    mConstantInUSourceTerm += phiI * this->mpEllipticPde->ComputeConstantInUSourceTermAtNode(*pNode);
    mLinearInUCoeffInSourceTerm += phiI * this->mpEllipticPde->ComputeLinearInUCoeffInSourceTermAtNode(*pNode);
}

template<unsigned DIM>
void TissueSimulationWithNutrientsAssembler<DIM>::InitialiseForSolve(Vec initialSolution)
{
    // linear system created here
    BaseClassType::InitialiseForSolve(initialSolution);
    
    this->mpLinearSystem->SetMatrixIsSymmetric(true);
}
    


/**
 * Helper structure that defines typedefs specifying where in the
 * hierarchy of assembler classes various methods are defined, so
 * that we can avoid virtual method overhead by setting which method
 * is called at compile time.
 */
template<unsigned DIM>
struct AssemblerTraits<TissueSimulationWithNutrientsAssembler<DIM> >
{
    /** The class in which ComputeVectorTerm is defined. */
    typedef TissueSimulationWithNutrientsAssembler<DIM> CVT_CLASS;
    /** The class in which ComputeMatrixTerm is defined. */
    typedef TissueSimulationWithNutrientsAssembler<DIM> CMT_CLASS;
    /**  The class in which IncrementInterpolatedQuantities and ResetInterpolatedQuantities are defined. */
    typedef TissueSimulationWithNutrientsAssembler<DIM> INTERPOLATE_CLASS;
};

//////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////

template class TissueSimulationWithNutrientsAssembler<1>;
template class TissueSimulationWithNutrientsAssembler<2>;
template class TissueSimulationWithNutrientsAssembler<3>;
