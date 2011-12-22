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

#include "ExtendedBidomainAssembler.hpp"
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "Exception.hpp"
#include "DistributedVector.hpp"
#include "PdeSimulationTime.hpp"
#include "ConstBoundaryCondition.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,3*(ELEMENT_DIM+1),3*(ELEMENT_DIM+1)>
    ExtendedBidomainAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(
            c_vector<double, ELEMENT_DIM+1> &rPhi,
            c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
            ChastePoint<SPACE_DIM> &rX,
            c_vector<double,3> &rU,
            c_matrix<double, 3, SPACE_DIM> &rGradU /* not used */,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    // get bidomain parameters
    double Am1 = mpExtendedBidomainTissue->GetAmFirstCell();
    double Am2 = mpExtendedBidomainTissue->GetAmSecondCell();
    double Cm1 = mpExtendedBidomainTissue->GetCmFirstCell();
    double Cm2 = mpExtendedBidomainTissue->GetCmSecondCell();

    const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_i_first_cell = mpExtendedBidomainTissue->rGetIntracellularConductivityTensor(pElement->GetIndex());
    const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_i_second_cell = mpExtendedBidomainTissue->rGetIntracellularConductivityTensorSecondCell(pElement->GetIndex());
    const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_e = mpExtendedBidomainTissue->rGetExtracellularConductivityTensor(pElement->GetIndex());

    double delta_t = PdeSimulationTime::GetPdeTimeStep();

    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp_1 = prod(sigma_i_first_cell, rGradPhi);
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_first_cell_grad_phi =
        prod(trans(rGradPhi), temp_1);

    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp_2 = prod(sigma_i_second_cell, rGradPhi);
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_second_cell_grad_phi =
        prod(trans(rGradPhi), temp_2);

    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> basis_outer_prod =
        outer_prod(rPhi, rPhi);

    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp_ext = prod(sigma_e, rGradPhi);
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_e_grad_phi =
        prod(trans(rGradPhi), temp_ext);


    c_matrix<double,3*(ELEMENT_DIM+1),3*(ELEMENT_DIM+1)> ret;

    // first equation, first unknown
    matrix_slice<c_matrix<double, 3*ELEMENT_DIM+3, 3*ELEMENT_DIM+3> >
    slice100(ret, slice (0, 3, ELEMENT_DIM+1), slice (0, 3, ELEMENT_DIM+1));
    slice100 = (Am1*Cm1/delta_t)*basis_outer_prod + grad_phi_sigma_i_first_cell_grad_phi;

    // first equation, second unknown
    matrix_slice<c_matrix<double, 3*ELEMENT_DIM+3, 3*ELEMENT_DIM+3> >
    slice200(ret, slice (0, 3, ELEMENT_DIM+1), slice (1, 3, ELEMENT_DIM+1));
    slice200 = zero_matrix<double>(ELEMENT_DIM+1, ELEMENT_DIM+1);

    // first equation, third unknown
    matrix_slice<c_matrix<double, 3*ELEMENT_DIM+3, 3*ELEMENT_DIM+3> >
    slice300(ret, slice (0, 3, ELEMENT_DIM+1), slice (2, 3, ELEMENT_DIM+1));
    slice300 = - (Am1*Cm1/delta_t)*basis_outer_prod;

    // second equation, first unknown
    matrix_slice<c_matrix<double, 3*ELEMENT_DIM+3, 3*ELEMENT_DIM+3> >
    slice010(ret, slice (1, 3, ELEMENT_DIM+1), slice (0, 3, ELEMENT_DIM+1));
    slice010 = zero_matrix<double>(ELEMENT_DIM+1, ELEMENT_DIM+1);

    // second equation, second unknown
    matrix_slice<c_matrix<double, 3*ELEMENT_DIM+3, 3*ELEMENT_DIM+3> >
    slice020(ret, slice (1, 3, ELEMENT_DIM+1), slice (1, 3, ELEMENT_DIM+1));
    slice020 = (Am2*Cm2/delta_t)*basis_outer_prod + grad_phi_sigma_i_second_cell_grad_phi;

    // second equation, third unknown
    matrix_slice<c_matrix<double, 3*ELEMENT_DIM+3, 3*ELEMENT_DIM+3> >
    slice030(ret, slice (1, 3, ELEMENT_DIM+1), slice (2, 3, ELEMENT_DIM+1));
    slice030 = - (Am2*Cm2/delta_t)*basis_outer_prod;

    // third equation, first unknown
    matrix_slice<c_matrix<double, 3*ELEMENT_DIM+3, 3*ELEMENT_DIM+3> >
    slice001(ret, slice (2, 3, ELEMENT_DIM+1), slice (0, 3, ELEMENT_DIM+1));
    slice001 =   - grad_phi_sigma_i_first_cell_grad_phi;

    // third equation, second unknown
    matrix_slice<c_matrix<double, 3*ELEMENT_DIM+3, 3*ELEMENT_DIM+3> >
    slice002(ret, slice (2, 3, ELEMENT_DIM+1), slice (1, 3, ELEMENT_DIM+1));
    slice002 =  - grad_phi_sigma_i_second_cell_grad_phi;

    // third equation, third unknown
    matrix_slice<c_matrix<double, 3*ELEMENT_DIM+3, 3*ELEMENT_DIM+3> >
    slice003(ret, slice (2, 3, ELEMENT_DIM+1), slice (2, 3, ELEMENT_DIM+1));
    slice003 =  - grad_phi_sigma_e_grad_phi;

    return ret;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ExtendedBidomainAssembler<ELEMENT_DIM,SPACE_DIM>::ExtendedBidomainAssembler(
                                AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                ExtendedBidomainTissue<SPACE_DIM>* pTissue,
                                unsigned numQuadPoints)
    : AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,3,true,true,NORMAL>(pMesh,pTissue,numQuadPoints),
              mpExtendedBidomainTissue(pTissue)
{
    assert(pTissue != NULL);
    mpConfig = HeartConfig::Instance();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ExtendedBidomainAssembler<ELEMENT_DIM,SPACE_DIM>::~ExtendedBidomainAssembler()
{
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class ExtendedBidomainAssembler<1,1>;
template class ExtendedBidomainAssembler<2,2>;
template class ExtendedBidomainAssembler<3,3>;
