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

#include "ExtendedBidomainAssembler.hpp"

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
    slice300 = grad_phi_sigma_i_first_cell_grad_phi;

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
    slice030 = grad_phi_sigma_i_second_cell_grad_phi;

    // third equation, first unknown
    matrix_slice<c_matrix<double, 3*ELEMENT_DIM+3, 3*ELEMENT_DIM+3> >
    slice001(ret, slice (2, 3, ELEMENT_DIM+1), slice (0, 3, ELEMENT_DIM+1));
    slice001 =   grad_phi_sigma_i_first_cell_grad_phi;

    // third equation, second unknown
    matrix_slice<c_matrix<double, 3*ELEMENT_DIM+3, 3*ELEMENT_DIM+3> >
    slice002(ret, slice (2, 3, ELEMENT_DIM+1), slice (1, 3, ELEMENT_DIM+1));
    slice002 =  grad_phi_sigma_i_second_cell_grad_phi;

    // third equation, third unknown
    matrix_slice<c_matrix<double, 3*ELEMENT_DIM+3, 3*ELEMENT_DIM+3> >
    slice003(ret, slice (2, 3, ELEMENT_DIM+1), slice (2, 3, ELEMENT_DIM+1));
    slice003 =  grad_phi_sigma_e_grad_phi + grad_phi_sigma_i_first_cell_grad_phi + grad_phi_sigma_i_second_cell_grad_phi;

    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ExtendedBidomainAssembler<ELEMENT_DIM,SPACE_DIM>::ExtendedBidomainAssembler(
                                AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                ExtendedBidomainTissue<SPACE_DIM>* pTissue)
    : AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,3,true,true,NORMAL>(pMesh,pTissue),
              mpExtendedBidomainTissue(pTissue)
{
    assert(pTissue != NULL);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ExtendedBidomainAssembler<ELEMENT_DIM,SPACE_DIM>::~ExtendedBidomainAssembler()
{
}

// Explicit instantiation
template class ExtendedBidomainAssembler<1,1>;
template class ExtendedBidomainAssembler<2,2>;
template class ExtendedBidomainAssembler<3,3>;
