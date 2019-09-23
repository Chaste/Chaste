
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

#include "BidomainAssembler.hpp"
#include "PdeSimulationTime.hpp"
#include "UblasIncludes.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)>
    BidomainAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(
            c_vector<double, ELEMENT_DIM+1> &rPhi,
            c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
            ChastePoint<SPACE_DIM> &rX,
            c_vector<double,2> &rU,
            c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    // get bidomain parameters
    double Am = this->mpConfig->GetSurfaceAreaToVolumeRatio();
    double Cm = this->mpConfig->GetCapacitance();

    const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_i = this->mpCardiacTissue->rGetIntracellularConductivityTensor(pElement->GetIndex());
    const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_e = this->mpCardiacTissue->rGetExtracellularConductivityTensor(pElement->GetIndex());


    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp = prod(sigma_i, rGradPhi);
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_grad_phi =
        prod(trans(rGradPhi), temp);

    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> basis_outer_prod =
        outer_prod(rPhi, rPhi);

    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp2 = prod(sigma_e, rGradPhi);
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_e_grad_phi =
        prod(trans(rGradPhi), temp2);


    c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ret;

    // even rows, even columns
    matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
    slice00(ret, slice (0, 2, ELEMENT_DIM+1), slice (0, 2, ELEMENT_DIM+1));
    slice00 = (Am*Cm*PdeSimulationTime::GetPdeTimeStepInverse())*basis_outer_prod + grad_phi_sigma_i_grad_phi;

    // odd rows, even columns
    matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
    slice10(ret, slice (1, 2, ELEMENT_DIM+1), slice (0, 2, ELEMENT_DIM+1));
    slice10 = grad_phi_sigma_i_grad_phi;

    // even rows, odd columns
    matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
    slice01(ret, slice (0, 2, ELEMENT_DIM+1), slice (1, 2, ELEMENT_DIM+1));
    slice01 = grad_phi_sigma_i_grad_phi;

    // odd rows, odd columns
    matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
    slice11(ret, slice (1, 2, ELEMENT_DIM+1), slice (1, 2, ELEMENT_DIM+1));
    slice11 = grad_phi_sigma_i_grad_phi + grad_phi_sigma_e_grad_phi;

    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BidomainAssembler<ELEMENT_DIM,SPACE_DIM>::BidomainAssembler(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            BidomainTissue<SPACE_DIM>* pTissue)
    : AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,CARDIAC>(pMesh,pTissue)
{
    assert(pTissue != NULL);
}

// Explicit instantiation
template class BidomainAssembler<1,1>;
template class BidomainAssembler<2,2>;
template class BidomainAssembler<3,3>;
