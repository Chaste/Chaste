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



#ifndef MONODOMAINPURKINJEVOLUMEASSEMBLER_CPP_
#define MONODOMAINPURKINJEVOLUMEASSEMBLER_CPP_

#include "MonodomainPurkinjeVolumeAssembler.hpp"
#include "PdeSimulationTime.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> MonodomainPurkinjeVolumeAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(
                c_vector<double, ELEMENT_DIM+1> &rPhi,
                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                ChastePoint<SPACE_DIM> &rX,
                c_vector<double,2> &rU,
                c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
                Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    c_vector<double,1> empty_u; /* should pass in rU(1) here, but it won't be used */
    c_matrix<double,1,SPACE_DIM> empty_grad_u; /* ditto above */

    c_matrix<double,ELEMENT_DIM+1,ELEMENT_DIM+1> normal_monodomain_mat
        = mMonodomainAssembler.ComputeMatrixTerm(rPhi,rGradPhi,rX,empty_u,empty_grad_u,pElement);

    c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ret
        = zero_matrix<double>(2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1));

    // even rows, even columns
    matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
    slice00(ret, slice(0, 2, ELEMENT_DIM+1), slice(0, 2, ELEMENT_DIM+1));
    slice00 = normal_monodomain_mat;

    return ret;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainPurkinjeVolumeAssembler<ELEMENT_DIM,SPACE_DIM>::MonodomainPurkinjeVolumeAssembler(
                        AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                        MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue)
    : AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,CARDIAC>(pMesh,pTissue),
      mMonodomainAssembler(pMesh,pTissue)
{
}



///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class MonodomainPurkinjeVolumeAssembler<2,2>;
template class MonodomainPurkinjeVolumeAssembler<3,3>;

#endif /*MONODOMAINPURKINJEVOLUMEASSEMBLER_CPP_*/
