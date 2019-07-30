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


#include "BidomainCorrectionTermAssembler.hpp"
#include "UblasIncludes.hpp"
#include <boost/numeric/ublas/vector_proxy.hpp>

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
BidomainCorrectionTermAssembler<ELEM_DIM,SPACE_DIM>::BidomainCorrectionTermAssembler(
     AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
     BidomainTissue<SPACE_DIM>* pTissue)
        : AbstractCorrectionTermAssembler<ELEM_DIM,SPACE_DIM,2>(pMesh,pTissue)
{
    mpConfig = HeartConfig::Instance();
    assert(mpConfig->GetUseStateVariableInterpolation());
}

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
c_vector<double,2*(ELEM_DIM+1)> BidomainCorrectionTermAssembler<ELEM_DIM,SPACE_DIM>::ComputeVectorTerm(
    c_vector<double, ELEM_DIM+1> &rPhi,
    c_matrix<double, SPACE_DIM, ELEM_DIM+1> &rGradPhi /* not used */,
    ChastePoint<SPACE_DIM> &rX /* not used */,
    c_vector<double,2> &rU,
    c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
    Element<ELEM_DIM,SPACE_DIM>* pElement /* not used */)
{
    double Am = mpConfig->GetSurfaceAreaToVolumeRatio();

    // compute the ionic current at this quadrature point using the
    // interpolated state variables, and a random choice of cell (all
    // should be the same)
    unsigned node_global_index = pElement->GetNodeGlobalIndex(0);
    AbstractCardiacCellInterface* p_any_cell = this->mpCardiacTissue->GetCardiacCellOrHaloCell(node_global_index);
    double ionic_sv_interp = p_any_cell->GetIIonic(&(this->mStateVariablesAtQuadPoint));

    c_vector<double,2*(ELEM_DIM+1)> ret;

    ublas::vector_slice<c_vector<double, 2*(ELEM_DIM+1)> > slice_V  (ret, slice (0, 2, ELEM_DIM+1));
    ublas::vector_slice<c_vector<double, 2*(ELEM_DIM+1)> > slice_Phi(ret, slice (1, 2, ELEM_DIM+1));

    // add on the SVI ionic current, and take away the original NCI (linearly
    // interpolated ionic current) that would have been added as part of
    // the matrix-based assembly stage.
    noalias(slice_V)   = rPhi * (-Am) * ( ionic_sv_interp - this->mIionicInterp );
    // no correction needed for elliptic equation
    noalias(slice_Phi) = zero_vector<double>(ELEM_DIM+1);

    return ret;
}

// Explicit instantiation
template class BidomainCorrectionTermAssembler<1,1>;
template class BidomainCorrectionTermAssembler<2,2>;
template class BidomainCorrectionTermAssembler<3,3>;
