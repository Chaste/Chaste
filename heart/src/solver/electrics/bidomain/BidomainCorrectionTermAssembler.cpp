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


#include "BidomainCorrectionTermAssembler.hpp"
#include "UblasIncludes.hpp"
#include <boost/numeric/ublas/vector_proxy.hpp>

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
BidomainCorrectionTermAssembler<ELEM_DIM,SPACE_DIM>::BidomainCorrectionTermAssembler(
     AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
     BidomainTissue<SPACE_DIM>* pTissue,
     unsigned numQuadPoints)
        : AbstractCorrectionTermAssembler<ELEM_DIM,SPACE_DIM,2>(pMesh,pTissue,numQuadPoints)
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
    AbstractCardiacCell* p_any_cell = this->mpCardiacTissue->GetCardiacCellOrHaloCell(node_global_index);
    double ionic_sv_interp = p_any_cell->GetIIonic(&(this->mStateVariablesAtQuadPoint));

    c_vector<double,2*(ELEM_DIM+1)> ret;

    vector_slice<c_vector<double, 2*(ELEM_DIM+1)> > slice_V  (ret, slice (0, 2, ELEM_DIM+1));
    vector_slice<c_vector<double, 2*(ELEM_DIM+1)> > slice_Phi(ret, slice (1, 2, ELEM_DIM+1));

    // add on the SVI ionic current, and take away the original NCI (linearly
    // interpolated ionic current) that would have been added as part of
    // the matrix-based assembly stage.
    noalias(slice_V)   = rPhi * (-Am) * ( ionic_sv_interp - this->mIionicInterp );
    // no correction needed for elliptic equation
    noalias(slice_Phi) = zero_vector<double>(ELEM_DIM+1);

    return ret;
}


///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class BidomainCorrectionTermAssembler<1,1>;
template class BidomainCorrectionTermAssembler<2,2>;
template class BidomainCorrectionTermAssembler<3,3>;
