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



#include "BidomainMassMatrixAssembler.hpp"
#include "HeartRegionCodes.hpp"

template<unsigned DIM>
c_matrix<double,2*(DIM+1),2*(DIM+1)> BidomainMassMatrixAssembler<DIM>::ComputeMatrixTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,2> &rU,
        c_matrix<double,2,DIM> &rGradU /* not used */,
        Element<DIM,DIM>* pElement)
{
    c_matrix<double,2*(DIM+1),2*(DIM+1)> ret = zero_matrix<double>(2*(DIM+1), 2*(DIM+1));

    if (!HeartRegionCode::IsRegionBath( pElement->GetRegion() ))
    {
        c_matrix<double, DIM+1, DIM+1> basis_outer_prod = outer_prod(rPhi, rPhi);

        // even rows, even columns
        matrix_slice<c_matrix<double, 2*DIM+2, 2*DIM+2> >
        slice00(ret, slice (0, 2, DIM+1), slice (0, 2, DIM+1));
        slice00 =  basis_outer_prod;

        // odd rows, even columns: are zero
        // even rows, odd columns: are zero
        // odd rows, odd columns: are zero
    }

    return ret;
}



///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class BidomainMassMatrixAssembler<1>;
template class BidomainMassMatrixAssembler<2>;
template class BidomainMassMatrixAssembler<3>;
