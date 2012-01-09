/*

Copyright (C) University of Oxford, 2005-2012

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

#include "ExtendedBidomainMassMatrixAssembler.hpp"
#include "HeartRegionCodes.hpp"


template<unsigned DIM>
c_matrix<double,3*(DIM+1),3*(DIM+1)> ExtendedBidomainMassMatrixAssembler<DIM>::ComputeMatrixTerm(
            c_vector<double, DIM+1> &rPhi,
            c_matrix<double, DIM, DIM+1> &rGradPhi,
            ChastePoint<DIM> &rX,
            c_vector<double,3> &rU,
            c_matrix<double,3,DIM> &rGradU /* not used */,
            Element<DIM,DIM>* pElement)
{

    c_matrix<double,3*(DIM+1),3*(DIM+1)> ret = zero_matrix<double>(3*(DIM+1),3*(DIM+1));

    if (!HeartRegionCode::IsRegionBath( pElement->GetRegion() )) // ie if a tissue element
    {
        c_matrix<double, DIM+1, DIM+1> basis_outer_prod = outer_prod(rPhi, rPhi);

        // first row,  first column
        matrix_slice<c_matrix<double, 3*DIM+3, 3*DIM+3> >
        slice00(ret, slice (0, 3, DIM+1), slice (0, 3, DIM+1));
        slice00 =  basis_outer_prod;

        // second row, second column
        matrix_slice<c_matrix<double, 3*DIM+3, 3*DIM+3> >
        slice11(ret, slice (1, 3, DIM+1), slice (1, 3, DIM+1));
        slice11 = basis_outer_prod;

        // third row, third column
        matrix_slice<c_matrix<double, 3*DIM+3, 3*DIM+3> >
        slice22(ret, slice (2, 3, DIM+1), slice (2, 3, DIM+1));
        slice22 = basis_outer_prod;
    }
    return ret;
}



///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class ExtendedBidomainMassMatrixAssembler<1>;
template class ExtendedBidomainMassMatrixAssembler<2>;
template class ExtendedBidomainMassMatrixAssembler<3>;
