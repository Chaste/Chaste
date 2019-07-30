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

    if (!HeartRegionCode::IsRegionBath( pElement->GetUnsignedAttribute() )) // ie if a tissue element
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

// Explicit instantiation
template class ExtendedBidomainMassMatrixAssembler<1>;
template class ExtendedBidomainMassMatrixAssembler<2>;
template class ExtendedBidomainMassMatrixAssembler<3>;
