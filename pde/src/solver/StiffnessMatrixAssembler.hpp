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

#ifndef STIFFNESSMATRIXASSEMBLER_HPP_
#define STIFFNESSMATRIXASSEMBLER_HPP_

#include "AbstractFeVolumeIntegralAssembler.hpp"

/**
 * Simple implementation of AbstractFeVolumeIntegralAssembler which provides the stiffness matrix
 *
 * K_{ij} =  integral_{domain} grad phi_i(x) . grad phi_j(x) dV
 *
 * where phi_i is the i-th (linear) basis function.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class StiffnessMatrixAssembler
    : public AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, 1, false /*no vectors*/, true/*assembles matrices*/, NORMAL>
{
public:

    /**
     * Implemented ComputeMatrixTerm(), defined in AbstractFeVolumeIntegralAssembler.
     * See documentation in that class.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases.
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i).
     * @param rX The point in space.
     * @param rU The unknown as a vector, u(i) = u_i.
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j).
     * @param pElement Pointer to the element.
     * @return computed stencil matrix
     */
    c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)>
        ComputeMatrixTerm(
                c_vector<double, ELEMENT_DIM+1> &rPhi,
                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                ChastePoint<SPACE_DIM> &rX,
                c_vector<double,1> &rU,
                c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
                Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        return prod( trans(rGradPhi), rGradPhi );
    }

    /**
     * Constructor.
     *
     * @param pMesh the mesh
     */
    StiffnessMatrixAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
        : AbstractFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,1,false,true,NORMAL>(pMesh)
    {
    }
};

#endif /*STIFFNESSMATRIXASSEMBLER_HPP_*/
