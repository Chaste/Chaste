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


#ifndef MONODOMAINCORRECTIONTERM_HPP_
#define MONODOMAINCORRECTIONTERM_HPP_

#include "AbstractCorrectionTermAssembler.hpp"
#include "MonodomainTissue.hpp"


/**
 * An assembler which computes the correction term to add to the
 * RHS vector if using state variable interpolation (SVI), as well as determining
 * which elements should be corrected on. The formula to determine which
 * elements SVI is used is delta Iionic > TOL, where delta Iionic is the max
 * difference between nodal ionic values, and TOL is chosen conservatively
 * to be 1uA/cm^2^. See wiki page ChasteGuides/StateVariableInterpolation
 * for more details.
 */
template<unsigned ELEM_DIM,unsigned SPACE_DIM>
class MonodomainCorrectionTermAssembler
    : public AbstractCorrectionTermAssembler<ELEM_DIM,SPACE_DIM,1>
{
protected:
    /**
     * This method is called by AssembleOnElement and tells the assembler
     * the contribution to add to the element stiffness vector.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, rU(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     * @return stencil matrix
     */
    c_vector<double,1*(ELEM_DIM+1)> ComputeVectorTerm(
                c_vector<double, ELEM_DIM+1> &rPhi,
                c_matrix<double, SPACE_DIM, ELEM_DIM+1> &rGradPhi /* not used */,
                ChastePoint<SPACE_DIM> &rX /* not used */,
                c_vector<double,1> &rU,
                c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
                Element<ELEM_DIM,SPACE_DIM>* pElement);
public:

    /**
     * Constructor.
     *
     * @param pMesh  pointer to the mesh
     * @param pTissue  pointer to the cardiac tissue
     */
    MonodomainCorrectionTermAssembler(AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
                                      MonodomainTissue<ELEM_DIM,SPACE_DIM>* pTissue);
};

#endif /*MONODOMAINCORRECTIONTERM_HPP_*/
