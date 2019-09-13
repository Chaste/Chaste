
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

#ifndef BIDOMAINMASSMATRIXASSEMBLER_HPP_
#define BIDOMAINMASSMATRIXASSEMBLER_HPP_

#include "AbstractFeVolumeIntegralAssembler.hpp"

/**
 *  Constructs a matrix with the mass matrix in the voltage-voltage block.
 *
 *  Ie. IF the bidomain unknowns were ordered [V1,..,Vn,phie_1,..,phie_n], the
 *  matrix would be, in block form
 *
 *  [ M 0 ]
 *  [ 0 0 ]
 *
 *  where M is the standard nxn mass matrix.
 *
 *  Since the bidomain ordering is not [V1,..,Vn,phie_1,..,phie_n]
 *  but [V1,phie1,..,Vn,phie_n], the matrix has a different form.
 *
 *  WORKS FOR BATH PROBLEMS AS WELL AS NON-BATH PROBLEMS
 *  (sets zeros in the voltage-voltage block for bath nodes, ie
 *  [ M 0 0 0 ]
 *  [ 0 0 0 0 ]
 *  [ 0 0 0 0 ]
 *  [ 0 0 0 0 ]
 *  where the ordering is Vtissue, Vbath, phi_tissue, phi_bath
 *
 */
template<unsigned DIM>
class BidomainMassMatrixAssembler : public AbstractFeVolumeIntegralAssembler<DIM,DIM,2,false,true,CARDIAC>
{
protected:
    /**
     * ComputeMatrixTerm()
     *
     * This method is called by AssembleOnElement() and tells the assembler
     * the contribution to add to the element stiffness matrix.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     * @return stencil matrix
     *
     */

    c_matrix<double,2*(DIM+1),2*(DIM+1)> ComputeMatrixTerm(
            c_vector<double, DIM+1> &rPhi,
            c_matrix<double, DIM, DIM+1> &rGradPhi,
            ChastePoint<DIM> &rX,
            c_vector<double,2> &rU,
            c_matrix<double,2,DIM> &rGradU /* not used */,
            Element<DIM,DIM>* pElement);

public:

    /**
     * Constructor
     *
     * @param pMesh pointer to the mesh
     */
    BidomainMassMatrixAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh)
        : AbstractFeVolumeIntegralAssembler<DIM,DIM,2,false,true,CARDIAC>(pMesh)
    {
    }

    /**
     * Destructor.
     */
    ~BidomainMassMatrixAssembler()
    {
    }
};


#endif /*BIDOMAINMASSMATRIXASSEMBLER_HPP_*/
