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

#ifndef MONODOMAINASSEMBLER_HPP_
#define MONODOMAINASSEMBLER_HPP_


#include "AbstractCardiacFeVolumeIntegralAssembler.hpp"
#include "MonodomainTissue.hpp"
#include "MassMatrixAssembler.hpp"
#include "MonodomainStiffnessMatrixAssembler.hpp"

/**
 *  Assembler, mainly used for assembling the LHS matrix of the linear system
 *  that arises when the monodomain equations are discretised.
 *
 *  The discretised monodomain equation leads to the linear system (see FEM
 *  implementations document)
 *
 *  ( (chi*C/dt) M  + K ) V^{n+1} = (chi*C/dt) M V^{n} + M F^{n} + c_surf
 *
 *  where chi is the surface-area to volume ratio, C the capacitance, dt the timestep
 *  M the mass matrix, K the stiffness matrix, V^{n} the vector of voltages at time n,
 *  F^{n} the vector of (chi*Iionic + Istim) at each node, and c_surf a vector
 *  arising from any surface stimuli (usually zero).
 *
 *  This assembler is used for assembling the matrix A :=(chi*C/dt) M  + K.
 *  Hence, this class inherits from AbstractCardiacFeVolumeIntegralAssembler and implements the
 *  method ComputeMatrixTerm().
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainAssembler
   : public AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,1,false,true,CARDIAC>
{
protected:
    /** This assembler uses another assembler, though just for calling the
     *  ComputeMatrixTerm() method. */
    MassMatrixAssembler<ELEMENT_DIM, SPACE_DIM> mMassMatrixAssembler;

    /** This assembler uses another assembler, though just for calling the
     *  ComputeMatrixTerm() method. */
    MonodomainStiffnessMatrixAssembler<ELEMENT_DIM, SPACE_DIM> mStiffnessMatrixAssembler;

public:

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
     */
    c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> ComputeMatrixTerm(
                c_vector<double, ELEMENT_DIM+1> &rPhi,
                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                ChastePoint<SPACE_DIM> &rX,
                c_vector<double,1> &rU,
                c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
                Element<ELEMENT_DIM,SPACE_DIM>* pElement);


    /**
     * Constructor
     *
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     */
    MonodomainAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                        MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue);
};

#endif /*MONODOMAINASSEMBLER_HPP_*/
