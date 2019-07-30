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

#ifndef MONODOMAINPURKINJECABLEMASSMATRIXASSEMBLER_HPP_
#define MONODOMAINPURKINJECABLEMASSMATRIXASSEMBLER_HPP_

#include "AbstractFeCableIntegralAssembler.hpp"
#include "HeartConfig.hpp"

/**
 * Simple implementation of AbstractFeCableIntegralAssembler which provides the cable part of the Monodomain mass matrix
 * for a given MixedDimesionMesh, multiplied by a scale factor if required. In other words, the matrix
 * If N is the space dimension, we compute the Matrix M (2N,2N) where
 *
 * M_{ij} = k integral_{domain}  phi_i(x) phi_j(x) dV, if i>=N and j>=N, and {domain} is a cable element.
 *
 * where phi_i is the i-th (linear) basis function and k the scale factor (constant
 * throughout the mesh).
 *
 * M_{i,j}= 0, if i<N or j<N;
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainPurkinjeCableMassMatrixAssembler : public AbstractFeCableIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,NORMAL>
{
private:
    /** Cable element dimension. */
    static const unsigned CABLE_ELEMENT_DIM = 1;

    /** Number of nodes in a cable element. */
    static const unsigned NUM_CABLE_ELEMENT_NODES = 2;

    /** Whether to use mass lumping or not. */
    bool mUseMassLumping;

public:
    /**
     * Implemented ComputeMatrixTerm(), defined in AbstractFeCableIntegralAssembler.
     * See documentation in that class.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases.
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i).
     * @param rX The point in space.
     * @param rU The unknown as a vector, u(i) = u_i.
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j).
     * @param pElement Pointer to the element.
     * @return the stencil matrix
     */
    c_matrix<double,2*NUM_CABLE_ELEMENT_NODES,2*NUM_CABLE_ELEMENT_NODES/*2=PROBLEM_DIM here*/>
        ComputeCableMatrixTerm(
            c_vector<double, NUM_CABLE_ELEMENT_NODES>& rPhi,
            c_matrix<double, SPACE_DIM, NUM_CABLE_ELEMENT_NODES>& rGradPhi,
            ChastePoint<SPACE_DIM>& rX,
            c_vector<double,2>& rU,
            c_matrix<double,2, SPACE_DIM>& rGradU,
            Element<CABLE_ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        c_matrix<double,4, 4> ret = zero_matrix<double>(4, 4);

        //We have to scale the assembled matrix by the cross sectional area of the Purkinje
        //fibre to ensure conservation of current at branch points. See #1899.
        const double fibre_cross_section_area = M_PI*pElement->GetAttribute()*pElement->GetAttribute();

        for(unsigned i=0; i<2; i++) // 2 = number of basis functions per cable element
        {
            for(unsigned j=0; j<2; j++)  // 2 = number of basis functions per cable element
            {
                ret(2*i,  2*j)   = 0;  // [V,V] block
                ret(2*i+1,2*j)   = 0;  // [Vpurkinje,V] block
                ret(2*i,  2*j+1) = 0;  // [V,Vpurkinje] block
                ret(2*i+1,2*j+1) = rPhi(i)*rPhi(j);

                ret(2*i+1,2*j+1) *= fibre_cross_section_area;
            }
        }

        if (mUseMassLumping)
        {
            // not implemented yet (and also not implemented yet in MonodomainPurkinjeVolumeMassMatrixAssembler)
            NEVER_REACHED;
        }

        return ret;
    }

   /**
    * Constructor.
    *
    * @param pMesh the mesh
    * @param useMassLumping whether to use mass matrix lumping or not
    */
    MonodomainPurkinjeCableMassMatrixAssembler(MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* pMesh, bool useMassLumping=false)
       : AbstractFeCableIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,NORMAL>(pMesh),
         mUseMassLumping(useMassLumping)
    {
    }
};

#endif /*MONODOMAINPURKINJECABLEMASSMATRIXASSEMBLER_HPP_*/
