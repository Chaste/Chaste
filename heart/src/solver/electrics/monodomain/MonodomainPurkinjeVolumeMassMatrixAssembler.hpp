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

#ifndef MONODOMAINPURKINJEVOLUMEMASSMATRIXASSEMBLER_HPP_
#define MONODOMAINPURKINJEVOLUMEMASSMATRIXASSEMBLER_HPP_

#include "AbstractFeVolumeIntegralAssembler.hpp"
#include "HeartConfig.hpp"

/**
 * Simple implementation of AbstractFeVolumeIntegralAssembler which provides the volume part of the Monodomain mass matrix
 * for a given mesh, multiplied by a scale factor if required. In other words, the matrix
 * If N is the space dimension, we compute the Matrix M (2N,2N) where
 *
 * M_{ij} = k integral_{domain}  phi_i(x) phi_j(x) dV, if i<N and j<N,
 *
 * where phi_i is the i-th (linear) basis function and k the scale factor (constant
 * throughout the mesh).
 *
 * M_{i,j}= 0, if i>=N or j>=N;
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainPurkinjeVolumeMassMatrixAssembler
    : public AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, 2, false /*no vectors*/, true/*assembles matrices*/, NORMAL>
{
private:

    /** The scale factor. */
    double mScaleFactor;

    /** Whether to use mass lumping or not. */
    bool mUseMassLumping;

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
     * @return the stencil matrix
     */
    c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)>
        ComputeMatrixTerm(
                c_vector<double, ELEMENT_DIM+1> &rPhi,
                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                ChastePoint<SPACE_DIM> &rX,
                c_vector<double,2> &rU,
                c_matrix<double,2,SPACE_DIM> &rGradU /* not used */,
                Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> mass_matrix = outer_prod(rPhi, rPhi);

        if (mUseMassLumping)
        {
            // not implemented yet (and also not implemented yet in MonodomainPurkinjeCableMassMatrixAssembler)
            NEVER_REACHED;
//            for (unsigned row=0; row<ELEMENT_DIM+1; row++)
//            {
//                for (unsigned column=0; column<ELEMENT_DIM+1; column++)
//                {
//                    if (row != column)
//                    {
//                        mass_matrix(row,row) += mass_matrix(row,column);
//                        mass_matrix(row,column) = 0.0;
//                    }
//                }
//            }
        }

        c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ret
                = zero_matrix<double>(2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1));
        // even rows, even columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice00(ret, slice(0, 2, ELEMENT_DIM+1), slice(0, 2, ELEMENT_DIM+1));
        slice00 = mScaleFactor*mass_matrix;
        return ret;
    }

    /**
    * Constructor.
    *
    * @param pMesh the mesh
    * @param scaleFactor the factor with which the multiply the mass matrix. Defaults to 1.0
    * @param useMassLumping whether to use mass matrix lumping or not
    */
    MonodomainPurkinjeVolumeMassMatrixAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh, bool useMassLumping=false, double scaleFactor=1.0)
        : AbstractFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,NORMAL>(pMesh),
          mScaleFactor(scaleFactor),
          mUseMassLumping(useMassLumping)
    {
    }
};

#endif /*MONODOMAINPURKINJEVOLUMEMASSMATRIXASSEMBLER_HPP_*/
