
/*

Copyright (c) 2005-2016, University of Oxford.
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


#ifndef BIDOMAINASSEMBLER_HPP_
#define BIDOMAINASSEMBLER_HPP_

#include "AbstractCardiacFeVolumeIntegralAssembler.hpp"
#include "BidomainTissue.hpp"
#include "HeartConfig.hpp"

/**
 *  Assembler, used for assembling the LHS matrix of the linear system
 *  that arises when the bidomain equations are discretised, and for assembling
 *  the contribution to the RHS vector that comes from a surface integral.
 *
 *  The discretised bidomain equation leads to the linear system for the default solver
 *  (see FEM implementations document) 
 *
 *  [ (chi*C/dt) M + K1    K1   ] [ V^{n+1}   ]  =  [  (chi*C/dt) M V^{n} + M F^{n} + c1_surf ]
 *  [        K1            K2   ] [ PhiE^{n+1}]     [              c2_surf                    ]
 *
 *  where chi is the surface-area to volume ratio, C the capacitance, dt the timestep
 *  M the mass matrix, K1 and K2 stiffness matrices, V^{n} and PhiE^{n} the vector of
 *  voltages and phi_e at time n, F^{n} the vector of (chi*Iionic + Istim) at each node,
 *  and c1_surf and c2_surf vectors arising from any surface stimuli (usually zero).
 *
 *  Note: The LHS matrix and RHS vector are different for other solvers,
 *  i.e., Backward Euler (Godunov), CN, and SDIRK2O2  
 *
 *  This assembler is used for assembling the coefficient matrix
 *
 *  [ (chi*C/dt) M + mAlpha*K1    mBeta*K1   ]
 *  [        mGamma*K1            mRho*K2    ]
 *
 *  Hence, this class inherits from AbstractCardiacFeVolumeIntegralAssembler and implements the
 *  methods ComputeMatrixTerm()
 * 
 *  mAlpha=mBeta=mGamma=mRho=1.0 for Backward Euler (Godunov) method
 *  mAlpha=0.5 mBeta=mGamma = 1.0 mRho=2.0 for CN method
 *  mAlpha=mBeta=mGamma=mRho=(2.0-sqrt(2.0))/2.0 for SDIRK2O2 method
 *  mAlpha = mBeta = mGamma = mRho = (3.0-sqrt(3.0)/6.0) for the sdirk2O3 method
 */


  /**
   * Enumeration determining which matrices to create.
   */
  namespace BidomainAssembler_helper
  {
      enum MatrixType{BACKWARDEULER, CRANKNICOLSON, SDIRK, SDIRK2O3,  DIRK3O3_stage1, DIRK3O3_stage2, DIRK3O3_stage3, SDIRK3O4};
  }

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainAssembler : public AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,CARDIAC>
{
protected:
    /** Local cache of the configuration singleton instance*/
    HeartConfig* mpConfig;

    /** The type of matrix that needs to be assembled, defaults to Backward Euler*/
    BidomainAssembler_helper::MatrixType mMatrixType;

    /** Factors for the matrix, dependent on which method is used*/
    double mAlpha, mBeta, mGamma, mRho;

    /**
     * ComputeMatrixTerm()
     *
     * This method is called by AssembleOnElement() and tells the assembler
     * the contribution to add to the element LHS matrix.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     * @return stencil matrix
     */
    virtual c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,2> &rU,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);

public:

    /**
     * Constructor stores the mesh and pde and sets up boundary conditions.
     *
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     */
    BidomainAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                      BidomainTissue<SPACE_DIM>* pTissue,
                      BidomainAssembler_helper::MatrixType matrixType = BidomainAssembler_helper::MatrixType());

    /**
     * Destructor.
     */
    virtual ~BidomainAssembler()
    {
    }
};

#endif /*BIDOMAINASSEMBLER_HPP_*/
