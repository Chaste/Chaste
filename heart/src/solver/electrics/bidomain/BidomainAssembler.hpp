
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
 *  The discretised bidomain equation leads to the linear system (see FEM
 *  implementations document)
 *
 *  [ (chi*C/dt) M + K1    K1   ] [ V^{n+1}   ]  =  [  (chi*C/dt) M V^{n} + M F^{n} + c1_surf ]
 *  [        K1            K2   ] [ PhiE^{n+1}]     [              c2_surf                    ]
 *
 *  where chi is the surface-area to volume ratio, C the capacitance, dt the timestep
 *  M the mass matrix, K1 and K2 stiffness matrices, V^{n} and PhiE^{n} the vector of
 *  voltages and phi_e at time n, F^{n} the vector of (chi*Iionic + Istim) at each node,
 *  and c1_surf and c2_surf vectors arising from any surface stimuli (usually zero).
 *
 *  This assembler is used to assemble the LHS matrix, ie
 *
 *  [ (chi*C/dt) M + K1    K1   ]
 *  [        K1            K2   ]
 *
 *  Hence, this class inherits from AbstractCardiacFeVolumeIntegralAssembler and implements the
 *  methods ComputeMatrixTerm()
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainAssembler : public AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,CARDIAC>
{
protected:
    /** Local cache of the configuration singleton instance*/
    HeartConfig* mpConfig;

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
     * @param numQuadPoints number of quadrature points in each dimension
     */
    BidomainAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                      BidomainTissue<SPACE_DIM>* pTissue,
                      unsigned numQuadPoints = 2);

    /**
     * Destructor.
     */
    ~BidomainAssembler()
    {
    }
};

#endif /*BIDOMAINASSEMBLER_HPP_*/
