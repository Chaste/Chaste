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
    /** Local cache of the configuration singleton instance*/
    HeartConfig* mpConfig;

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
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    MonodomainAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                        MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                        unsigned numQuadPoints = 2);
};

#endif /*MONODOMAINASSEMBLER_HPP_*/
