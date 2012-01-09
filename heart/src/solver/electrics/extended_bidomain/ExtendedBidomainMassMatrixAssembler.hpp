
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



#ifndef EXTENDEDBIDOMAINMASSMATRIXASSEMBLER_HPP_
#define EXTENDEDBIDOMAINMASSMATRIXASSEMBLER_HPP_



#include "AbstractCardiacFeVolumeIntegralAssembler.hpp"

/**
 *  Constructs a matrix with the mass matrix in the voltage-voltage block.
 *
 *  Ie. IF the extended bidomain unknowns were ordered [phi1_1,..,phi1_n, phi2_1, ..., phi2_n, phie_1,..,phie_n], the
 *  matrix would be, in block form
 *
 *  [ M 0 0]
 *  [ 0 M 0]
 *  [ 0 0 M]
 *
 *  where M is the standard nxn mass matrix.
 *
 *  Since the bidomain ordering is not [phi1_1,..,phi1_n,phi2_1,..,phi2_n, phie_1,...phie_n]
 *  but [phi1_1,phi2_1,phie_1,..,phi1_n,phi2_n,phie_n], the matrix has a different form.
 *
 *
 */
template<unsigned DIM>
class ExtendedBidomainMassMatrixAssembler : public AbstractFeVolumeIntegralAssembler<DIM,DIM,3,false,true,NORMAL>
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
     */

    virtual c_matrix<double,3*(DIM+1),3*(DIM+1)> ComputeMatrixTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,3> &rU,
        c_matrix<double,3,DIM> &rGradU /* not used */,
        Element<DIM,DIM>* pElement);

public:

    /**
     * Constructor
     *
     * @param pMesh pointer to the mesh
     */
    ExtendedBidomainMassMatrixAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh)
        : AbstractFeVolumeIntegralAssembler<DIM,DIM,3,false,true,NORMAL>(pMesh)
    {
    }

    /**
     * Destructor.
     */
    ~ExtendedBidomainMassMatrixAssembler()
    {
    }
};


#endif /*EXTENDEDBIDOMAINMASSMATRIXASSEMBLER_HPP_*/
