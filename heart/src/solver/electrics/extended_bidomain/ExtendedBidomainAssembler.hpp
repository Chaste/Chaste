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


#ifndef EXTENDEDBIDOMAINASSEMBLER_HPP_
#define EXTENDEDBIDOMAINASSEMBLER_HPP_

#include "UblasIncludes.hpp"

//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "AbstractCardiacFeVolumeIntegralAssembler.hpp"

#include "ExtendedBidomainTissue.hpp"
#include "HeartConfig.hpp"
#include "Element.hpp"
#include "BoundaryElement.hpp"
#include "ChastePoint.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"

/**
 * This is a derived class of AbstractCardiacFeVolumeIntegralAssembler and takes care of
 * assembling LHS for solution of the extended bidomain euqations
 *
 * See Buist ML, Poh YC.
 * An Extended Bidomain Framework Incorporating Multiple Cell Types.
 * Biophysical Journal, Volume 99, Issue 1, 13-18, 7 July 2010.
 *
 * Major differences with bidomain assembler:
 *
 * - there are two cells (instead of one)
 * - extracellular stimulus is present
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ExtendedBidomainAssembler
    : public AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,3,true,true,NORMAL>
{

protected:

    /** The tissue for which we assemble the matrix */
    ExtendedBidomainTissue<SPACE_DIM>* mpExtendedBidomainTissue;

    /** Local cache of the configuration singleton instance*/
    HeartConfig* mpConfig;


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
    virtual c_matrix<double,3*(ELEMENT_DIM+1),3*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,3> &rU,
        c_matrix<double, 3, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);


public:

    /**
     * Constructor stores the mesh and pde and sets up boundary conditions.
     *
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    ExtendedBidomainAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                            ExtendedBidomainTissue<SPACE_DIM>* pTissue,
                            unsigned numQuadPoints = 2);

    /**
     * Destructor.
     */
    ~ExtendedBidomainAssembler();

};


#endif /*EXTENDEDBIDOMAINASSEMBLER_HPP_*/
