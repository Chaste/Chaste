
/*

Copyright (C) University of Oxford, 2005-2011

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


#ifndef BIDOMAINWITHBATHASSEMBLER_HPP_
#define BIDOMAINWITHBATHASSEMBLER_HPP_

#include "BidomainAssembler.hpp"

/**
 *  Assembler for assembling the LHS matrix system solved in bidomain
 *  problems with a perfusing bath. See FEM implementations document for
 *  the exact definition of this.
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainWithBathAssembler : public BidomainAssembler<ELEMENT_DIM,SPACE_DIM>
{
protected:
    /**
     * Overloaded ComputeMatrixTerm() - calls the base class version for
     * tissue elements.
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
     * Constructor.
     *
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    BidomainWithBathAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                              BidomainTissue<SPACE_DIM>* pTissue,
                              unsigned numQuadPoints = 2)
      : BidomainAssembler<ELEMENT_DIM,SPACE_DIM>(pMesh,pTissue,numQuadPoints)
    {
    }

    /**
     * Destructor.
     */
    ~BidomainWithBathAssembler()
    {
    }
};

#endif /*BIDOMAINWITHBATHASSEMBLER_HPP_*/
