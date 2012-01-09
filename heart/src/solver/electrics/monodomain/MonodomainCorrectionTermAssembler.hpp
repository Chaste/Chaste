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


#ifndef MONODOMAINCORRECTIONTERM_HPP_
#define MONODOMAINCORRECTIONTERM_HPP_

#include "AbstractCorrectionTermAssembler.hpp"
#include "MonodomainTissue.hpp"


/**
 * An assembler which computes the correction term to add to the
 * RHS vector if using state variable interpolation (SVI), as well as determining
 * which elements should be corrected on. The formula to determine which
 * elements SVI is used is delta Iionic > TOL, where delta Iionic is the max
 * difference between nodal ionic values, and TOL is chosen conservatively
 * to be 1uA/cm^2^. See wiki page ChasteGuides/StateVariableInterpolation
 * for more details.
 */
template<unsigned ELEM_DIM,unsigned SPACE_DIM>
class MonodomainCorrectionTermAssembler
    : public AbstractCorrectionTermAssembler<ELEM_DIM,SPACE_DIM,1>
{
protected:
    /** Local cache of the configuration singleton instance */
    HeartConfig* mpConfig;

    /**
     * This method is called by AssembleOnElement and tells the assembler
     * the contribution to add to the element stiffness vector.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, rU(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    c_vector<double,1*(ELEM_DIM+1)> ComputeVectorTerm(
                c_vector<double, ELEM_DIM+1> &rPhi,
                c_matrix<double, SPACE_DIM, ELEM_DIM+1> &rGradPhi /* not used */,
                ChastePoint<SPACE_DIM> &rX /* not used */,
                c_vector<double,1> &rU,
                c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
                Element<ELEM_DIM,SPACE_DIM>* pElement);
public:

    /**
     * Constructor.
     *
     * @param pMesh  pointer to the mesh
     * @param pTissue  pointer to the cardiac tissue
     * @param numQuadPoints  number of quadrature points
     */
    MonodomainCorrectionTermAssembler(AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
                                      MonodomainTissue<ELEM_DIM,SPACE_DIM>* pTissue,
                                      unsigned numQuadPoints = 2);
};



#endif /*MONODOMAINCORRECTIONTERM_HPP_*/
