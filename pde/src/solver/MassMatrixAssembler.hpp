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

#ifndef MASSMATRIXASSEMBLER_HPP_
#define MASSMATRIXASSEMBLER_HPP_

#include "AbstractFeVolumeIntegralAssembler.hpp"
#include "HeartConfig.hpp"

/**
 * Simple implementation of AbstractFeVolumeIntegralAssembler which provides mass matrices
 * for a given mesh, multiplied by a scale factor if required. In other words, the matrix
 *
 * M_{ij} = k integral_{domain}  phi_i(x) phi_j(x) dV
 *
 * where phi_i is the i-th (linear) basis function and k the scale factor (constant
 * throughout the mesh).
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MassMatrixAssembler
    : public AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, 1, false /*no vectors*/, true/*assembles matrices*/, NORMAL>
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
     */
    c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)>
        ComputeMatrixTerm(
                c_vector<double, ELEMENT_DIM+1> &rPhi,
                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                ChastePoint<SPACE_DIM> &rX,
                c_vector<double,1> &rU,
                c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
                Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> mass_matrix = outer_prod(rPhi, rPhi);

        if (mUseMassLumping)
        {
            for (unsigned row=0; row<ELEMENT_DIM+1; row++)
            {
                for (unsigned column=0; column<ELEMENT_DIM+1; column++)
                {
                    if (row != column)
                    {
                        mass_matrix(row,row) += mass_matrix(row,column);
                        mass_matrix(row,column) = 0.0;
                    }
                }
            }
        }

        return mScaleFactor*mass_matrix;
    }

    /**
     * Constructor.
     *
     * @param pMesh the mesh
     * @param scaleFactor the factor with which the multiply the mass matrix. Defaults to 1.0
     * @param useMassLumping whether to use mass matrix lumping or not
     */
    MassMatrixAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh, bool useMassLumping=false, double scaleFactor=1.0)
        : AbstractFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,1,false,true,NORMAL>(pMesh),
          mScaleFactor(scaleFactor),
          mUseMassLumping(useMassLumping)
    {
    }
};

#endif /*MASSMATRIXASSEMBLER_HPP_*/
