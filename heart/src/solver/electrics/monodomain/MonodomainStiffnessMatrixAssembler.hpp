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

#ifndef MONODOMAINSTIFFNESSMATRIXASSEMBLER_HPP_
#define MONODOMAINSTIFFNESSMATRIXASSEMBLER_HPP_


#include "AbstractCardiacFeVolumeIntegralAssembler.hpp"
#include "HeartConfig.hpp"
#include "AbstractCardiacTissue.hpp"

/**
 *  Implementation of AbstractFeVolumeIntegralAssembler which provides stiffness matrices
 *  required in monodomain problems:
 *
 *  K_{ij} = integral_{domain}  grad_phi_i(x)^T (sigma * grad_phi_j(x)) dV
 *
 *  where phi_i is the i-th (linear) basis function
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainStiffnessMatrixAssembler
   : public AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, 1, false /*no vectors*/, true/*assembles matrices*/, NORMAL>
{
public:
    /** Implemented ComputeMatrixTerm(), defined in AbstractFeVolumeIntegralAssembler. See
     *  documentation in that class.
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
        const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_i = this->mpCardiacTissue->rGetIntracellularConductivityTensor(pElement->GetIndex());

        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp = prod(sigma_i, rGradPhi);
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_grad_phi =
            prod(trans(rGradPhi), temp);

        return grad_phi_sigma_i_grad_phi;
    }

    /**
     *  Constructor
     *  @param pMesh the mesh
     *  @param pTissue  pointer to the tissue used for getting conductivity values
     */
    MonodomainStiffnessMatrixAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                             AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* pTissue)
        : AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,1,false,true,NORMAL>(pMesh,pTissue)
    {
    }
};

#endif /*MONODOMAINSTIFFNESSMATRIXASSEMBLER_HPP_*/
