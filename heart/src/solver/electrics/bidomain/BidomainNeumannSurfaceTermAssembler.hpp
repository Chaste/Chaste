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

#ifndef BIDOMAINNEUMANNSURFACETERMASSEMBLER_HPP_
#define BIDOMAINNEUMANNSURFACETERMASSEMBLER_HPP_

#include "AbstractFeSurfaceIntegralAssembler.hpp"



/**
 *  Assembler which sets up the surface integral integrals for the bidomain equations, assuming
 *  that the boundary conditions are written:  div(sigma_i grad phi_i) . n = g1  and
 *  div(sigma_e grad phi_e) dot n = g2.
 *
 *  These are not 'natural' boundary conditions for the para-elliptic bidomain equations (natural BCs for the second
 *
 *  Hence we don't use the NaturalNeumannSurfaceTermAssembler and have a special class here. It means that
 *  any BCs specified for bidomain and put in a BoundaryConditionsContainer should be for
 *  div(sigma_i grad phi_i) . n and div(sigma_e grad phi_e) . n.
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainNeumannSurfaceTermAssembler : public AbstractFeSurfaceIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2>
{
protected:
    /**
     * This method returns the vector to be added to full vector
     * for a given Gauss point in BoundaryElement, ie, essentially the
     * INTEGRAND in the boundary integral part of the definition of the vector.
     * The arguments are the bases, x and current solution computed at the
     * Gauss point.
     *
     * @param rSurfaceElement the element which is being considered.
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rX The point in space
     */
    virtual c_vector<double, 2*ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
        c_vector<double, ELEMENT_DIM>& rPhi,
        ChastePoint<SPACE_DIM>& rX);

public:
    /**
     * Constructor
     *
     * @param pMesh The mesh
     * @param pBoundaryConditions The boundary conditions container
     * @param numQuadPoints Number of quad points (per dimension) to use
     */
    BidomainNeumannSurfaceTermAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBoundaryConditions,
                                        unsigned numQuadPoints = 2)
        : AbstractFeSurfaceIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2>(pMesh, pBoundaryConditions, numQuadPoints)
    {
    }
};



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 2*ELEMENT_DIM> BidomainNeumannSurfaceTermAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
        c_vector<double, ELEMENT_DIM>& rPhi,
        ChastePoint<SPACE_DIM>& rX)
{
    // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
    double sigma_i_times_grad_phi_i_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 0);
    double sigma_e_times_grad_phi_e_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 1);

    c_vector<double, 2*ELEMENT_DIM> ret;
    for (unsigned i=0; i<ELEMENT_DIM; i++)
    {
        ret(2*i)   = rPhi(i)*sigma_i_times_grad_phi_i_dot_n;
        ret(2*i+1) = rPhi(i)*(sigma_i_times_grad_phi_i_dot_n + sigma_e_times_grad_phi_e_dot_n);
    }

    return ret;
}

#endif // BIDOMAINNEUMANNSURFACETERMASSEMBLER_HPP_
