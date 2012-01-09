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

#ifndef EXTENDEDBIDOMAINNEUMANNSURFACETERMASSEMBLER_HPP_
#define EXTENDEDBIDOMAINNEUMANNSURFACETERMASSEMBLER_HPP_

#include "AbstractFeSurfaceIntegralAssembler.hpp"


/**
 *  Assembler which sets up the surface integral integrals for the extended bidomain equations, assuming
 *  that the boundary conditions are written:  div(sigma_i_1 grad phi_i_1) . n = g1,
 *  div(sigma_i_2 grad phi_i_2) . n = g2   and
 *  div(sigma_e grad phi_e) dot n = g3.
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ExtendedBidomainNeumannSurfaceTermAssembler : public AbstractFeSurfaceIntegralAssembler<ELEMENT_DIM,SPACE_DIM,3>
{
protected:
    /**
     * ComputeVectorSurfaceTerm()
     *
     * This method is called by AssembleOnSurfaceElement() and tells the
     * assembler what to add to the element stiffness matrix arising
     * from surface element contributions.
     *
     * NOTE: this method has to be implemented but shouldn't ever be called -
     * because all bidomain problems (currently) just have zero Neumann boundary
     * conditions and the AbstractLinearAssmebler::AssembleSystem() method
     * will realise this and not loop over surface elements.
     *
     * @param rSurfaceElement the element which is being considered.
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rX The point in space
     */

    virtual c_vector<double, 3*ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
        c_vector<double,ELEMENT_DIM> &rPhi,
        ChastePoint<SPACE_DIM> &rX);

public:
    /**
     * Constructor
     *
     * @param pMesh The mesh
     * @param pBoundaryConditions The boundary conditions container
     * @param numQuadPoints Number of quad points (per dimension) to use
     */
    ExtendedBidomainNeumannSurfaceTermAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,3>* pBoundaryConditions,
                                        unsigned numQuadPoints = 2)
        : AbstractFeSurfaceIntegralAssembler<ELEMENT_DIM,SPACE_DIM,3>(pMesh, pBoundaryConditions, numQuadPoints)
    {
    }
};


#define COVERAGE_IGNORE //no non-zero Neumann BC allowed at the moment in extended bidomain problems
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3*ELEMENT_DIM> ExtendedBidomainNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorSurfaceTerm(
    const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
    c_vector<double,ELEMENT_DIM> &rPhi,
    ChastePoint<SPACE_DIM> &rX)
{
    // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
    double sigma_i_times_grad_phi_i_first_cell_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 0);
    double sigma_i_times_grad_phi_i_second_cell_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 1);
    double sigma_e_times_grad_phi_e_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 2);

    c_vector<double, 3*ELEMENT_DIM> ret;

    for (unsigned i=0; i<3*ELEMENT_DIM; i = i + 3)
    {
        ret(i) = rPhi(i)*sigma_i_times_grad_phi_i_first_cell_dot_n;
        ret(i+1) = rPhi(i)*sigma_i_times_grad_phi_i_second_cell_dot_n;
        ret(i+2) = rPhi(i)*(sigma_i_times_grad_phi_i_first_cell_dot_n + sigma_i_times_grad_phi_i_second_cell_dot_n + sigma_e_times_grad_phi_e_dot_n);
    }
    return ret;
}
#undef COVERAGE_IGNORE

#endif // EXTENDEDBIDOMAINNEUMANNSURFACETERMASSEMBLER_HPP_
