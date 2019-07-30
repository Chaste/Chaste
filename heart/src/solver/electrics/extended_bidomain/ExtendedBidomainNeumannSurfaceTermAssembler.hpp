/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
     * @return stencil vector
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
     */
    ExtendedBidomainNeumannSurfaceTermAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,3>* pBoundaryConditions)
        : AbstractFeSurfaceIntegralAssembler<ELEMENT_DIM,SPACE_DIM,3>(pMesh, pBoundaryConditions)
    {
    }
};


// LCOV_EXCL_START //no non-zero Neumann BC allowed at the moment in extended bidomain problems
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

//    for (unsigned i=0; i<3*ELEMENT_DIM; i = i + 3)
//    {
//        ret(i) = rPhi(i)*sigma_i_times_grad_phi_i_first_cell_dot_n;
//        ret(i+1) = rPhi(i)*sigma_i_times_grad_phi_i_second_cell_dot_n;
//        ret(i+2) = rPhi(i)*(sigma_i_times_grad_phi_i_first_cell_dot_n + sigma_i_times_grad_phi_i_second_cell_dot_n + sigma_e_times_grad_phi_e_dot_n);
//    }
    for (unsigned i=0; i<ELEMENT_DIM; i++)
    {
        ret(3*i) = rPhi(i)*sigma_i_times_grad_phi_i_first_cell_dot_n;
        ret(3*i+1) = rPhi(i)*sigma_i_times_grad_phi_i_second_cell_dot_n;
        ret(3*i+2) = rPhi(i)*(sigma_i_times_grad_phi_i_first_cell_dot_n + sigma_i_times_grad_phi_i_second_cell_dot_n + sigma_e_times_grad_phi_e_dot_n);
    }
    return ret;
}
// LCOV_EXCL_STOP

#endif // EXTENDEDBIDOMAINNEUMANNSURFACETERMASSEMBLER_HPP_
