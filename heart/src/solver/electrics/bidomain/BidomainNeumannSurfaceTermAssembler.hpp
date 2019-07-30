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
     * @return stencil vector
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
     */
    BidomainNeumannSurfaceTermAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBoundaryConditions)
        : AbstractFeSurfaceIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2>(pMesh, pBoundaryConditions)
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
