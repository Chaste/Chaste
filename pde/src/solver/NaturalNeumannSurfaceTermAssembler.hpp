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

#ifndef NATURALNEUMANNSURFACETERMASSEMBLER_HPP
#define NATURALNEUMANNSURFACETERMASSEMBLER_HPP

#include "AbstractFeSurfaceIntegralAssembler.hpp"


/**
 *   An assembler for assembling surface element contributions to a vector coming from
 *   Neumann boundary conditions, assuming the prescribed BCs are NATURAL boundary conditions.
 *
 *   Examples of natural BCs:
 *   u_t = u_{xx} +f                ---> natural BCs are specification of (u_x . n) = g
 *   u_t = Laplacian(u) +f          ---> natural BCs are specification of (grad(u) . n) = g
 *   u_t = Div (D grad(u))          ---> natural BCs are specification of ( (D grad(u)) . n) = g
 *   0   = Div (D grad(u))          ---> natural BCs are specification of ( (D grad(u)) . n) = g
 *   0   = Div (D1 grad(u1)) + Div (D2 grad(u2)) (one equation of a 2-unknown problem)
 *        ---> natural BCs are specification of ( (D1 grad(u1)) . n + (D2 grad(u2)) . n) = g1
 *
 *   In all cases, when in weak form, the surface integral is integral(gv dS), where v is the test function
 *   (or integral_over_Gamma(g1*v1+g2*v2 dS) in the 2-unknown case).
 *
 *   The contribution of the finite element vector is therefore always: the STRIPED version of the following
 *   BLOCK vector: [c1 c2 .. cP], where P = PROBLEM_DIM, each c_i is of dimension N (number of nodes/bases), and
 *   c1 has entries `integral_over_Gamma (g1*phi_i dS)`, etc.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class NaturalNeumannSurfaceTermAssembler : public AbstractFeSurfaceIntegralAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
protected:
    /** Scale factor to multiply the integrals */
    double mScaleFactor;


    /**
     * This method returns the vector to be added to full vector
     * for a given Gauss point in BoundaryElement, ie, essentially the
     * INTEGRAND in the boundary integral part of the definition of the vector.
     * The arguments are the bases, x and current solution computed at the
     * Gauss point.
     *
     * Returns the integrand corresponding to natural Neumann BCs being
     * specified.
     *
     * @param rSurfaceElement the element which is being considered.
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rX The point in space
     */
    virtual c_vector<double, PROBLEM_DIM*ELEMENT_DIM> ComputeVectorSurfaceTerm(
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
    NaturalNeumannSurfaceTermAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                       BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* pBoundaryConditions,
                                       unsigned numQuadPoints = 2)
        : AbstractFeSurfaceIntegralAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(pMesh, pBoundaryConditions, numQuadPoints),
          mScaleFactor(1.0)
    {
    }

    /**
     * Set a scale factor to multiply the contribution to the vector that is assembled by this class
     * @param scaleFactor
     */
    void SetScaleFactor(double scaleFactor)
    {
        mScaleFactor = scaleFactor;
    }
};



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
c_vector<double, PROBLEM_DIM*ELEMENT_DIM> NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
        c_vector<double, ELEMENT_DIM>& rPhi,
        ChastePoint<SPACE_DIM>& rX)
{
    c_vector<double, ELEMENT_DIM*PROBLEM_DIM> ret;
    c_vector<double, PROBLEM_DIM> neumann_bc_values;

    for (unsigned i=0; i<ELEMENT_DIM; i++)
    {
        for(unsigned problem_dim = 0; problem_dim<PROBLEM_DIM; problem_dim++)
        {
            double neumann_bc_value = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, problem_dim);

            ret(PROBLEM_DIM*i + problem_dim)  =  mScaleFactor * rPhi(i) * neumann_bc_value;
        }
    }

    return ret;
}

#endif /* NATURALNEUMANNSURFACETERMASSEMBLER_HPP */
