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

#ifndef ABSTRACTFESURFACENTEGRALASSEMBLER_HPP_
#define ABSTRACTFESURFACENTEGRALASSEMBLER_HPP_

#include "AbstractFeAssemblerCommon.hpp"
#include "GaussianQuadratureRule.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "PetscVecTools.hpp"
#include "PetscMatTools.hpp"


/**
 *  Similar to AbstractFeVolumeIntegralAssembler but is used for constructing finite element objects
 *  that are based on SURFACE INTEGRALS, as opposed to volume integrals.
 *
 *  This class assumes that the concrete class only needs to assemble a vector, not a matrix.
 *  (Can be extended in the future if needed).
 *
 *  Hence, the (effectively) pure method, that needs to be implemented, is ComputeVectorSurfaceTerm().
 *
 *  The surface terms is assumed to come from Neumann BCs, so only the surface elements containing
 *  non-zero Neumann BCs (from the BoundaryConditionsContainer given) are assembled on.
 *
 *  The interface is the same the volume assemblers.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractFeSurfaceIntegralAssembler : public AbstractFeAssemblerCommon<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,true,false,NORMAL>
{
protected:
    /** Mesh to be solved on. */
    AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* mpMesh;

    /** Boundary conditions container */
    BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpBoundaryConditions;

    /** Quadrature rule for use on boundary elements. */
    GaussianQuadratureRule<ELEMENT_DIM-1>* mpSurfaceQuadRule;

    /** Basis function for use with boundary elements. */
    typedef LinearBasisFunction<ELEMENT_DIM-1> SurfaceBasisFunction;

    /**
     * This method returns the vector to be added to full vector
     * for a given Gauss point in BoundaryElement, ie, essentially the
     * INTEGRAND in the boundary integral part of the definition of the vector.
     * The arguments are the bases, x and current solution computed at the
     * Gauss point.
     *
     *  ** This method needs to be overloaded in the concrete class **
     *
     * @param rSurfaceElement the element which is being considered.
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rX The point in space
     */
    virtual c_vector<double, PROBLEM_DIM*ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
        c_vector<double, ELEMENT_DIM>& rPhi,
        ChastePoint<SPACE_DIM>& rX)
    {
        // If this line is reached this means this method probably hasn't been over-ridden correctly in
        // the concrete class
        NEVER_REACHED;
        return zero_vector<double>(ELEMENT_DIM*PROBLEM_DIM);
    }

    /**
     * Calculate the contribution of a single surface element with Neumann
     * boundary condition to the linear system.
     *
     * @param rSurfaceElement The element to assemble on.
     * @param rBSurfElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     */
    virtual void AssembleOnSurfaceElement(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
                                          c_vector<double, PROBLEM_DIM*ELEMENT_DIM>& rBSurfElem);


    /**
     * Main assemble method. Users should call Assemble() however
     */
    void DoAssemble();


public:
    /**
     * Constructor
     *
     * @param pMesh The mesh
     * @param pBoundaryConditions The boundary conditions container
     * @param numQuadPoints Number of quad points (per dimension) to use
     */
    AbstractFeSurfaceIntegralAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                       BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* pBoundaryConditions,
                                       unsigned numQuadPoints = 2);

    /**
     * Destructor
     */
    ~AbstractFeSurfaceIntegralAssembler();

    /**
     * Reset the internal boundary conditions container pointer
     * @param pBoundaryConditions
     */
    void ResetBoundaryConditionsContainer(BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* pBoundaryConditions)
    {
        assert(pBoundaryConditions);
        this->mpBoundaryConditions = pBoundaryConditions;
    }
};


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractFeSurfaceIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AbstractFeSurfaceIntegralAssembler(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* pBoundaryConditions,
            unsigned numQuadPoints)
    : AbstractFeAssemblerCommon<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,true,false,NORMAL>(),
      mpMesh(pMesh),
      mpBoundaryConditions(pBoundaryConditions)
{
    assert(pMesh);
    assert(pBoundaryConditions);
    assert(numQuadPoints > 0);

    mpSurfaceQuadRule = new GaussianQuadratureRule<ELEMENT_DIM-1>(numQuadPoints);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractFeSurfaceIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::~AbstractFeSurfaceIntegralAssembler()
{
    delete mpSurfaceQuadRule;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractFeSurfaceIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::DoAssemble()
{
    assert(this->mAssembleVector);

    HeartEventHandler::BeginEvent(HeartEventHandler::NEUMANN_BCS);

    // Loop over surface elements with non-zero Neumann boundary conditions
    if (mpBoundaryConditions->AnyNonZeroNeumannConditions())
    {
        typename BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::NeumannMapIterator
            neumann_iterator = mpBoundaryConditions->BeginNeumann();

        c_vector<double, PROBLEM_DIM*ELEMENT_DIM> b_surf_elem;

        // Iterate over defined conditions
        while (neumann_iterator != mpBoundaryConditions->EndNeumann())
        {
            const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& r_surf_element = *(neumann_iterator->first);
            AssembleOnSurfaceElement(r_surf_element, b_surf_elem);

            const size_t STENCIL_SIZE=PROBLEM_DIM*ELEMENT_DIM; // problem_dim*num_nodes_on_surface_element
            unsigned p_indices[STENCIL_SIZE];
            r_surf_element.GetStiffnessMatrixGlobalIndices(PROBLEM_DIM, p_indices);
            PetscVecTools::AddMultipleValues<STENCIL_SIZE>(this->mVectorToAssemble, p_indices, b_surf_elem);
            ++neumann_iterator;
        }
    }

    HeartEventHandler::EndEvent(HeartEventHandler::NEUMANN_BCS);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractFeSurfaceIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AssembleOnSurfaceElement(
            const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
            c_vector<double, PROBLEM_DIM*ELEMENT_DIM>& rBSurfElem)
{
    c_vector<double, SPACE_DIM> weighted_direction;
    double jacobian_determinant;
    mpMesh->GetWeightedDirectionForBoundaryElement(rSurfaceElement.GetIndex(), weighted_direction, jacobian_determinant);

    rBSurfElem.clear();

    // Allocate memory for the basis function values
    c_vector<double, ELEMENT_DIM>  phi;

    // Loop over Gauss points
    for (unsigned quad_index=0; quad_index<mpSurfaceQuadRule->GetNumQuadPoints(); quad_index++)
    {
        const ChastePoint<ELEMENT_DIM-1>& quad_point = mpSurfaceQuadRule->rGetQuadPoint(quad_index);

        SurfaceBasisFunction::ComputeBasisFunctions(quad_point, phi);

        //////////////////////////////
        // Interpolation: X only
        //////////////////////////////

        // The location of the Gauss point in the original element will be stored in x
        ChastePoint<SPACE_DIM> x(0,0,0);

        this->ResetInterpolatedQuantities();
        for (unsigned i=0; i<rSurfaceElement.GetNumNodes(); i++)
        {
            const c_vector<double, SPACE_DIM> node_loc = rSurfaceElement.GetNode(i)->rGetLocation();
            x.rGetLocation() += phi(i)*node_loc;

            // Allow the concrete version of the assembler to interpolate any desired quantities
            this->IncrementInterpolatedQuantities(phi(i), rSurfaceElement.GetNode(i));
        }


        //////////////////////////////////
        // Create elemental contribution
        //////////////////////////////////

        double wJ = jacobian_determinant * mpSurfaceQuadRule->GetWeight(quad_index);
        ///\todo #1321 Improve efficiency of Neumann BC implementation
        noalias(rBSurfElem) += ComputeVectorSurfaceTerm(rSurfaceElement, phi, x) * wJ;
    }
};

#endif // ABSTRACTFESURFACENTEGRALASSEMBLER_HPP_
