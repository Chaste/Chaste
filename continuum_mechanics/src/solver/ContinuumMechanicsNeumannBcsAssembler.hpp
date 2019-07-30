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

#ifndef CONTINUUMMECHANICSNEUMANNBCSASSEMBLER_HPP_
#define CONTINUUMMECHANICSNEUMANNBCSASSEMBLER_HPP_

#include "AbstractFeAssemblerInterface.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "QuadraticMesh.hpp"
#include "DistributedQuadraticMesh.hpp"
#include "LinearBasisFunction.hpp"
#include "QuadraticBasisFunction.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "PetscTools.hpp"
#include "PetscVecTools.hpp"
#include "PetscMatTools.hpp"
#include "GaussianQuadratureRule.hpp"
#include "ContinuumMechanicsProblemDefinition.hpp"


/**
 *  Abstract class for assembling surface-integral parts of vectors in continuum
 *  mechanics problems.
 *
 *  For such problems, the matrix has the form
 *  [A     B1]
 *  [B2^T  C ]
 *  (where often B1=B2 and C=0) and the vector has the form
 *  [b1]
 *  [b2]
 *
 *  This class adds surface integral components, arising from natural Neumann boundary conditions,
 *  to the spatial part of the RHS vector, ie to b1.
 *
 */
template<unsigned DIM>
class ContinuumMechanicsNeumannBcsAssembler : public AbstractFeAssemblerInterface<true,false>
{
    /** Number of vertices per (boundary) element. */
    static const unsigned NUM_VERTICES_PER_ELEMENT = DIM;

    /** Number of nodes per (boundary) element. */
    static const unsigned NUM_NODES_PER_ELEMENT = DIM*(DIM+1)/2; // assuming quadratic

    /** Size of the spatial-block (the number or rows or columns in the submatrix A), restricted to one element */
    static const unsigned SPATIAL_BLOCK_SIZE_ELEMENTAL = DIM*NUM_NODES_PER_ELEMENT;
    /** Size of the pressure-block (the number of rows or columns in the submatrix C), restricted to one element */
    static const unsigned PRESSURE_BLOCK_SIZE_ELEMENTAL = NUM_VERTICES_PER_ELEMENT;

    /** Stencil size. */
    static const unsigned STENCIL_SIZE = DIM*NUM_NODES_PER_ELEMENT + NUM_VERTICES_PER_ELEMENT;

protected:
    /** The quadratic mesh */
    AbstractTetrahedralMesh<DIM,DIM>* mpMesh;

    /** Problem definition, containing the boundary conditions */
    ContinuumMechanicsProblemDefinition<DIM>* mpProblemDefinition;

    /** Quadrature rule for surface integrals */
    GaussianQuadratureRule<DIM-1>* mpQuadRule;

    /**
     * The main assembly method. Protected, should only be called through Assemble(),
     * or AssembleVector() in the parent class.
     */
    void DoAssemble();


    /**
     * Calculate the contribution of a single element to the linear system.
     *
     * @param rElement The element to assemble on.
     * @param rBElem The element's contribution to the RHS vector is returned in this
     *    vector of length n, the no. of nodes in this element. There is no
     *    need to zero this vector before calling.
     * @param boundaryConditionIndex into the traction BCs structures in the problem
     *    definition
     */
    void AssembleOnBoundaryElement(BoundaryElement<DIM-1,DIM>& rElement,
                                   c_vector<double, STENCIL_SIZE>& rBElem,
                                   unsigned boundaryConditionIndex);

public:
    /** Constructor
     *  @param pMesh Pointer to the mesh
     *  @param pProblemDefinition Pointer to the problem definition object
     */
    ContinuumMechanicsNeumannBcsAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh,
                                          ContinuumMechanicsProblemDefinition<DIM>* pProblemDefinition)
        : AbstractFeAssemblerInterface<true,false>(),
          mpMesh(pMesh),
          mpProblemDefinition(pProblemDefinition)
    {
        assert(pMesh);
        assert(pProblemDefinition);

        //Check that the mesh is Quadratic
        QuadraticMesh<DIM>* p_quad_mesh = dynamic_cast<QuadraticMesh<DIM>* >(pMesh);
        DistributedQuadraticMesh<DIM>* p_distributed_quad_mesh = dynamic_cast<DistributedQuadraticMesh<DIM>* >(pMesh);

        if ((p_quad_mesh == NULL) && (p_distributed_quad_mesh == NULL))
        {
            EXCEPTION("Continuum mechanics solvers require a quadratic mesh");
        }
        // In general a mechanics problem is non-polynomial.
        // We therefore use the highest order integration rule available.
        mpQuadRule = new GaussianQuadratureRule<DIM-1>(3);
    }

    /**
     * Destructor.
     */
    virtual ~ContinuumMechanicsNeumannBcsAssembler()
    {
        delete mpQuadRule;
    }
};


template<unsigned DIM>
void ContinuumMechanicsNeumannBcsAssembler<DIM>::DoAssemble()
{
    if (this->mVectorToAssemble==NULL)
    {
        EXCEPTION("Vector to be assembled has not been set");
    }

    if (PetscVecTools::GetSize(this->mVectorToAssemble) != (DIM+1)*mpMesh->GetNumNodes() )
    {
        EXCEPTION("Vector provided to be assembled has size " << PetscVecTools::GetSize(this->mVectorToAssemble) << ", not expected size of " << (DIM+1)*mpMesh->GetNumNodes() << " ((dim+1)*num_nodes)");
    }

    // Zero the matrix/vector if it is to be assembled
    if (this->mZeroVectorBeforeAssembly)
    {
        PetscVecTools::Zero(this->mVectorToAssemble);
    }


    if (mpProblemDefinition->GetTractionBoundaryConditionType() != NO_TRACTIONS)
    {
        c_vector<double, STENCIL_SIZE> b_elem = zero_vector<double>(STENCIL_SIZE);

        for (unsigned bc_index=0; bc_index<mpProblemDefinition->rGetTractionBoundaryElements().size(); bc_index++)
        {
            BoundaryElement<DIM-1,DIM>& r_boundary_element = *(mpProblemDefinition->rGetTractionBoundaryElements()[bc_index]);
            AssembleOnBoundaryElement(r_boundary_element, b_elem, bc_index);

            unsigned p_indices[STENCIL_SIZE];
            for (unsigned i=0; i<NUM_NODES_PER_ELEMENT; i++)
            {
                for (unsigned j=0; j<DIM; j++)
                {
                    p_indices[DIM*i+j] = (DIM+1)*r_boundary_element.GetNodeGlobalIndex(i) + j;
                }
            }
            // Note: The pressure block of b_elem will be zero, but this bit still needs to be
            // set to avoid memory leaks.
            for (unsigned i=0; i<DIM /*vertices per boundary elem */; i++)
            {
                p_indices[DIM*NUM_NODES_PER_ELEMENT + i] = (DIM+1)*r_boundary_element.GetNodeGlobalIndex(i)+DIM;
            }

            PetscVecTools::AddMultipleValues<STENCIL_SIZE>(this->mVectorToAssemble, p_indices, b_elem);
        }
    }
}

template<unsigned DIM>
void ContinuumMechanicsNeumannBcsAssembler<DIM>::AssembleOnBoundaryElement(BoundaryElement<DIM-1,DIM>& rBoundaryElement,
                                                                           c_vector<double,STENCIL_SIZE>& rBelem,
                                                                           unsigned boundaryConditionIndex)
{
    rBelem.clear();

    c_vector<double, DIM> weighted_direction;
    double jacobian_determinant;
    mpMesh->GetWeightedDirectionForBoundaryElement(rBoundaryElement.GetIndex(), weighted_direction, jacobian_determinant);

    c_vector<double,NUM_NODES_PER_ELEMENT> phi;

    for (unsigned quad_index=0; quad_index<mpQuadRule->GetNumQuadPoints(); quad_index++)
    {
        double wJ = jacobian_determinant * mpQuadRule->GetWeight(quad_index);
        const ChastePoint<DIM-1>& quad_point = mpQuadRule->rGetQuadPoint(quad_index);
        QuadraticBasisFunction<DIM-1>::ComputeBasisFunctions(quad_point, phi);

        c_vector<double,DIM> traction = zero_vector<double>(DIM);
        switch (mpProblemDefinition->GetTractionBoundaryConditionType())
        {
            case ELEMENTWISE_TRACTION:
            {
                traction = mpProblemDefinition->rGetElementwiseTractions()[boundaryConditionIndex];
                break;
            }
            default:
                // Functional traction not implemented yet..
                NEVER_REACHED;
        }

        for (unsigned index=0; index<NUM_NODES_PER_ELEMENT*DIM; index++)
        {
            unsigned spatial_dim = index%DIM;
            unsigned node_index = (index-spatial_dim)/DIM;

            assert(node_index < NUM_NODES_PER_ELEMENT);

            rBelem(index) += traction(spatial_dim) * phi(node_index) * wJ;
        }
    }
}


#endif // CONTINUUMMECHANICSNEUMANNBCSASSEMBLER_HPP_
