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

#ifndef ABSTRACTCONTINUUMMECHANICSASSEMBLER_HPP_
#define ABSTRACTCONTINUUMMECHANICSASSEMBLER_HPP_

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


/**
 *  Abstract class for assembling volume-integral parts of matrices and vectors in continuum
 *  mechanics problems.
 *
 *  For such problems, the matrix has, essentially, the form
 *  [A     B1]
 *  [B2^T  C ]
 *  (where often B1=B2 and C=0) and the vector has the form
 *  [b1]
 *  [b2]
 *  A is the spatial-spatial part, B1 the spatial-pressure part, etc.
 *
 *  Currently B1=B2 is assumed, this can be changed in the future.
 *
 *  This class works in the same way as the volume assembler in pde (AbstractFeVolumeIntegralAssembler),
 *  except the concrete class has to provide up to 6 methods, for each of the blocks A,B1,B2 and for b1
 *  and b2.
 *
 *  The assembler is main used for fluids - currently it is used for assembling the matrix for Stokes
 *  flow [A B ; B^T 0], and the preconditioner for Stokes flow [A, B; B^T M], and will be useful for
 *  further fluids assemblies, and could be used for mixed-problem incompressible linear elasticity.
 *  We DON'T use this assembler in the incompressible nonlinear elasticity solver, because even though
 *  the matrix is of the right form, things like stress and stress-derivative need to be computed at
 *  the AssembleOnElement level and the assembly is in general too complex for this class to be used.
 *
 *  NOTE: The elemental matrix and vector is as above. The full matrix and vector uses a completely
 *  different ordering: for parallelisation reasons the pressure variables are interleaved with the
 *  spatial variables and dummy pressure variables are used for internal nodes. For example, in 2d,
 *  the ordering is
 *  [U1 V1 P1 , .. , Un Vn, Pn]
 *  where n is the total number of nodes.
 *
 */
template<unsigned DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX>
class AbstractContinuumMechanicsAssembler : public AbstractFeAssemblerInterface<CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX>
{
protected:
    /** Whether the matrix is block symmetric (B1=B2). Currently fixed to true, in
     *  the future this may become a template.
     */
    static const bool BLOCK_SYMMETRIC_MATRIX = true; //generalise to non-block symmetric matrices later (when needed maybe)

    /** Number of vertices per element. */
    static const unsigned NUM_VERTICES_PER_ELEMENT = DIM+1;

    /** Number of nodes per element. */
    static const unsigned NUM_NODES_PER_ELEMENT = (DIM+1)*(DIM+2)/2; // assuming quadratic

    /** Size of the spatial-block (the number or rows or columns in the submatrix A), restricted to one element */
    static const unsigned SPATIAL_BLOCK_SIZE_ELEMENTAL = DIM*NUM_NODES_PER_ELEMENT;
    /** Size of the pressure-block (the number of rows or columns in the submatrix C), restricted to one element */
    static const unsigned PRESSURE_BLOCK_SIZE_ELEMENTAL = NUM_VERTICES_PER_ELEMENT;

    /** Stencil size. */
    static const unsigned STENCIL_SIZE = DIM*NUM_NODES_PER_ELEMENT + NUM_VERTICES_PER_ELEMENT;

    /** The quadratic mesh */
    AbstractTetrahedralMesh<DIM,DIM>* mpMesh;

    /** Quadrature rule for volume integrals */
    GaussianQuadratureRule<DIM>* mpQuadRule;

    /**
     * The main assembly method. Protected, should only be called through Assemble(),
     * AssembleMatrix() or AssembleVector() which set mAssembleMatrix, mAssembleVector
     * accordingly. Involves looping over elements, and computing
     * integrals and adding them to the vector or matrix
     */
    void DoAssemble();


    /**
     *  For a continuum mechanics problem in mixed form (displacement-pressure or velocity-pressure), the matrix
     *  has the form (except see comments about ordering above)
     *  [A     B1]
     *  [B2^T  C ]
     *  (where often B1=B2 and C=0). The function is related to the spatial-spatial block, ie matrix A.
     *
     *  For the contribution to A from a given element, this method should return the INTEGRAND in the definition of A.
     *  See concrete classes for examples. Needed to be implemented (overridden) if the concrete class
     *  is going to assemble matrices (ie if CAN_ASSEMBLE_MATRIX is true).
     *
     *  Default implementation returns a zero matrix - ie the block will be zero if this is not over-ridden
     *
     *  @param rQuadPhi  All the quadratic basis functions on this element, evaluated at the current quad point
     *  @param rGradQuadPhi  Gradients of all the quadratic basis functions on this element, evaluated at the current quad point
     *  @param rX Current location (physical position)
     *  @param pElement Current element
     *  @return stencil matrix
     */
    virtual c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,SPATIAL_BLOCK_SIZE_ELEMENTAL> ComputeSpatialSpatialMatrixTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        return zero_matrix<double>(SPATIAL_BLOCK_SIZE_ELEMENTAL,SPATIAL_BLOCK_SIZE_ELEMENTAL);
    }

    /**
     *  For a continuum mechanics problem in mixed form (displacement-pressure or velocity-pressure), the matrix
     *  has the form (except see comments about ordering above)
     *  [A     B1]
     *  [B2^T  C ]
     *  (where often B1=B2 and C=0). The function is related to the spatial-pressure block, ie matrix B1. If
     *  BLOCK_SYMMETRIC_MATRIX is true, B1=B2 is assumed, so it also relates to B2.
     *
     *  For the contribution to A from a given element, this method should return the INTEGRAND in the definition of B.
     *  See concrete classes for examples. Needed to be implemented (overridden) if the concrete class
     *  is going to assemble matrices (ie if CAN_ASSEMBLE_MATRIX is true).
     *
     *  Default implementation returns a zero matrix - ie the block will be zero if this is not over-ridden
     *
     *  @param rQuadPhi  All the quadratic basis functions on this element, evaluated at the current quad point
     *  @param rGradQuadPhi  Gradients of all the quadratic basis functions on this element, evaluated at the current quad point
     *  @param rLinearPhi  All the linear basis functions on this element, evaluated at the current quad point
     *  @param rGradLinearPhi  Gradients of all the linear basis functions on this element, evaluated at the current quad point
     *  @param rX Current location (physical position corresponding to quad point)
     *  @param pElement Current element
     *  @return stencil matrix
     */
    virtual c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputeSpatialPressureMatrixTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
        c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        return zero_matrix<double>(SPATIAL_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL);
    }


    /**
     *  For a continuum mechanics problem in mixed form (displacement-pressure or velocity-pressure), the matrix
     *  has the form (except see comments about ordering above)
     *  [A     B1]
     *  [B2^T  C ]
     *  (where often B1=B2 and C=0). The function is related to the pressure-pressure block, ie matrix C.
     *
     *  For the contribution to A from a given element, this method should return the INTEGRAND in the definition of C.
     *  See concrete classes for examples. Needed to be implemented (overridden) if the concrete class
     *  is going to assemble matrices (ie if CAN_ASSEMBLE_MATRIX is true).
     *
     *  Default implementation returns a zero matrix - ie the block will be zero if this is not over-ridden
     *
     *  @param rLinearPhi  All the linear basis functions on this element, evaluated at the current quad point
     *  @param rGradLinearPhi  Gradients of all the linear basis functions on this element, evaluated at the current quad point
     *  @param rX Current location (physical position)
     *  @param pElement Current element
     *  @return stencil matrix
     */
    virtual c_matrix<double,PRESSURE_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputePressurePressureMatrixTerm(
        c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
        c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        return zero_matrix<double>(PRESSURE_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL);
    }


    /**
     *  For a continuum mechanics problem in mixed form (displacement-pressure or velocity-pressure), the matrix
     *  has the form (except see comments about ordering above)
     *  [A     B1]
     *  [B2^T  C ]
     *  (where often B1=B2 and C=0) and the vector has the form
     *  [b1]
     *  [b2]
     *  The function is related to the spatial-block in the vector, ie b1.
     *
     *  For the contribution to b1 from a given element, this method should return the INTEGRAND in the definition of b1.
     *  See concrete classes for examples. Needed to be implemented (overridden) if the concrete class
     *  is going to assemble vectors (ie if CAN_ASSEMBLE_VECTOR is true).
     *
     *  No default implementation - this method must be over-ridden
     *
     *  @param rQuadPhi  All the quadratic basis functions on this element, evaluated at the current quad point
     *  @param rGradQuadPhi  Gradients of all the quadratic basis functions on this element, evaluated at the current quad point
     *  @param rX Current location (physical position)
     *  @param pElement Current element
     *  @return stencil vector
     */
    virtual c_vector<double,SPATIAL_BLOCK_SIZE_ELEMENTAL> ComputeSpatialVectorTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement) = 0;


    /**
     *  For a continuum mechanics problem in mixed form (displacement-pressure or velocity-pressure), the matrix
     *  has the form (except see comments about ordering above)
     *  [A     B1]
     *  [B2^T  C ]
     *  (where often B1=B2 and C=0) and the vector has the form
     *  [b1]
     *  [b2]
     *  The function is related to the pressure-block in the vector, ie b2.
     *
     *  For the contribution to b1 from a given element, this method should return the INTEGRAND in the definition of b2.
     *  See concrete classes for examples. Needed to be implemented (overridden) if the concrete class
     *  is going to assemble vectors (ie if CAN_ASSEMBLE_VECTOR is true).
     *
     *  Default implementation returns a zero vector - ie the block will be zero if this is not over-ridden
     *
     *  @param rLinearPhi  All the linear basis functions on this element, evaluated at the current quad point
     *  @param rGradLinearPhi  Gradients of all the linear basis functions on this element, evaluated at the current quad point
     *  @param rX Current location (physical position)
     *  @param pElement Current element
     *  @return stencil vector
     */
    virtual c_vector<double,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputePressureVectorTerm(
            c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
            c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
            c_vector<double,DIM>& rX,
            Element<DIM,DIM>* pElement)
    {
        return zero_vector<double>(PRESSURE_BLOCK_SIZE_ELEMENTAL);
    }

    /**
     * Calculate the contribution of a single element to the linear system.
     *
     * @param rElement The element to assemble on.
     * @param rAElem The element's contribution to the LHS matrix is returned in this
     *    n by n matrix, where n is the no. of nodes in this element. There is no
     *    need to zero this matrix before calling.
     * @param rBElem The element's contribution to the RHS vector is returned in this
     *    vector of length n, the no. of nodes in this element. There is no
     *    need to zero this vector before calling.
     */
    void AssembleOnElement(Element<DIM, DIM>& rElement,
                           c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElem,
                           c_vector<double, STENCIL_SIZE>& rBElem);

public:
    /** Constructor
     *  @param pMesh Pointer to the mesh
     */
    AbstractContinuumMechanicsAssembler(AbstractTetrahedralMesh<DIM, DIM>* pMesh)
        : AbstractFeAssemblerInterface<CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX>(),
          mpMesh(pMesh)
    {
        assert(pMesh);

        // Check that the mesh is quadratic
        QuadraticMesh<DIM>* p_quad_mesh = dynamic_cast<QuadraticMesh<DIM>* >(pMesh);
        DistributedQuadraticMesh<DIM>* p_distributed_quad_mesh = dynamic_cast<DistributedQuadraticMesh<DIM>* >(pMesh);

        if ((p_quad_mesh == NULL) && (p_distributed_quad_mesh == NULL))
        {
            EXCEPTION("Continuum mechanics assemblers require a quadratic mesh");
        }

        // In general the Jacobian for a mechanics problem is non-polynomial.
        // We therefore use the highest order integration rule available
        mpQuadRule = new GaussianQuadratureRule<DIM>(3);
    }

//    void SetCurrentSolution(Vec currentSolution);

    /**
     * Destructor.
     */
    virtual ~AbstractContinuumMechanicsAssembler()
    {
        delete mpQuadRule;
    }
};


//// add this method when needed..
//template<unsigned DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX>
//void AbstractContinuumMechanicsAssembler<DIM,CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX>::SetCurrentSolution(Vec currentSolution)
//{
//    assert(currentSolution != NULL);
//
//    // Replicate the current solution and store so can be used in AssembleOnElement
//    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
//    mCurrentSolutionOrGuessReplicated.ReplicatePetscVector(currentSolution);
//    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
//
//    // The AssembleOnElement type methods will determine if a current solution or
//    // current guess exists by looking at the size of the replicated vector, so
//    // check the size is zero if there isn't a current solution.
//    assert(mCurrentSolutionOrGuessReplicated.GetSize() > 0);
//}

template<unsigned DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX>
void AbstractContinuumMechanicsAssembler<DIM,CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX>::DoAssemble()
{
    assert(this->mAssembleMatrix || this->mAssembleVector);
    if (this->mAssembleMatrix)
    {
        if (this->mMatrixToAssemble == NULL)
        {
            EXCEPTION("Matrix to be assembled has not been set");
        }
        if (PetscMatTools::GetSize(this->mMatrixToAssemble) != (DIM+1)*mpMesh->GetNumNodes())
        {
            EXCEPTION("Matrix provided to be assembled has size " << PetscMatTools::GetSize(this->mMatrixToAssemble) << ", not expected size of " << (DIM+1)*mpMesh->GetNumNodes() << " ((dim+1)*num_nodes)");
        }
    }

    if (this->mAssembleVector)
    {
        if (this->mVectorToAssemble == NULL)
        {
            EXCEPTION("Vector to be assembled has not been set");
        }
        if (PetscVecTools::GetSize(this->mVectorToAssemble) != (DIM+1)*mpMesh->GetNumNodes())
        {
            EXCEPTION("Vector provided to be assembled has size " << PetscVecTools::GetSize(this->mVectorToAssemble) << ", not expected size of " << (DIM+1)*mpMesh->GetNumNodes() << " ((dim+1)*num_nodes)");
        }
    }

    // Zero the matrix/vector if it is to be assembled
    if (this->mAssembleVector && this->mZeroVectorBeforeAssembly)
    {
        // Note PetscVecTools::Finalise(this->mVectorToAssemble); on an unused matrix
        // would "compress" data and make any pre-allocated entries redundant.
        PetscVecTools::Zero(this->mVectorToAssemble);
    }
    if (this->mAssembleMatrix && this->mZeroMatrixBeforeAssembly)
    {
        // Note PetscMatTools::Finalise(this->mMatrixToAssemble); on an unused matrix
        // would "compress" data and make any pre-allocated entries redundant.
        PetscMatTools::Zero(this->mMatrixToAssemble);
    }

    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem = zero_matrix<double>(STENCIL_SIZE,STENCIL_SIZE);
    c_vector<double, STENCIL_SIZE> b_elem = zero_vector<double>(STENCIL_SIZE);


    // Loop over elements
    for (typename AbstractTetrahedralMesh<DIM, DIM>::ElementIterator iter = mpMesh->GetElementIteratorBegin();
         iter != mpMesh->GetElementIteratorEnd();
         ++iter)
    {
        Element<DIM, DIM>& r_element = *iter;

        // Test for ownership first, since it's pointless to test the criterion on something which we might know nothing about.
        if (r_element.GetOwnership() == true  /*&& ElementAssemblyCriterion(r_element)==true*/)
        {
            // LCOV_EXCL_START
            // note: if assemble matrix only
            if (CommandLineArguments::Instance()->OptionExists("-mech_very_verbose") && this->mAssembleMatrix)
            {
                std::cout << "\r[" << PetscTools::GetMyRank() << "]: Element " << r_element.GetIndex() << " of " << mpMesh->GetNumElements() << std::flush;
            }
            // LCOV_EXCL_STOP

            AssembleOnElement(r_element, a_elem, b_elem);

            // Note that a different ordering is used for the elemental matrix compared to the global matrix.
            // See comments about ordering above.
            unsigned p_indices[STENCIL_SIZE];
            // Work out the mapping for spatial terms
            for (unsigned i=0; i<NUM_NODES_PER_ELEMENT; i++)
            {
                for (unsigned j=0; j<DIM; j++)
                {
                    // DIM+1 on the right-hand side here is the problem dimension
                    p_indices[DIM*i+j] = (DIM+1)*r_element.GetNodeGlobalIndex(i) + j;
                }
            }
            // Work out the mapping for pressure terms
            for (unsigned i=0; i<NUM_VERTICES_PER_ELEMENT; i++)
            {
                p_indices[DIM*NUM_NODES_PER_ELEMENT + i] = (DIM+1)*r_element.GetNodeGlobalIndex(i)+DIM;
            }


            if (this->mAssembleMatrix)
            {
                PetscMatTools::AddMultipleValues<STENCIL_SIZE>(this->mMatrixToAssemble, p_indices, a_elem);
            }

            if (this->mAssembleVector)
            {
                PetscVecTools::AddMultipleValues<STENCIL_SIZE>(this->mVectorToAssemble, p_indices, b_elem);
            }
        }
    }
}

template<unsigned DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX>
void AbstractContinuumMechanicsAssembler<DIM,CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX>::AssembleOnElement(Element<DIM, DIM>& rElement,
                                                                                                         c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElem,
                                                                                                         c_vector<double, STENCIL_SIZE>& rBElem)
{
    static c_matrix<double,DIM,DIM> jacobian;
    static c_matrix<double,DIM,DIM> inverse_jacobian;
    double jacobian_determinant;

    mpMesh->GetInverseJacobianForElement(rElement.GetIndex(), jacobian, jacobian_determinant, inverse_jacobian);

    if (this->mAssembleMatrix)
    {
        rAElem.clear();
    }

    if (this->mAssembleVector)
    {
        rBElem.clear();
    }

    // Allocate memory for the basis functions values and derivative values
    static c_vector<double, NUM_VERTICES_PER_ELEMENT> linear_phi;
    static c_vector<double, NUM_NODES_PER_ELEMENT> quad_phi;
    static c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> grad_quad_phi;
    static c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT> grad_linear_phi;

    c_vector<double,DIM> body_force;

    // Loop over Gauss points
    for (unsigned quadrature_index=0; quadrature_index < mpQuadRule->GetNumQuadPoints(); quadrature_index++)
    {
        double wJ = jacobian_determinant * mpQuadRule->GetWeight(quadrature_index);
        const ChastePoint<DIM>& quadrature_point = mpQuadRule->rGetQuadPoint(quadrature_index);

        // Set up basis function info
        LinearBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, linear_phi);
        QuadraticBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, quad_phi);
        QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, inverse_jacobian, grad_quad_phi);
        LinearBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, inverse_jacobian, grad_linear_phi);

        // interpolate X (ie physical location of this quad point).
        c_vector<double,DIM> X = zero_vector<double>(DIM);
        for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
        {
            for (unsigned j=0; j<DIM; j++)
            {
                X(j) += linear_phi(vertex_index)*rElement.GetNode(vertex_index)->rGetLocation()(j);
            }
        }

        if (this->mAssembleVector)
        {
            c_vector<double,SPATIAL_BLOCK_SIZE_ELEMENTAL> b_spatial
                = ComputeSpatialVectorTerm(quad_phi, grad_quad_phi, X, &rElement);
            c_vector<double,PRESSURE_BLOCK_SIZE_ELEMENTAL> b_pressure = ComputePressureVectorTerm(linear_phi, grad_linear_phi, X, &rElement);

            for (unsigned i=0; i<SPATIAL_BLOCK_SIZE_ELEMENTAL; i++)
            {
                rBElem(i) += b_spatial(i)*wJ;
            }

            for (unsigned i=0; i<PRESSURE_BLOCK_SIZE_ELEMENTAL; i++)
            {
                rBElem(SPATIAL_BLOCK_SIZE_ELEMENTAL + i) += b_pressure(i)*wJ;
            }
        }

        if (this->mAssembleMatrix)
        {
            c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,SPATIAL_BLOCK_SIZE_ELEMENTAL> a_spatial_spatial
                = ComputeSpatialSpatialMatrixTerm(quad_phi, grad_quad_phi, X, &rElement);

            c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> a_spatial_pressure
                = ComputeSpatialPressureMatrixTerm(quad_phi, grad_quad_phi, linear_phi, grad_linear_phi, X, &rElement);

            c_matrix<double,PRESSURE_BLOCK_SIZE_ELEMENTAL,SPATIAL_BLOCK_SIZE_ELEMENTAL> a_pressure_spatial;
            if (!BLOCK_SYMMETRIC_MATRIX)
            {
                NEVER_REACHED; // to-come: non-mixed problems
                //a_pressure_spatial = ComputeSpatialPressureMatrixTerm(quad_phi, grad_quad_phi, lin_phi, grad_lin_phi, x, &rElement);
            }

            c_matrix<double,PRESSURE_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> a_pressure_pressure
                = ComputePressurePressureMatrixTerm(linear_phi, grad_linear_phi, X, &rElement);

            for (unsigned i=0; i<SPATIAL_BLOCK_SIZE_ELEMENTAL; i++)
            {
                for (unsigned j=0; j<SPATIAL_BLOCK_SIZE_ELEMENTAL; j++)
                {
                    rAElem(i,j) += a_spatial_spatial(i,j)*wJ;
                }

                for (unsigned j=0; j<PRESSURE_BLOCK_SIZE_ELEMENTAL; j++)
                {
                    rAElem(i, SPATIAL_BLOCK_SIZE_ELEMENTAL + j) += a_spatial_pressure(i,j)*wJ;
                }
            }

            for (unsigned i=0; i<PRESSURE_BLOCK_SIZE_ELEMENTAL; i++)
            {
                if (BLOCK_SYMMETRIC_MATRIX)
                {
                    for (unsigned j=0; j<SPATIAL_BLOCK_SIZE_ELEMENTAL; j++)
                    {
                        rAElem(SPATIAL_BLOCK_SIZE_ELEMENTAL + i, j) += a_spatial_pressure(j,i)*wJ;
                    }
                }
                else
                {
                    NEVER_REACHED; // to-come: non-mixed problems
                }

                for (unsigned j=0; j<PRESSURE_BLOCK_SIZE_ELEMENTAL; j++)
                {
                    rAElem(SPATIAL_BLOCK_SIZE_ELEMENTAL + i, SPATIAL_BLOCK_SIZE_ELEMENTAL + j) += a_pressure_pressure(i,j)*wJ;
                }
            }
        }
    }
}

#endif // ABSTRACTCONTINUUMMECHANICSASSEMBLER_HPP_
