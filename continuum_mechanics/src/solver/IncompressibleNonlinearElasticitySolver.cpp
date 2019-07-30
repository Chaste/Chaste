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

/*
 * NOTE ON COMPILATION ERRORS:
 *
 * This file won't compile with Intel icpc version 9.1.039, with error message:
 * "Terminate with:
  (0): internal error: backend signals"
 *
 * Try recompiling with icpc version 10.0.025.
 */

#include "IncompressibleNonlinearElasticitySolver.hpp"
#include "LinearBasisFunction.hpp"
#include "QuadraticBasisFunction.hpp"
#include <algorithm>

template<size_t DIM>
void IncompressibleNonlinearElasticitySolver<DIM>::AssembleSystem(bool assembleResidual,
                                                                  bool assembleJacobian)
{
    // Check we've actually been asked to do something!
    assert(assembleResidual || assembleJacobian);
    assert(this->mCurrentSolution.size()==this->mNumDofs);

    // Zero the matrix/vector if it is to be assembled
    if (assembleResidual)
    {
        PetscVecTools::Finalise(this->mResidualVector);
        PetscVecTools::Zero(this->mResidualVector);
    }
    if (assembleJacobian)
    {
        PetscMatTools::Zero(this->mrJacobianMatrix);
        PetscMatTools::Zero(this->mPreconditionMatrix);
    }

    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem;
    // The (element) preconditioner matrix: this is the same as the jacobian, but
    // with the mass matrix (ie \intgl phi_i phi_j) in the pressure-pressure block.
    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem_precond;
    c_vector<double, STENCIL_SIZE> b_elem;

    // Loop over elements
    for (typename AbstractTetrahedralMesh<DIM, DIM>::ElementIterator iter = this->mrQuadMesh.GetElementIteratorBegin();
         iter != this->mrQuadMesh.GetElementIteratorEnd();
         ++iter)
    {
        // LCOV_EXCL_START
        // Note: if assembleJacobian only
        if (CommandLineArguments::Instance()->OptionExists("-mech_very_verbose") && assembleJacobian)
        {
            std::cout << "\r[" << PetscTools::GetMyRank() << "]: Element " << (*iter).GetIndex() << " of " << this->mrQuadMesh.GetNumElements() << std::flush;
        }
        // LCOV_EXCL_STOP

        Element<DIM, DIM>& element = *iter;

        if (element.GetOwnership() == true)
        {
            AssembleOnElement(element, a_elem, a_elem_precond, b_elem, assembleResidual, assembleJacobian);

            //// todo: assemble quickly by commenting the AssembleOnElement() and doing
            //// the following, to determine exact non-zeroes per row, and reallocate
            //// with correct nnz (by destroying old matrix and creating a new one)
            //for (unsigned i=0; i<STENCIL_SIZE; i++)
            //{
            //    for (unsigned j=0; j<STENCIL_SIZE; j++)
            //    {
            //        a_elem(i,j)=1.0;
            //    }
            //}


            /////////////////////////////////////////////////////////////////////////////////////////
            // See comments about ordering at the elemental level vs ordering of the global mat/vec
            // in eg AbstractContinuumMechanicsAssembler
            /////////////////////////////////////////////////////////////////////////////////////////

            unsigned p_indices[STENCIL_SIZE];
            for (unsigned i=0; i<NUM_NODES_PER_ELEMENT; i++)
            {
                for (unsigned j=0; j<DIM; j++)
                {
                    // note: DIM+1 is the problem dimension (= this->mProblemDimension)
                    p_indices[DIM*i+j] = (DIM+1)*element.GetNodeGlobalIndex(i) + j;
                }
            }

            for (unsigned i=0; i<NUM_VERTICES_PER_ELEMENT; i++)
            {
                // We assume the vertices are the first num_vertices nodes in the list of nodes
                // in the element. Hence:
                unsigned vertex_index = element.GetNodeGlobalIndex(i);
                // note: DIM+1 is the problem dimension (= this->mProblemDimension)
                p_indices[DIM*NUM_NODES_PER_ELEMENT + i] = (DIM+1)*vertex_index + DIM;
            }

            if (assembleJacobian)
            {
                PetscMatTools::AddMultipleValues<STENCIL_SIZE>(this->mrJacobianMatrix, p_indices, a_elem);
                PetscMatTools::AddMultipleValues<STENCIL_SIZE>(this->mPreconditionMatrix, p_indices, a_elem_precond);
            }

            if (assembleResidual)
            {
                PetscVecTools::AddMultipleValues<STENCIL_SIZE>(this->mResidualVector, p_indices, b_elem);
            }
        }
    }

    // Loop over specified boundary elements and compute surface traction terms
    c_vector<double, BOUNDARY_STENCIL_SIZE> b_boundary_elem; // note BOUNDARY_STENCIL_SIZE = DIM*NUM_BOUNDARY_NODES, as all pressure block is zero
    c_matrix<double, BOUNDARY_STENCIL_SIZE, BOUNDARY_STENCIL_SIZE> a_boundary_elem;

    if (this->mrProblemDefinition.GetTractionBoundaryConditionType() != NO_TRACTIONS)
    {
        for (unsigned bc_index=0; bc_index<this->mrProblemDefinition.rGetTractionBoundaryElements().size(); bc_index++)
        {
            BoundaryElement<DIM-1,DIM>& r_boundary_element = *(this->mrProblemDefinition.rGetTractionBoundaryElements()[bc_index]);

            // If the BCs are tractions applied on a given surface, the boundary integral is independent of u,
            // so a_boundary_elem will be zero (no contribution to jacobian).
            // If the BCs are normal pressure applied to the deformed body, the boundary depends on the deformation,
            // so there is a contribution to the jacobian, and a_boundary_elem is non-zero. Note however that
            // the AssembleOnBoundaryElement() method might decide not to include this, as it can actually
            // cause divergence if the current guess is not close to the true solution
            this->AssembleOnBoundaryElement(r_boundary_element, a_boundary_elem, b_boundary_elem, assembleResidual, assembleJacobian, bc_index);

            unsigned p_indices[BOUNDARY_STENCIL_SIZE];
            for (unsigned i=0; i<NUM_NODES_PER_BOUNDARY_ELEMENT; i++)
            {
                for (unsigned j=0; j<DIM; j++)
                {
                    // note: DIM+1, on the right hand side of the below, is the problem dimension (= this->mProblemDimension)
                    p_indices[DIM*i+j] = (DIM+1)*r_boundary_element.GetNodeGlobalIndex(i) + j;
                }
            }

            if (assembleJacobian)
            {
                PetscMatTools::AddMultipleValues<BOUNDARY_STENCIL_SIZE>(this->mrJacobianMatrix, p_indices, a_boundary_elem);
                PetscMatTools::AddMultipleValues<BOUNDARY_STENCIL_SIZE>(this->mPreconditionMatrix, p_indices, a_boundary_elem);
            }

            if (assembleResidual)
            {
                PetscVecTools::AddMultipleValues<BOUNDARY_STENCIL_SIZE>(this->mResidualVector, p_indices, b_boundary_elem);
            }
        }
    }


    if (assembleResidual)
    {
        PetscVecTools::Finalise(this->mResidualVector);
    }
    if (assembleJacobian)
    {
        PetscMatTools::SwitchWriteMode(this->mrJacobianMatrix);
        PetscMatTools::SwitchWriteMode(this->mPreconditionMatrix);
    }

    if (assembleJacobian)
    {
       this->AddIdentityBlockForDummyPressureVariables(NONLINEAR_PROBLEM_APPLY_TO_EVERYTHING);
    }
    else if (assembleResidual)
    {
       this->AddIdentityBlockForDummyPressureVariables(NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY);
    }

    this->FinishAssembleSystem(assembleResidual, assembleJacobian);
}

template<size_t DIM>
void IncompressibleNonlinearElasticitySolver<DIM>::AssembleOnElement(
            Element<DIM, DIM>& rElement,
            c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElem,
            c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElemPrecond,
            c_vector<double, STENCIL_SIZE>& rBElem,
            bool assembleResidual,
            bool assembleJacobian)
{
    static c_matrix<double,DIM,DIM> jacobian;
    static c_matrix<double,DIM,DIM> inverse_jacobian;
    double jacobian_determinant;

    this->mrQuadMesh.GetInverseJacobianForElement(rElement.GetIndex(), jacobian, jacobian_determinant, inverse_jacobian);

    if (assembleJacobian)
    {
        rAElem.clear();
        rAElemPrecond.clear();
    }

    if (assembleResidual)
    {
        rBElem.clear();
    }

    // Get the current displacement at the nodes
    static c_matrix<double,DIM,NUM_NODES_PER_ELEMENT> element_current_displacements;
    static c_vector<double,NUM_VERTICES_PER_ELEMENT> element_current_pressures;
    for (unsigned II=0; II<NUM_NODES_PER_ELEMENT; II++)
    {
        for (unsigned JJ=0; JJ<DIM; JJ++)
        {
            // note: DIM+1, on the right hand side of the below, is the problem dimension (= this->mProblemDimension)
            element_current_displacements(JJ,II) = this->mCurrentSolution[(DIM+1)*rElement.GetNodeGlobalIndex(II) + JJ];
        }
    }

    // Get the current pressure at the vertices
    for (unsigned II=0; II<NUM_VERTICES_PER_ELEMENT; II++)
    {
        // At the moment we assume the vertices are the first num_vertices nodes in the list of nodes
        // in the mesh. Hence:
        unsigned vertex_index = rElement.GetNodeGlobalIndex(II);

        // note: DIM+1, on the right hand side of the below, is the problem dimension (= this->mProblemDimension)
        element_current_pressures(II) = this->mCurrentSolution[(DIM+1)*vertex_index + DIM];
    }

    // Allocate memory for the basis functions values and derivative values
    static c_vector<double, NUM_VERTICES_PER_ELEMENT> linear_phi;
    static c_vector<double, NUM_NODES_PER_ELEMENT> quad_phi;
    static c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> grad_quad_phi;
    static c_matrix<double, NUM_NODES_PER_ELEMENT, DIM> trans_grad_quad_phi;

    // Get the material law
    AbstractIncompressibleMaterialLaw<DIM>* p_material_law
        = this->mrProblemDefinition.GetIncompressibleMaterialLaw(rElement.GetIndex());

    static c_matrix<double,DIM,DIM> grad_u; // grad_u = (du_i/dX_M)

    static c_matrix<double,DIM,DIM> F;      // the deformation gradient, F = dx/dX, F_{iM} = dx_i/dX_M
    static c_matrix<double,DIM,DIM> C;      // Green deformation tensor, C = F^T F
    static c_matrix<double,DIM,DIM> inv_C;  // inverse(C)
    static c_matrix<double,DIM,DIM> inv_F;  // inverse(F)
    static c_matrix<double,DIM,DIM> T;      // Second Piola-Kirchoff stress tensor (= dW/dE = 2dW/dC)

    static c_matrix<double,DIM,DIM> F_T;    // F*T
    static c_matrix<double,DIM,NUM_NODES_PER_ELEMENT> F_T_grad_quad_phi; // F*T*grad_quad_phi

    c_vector<double,DIM> body_force;

    static FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE;    // dTdE(M,N,P,Q) = dT_{MN}/dE_{PQ}
    static FourthOrderTensor<DIM,DIM,DIM,DIM> dSdF;    // dSdF(M,i,N,j) = dS_{Mi}/dF_{jN}

    static FourthOrderTensor<NUM_NODES_PER_ELEMENT,DIM,DIM,DIM> temp_tensor;
    static FourthOrderTensor<NUM_NODES_PER_ELEMENT,DIM,NUM_NODES_PER_ELEMENT,DIM> dSdF_quad_quad;

    static c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> temp_matrix;
    static c_matrix<double,NUM_NODES_PER_ELEMENT,DIM> grad_quad_phi_times_invF;


    if (this->mSetComputeAverageStressPerElement)
    {
        this->mAverageStressesPerElement[rElement.GetIndex()] = zero_vector<double>(DIM*(DIM+1)/2);
    }

    // Loop over Gauss points
    for (unsigned quadrature_index=0; quadrature_index < this->mpQuadratureRule->GetNumQuadPoints(); quadrature_index++)
    {
        // This is needed by the cardiac mechanics solver
        unsigned current_quad_point_global_index =   rElement.GetIndex()*this->mpQuadratureRule->GetNumQuadPoints()
                                                   + quadrature_index;

        double wJ = jacobian_determinant * this->mpQuadratureRule->GetWeight(quadrature_index);

        const ChastePoint<DIM>& quadrature_point = this->mpQuadratureRule->rGetQuadPoint(quadrature_index);

        // Set up basis function information
        LinearBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, linear_phi);
        QuadraticBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, quad_phi);
        QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, inverse_jacobian, grad_quad_phi);
        trans_grad_quad_phi = trans(grad_quad_phi);

        // Get the body force, interpolating X if necessary
        if (assembleResidual)
        {
            switch (this->mrProblemDefinition.GetBodyForceType())
            {
                case FUNCTIONAL_BODY_FORCE:
                {
                    c_vector<double,DIM> X = zero_vector<double>(DIM);
                    // interpolate X (using the vertices and the /linear/ bases, as no curvilinear elements
                    for (unsigned node_index=0; node_index<NUM_VERTICES_PER_ELEMENT; node_index++)
                    {
                        X += linear_phi(node_index)*this->mrQuadMesh.GetNode( rElement.GetNodeGlobalIndex(node_index) )->rGetLocation();
                    }
                    body_force = this->mrProblemDefinition.EvaluateBodyForceFunction(X, this->mCurrentTime);
                    break;
                }
                case CONSTANT_BODY_FORCE:
                {
                    body_force = this->mrProblemDefinition.GetConstantBodyForce();
                    break;
                }
                default:
                    NEVER_REACHED;
            }
        }

        // Interpolate grad_u and p
        grad_u = zero_matrix<double>(DIM,DIM);

        for (unsigned node_index=0; node_index<NUM_NODES_PER_ELEMENT; node_index++)
        {
            for (unsigned i=0; i<DIM; i++)
            {
                for (unsigned M=0; M<DIM; M++)
                {
                    grad_u(i,M) += grad_quad_phi(M,node_index)*element_current_displacements(i,node_index);
                }
            }
        }

        double pressure = 0;
        for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
        {
            pressure += linear_phi(vertex_index)*element_current_pressures(vertex_index);
        }

        // Calculate C, inv(C) and T
        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned M=0; M<DIM; M++)
            {
                F(i,M) = (i==M?1:0) + grad_u(i,M);
            }
        }

        C = prod(trans(F),F);
        inv_C = Inverse(C);
        inv_F = Inverse(F);

        double detF = Determinant(F);

        // Compute the passive stress, and dTdE corresponding to passive stress
        this->SetupChangeOfBasisMatrix(rElement.GetIndex(), current_quad_point_global_index);
        p_material_law->SetChangeOfBasisMatrix(this->mChangeOfBasisMatrix);
        p_material_law->ComputeStressAndStressDerivative(C, inv_C, pressure, T, dTdE, assembleJacobian);

        if (this->mIncludeActiveTension)
        {
            // Add any active stresses, if there are any. Requires subclasses to overload this method,
            // see for example the cardiac mechanics assemblers.
            this->AddActiveStressAndStressDerivative(C, rElement.GetIndex(), current_quad_point_global_index,
                                                     T, dTdE, assembleJacobian);
        }

        if (this->mSetComputeAverageStressPerElement)
        {
            this->AddStressToAverageStressPerElement(T,rElement.GetIndex());
        }

        // Residual vector
        if (assembleResidual)
        {
            F_T = prod(F,T);
            F_T_grad_quad_phi = prod(F_T, grad_quad_phi);

            for (unsigned index=0; index<NUM_NODES_PER_ELEMENT*DIM; index++)
            {
                unsigned spatial_dim = index%DIM;
                unsigned node_index = (index-spatial_dim)/DIM;

                rBElem(index) +=  - this->mrProblemDefinition.GetDensity()
                                  * body_force(spatial_dim)
                                  * quad_phi(node_index)
                                  * wJ;

                // The T(M,N)*F(spatial_dim,M)*grad_quad_phi(N,node_index) term
                rBElem(index) +=   F_T_grad_quad_phi(spatial_dim,node_index)
                                 * wJ;
            }

            for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
            {
                rBElem( NUM_NODES_PER_ELEMENT*DIM + vertex_index ) +=   linear_phi(vertex_index)
                                                                      * (detF - 1)
                                                                      * wJ;
            }
        }

        // Jacobian matrix
        if (assembleJacobian)
        {
            // Save trans(grad_quad_phi) * invF
            grad_quad_phi_times_invF = prod(trans_grad_quad_phi, inv_F);

            /////////////////////////////////////////////////////////////////////////////////////////////
            // Set up the tensor dSdF
            //
            // dSdF as a function of T and dTdE (which is what the material law returns) is given by:
            //
            // dS_{Mi}/dF_{jN} = (dT_{MN}/dC_{PQ}+dT_{MN}/dC_{PQ}) F{iP} F_{jQ}  + T_{MN} delta_{ij}
            //
            // todo1: this should probably move into the material law (but need to make sure
            // memory is handled efficiently
            // todo2: get material law to return this immediately, not dTdE
            /////////////////////////////////////////////////////////////////////////////////////////////

            // Set up the tensor 0.5(dTdE(M,N,P,Q) + dTdE(M,N,Q,P))
            for (unsigned M=0; M<DIM; M++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    for (unsigned P=0; P<DIM; P++)
                    {
                        for (unsigned Q=0; Q<DIM; Q++)
                        {
                            // This is NOT dSdF, just using this as storage space
                            dSdF(M,N,P,Q) = 0.5*(dTdE(M,N,P,Q) + dTdE(M,N,Q,P));
                        }
                    }
                }
            }

            // This is NOT dTdE, just reusing memory. A^{MdPQ}  = F^d_N * dTdE_sym^{MNPQ}
            dTdE.template SetAsContractionOnSecondDimension<DIM>(F, dSdF);

            // dSdF{MdPe} := F^d_N * F^e_Q * dTdE_sym^{MNPQ}
            dSdF.template SetAsContractionOnFourthDimension<DIM>(F, dTdE);

            // Now add the T_{MN} delta_{ij} term
            for (unsigned M=0; M<DIM; M++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    for (unsigned i=0; i<DIM; i++)
                    {
                        dSdF(M,i,N,i) += T(M,N);
                    }
                }
            }

            ///////////////////////////////////////////////////////
            // Set up the tensor
            //   dSdF_quad_quad(node_index1, spatial_dim1, node_index2, spatial_dim2)
            //            =    dS_{M,spatial_dim1}/d_F{spatial_dim2,N}
            //               * grad_quad_phi(M,node_index1)
            //               * grad_quad_phi(P,node_index2)
            //
            //            =    dSdF(M,spatial_index1,N,spatial_index2)
            //               * grad_quad_phi(M,node_index1)
            //               * grad_quad_phi(P,node_index2)
            //
            ///////////////////////////////////////////////////////
            temp_tensor.template SetAsContractionOnFirstDimension<DIM>( trans_grad_quad_phi, dSdF );
            dSdF_quad_quad.template SetAsContractionOnThirdDimension<DIM>( trans_grad_quad_phi, temp_tensor );

            for (unsigned index1=0; index1<NUM_NODES_PER_ELEMENT*DIM; index1++)
            {
                unsigned spatial_dim1 = index1%DIM;
                unsigned node_index1 = (index1-spatial_dim1)/DIM;


                for (unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
                {
                    unsigned spatial_dim2 = index2%DIM;
                    unsigned node_index2 = (index2-spatial_dim2)/DIM;

                    // The dSdF*grad_quad_phi*grad_quad_phi term
                    rAElem(index1,index2)  +=   dSdF_quad_quad(node_index1,spatial_dim1,node_index2,spatial_dim2)
                                              * wJ;
                }

                for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
                {
                    unsigned index2 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;

                    // The -invF(M,spatial_dim1)*grad_quad_phi(M,node_index1)*linear_phi(vertex_index) term
                    rAElem(index1,index2)  +=  - grad_quad_phi_times_invF(node_index1,spatial_dim1)
                                               * linear_phi(vertex_index)
                                               * wJ;
                }
            }

            for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
            {
                unsigned index1 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;

                for (unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
                {
                    unsigned spatial_dim2 = index2%DIM;
                    unsigned node_index2 = (index2-spatial_dim2)/DIM;

                    // Same as (negative of) the opposite block (ie a few lines up), except for detF
                    rAElem(index1,index2) +=   detF
                                             * grad_quad_phi_times_invF(node_index2,spatial_dim2)
                                             * linear_phi(vertex_index)
                                             * wJ;
                }

                /////////////////////////////////////////////////////
                // Preconditioner matrix
                // Fill the mass matrix (ie \intgl phi_i phi_j) in the
                // pressure-pressure block. Note, the rest of the
                // entries are filled in at the end
                /////////////////////////////////////////////////////
                for (unsigned vertex_index2=0; vertex_index2<NUM_VERTICES_PER_ELEMENT; vertex_index2++)
                {
                    unsigned index2 = NUM_NODES_PER_ELEMENT*DIM + vertex_index2;
                    rAElemPrecond(index1,index2) +=   linear_phi(vertex_index)
                                                    * linear_phi(vertex_index2)
                                                    * wJ;
                }
            }
        }
    }

    if (assembleJacobian)
    {
        if (this->mPetscDirectSolve)
        {
            // Petsc will do an LU factorisation of the preconditioner, which we
            // set equal to  [ A  B1^T ]
            //               [ B2  M   ]
            // The reason for the mass matrix is to avoid zeros on the diagonal
            rAElemPrecond = rAElemPrecond + rAElem;
        }
        else
        {
            // Fill in the other blocks of the preconditioner matrix, by adding
            // the Jacobian matrix (this doesn't effect the pressure-pressure block
            // of rAElemPrecond as the pressure-pressure block of rAElem is zero),
            // and the zero a block.
            //
            // The following altogether gives the preconditioner  [ A  B1^T ]
            //                                                    [ 0  M    ]
            rAElemPrecond = rAElemPrecond + rAElem;

            for (unsigned i=NUM_NODES_PER_ELEMENT*DIM; i<STENCIL_SIZE; i++)
            {
                for (unsigned j=0; j<NUM_NODES_PER_ELEMENT*DIM; j++)
                {
                    rAElemPrecond(i,j) = 0.0;
                }
            }
        }
    }

    if (this->mSetComputeAverageStressPerElement)
    {
        for (unsigned i=0; i<DIM*(DIM+1)/2; i++)
        {
            this->mAverageStressesPerElement[rElement.GetIndex()](i) /= this->mpQuadratureRule->GetNumQuadPoints();
        }
    }
}

template<size_t DIM>
void IncompressibleNonlinearElasticitySolver<DIM>::FormInitialGuess()
{
    this->mCurrentSolution.resize(this->mNumDofs, 0.0);

    for (typename AbstractTetrahedralMesh<DIM, DIM>::ElementIterator iter = this->mrQuadMesh.GetElementIteratorBegin();
         iter != this->mrQuadMesh.GetElementIteratorEnd();
         ++iter)
    {
        ///\todo #2223 This is unlikely to work using DistributedQuadraticMesh with Non-homogeneous material laws
        double zero_strain_pressure
           = this->mrProblemDefinition.GetIncompressibleMaterialLaw(iter->GetIndex())->GetZeroStrainPressure();


        // Loop over vertices and set pressure solution to be zero-strain-pressure
        for (unsigned j=0; j<NUM_VERTICES_PER_ELEMENT; j++)
        {
            // We assume the vertices are the first num_vertices nodes in the list of nodes
            // in the element. Hence:
            unsigned vertex_index = iter->GetNodeGlobalIndex(j);
            // note: DIM+1 is the problem dimension (= this->mProblemDimension)
            this->mCurrentSolution[ (DIM+1)*vertex_index + DIM ] = zero_strain_pressure;
        }
    }
}

template<size_t DIM>
IncompressibleNonlinearElasticitySolver<DIM>::IncompressibleNonlinearElasticitySolver(
        AbstractTetrahedralMesh<DIM,DIM>& rQuadMesh,
        SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
        std::string outputDirectory)
    : AbstractNonlinearElasticitySolver<DIM>(rQuadMesh,
                                             rProblemDefinition,
                                             outputDirectory,
                                             INCOMPRESSIBLE)
{
    if (rProblemDefinition.GetCompressibilityType() != INCOMPRESSIBLE)
    {
        EXCEPTION("SolidMechanicsProblemDefinition object contains compressible material laws");
    }

    FormInitialGuess();
}

// Explicit instantiation
template class IncompressibleNonlinearElasticitySolver<2>;
template class IncompressibleNonlinearElasticitySolver<3>;
