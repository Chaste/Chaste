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

#ifndef ABSTRACTFECABLEINTEGRALASSEMBLER_HPP_
#define ABSTRACTFECABLEINTEGRALASSEMBLER_HPP_

#include "AbstractFeAssemblerCommon.hpp"
#include "GaussianQuadratureRule.hpp"
#include "PetscVecTools.hpp"
#include "PetscMatTools.hpp"

/**
 * The class in similar to AbstractFeVolumeIntegralAssembler (see documentation for this), but is for
 * creating a finite element matrices or vectors that involve integrals over CABLES, ie 1d regions
 * in a 2d/3d mesh. Required for cardiac simulations with a Purkinje network. Uses
 * a MixedDimensionMesh, which is composed of the normal mesh plus cables.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
class AbstractFeCableIntegralAssembler : public AbstractFeAssemblerCommon<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX,INTERPOLATION_LEVEL>
{
protected:
    /** Cable element dimension. */
    static const unsigned CABLE_ELEMENT_DIM = 1;

    /** Number of nodes in a cable element. */
    static const unsigned NUM_CABLE_ELEMENT_NODES = 2;

    /** Mesh to be solved on. */
    MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>* mpMesh;

    /** Quadrature rule for use on cable elements. */
    GaussianQuadratureRule<1>* mpCableQuadRule;

    /** Basis function for use with normal elements. */
    typedef LinearBasisFunction<1> CableBasisFunction;

    /**
     * Compute the derivatives of all basis functions at a point within a cable element.
     * This method will transform the results, for use within Gaussian quadrature
     * for example.
     *
     * @param rPoint The point at which to compute the basis functions. The
     *     results are undefined if this is not within the canonical element.
     * @param rInverseJacobian The inverse of the Jacobian matrix mapping the real
     *     element into the canonical element.
     * @param rReturnValue A reference to a vector, to be filled in
     * @return The derivatives of the basis functions, in local index order. Each
     *     entry is a vector (c_vector<double, SPACE_DIM> instance) giving the
     *     derivative along each axis.
     */
    void ComputeTransformedBasisFunctionDerivatives(const ChastePoint<CABLE_ELEMENT_DIM>& rPoint,
                                                    const c_matrix<double, CABLE_ELEMENT_DIM, SPACE_DIM>& rInverseJacobian,
                                                    c_matrix<double, SPACE_DIM, NUM_CABLE_ELEMENT_NODES>& rReturnValue);

    /**
     * The main assembly method. Protected, should only be called through Assemble(),
     * AssembleMatrix() or AssembleVector() which set mAssembleMatrix, mAssembleVector
     * accordingly.
     */
    void DoAssemble();

    /**
     * @return the matrix to be added to element stiffness matrix
     * for a given Gauss point, ie, essentially the INTEGRAND in the integral
     * definition of the matrix (integral over cable region).
     *
     * The arguments are the bases, bases gradients,
     * x and current solution computed at the Gauss point. The returned matrix
     * will be multiplied by the Gauss weight and Jacobian determinant and
     * added to the element stiffness matrix (see AssembleOnElement()).
     *
     *  ** This method has to be implemented in the concrete class if CAN_ASSEMBLE_MATRIX is true. **
     *
     * NOTE: for linear problems rGradU is NOT set up correctly because it should
     * not be needed.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases.
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i).
     * @param rX The point in space.
     * @param rU The unknown as a vector, u(i) = u_i.
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j).
     * @param pElement Pointer to the element.
     */
    // LCOV_EXCL_START
    virtual c_matrix<double,PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES,PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES> ComputeCableMatrixTerm(
        c_vector<double, NUM_CABLE_ELEMENT_NODES>& rPhi,
        c_matrix<double, SPACE_DIM, NUM_CABLE_ELEMENT_NODES>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
        Element<CABLE_ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        // If this line is reached this means this method probably hasn't been over-ridden correctly in
        // the concrete class
        NEVER_REACHED;
        return zero_matrix<double>(PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES,PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES);
    }
    // LCOV_EXCL_STOP

    /**
     * @return the vector to be added to element stiffness vector
     * for a given Gauss point, ie, essentially the INTEGRAND in the integral
     * definition of the vector. The arguments are the bases,
     * x and current solution computed at the Gauss point. The returned vector
     * will be multiplied by the Gauss weight and Jacobian determinant and
     * added to the element stiffness matrix (see AssembleOnElement()).
     *
     * ** This method has to be implemented in the concrete class if CAN_ASSEMBLE_VECTOR is true. **
     *
     * NOTE: for linear problems rGradPhi and rGradU are NOT set up correctly because
     * they should not be needed.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    // LCOV_EXCL_START
    virtual c_vector<double,PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES> ComputeCableVectorTerm(
        c_vector<double, NUM_CABLE_ELEMENT_NODES>& rPhi,
        c_matrix<double, SPACE_DIM, NUM_CABLE_ELEMENT_NODES>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
        Element<CABLE_ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        // If this line is reached this means this method probably hasn't been over-ridden correctly in
        // the concrete class
        NEVER_REACHED;
        return zero_vector<double>(PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES);
    }
    // LCOV_EXCL_STOP

    /**
     * Calculate the contribution of a single cable element to the linear system.
     *
     * @param rElement The element to assemble on.
     * @param rAElem The element's contribution to the LHS matrix is returned in this
     *    n by n matrix, where n is the no. of nodes in this element. There is no
     *    need to zero this matrix before calling.
     * @param rBElem The element's contribution to the RHS vector is returned in this
     *    vector of length n, the no. of nodes in this element. There is no
     *    need to zero this vector before calling.
     *
     * Called by AssembleSystem().
     * Calls ComputeCableMatrixTerm() etc.
     */
    virtual void AssembleOnCableElement(Element<CABLE_ELEMENT_DIM,SPACE_DIM>& rElement,
                                        c_matrix<double, PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES, PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES >& rAElem,
                                        c_vector<double, PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES>& rBElem);


    /**
     * @return true if we should include this (cable) element when assembling. Returns true
     * here but can be overridden by the concrete assembler if not all
     * elements should be included.
     *
     * @param rElement the element
     */
    virtual bool ElementAssemblyCriterion(Element<CABLE_ELEMENT_DIM,SPACE_DIM>& rElement)
    {
        return true;
    }

public:

    /**
     * Constructor.
     *
     * @param pMesh The mesh
     */
    AbstractFeCableIntegralAssembler(MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /**
     * Destructor.
     */
    virtual ~AbstractFeCableIntegralAssembler()
    {
        delete mpCableQuadRule;
    }
};

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
AbstractFeCableIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::AbstractFeCableIntegralAssembler(
            MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    : AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>(),
      mpMesh(pMesh)
{
    assert(pMesh);
    assert(CAN_ASSEMBLE_VECTOR || CAN_ASSEMBLE_MATRIX);
    // Default to 2nd order quadrature.  Our default basis functions are piecewise linear
    // which means that we are integrating functions which in the worst case (mass matrix)
    // are quadratic.
    mpCableQuadRule = new GaussianQuadratureRule<CABLE_ELEMENT_DIM>(2);

    // Not supporting this yet - if a nonlinear assembler on cable elements is required, uncomment code
    // in AssembleOnCableElement below (search for NONLINEAR)
    assert(INTERPOLATION_LEVEL!=NONLINEAR);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractFeCableIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::DoAssemble()
{
    HeartEventHandler::EventType assemble_event;
    if (this->mAssembleMatrix)
    {
        assemble_event = HeartEventHandler::ASSEMBLE_SYSTEM;
    }
    else
    {
        assemble_event = HeartEventHandler::ASSEMBLE_RHS;
    }

    if (this->mAssembleMatrix && this->mMatrixToAssemble==NULL)
    {
        EXCEPTION("Matrix to be assembled has not been set");
    }
    if (this->mAssembleVector && this->mVectorToAssemble==NULL)
    {
        EXCEPTION("Vector to be assembled has not been set");
    }

    // This has to be below PrepareForAssembleSystem as in that method the ODEs are solved in cardiac problems
    HeartEventHandler::BeginEvent(assemble_event);

    // Zero the matrix/vector if it is to be assembled
    if (this->mAssembleVector && this->mZeroVectorBeforeAssembly)
    {
        PetscVecTools::Zero(this->mVectorToAssemble);
    }
    if (this->mAssembleMatrix && this->mZeroMatrixBeforeAssembly)
    {
        PetscMatTools::Zero(this->mMatrixToAssemble);
    }

    const size_t STENCIL_SIZE=PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES;
    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem;
    c_vector<double, STENCIL_SIZE> b_elem;

    // Loop over elements
    if (this->mAssembleMatrix || this->mAssembleVector)
    {
        for (typename MixedDimensionMesh<CABLE_ELEMENT_DIM, SPACE_DIM>::CableElementIterator iter = mpMesh->GetCableElementIteratorBegin();
             iter != mpMesh->GetCableElementIteratorEnd();
             ++iter)
        {
            Element<CABLE_ELEMENT_DIM, SPACE_DIM>& r_element = *(*iter);

            // Test for ownership first, since it's pointless to test the criterion on something which we might know nothing about.
            if (r_element.GetOwnership() == true && ElementAssemblyCriterion(r_element)==true)
            {
                AssembleOnCableElement(r_element, a_elem, b_elem);

                unsigned p_indices[STENCIL_SIZE];
                r_element.GetStiffnessMatrixGlobalIndices(PROBLEM_DIM, p_indices);

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

    HeartEventHandler::EndEvent(assemble_event);
}

///////////////////////////////////////////////////////////////////////////////////
// Implementation - AssembleOnCableElement and smaller
///////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractFeCableIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::ComputeTransformedBasisFunctionDerivatives(
        const ChastePoint<CABLE_ELEMENT_DIM>& rPoint,
        const c_matrix<double, CABLE_ELEMENT_DIM, SPACE_DIM>& rInverseJacobian,
        c_matrix<double, SPACE_DIM, NUM_CABLE_ELEMENT_NODES>& rReturnValue)
{
    static c_matrix<double, CABLE_ELEMENT_DIM, NUM_CABLE_ELEMENT_NODES> grad_phi;

    LinearBasisFunction<CABLE_ELEMENT_DIM>::ComputeBasisFunctionDerivatives(rPoint, grad_phi);
    rReturnValue = prod(trans(rInverseJacobian), grad_phi);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractFeCableIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::AssembleOnCableElement(
    Element<CABLE_ELEMENT_DIM,SPACE_DIM>& rElement,
    c_matrix<double, PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES, PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES >& rAElem,
    c_vector<double, PROBLEM_DIM*NUM_CABLE_ELEMENT_NODES>& rBElem)
{
    /**
     * \todo #1320 This assumes that the Jacobian is constant on an element.
     * This is true for linear basis functions, but not for any other type of
     * basis function.
     */
    c_matrix<double, SPACE_DIM, CABLE_ELEMENT_DIM> jacobian;
    c_matrix<double, CABLE_ELEMENT_DIM, SPACE_DIM> inverse_jacobian;
    double jacobian_determinant;

    ////NOTE: the normal assembler calls this, but this wouldn't work here - it would end up looking at the volume
    //// element with the same index as this cable element.
    // mpMesh->GetInverseJacobianForElement(rElement.GetIndex(), jacobian, jacobian_determinant, inverse_jacobian);
    //// So call this instead
    rElement.CalculateInverseJacobian(jacobian, jacobian_determinant, inverse_jacobian);

    if (this->mAssembleMatrix)
    {
        rAElem.clear();
    }

    if (this->mAssembleVector)
    {
        rBElem.clear();
    }

    const unsigned num_nodes = rElement.GetNumNodes();

    // Allocate memory for the basis functions values and derivative values
    c_vector<double, NUM_CABLE_ELEMENT_NODES> phi;
    c_matrix<double, SPACE_DIM, NUM_CABLE_ELEMENT_NODES> grad_phi;

    // Loop over Gauss points
    for (unsigned quad_index=0; quad_index < mpCableQuadRule->GetNumQuadPoints(); quad_index++)
    {
        const ChastePoint<CABLE_ELEMENT_DIM>& quad_point = mpCableQuadRule->rGetQuadPoint(quad_index);

        CableBasisFunction::ComputeBasisFunctions(quad_point, phi);

        if (this->mAssembleMatrix || INTERPOLATION_LEVEL==NONLINEAR)
        {
            ComputeTransformedBasisFunctionDerivatives(quad_point, inverse_jacobian, grad_phi);
        }

        // Location of the Gauss point in the original element will be stored in x
        // Where applicable, u will be set to the value of the current solution at x
        ChastePoint<SPACE_DIM> x(0,0,0);

        c_vector<double,PROBLEM_DIM> u = zero_vector<double>(PROBLEM_DIM);
        c_matrix<double,PROBLEM_DIM,SPACE_DIM> grad_u = zero_matrix<double>(PROBLEM_DIM,SPACE_DIM);

        // Allow the concrete version of the assembler to interpolate any desired quantities
        this->ResetInterpolatedQuantities();

        // Interpolation
        for (unsigned i=0; i<num_nodes; i++)
        {
            const Node<SPACE_DIM>* p_node = rElement.GetNode(i);

            if (INTERPOLATION_LEVEL != CARDIAC) // don't interpolate X if cardiac problem
            {
                const c_vector<double, SPACE_DIM>& r_node_loc = p_node->rGetLocation();
                // interpolate x
                x.rGetLocation() += phi(i)*r_node_loc;
            }

            // Interpolate u and grad u if a current solution or guess exists
            unsigned node_global_index = rElement.GetNodeGlobalIndex(i);
            if (this->mCurrentSolutionOrGuessReplicated.GetSize() > 0)
            {
                for (unsigned index_of_unknown=0; index_of_unknown<(INTERPOLATION_LEVEL!=CARDIAC ? PROBLEM_DIM : 1); index_of_unknown++)
                {
                    /*
                     * If we have a solution (e.g. this is a dynamic problem) then
                     * interpolate the value at this quadrature point.
                     *
                     * NOTE: the following assumes that if, say, there are two unknowns
                     * u and v, they are stored in the current solution vector as
                     * [U1 V1 U2 V2 ... U_n V_n].
                     */
                    double u_at_node = this->GetCurrentSolutionOrGuessValue(node_global_index, index_of_unknown);
                    u(index_of_unknown) += phi(i)*u_at_node;

//// See assertion in constructor
//                    if (INTERPOLATION_LEVEL==NONLINEAR) // don't need to construct grad_phi or grad_u in other cases
//                    {
//                        for (unsigned j=0; j<SPACE_DIM; j++)
//                        {
//                            grad_u(index_of_unknown,j) += grad_phi(j,i)*u_at_node;
//                        }
//                    }
                }
            }

            // Allow the concrete version of the assembler to interpolate any desired quantities
            this->IncrementInterpolatedQuantities(phi(i), p_node);
        }

        double wJ = jacobian_determinant * mpCableQuadRule->GetWeight(quad_index);

        // Create rAElem and rBElem
        if (this->mAssembleMatrix)
        {
            noalias(rAElem) += ComputeCableMatrixTerm(phi, grad_phi, x, u, grad_u, &rElement) * wJ;
        }

        if (this->mAssembleVector)
        {
            noalias(rBElem) += ComputeCableVectorTerm(phi, grad_phi, x, u, grad_u, &rElement) * wJ;
        }
    }
}

#endif /*ABSTRACTFECABLEINTEGRALASSEMBLER_HPP_*/
