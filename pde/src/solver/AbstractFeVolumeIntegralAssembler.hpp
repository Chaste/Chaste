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

#ifndef ABSTRACTFEVOLUMEINTEGRALASSEMBLER_HPP_
#define ABSTRACTFEVOLUMEINTEGRALASSEMBLER_HPP_

#include "AbstractFeAssemblerCommon.hpp"
#include "GaussianQuadratureRule.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "PetscVecTools.hpp"
#include "PetscMatTools.hpp"

/**
 *
 * An abstract class for creating finite element vectors or matrices that are defined
 * by integrals over the computational domain of functions of basis functions (for
 * example, stiffness or mass matrices), that require assembly by looping over
 * each element in the mesh and computing the element-wise integrals and adding it to
 * the full matrix or vector.
 *
 * This class is used for VOLUME integrals. For surface integrals there is a
 * similar class, AbstractFeSurfaceIntegralAssembler.
 *
 * This class can be used to assemble a matrix OR a vector OR one of each. The
 * template booleans CAN_ASSEMBLE_VECTOR and CAN_ASSEMBLE_MATRIX should be chosen
 * accordingly.
 *
 * The class provides the functionality to loop over elements, perform element-wise
 * integration (using Gaussian quadrature and linear basis functions), and add the
 * results to the final matrix or vector. The concrete class which inherits from this
 * must implement either COMPUTE_MATRIX_TERM or COMPUTE_VECTOR_TERM or both, which
 * should return the INTEGRAND, as a function of the basis functions.
 *
 * The final template parameter defines how much interpolation (onto quadrature points)
 * is required by the concrete class.
 *
 * CARDIAC: only interpolates the first component of the unknown (ie the voltage)
 * NORMAL: interpolates the position X and all components of the unknown u
 * NONLINEAR: interpolates X, u and grad(u). Also computes the gradient of the
 *  basis functions when assembling vectors.
 *
 * This class inherits from AbstractFeAssemblerCommon which is where some member variables
 * (the matrix/vector to be created, for example) are defined.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
class AbstractFeVolumeIntegralAssembler :
     public AbstractFeAssemblerCommon<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX,INTERPOLATION_LEVEL>
{
protected:
    /** Mesh to be solved on. */
    AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* mpMesh;

    /** Quadrature rule for use on normal elements. */
    GaussianQuadratureRule<ELEMENT_DIM>* mpQuadRule;

    /** Basis function for use with normal elements. */
    typedef LinearBasisFunction<ELEMENT_DIM> BasisFunction;

    /**
     * Compute the derivatives of all basis functions at a point within an element.
     * This method will transform the results, for use within Gaussian quadrature
     * for example.
     *
     * This is almost identical to LinearBasisFunction::ComputeTransformedBasisFunctionDerivatives,
     * except that it is also templated over SPACE_DIM and can handle cases such as 1d in 3d space.
     *
     * \todo #1319 Template LinearBasisFunction over SPACE_DIM and remove this method?
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
    void ComputeTransformedBasisFunctionDerivatives(const ChastePoint<ELEMENT_DIM>& rPoint,
                                                    const c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian,
                                                    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rReturnValue);

    /**
     * The main assembly method. Should only be called through Assemble(),
     * AssembleMatrix() or AssembleVector() which set mAssembleMatrix, mAssembleVector
     * accordingly.
     */
    void DoAssemble();

protected:

    /**
     * @return the matrix to be added to element stiffness matrix
     * for a given Gauss point, ie, essentially the INTEGRAND in the integral
     * definition of the matrix. The arguments are the bases, bases gradients,
     * x and current solution computed at the Gauss point. The returned matrix
     * will be multiplied by the Gauss weight and Jacobian determinant and
     * added to the element stiffness matrix (see AssembleOnElement()).
     *
     *  ** This method has to be implemented in the concrete class if CAN_ASSEMBLE_MATRIX is true. **
     *
     * NOTE: When INTERPOLATION_LEVEL==NORMAL, rGradU does not get set up and should not be used.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases.
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i).
     * @param rX The point in space.
     * @param rU The unknown as a vector, u(i) = u_i.
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j).
     * @param pElement Pointer to the element.
     */
    // LCOV_EXCL_START
    virtual c_matrix<double,PROBLEM_DIM*(ELEMENT_DIM+1),PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        // If this line is reached this means this method probably hasn't been over-ridden correctly in
        // the concrete class
        NEVER_REACHED;
        return zero_matrix<double>(PROBLEM_DIM*(ELEMENT_DIM+1),PROBLEM_DIM*(ELEMENT_DIM+1));
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
     * NOTE: When INTERPOLATION_LEVEL==NORMAL, rGradPhi and rGradU do not get set up and should not be used.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    // LCOV_EXCL_START
    virtual c_vector<double,PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        // If this line is reached this means this method probably hasn't been over-ridden correctly in
        // the concrete class
        NEVER_REACHED;
        return zero_vector<double>(PROBLEM_DIM*(ELEMENT_DIM+1));
    }
    // LCOV_EXCL_STOP


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
     *
     * Called by AssembleSystem().
     * Calls ComputeMatrixTerm() etc.
     */
    virtual void AssembleOnElement(Element<ELEMENT_DIM,SPACE_DIM>& rElement,
                                   c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1) >& rAElem,
                                   c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)>& rBElem);

    /**
     * @return true if we should include this (volume) element when assembling. Returns true
     * here but can be overridden by the concrete assembler if not all
     * elements should be included.
     *
     * @param rElement the element
     */
    virtual bool ElementAssemblyCriterion(Element<ELEMENT_DIM,SPACE_DIM>& rElement)
    {
        return true;
    }


public:

    /**
     * Constructor.
     *
     * @param pMesh The mesh
     */
    AbstractFeVolumeIntegralAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /**
     * Destructor.
     */
    virtual ~AbstractFeVolumeIntegralAssembler()
    {
        delete mpQuadRule;
    }
};

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::AbstractFeVolumeIntegralAssembler(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    : AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>(),
      mpMesh(pMesh)
{
    assert(pMesh);
    // Default to 2nd order quadrature.  Our default basis functions are piecewise linear
    // which means that we are integrating functions which in the worst case (mass matrix)
    // are quadratic.
    mpQuadRule = new GaussianQuadratureRule<ELEMENT_DIM>(2);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::DoAssemble()
{
    assert(this->mAssembleMatrix || this->mAssembleVector);

    HeartEventHandler::EventType assemble_event;
    if (this->mAssembleMatrix)
    {
        assemble_event = HeartEventHandler::ASSEMBLE_SYSTEM;
    }
    else
    {
        assemble_event = HeartEventHandler::ASSEMBLE_RHS;
    }

    if (this->mAssembleMatrix && this->mMatrixToAssemble==nullptr)
    {
        EXCEPTION("Matrix to be assembled has not been set");
    }
    if (this->mAssembleVector && this->mVectorToAssemble==nullptr)
    {
        EXCEPTION("Vector to be assembled has not been set");
    }

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

    const size_t STENCIL_SIZE=PROBLEM_DIM*(ELEMENT_DIM+1);
    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem;
    c_vector<double, STENCIL_SIZE> b_elem;

    // Loop over elements
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = mpMesh->GetElementIteratorBegin();
         iter != mpMesh->GetElementIteratorEnd();
         ++iter)
    {
        Element<ELEMENT_DIM, SPACE_DIM>& r_element = *iter;

        // Test for ownership first, since it's pointless to test the criterion on something which we might know nothing about.
        if (r_element.GetOwnership() == true && ElementAssemblyCriterion(r_element)==true)
        {
            AssembleOnElement(r_element, a_elem, b_elem);

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

    HeartEventHandler::EndEvent(assemble_event);
}


///////////////////////////////////////////////////////////////////////////////////
// Implementation - AssembleOnElement and smaller
///////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::ComputeTransformedBasisFunctionDerivatives(
        const ChastePoint<ELEMENT_DIM>& rPoint,
        const c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rReturnValue)
{
    assert(ELEMENT_DIM < 4 && ELEMENT_DIM > 0);
    static c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> grad_phi;

    LinearBasisFunction<ELEMENT_DIM>::ComputeBasisFunctionDerivatives(rPoint, grad_phi);
    rReturnValue = prod(trans(rInverseJacobian), grad_phi);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::AssembleOnElement(
    Element<ELEMENT_DIM,SPACE_DIM>& rElement,
    c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1) >& rAElem,
    c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)>& rBElem)
{
    /**
     * \todo #1320 This assumes that the Jacobian is constant on an element.
     * This is true for linear basis functions, but not for any other type of
     * basis function.
     */
    c_matrix<double, SPACE_DIM, ELEMENT_DIM> jacobian;
    c_matrix<double, ELEMENT_DIM, SPACE_DIM> inverse_jacobian;
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

    const unsigned num_nodes = rElement.GetNumNodes();

    // Allocate memory for the basis functions values and derivative values
    c_vector<double, ELEMENT_DIM+1> phi;
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> grad_phi;

    // Loop over Gauss points
    for (unsigned quad_index=0; quad_index < mpQuadRule->GetNumQuadPoints(); quad_index++)
    {
        const ChastePoint<ELEMENT_DIM>& quad_point = mpQuadRule->rGetQuadPoint(quad_index);

        BasisFunction::ComputeBasisFunctions(quad_point, phi);

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

            if (INTERPOLATION_LEVEL != CARDIAC) // don't even interpolate X if cardiac problem
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

                    if (INTERPOLATION_LEVEL==NONLINEAR) // don't need to construct grad_phi or grad_u in other cases
                    {
                        for (unsigned j=0; j<SPACE_DIM; j++)
                        {
                            grad_u(index_of_unknown,j) += grad_phi(j,i)*u_at_node;
                        }
                    }
                }
            }

            // Allow the concrete version of the assembler to interpolate any desired quantities
            this->IncrementInterpolatedQuantities(phi(i), p_node);
            if (this->mAssembleMatrix || INTERPOLATION_LEVEL==NONLINEAR)
            {
                this->IncrementInterpolatedGradientQuantities(grad_phi, i, p_node);
            }
        }

        double wJ = jacobian_determinant * mpQuadRule->GetWeight(quad_index);

        // Create rAElem and rBElem
        if (this->mAssembleMatrix)
        {
            noalias(rAElem) += ComputeMatrixTerm(phi, grad_phi, x, u, grad_u, &rElement) * wJ;
        }

        if (this->mAssembleVector)
        {
            noalias(rBElem) += ComputeVectorTerm(phi, grad_phi, x, u, grad_u, &rElement) * wJ;
        }
    }
}


#endif /*ABSTRACTFEVOLUMEINTEGRALASSEMBLER_HPP_*/
