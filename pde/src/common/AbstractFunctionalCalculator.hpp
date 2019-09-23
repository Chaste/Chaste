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

#ifndef ABSTRACTFUNCTIONALCALCULATOR_HPP_
#define ABSTRACTFUNCTIONALCALCULATOR_HPP_

#include "LinearBasisFunction.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "ReplicatableVector.hpp"

/**
 * This is an abstract class for computing user-defined integral-based
 * functionals of a solution on the mesh. The user needs to define
 * GetIntegrand() in the concrete class, and this class can then be used
 * to compute the integral of f(x,u,grad_u) over the mesh, where x is
 * position, u is the solution at x (possibly with multiple components, for
 * which PROBLEM_DIM>1), grad_u the gradient of u and f the integrand as
 * defined in the concrete class.
 *
 * Note linear basis functions are the default, but that the functional may be
 * non-polynomial.  The order of integration used is the highest degree available
 * (3rd order Gaussian quadrature)
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractFunctionalCalculator
{
private:

    /** Replicated store of the solution vector. */
    ReplicatableVector mSolutionReplicated;

    /**
     * @return the integrand. Must be defined by the user.
     *
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     */
    virtual double GetIntegrand(ChastePoint<SPACE_DIM>& rX,
                                c_vector<double,PROBLEM_DIM>& rU,
                                c_matrix<double,PROBLEM_DIM,SPACE_DIM>& rGradU)=0;

    /**
     * Whether we should not calculate the functional on this element for any reason
     *
     * @param rElement  the element of interest
     * @return  whether we should skip this element.
     */
    virtual bool ShouldSkipThisElement(Element<ELEMENT_DIM,SPACE_DIM>& rElement);

public:

    /**
     * Destructor.
     */
    virtual ~AbstractFunctionalCalculator()
    {
    }

    /**
     * @return calculated integral over the given mesh, using the given solution
     * vector on the mesh.
     *
     * Note that, in parallel, this method uses a collective reduction step and
     * should therefore always be called collectively.
     *
     * @param rMesh  The mesh
     * @param solution  The solution vector
     */
    double Calculate(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh, Vec solution);

    /**
     * @return computed the contribution to the integral from one element.
     *
     * @param rElement The element
     */
    double CalculateOnElement(Element<ELEMENT_DIM,SPACE_DIM>& rElement);
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractFunctionalCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::CalculateOnElement(Element<ELEMENT_DIM, SPACE_DIM>& rElement)
{
    double result_on_element = 0;

    // Third order quadrature.  Note that the functional may be non-polynomial (see documentation of class).
    GaussianQuadratureRule<ELEMENT_DIM> quad_rule(3);

    /// NOTE: This assumes that the Jacobian is constant on an element, ie
    /// no curvilinear bases were used for position
    double jacobian_determinant;
    c_matrix<double, SPACE_DIM, ELEMENT_DIM> jacobian;
    c_matrix<double, ELEMENT_DIM, SPACE_DIM> inverse_jacobian;
    rElement.CalculateInverseJacobian(jacobian, jacobian_determinant, inverse_jacobian);

    const unsigned num_nodes = rElement.GetNumNodes();

    // Loop over Gauss points
    for (unsigned quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
    {
        const ChastePoint<ELEMENT_DIM>& quad_point = quad_rule.rGetQuadPoint(quad_index);

        c_vector<double, ELEMENT_DIM+1> phi;
        LinearBasisFunction<ELEMENT_DIM>::ComputeBasisFunctions(quad_point, phi);
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> grad_phi;
        LinearBasisFunction<ELEMENT_DIM>::ComputeTransformedBasisFunctionDerivatives(quad_point, inverse_jacobian, grad_phi);

        // Location of the Gauss point in the original element will be stored in x
        ChastePoint<SPACE_DIM> x(0,0,0);
        c_vector<double,PROBLEM_DIM> u = zero_vector<double>(PROBLEM_DIM);
        c_matrix<double,PROBLEM_DIM,SPACE_DIM> grad_u = zero_matrix<double>(PROBLEM_DIM,SPACE_DIM);

        for (unsigned i=0; i<num_nodes; i++)
        {
            const c_vector<double, SPACE_DIM>& r_node_loc = rElement.GetNode(i)->rGetLocation();

            // Interpolate x
            x.rGetLocation() += phi(i)*r_node_loc;

            // Interpolate u and grad u
            unsigned node_global_index = rElement.GetNodeGlobalIndex(i);
            for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
            {
                // NOTE - following assumes that, if say there are two unknowns u and v, they
                // are stored in the current solution vector as
                // [U1 V1 U2 V2 ... U_n V_n]
                unsigned index_into_vec = PROBLEM_DIM*node_global_index + index_of_unknown;

                double u_at_node = mSolutionReplicated[index_into_vec];
                u(index_of_unknown) += phi(i)*u_at_node;
                // NB. grad_u is PROBLEM_DIM x SPACE_DIM but grad_phi is ELEMENT_DIM x (ELEMENT_DIM+1)
                // Assume here that SPACE_DIM == ELEMENT_DIM and assert it in calling function
                for (unsigned j=0; j<ELEMENT_DIM; j++)
                {
                    grad_u(index_of_unknown,j) += grad_phi(j,i)*u_at_node;
                }
            }
        }

        double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);
        result_on_element += GetIntegrand(x, u, grad_u) * wJ;
    }

    return result_on_element;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractFunctionalCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::Calculate(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh, Vec solution)
{
    assert(ELEMENT_DIM == SPACE_DIM);    // LCOV_EXCL_LINE
    assert(solution);
    mSolutionReplicated.ReplicatePetscVector(solution);
    if (mSolutionReplicated.GetSize() != rMesh.GetNumNodes() * PROBLEM_DIM)
    {
        EXCEPTION("The solution size does not match the mesh");
    }

    double local_result = 0;

    try
    {
        for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = rMesh.GetElementIteratorBegin();
             iter != rMesh.GetElementIteratorEnd();
             ++iter)
        {
            if (rMesh.CalculateDesignatedOwnershipOfElement((*iter).GetIndex()) == true && !ShouldSkipThisElement(*iter))
            {
                local_result += CalculateOnElement(*iter);
            }
        }
    }
    catch (Exception &exception_in_integral)
    {
        PetscTools::ReplicateException(true);
        throw exception_in_integral;
    }
    PetscTools::ReplicateException(false);

    double final_result;
    MPI_Allreduce(&local_result, &final_result, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    return final_result;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractFunctionalCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ShouldSkipThisElement(Element<ELEMENT_DIM,SPACE_DIM>& rElement)
{
    return false;
}

#endif /*ABSTRACTFUNCTIONALCALCULATOR_HPP_*/
