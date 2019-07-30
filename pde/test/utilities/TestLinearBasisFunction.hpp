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

#ifndef _TESTLINEARBASISFUNCTION_HPP_
#define _TESTLINEARBASISFUNCTION_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>

#include "LinearBasisFunction.hpp"
#include "GaussianQuadratureRule.hpp"
#include "BasisFunctionsCheckers.hpp"
#include "Element.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestLinearBasisFunction : public CxxTest::TestSuite
{
public:

    void TestLinearBasisFunction0d()
    {
        ChastePoint<0> zero;
        TS_ASSERT_DELTA(LinearBasisFunction<0>::ComputeBasisFunction(zero, 0), 1.0, 1e-12);

        c_vector<double, 1> basis_function_vector;
        LinearBasisFunction<0>::ComputeBasisFunctions(zero,basis_function_vector);
        TS_ASSERT_EQUALS(basis_function_vector.size(), 1u);
        TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);

        // Check link with 0d quad rule works ok
        GaussianQuadratureRule<0>  quad_rule(0);
        const ChastePoint<0>& quad_point = quad_rule.rGetQuadPoint(0);

        c_vector<double, 1> basis_function_vector2;
        LinearBasisFunction<0>::ComputeBasisFunctions(quad_point,basis_function_vector2);
        TS_ASSERT_EQUALS(basis_function_vector.size(), 1u);
        TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);
    }

    void TestLinearBasisFunction1d()
    {
        ChastePoint<1> zero(0);
        ChastePoint<1> one(1);

        std::vector<ChastePoint<1>*> evaluation_points;
        evaluation_points.push_back(&zero);
        evaluation_points.push_back(&one);

        BasisFunctionsCheckers<1> checker;
        checker.checkLinearBasisFunctions(evaluation_points);

        // Derivatives
        c_matrix<double, 1, 2> derivatives;
        LinearBasisFunction<1>::ComputeBasisFunctionDerivatives(one,derivatives);
        TS_ASSERT_DELTA(derivatives(0,0), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1),  1, 1e-12);
    }

    void TestLinearBasisFunction2d()
    {
        ChastePoint<2> zero(0,0);
        ChastePoint<2> onezero(1,0);
        ChastePoint<2> zeroone(0,1);

        std::vector<ChastePoint<2>*> evaluation_points;
        evaluation_points.push_back(&zero);
        evaluation_points.push_back(&onezero);
        evaluation_points.push_back(&zeroone);

        BasisFunctionsCheckers<2> checker;
        checker.checkLinearBasisFunctions(evaluation_points);

        // Derivatives
        c_matrix<double, 2, 3> derivatives;
        LinearBasisFunction<2>::ComputeBasisFunctionDerivatives(onezero,derivatives);
        TS_ASSERT_DELTA(derivatives(0,0), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1),  1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,2),  0, 1e-12);
    }

    void TestLinearBasisFunction3d()
    {
        ChastePoint<3> zero(0,0,0);
        ChastePoint<3> zerozeroone(0,0,1);
        ChastePoint<3> zeroonezero(0,1,0);
        ChastePoint<3> onezerozero(1,0,0);

        std::vector<ChastePoint<3>*> evaluation_points;
        evaluation_points.push_back(&zero);
        evaluation_points.push_back(&onezerozero);
        evaluation_points.push_back(&zeroonezero);
        evaluation_points.push_back(&zerozeroone);

        BasisFunctionsCheckers<3> checker;
        checker.checkLinearBasisFunctions(evaluation_points);

        // Derivatives
        c_matrix<double, 3, 4> derivatives;
        LinearBasisFunction<3>::ComputeBasisFunctionDerivatives(onezerozero,derivatives);
        TS_ASSERT_DELTA(derivatives(0,0), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1),  1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,2),  0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,3),  0, 1e-12);
    }

    void TestComputeTransformedBasisFunctionDerivatives()
    {
        // 1D
        ChastePoint<1> one(1);

        c_matrix<double, 1, 1> inv_J;
        inv_J(0,0) = 0.5;

        c_matrix<double, 1, 2> trans_deriv;
        LinearBasisFunction<1>::ComputeTransformedBasisFunctionDerivatives(one, inv_J, trans_deriv);

        TS_ASSERT_DELTA(trans_deriv(0,0), -0.5, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,1),  0.5, 1e-12);

        // 2D
        ChastePoint<2> oneone(1,1);

        c_matrix<double, 2, 2> inv_J2 = 0.5 * identity_matrix<double>(2);

        c_matrix<double, 2, 3> trans_deriv2;
        LinearBasisFunction<2>::ComputeTransformedBasisFunctionDerivatives(oneone, inv_J2, trans_deriv2);

        TS_ASSERT_DELTA(trans_deriv2(0,0), -0.5, 1e-12);
        TS_ASSERT_DELTA(trans_deriv2(0,1),  0.5, 1e-12);
        TS_ASSERT_DELTA(trans_deriv2(0,2),    0, 1e-12);

        // 3D
        ChastePoint<3> oneoneone(1,1,1);

        c_matrix<double, 3, 3> inv_J3 = 0.5 * identity_matrix<double>(3);

        c_matrix<double, 3, 4> trans_deriv3;
        LinearBasisFunction<3>::ComputeTransformedBasisFunctionDerivatives(oneoneone, inv_J3, trans_deriv3);

        TS_ASSERT_DELTA(trans_deriv3(0,0), -0.5, 1e-12);
        TS_ASSERT_DELTA(trans_deriv3(0,1),  0.5, 1e-12);
        TS_ASSERT_DELTA(trans_deriv3(0,2),    0, 1e-12);
        TS_ASSERT_DELTA(trans_deriv3(0,3),    0, 1e-12);
    }

    void TestComputeTransformedBasisFunction2()
    {
        // 2D - with better test data
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 4.0, 3.0));
        nodes.push_back(new Node<2>(1, false, 6.0, 4.0));
        nodes.push_back(new Node<2>(2, false, 3.0, 5.0));
        Element<2,2> element(INDEX_IS_NOT_USED, nodes);

        c_matrix<double, 2, 2> jacobian, inverse_jacobian;
        double determinant;
        element.CalculateInverseJacobian(jacobian, determinant, inverse_jacobian);
        ChastePoint<2> evaluation_point(1,1);
        c_matrix<double, 2, 3> trans_deriv;
        LinearBasisFunction<2>::ComputeTransformedBasisFunctionDerivatives(evaluation_point,
                                                                           inverse_jacobian,
                                                                           trans_deriv);

        TS_ASSERT_DELTA(trans_deriv(0,0),-0.2, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,0),-0.6, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,1),0.4, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,1),0.2, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,2),-0.2, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,2),0.4, 1e-12);

        delete nodes[0];
        delete nodes[1];
        delete nodes[2];
    }
};

#endif //_TESTLINEARBASISFUNCTION_HPP_
