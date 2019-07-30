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

#ifndef _TESTQUADRATICBASISFUNCTION_HPP_
#define _TESTQUADRATICBASISFUNCTION_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include "QuadraticBasisFunction.hpp"
#include "BasisFunctionsCheckers.hpp"
#include "Element.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestQuadraticBasisFunction : public CxxTest::TestSuite
{
public:

    void TestQuadraticBasisFunction0d()
    {
        ChastePoint<0> zero;
        QuadraticBasisFunction<0> basis_func;
        TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 0), 1.0, 1e-12);

        c_vector<double, 1> basis_functions;
        basis_func.ComputeBasisFunctions(zero, basis_functions);
        TS_ASSERT_DELTA(basis_functions(0), 1.0, 1e-12);
    }

    void TestQuadraticBasisFunction1d()
    {
        std::vector<ChastePoint<1>*> evaluation_points;
        ChastePoint<1> zero(0);
        ChastePoint<1> one(1);
        ChastePoint<1> half(0.5);
        evaluation_points.push_back(&zero);
        evaluation_points.push_back(&one);
        evaluation_points.push_back(&half);

        QuadraticBasisFunction<1> basis_func;

        BasisFunctionsCheckers<1> checker;
        checker.checkQuadraticBasisFunctions(evaluation_points);

        // Derivatives
        c_matrix<double, 1, 3> derivatives;
        basis_func.ComputeBasisFunctionDerivatives(zero, derivatives);
        TS_ASSERT_DELTA(derivatives(0,0), -3.0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1),  -1.0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,2),  4.0, 1e-12);

        basis_func.ComputeBasisFunctionDerivatives(one, derivatives);
        TS_ASSERT_DELTA(derivatives(0,0), 1.0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1),  3.0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,2),  -4.0, 1e-12);

        basis_func.ComputeBasisFunctionDerivatives(half, derivatives);
        TS_ASSERT_DELTA(derivatives(0,0), -1.0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1),  1.0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,2),  0.0, 1e-12);
    }

    void TestQuadraticBasisFunction2d()
    {
        ChastePoint<2> zero(0,0);
        ChastePoint<2> onezero(1,0);
        ChastePoint<2> zeroone(0,1);
        ChastePoint<2> halfzero(0.5,0);
        ChastePoint<2> zerohalf(0,0.5);
        ChastePoint<2> halfhalf(0.5,0.5);

        std::vector<ChastePoint<2>*> evaluation_points;
        evaluation_points.push_back(&zero);
        evaluation_points.push_back(&onezero);
        evaluation_points.push_back(&zeroone);
        evaluation_points.push_back(&halfhalf);
        evaluation_points.push_back(&zerohalf);
        evaluation_points.push_back(&halfzero);

        QuadraticBasisFunction<2> basis_func;

        BasisFunctionsCheckers<2> checker;
        checker.checkQuadraticBasisFunctions(evaluation_points);

        // Derivatives
        c_matrix<double, 2, 6> derivatives;
        basis_func.ComputeBasisFunctionDerivatives(onezero, derivatives);
        TS_ASSERT_DELTA(derivatives(0,0), 1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1), 3, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,2), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,3), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,4), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,5), -4, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,0), 1, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,1), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,2), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,3), 4, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,4), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,5), -4, 1e-12);

        basis_func.ComputeBasisFunctionDerivatives(zero, derivatives);
        TS_ASSERT_DELTA(derivatives(0,0), -3, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,2), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,3), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,4), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,5), 4, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,0), -3, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,1), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,2), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,3), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,4), 4, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,5), 0, 1e-12);

        basis_func.ComputeBasisFunctionDerivatives(zeroone, derivatives);
        TS_ASSERT_DELTA(derivatives(0,0), 1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,2), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,3), 4, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,4), -4, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,5), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,0), 1, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,1), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,2), 3, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,3), 0, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,4), -4, 1e-12);
        TS_ASSERT_DELTA(derivatives(1,5), 0, 1e-12);
    }

    void TestQuadraticBasisFunction3d()
    {
        ChastePoint<3> zero(0,0,0);
        ChastePoint<3> zerozeroone(0,0,1);
        ChastePoint<3> zeroonezero(0,1,0);
        ChastePoint<3> onezerozero(1,0,0);
        ChastePoint<3> halfhalfzero(0.5,0.5,0);
        ChastePoint<3> halfzerozero(0.5,0,0);
        ChastePoint<3> halfzerohalf(0.5,0.0,0.5);
        ChastePoint<3> zerohalfhalf(0.0,0.5,0.5);
        ChastePoint<3> zerohalfzero(0.0,0.5,0.0);
        ChastePoint<3> zerozerohalf(0.0,0.0,0.5);

        std::vector<ChastePoint<3>*> evaluation_points;
        evaluation_points.push_back(&zero);
        evaluation_points.push_back(&onezerozero);
        evaluation_points.push_back(&zeroonezero);
        evaluation_points.push_back(&zerozeroone);
        evaluation_points.push_back(&halfzerozero);
        evaluation_points.push_back(&halfhalfzero);
        evaluation_points.push_back(&zerohalfzero);
        evaluation_points.push_back(&zerozerohalf);
        evaluation_points.push_back(&halfzerohalf);
        evaluation_points.push_back(&zerohalfhalf);

        QuadraticBasisFunction<3> basis_func;

        BasisFunctionsCheckers<3> checker;
        checker.checkQuadraticBasisFunctions(evaluation_points);

        // Derivatives
        c_matrix<double, 3, 10> derivatives;
        basis_func.ComputeBasisFunctionDerivatives(onezerozero, derivatives);
        TS_ASSERT_DELTA(derivatives(0, 0),  1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0, 1),  3, 1e-12);
        TS_ASSERT_DELTA(derivatives(0, 2),  0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0, 3),  0, 1e-12);
    }

    void TestComputeTransformedQuadraticBasisFunctionDerivatives1d()
    {
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0, false, 3.0));
        nodes.push_back(new Node<1>(1, false, 5.0));

        Element<1,1> element(INDEX_IS_NOT_USED, nodes);
        QuadraticBasisFunction<1> basis_function;

        c_matrix<double,1,1> jacobian, inverse_jacobian;
        double determinant;
        element.CalculateInverseJacobian(jacobian, determinant, inverse_jacobian);
        ChastePoint<1> evaluation_point(0.2);
        c_matrix<double, 1, 3> trans_deriv;
        basis_function.ComputeTransformedBasisFunctionDerivatives(evaluation_point,
                                                                  inverse_jacobian,
                                                                  trans_deriv);

        TS_ASSERT_DELTA(trans_deriv(0,0), -1.1, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,1), -0.1, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,2),  1.2, 1e-12);

        // Free nodes
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestComputeTransformedQuadraticBasisFunction2d()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 4.0, 3.0));
        nodes.push_back(new Node<2>(1, false, 6.0, 4.0));
        nodes.push_back(new Node<2>(2, false, 3.0, 5.0));
//        nodes.push_back(new Node<2>(3, false, 5.0, 3.5)); // extra nodes
//        nodes.push_back(new Node<2>(4, false, 3.5, 4.0));
//        nodes.push_back(new Node<2>(5, false, 4.5, 4.5));
        Element<2,2> element(INDEX_IS_NOT_USED, nodes);
        QuadraticBasisFunction<2> basis_function;

        c_matrix<double,2,2> jacobian, inverse_jacobian;
        double determinant;
        element.CalculateInverseJacobian(jacobian, determinant, inverse_jacobian);
        ChastePoint<2> evaluation_point(0.3, 0.6);
        c_matrix<double,2,6> trans_deriv;
        basis_function.ComputeTransformedBasisFunctionDerivatives(evaluation_point,
                                                                  inverse_jacobian,
                                                                  trans_deriv);

        TS_ASSERT_DELTA(trans_deriv(0,0),0.12, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,0),0.36, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,1),0.08, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,1),0.04, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,2),-0.28, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,2),0.56, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,3),0.72, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,3),0.96, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,4),-0.56, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,4),-1.28, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,5),-0.08, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,5),-0.64, 1e-12);

        // Free nodes
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestComputeTransformedQuadraticBasisFunction3d()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 4.0, 3.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 6.0, 4.0, 1.0));
        nodes.push_back(new Node<3>(2, false, 3.0, 5.0, 2.0));
        nodes.push_back(new Node<3>(3, false, 5.0, 4.0, 3.0));
//        nodes.push_back(new Node<3>(4, false, 5.0, 3.5, 0.5));
//        nodes.push_back(new Node<3>(5, false, 3.5, 4.0, 1.0));
//        nodes.push_back(new Node<3>(6, false, 4.5, 3.5, 1.5));
//        nodes.push_back(new Node<3>(7, false, 4.5, 4.5, 1.5));
//        nodes.push_back(new Node<3>(8, false, 5.5, 4.0, 2.0));
//        nodes.push_back(new Node<3>(9, false, 4.0, 4.5, 2.5));
        Element<3,3> element(INDEX_IS_NOT_USED, nodes);
        QuadraticBasisFunction<3> basis_function;

        c_matrix<double,3,3> jacobian, inverse_jacobian;
        double determinant;
        element.CalculateInverseJacobian(jacobian, determinant, inverse_jacobian);
        ChastePoint<3> evaluation_point(0.3, 0.1, 0.2);
        c_matrix<double,3,10> trans_deriv;
        basis_function.ComputeTransformedBasisFunctionDerivatives(evaluation_point,
                                                                  inverse_jacobian,
                                                                  trans_deriv);

        TS_ASSERT_DELTA(trans_deriv(0,0),-0.12, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,0),-0.3, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(2,0),-0.06, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,1),0.08, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,1),0.1, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(2,1),-0.06, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,2),0.12, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,2),-0.3, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(2,2),0.06, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,3),0.0, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,3),0.1, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(2,3),-0.1, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,4),0.4, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,4),0.2, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(2,4),-0.6, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,5),-0.08, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,5),0.8, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(2,5),-0.24, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,6),-0.4, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,6),0.6, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(2,6),-0.2, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,7),-0.16, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,7),-1.2, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(2,7),0.72, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,8),0.32, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,8),-0.2, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(2,8),0.36, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,9),-0.16, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,9),0.2, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(2,9),0.12, 1e-12);

        // Free nodes
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif //_TESTQUADRATICBASISFUNCTION_HPP_
