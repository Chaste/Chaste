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

#ifndef _TESTGAUSSIANQUADRATURERULE_HPP_
#define _TESTGAUSSIANQUADRATURERULE_HPP_

#include <cxxtest/TestSuite.h>
#include "GaussianQuadratureRule.hpp"
#include "Node.hpp"
#include "Element.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestGaussianQuadratureRule : public CxxTest::TestSuite
{
public:

    /**
     * Check points and weights are in the right ranges.
     */
    void TestTheGaussianQuadratureRule()
    {
        // 0d
        GaussianQuadratureRule<0> quad_rule1(0);
        TS_ASSERT_EQUALS(quad_rule1.GetNumQuadPoints(),1U);

        GaussianQuadratureRule<0> quad_rule2(1);
        TS_ASSERT_EQUALS(quad_rule2.GetNumQuadPoints(),1U);

        TS_ASSERT_DELTA(quad_rule1.GetWeight(0), 1.0, DBL_EPSILON);

        // 1d
        for (unsigned order=0; order<3; order++)
        {
            GaussianQuadratureRule<1> quad_rule(order);

            for (unsigned i=0; i<quad_rule.GetNumQuadPoints(); i++)
            {
                TS_ASSERT_LESS_THAN_EQUALS( quad_rule.GetWeight(i),1);
                TS_ASSERT_LESS_THAN_EQUALS(0, quad_rule.GetWeight(i));

                TS_ASSERT_LESS_THAN(quad_rule.rGetQuadPoint(i)[0], 1); // x<1
                TS_ASSERT_LESS_THAN(0, quad_rule.rGetQuadPoint(i)[0]); // x>0
            }
        }

        // 2d
        for (unsigned order=0; order<3; order++)
        {
            GaussianQuadratureRule<2> quad_rule(order);

            for (unsigned i=0; i<quad_rule.GetNumQuadPoints(); i++)
            {
                TS_ASSERT_LESS_THAN_EQUALS( quad_rule.GetWeight(i),1);
                TS_ASSERT_LESS_THAN_EQUALS(-quad_rule.GetWeight(i),0);

                TS_ASSERT_LESS_THAN(-(1-quad_rule.rGetQuadPoint(i)[0]
                                      -quad_rule.rGetQuadPoint(i)[1]),0); // 1-x-y>0
                TS_ASSERT_LESS_THAN(-quad_rule.rGetQuadPoint(i)[0],0);  // x>0
                TS_ASSERT_LESS_THAN(-quad_rule.rGetQuadPoint(i)[1],0);  // y>0
            }
        }

        // 3d
        for (unsigned order=0; order<4; order++)
        {
            GaussianQuadratureRule<3> quad_rule(order);

            for (unsigned i=0; i<quad_rule.GetNumQuadPoints(); i++)
            {
                TS_ASSERT_LESS_THAN_EQUALS( quad_rule.GetWeight(i),1);
                TS_ASSERT_LESS_THAN_EQUALS(-quad_rule.GetWeight(i),0);

                TS_ASSERT_LESS_THAN(-(1-quad_rule.rGetQuadPoint(i)[0]
                                      -quad_rule.rGetQuadPoint(i)[1]
                                      -quad_rule.rGetQuadPoint(i)[2]),0); // 1-x-y-z>0
                TS_ASSERT_LESS_THAN(-quad_rule.rGetQuadPoint(i)[0],0);  // x>0
                TS_ASSERT_LESS_THAN(-quad_rule.rGetQuadPoint(i)[1],0);  // y>0
                TS_ASSERT_LESS_THAN(-quad_rule.rGetQuadPoint(i)[2],0);  // z>0
            }
        }

        // Exceptions (unsupported cases)
        TS_ASSERT_THROWS_THIS(GaussianQuadratureRule<1> quad_rule(6),
                "Gauss quadrature order not supported.");
        TS_ASSERT_THROWS_THIS(GaussianQuadratureRule<2> quad_rule(5),
                "Gauss quadrature order not supported.");
        TS_ASSERT_THROWS_THIS(GaussianQuadratureRule<3> quad_rule(4),
                "Gauss quadrature order not supported.");
        TS_ASSERT_THROWS_THIS(GaussianQuadratureRule<4> quad_rule(1),
                "Gauss quadrature rule not available for this dimension.");
    }

    /**
     * Test by integrating polynomials of the form x^n.
     * 1d case.
     */
    void TestGaussianQuadratureRuleIntegralOneD()
    {
        for (unsigned order=0; order<6; order++)
        {
            GaussianQuadratureRule<1> quad_rule(order);

            //Test that the rule can be used to integrate the polynomials 0, x, x^2, x^3 ... on the interval [1,3]
            for (unsigned poly_degree=0; poly_degree<=order; poly_degree++)
            {
                std::vector<Node<1>*> nodes2;
                nodes2.push_back(new Node<1>(0, false, 1.0));
                nodes2.push_back(new Node<1>(1, false, 3.0));
                Element<1,1> element(INDEX_IS_NOT_USED, nodes2);

                double integral=0;
                double jacobian_determinant;
                c_matrix<double, 1, 1> jacobian;
                element.CalculateJacobian(jacobian, jacobian_determinant);

                // This assumes linear basis functions in 1d
                double x1 = element.GetNodeLocation(0,0);
                double x2 = element.GetNodeLocation(1,0);

                for (unsigned quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
                {
                    const ChastePoint<1>& quad_point = quad_rule.rGetQuadPoint(quad_index);
                    const ChastePoint<1> transformed_quad_point =
                        ChastePoint<1>((1-quad_point[0])*x1 + quad_point[0]*x2);

                    double integrand_value = pow(transformed_quad_point[0],(double)poly_degree);

                    integral+= integrand_value*jacobian_determinant
                               *quad_rule.GetWeight(quad_index);
                }

                //The analytic solution is
                //
                // \int_0^3 x^p dx = [(1/(p+1))*x^(p+1)]_1^3
                //
                TS_ASSERT_DELTA(integral, 1.0/(poly_degree+1.0)*(pow(3.0,(double)poly_degree+1)-1), 1e-12); // Largest integral is about 100

                delete nodes2[0];
                delete nodes2[1];
            }
        }
    }

    /**
     * Test by integrating polynomials up to degree 2p-2, where p is the no. of
     * points in each dimension. This is the 2d case.
     *
     * We integrate things like x, y, x^2, xy, y^2, ...
     */
    void TestGaussianQuadratureRuleIntegralTwoD()
    {
        // Expected answers [degree_x][degree_y]
        double expected[5][5] = { {1.0/2,  1.0/6,   1.0/12, 1.0/20, 1.0/30},
                                  {1.0/6,  1.0/24,  1.0/60, 1.0/120, 0},
                                  {1.0/12, 1.0/60,  1.0/180, 0, 0},
                                  {1.0/20, 1.0/120, 0,      0, 0},
                                  {1.0/30, 0.0, 0.0, 0.0, 0.0} };

        for (unsigned order=0; order<4; order++)
        {
            GaussianQuadratureRule<2> quad_rule(order);

            for (unsigned poly_degree_x=0; poly_degree_x<=order; poly_degree_x++)
            {

                for (unsigned poly_degree_y=0; poly_degree_y<=order-poly_degree_x; poly_degree_y++)
                {
                    double integral = 0.0;

                    for (unsigned quad_index=0;
                         quad_index<quad_rule.GetNumQuadPoints();
                         quad_index++)
                    {
                        const ChastePoint<2>& quad_point = quad_rule.rGetQuadPoint(quad_index);

                        integral += pow(quad_point[0], (double)poly_degree_x)
                                    *pow(quad_point[1], (double)poly_degree_y)
                                    *quad_rule.GetWeight(quad_index);
                    }

                    TS_ASSERT_DELTA(integral,
                                    expected[poly_degree_x][poly_degree_y],
                                    1e-15);
                }
            }
        }
    }

    /**
     * Test by integrating polynomials up to degree p=3. This is the 3d case.
     *
     * We integrate things like x, y, z, x^2, y^2, z^2, xy, xz, yz...
     */
    void TestGaussianQuadratureRuleIntegralThreeD()
    {
        // Expected answers [degree_x][degree_y][degree_z]
        // Final row/column/thing is unused (we can't yet do order 4)
        double expected[5][5][5]
        =
            {
                {
                    {1.0/6,      1.0/24,      1.0/60,      1.0/120,      1.0/210},
                    {1.0/24,     1.0/120,     1.0/360,     1.0/840,     1.0/1680},
                    {1.0/60,     1.0/360,     1.0/1260,    1.0/3360,     1.0/7560},
                    {1.0/120,    1.0/840,     1.0/3360,    1.0/10080,    1.0/25200},
                    {1.0/210,    1.0/1680,    1.0/7560,    1.0/25200,    1.0/69300}
                },
                {
                    {1.0/24,    1.0/120,    1.0/360,    1.0/840,   1.0/1680},
                    {1.0/120,   1.0/720,   1.0/2520,   1.0/6720,  1.0/15120},
                    {1.0/360,  1.0/2520,  1.0/10080,  1.0/30240,  1.0/75600},
                    {1.0/840,  1.0/6720,  1.0/30240, 1.0/100800, 1.0/277200},
                    {1.0/1680, 1.0/15120, 1.0/75600, 1.0/277200, 1.0/831600}
                },
                {
                    {1.0/60.0 , 1.0/360.0 , 1.0/1260.0 , 1.0/3360.0 , 1.0/7560.0},
                    {1.0/360.0 , 1.0/2520.0 , 1.0/10080.0 , 1.0/30240.0 , 1.0/75600.0},
                    {1.0/1260.0 , 1.0/10080.0 , 1.0/45360.0 , 1.0/151200.0 , 1.0/415800.0},
                    {1.0/3360.0 , 1.0/30240.0 , 1.0/151200.0 , 1.0/554400.0 , 1.0/1663200.0},
                    {1.0/7560.0 , 1.0/75600.0 , 1.0/415800.0 , 1.0/1663200.0 , 1.0/5405400.0}
                },
                {
                    {1.0/120.0 , 1.0/840.0 , 1.0/3360.0 , 1.0/10080.0 , 1.0/25200.0},
                    {1.0/840.0 , 1.0/6720.0 , 1.0/30240.0 , 1.0/100800.0 , 1.0/277200.0},
                    {1.0/3360.0 , 1.0/30240.0 , 1.0/151200.0 , 1.0/554400.0 , 1.0/1663200.0},
                    {1.0/10080.0 , 1.0/100800.0 , 1.0/554400.0 , 1.0/2217600.0 , 1.0/7207200.0},
                    {1.0/25200.0 , 1.0/277200.0 , 1.0/1663200.0 , 1.0/7207200.0 , 1.0/25225200.0}
                },
                {
                    {1.0/210.0 , 1.0/1680.0 , 1.0/7560.0 , 1.0/2520.0 , 1.0/69300.0},
                    {1.0/1680.0 , 1.0/15120.0 , 1.0/75600.0 , 1.0/277200.0 , 1.0/831600.0},
                    {1.0/7560.0 , 1.0/75600.0 , 1.0/415800.0 , 1.0/1663200.0 , 1.0/5405400.0},
                    {1.0/25200.0 , 1.0/277200.0 , 1.0/1663200.0 ,1.0/7207200.0 ,1.0/25225200.0},
                    {1.0/69300.0 ,1.0/831600.0 ,1.0/5405400.0 ,1.0/25225200.0 ,1.0/94594500.0}
                }
            };

        // Test up to order 3 in 3D
        for (unsigned order=0; order<4; order++)
        {
            GaussianQuadratureRule<3> quad_rule(order);

            for (unsigned poly_degree_x=0;
                 poly_degree_x<=order;
                 poly_degree_x++)
            {

                for (unsigned poly_degree_y=0;
                     poly_degree_y<=order-poly_degree_x;
                     poly_degree_y++)
                {

                    for (unsigned poly_degree_z=0;
                         poly_degree_z<=order-poly_degree_x-poly_degree_y;
                         poly_degree_z++)
                    {
                        double integral = 0.0;

                        for (unsigned quad_index=0;
                             quad_index<quad_rule.GetNumQuadPoints();
                             quad_index++)
                        {
                            const ChastePoint<3>& quad_point = quad_rule.rGetQuadPoint(quad_index);

                            integral += pow(quad_point[0], (double)poly_degree_x)
                                        * pow(quad_point[1], (double)poly_degree_y)
                                        * pow(quad_point[2], (double)poly_degree_z)
                                        * quad_rule.GetWeight(quad_index);
                        }

                        TS_ASSERT_DELTA(integral, expected[poly_degree_x][poly_degree_y][poly_degree_z],
                                        1e-15);
                    }
                }
            }
        }
    }
};

#endif //_TESTGAUSSIANQUADRATURERULE_HPP_
