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

#ifndef _TESTGAUSSIANQUADRATURERULE_HPP_
#define _TESTGAUSSIANQUADRATURERULE_HPP_

#include <cxxtest/TestSuite.h>
#include "GaussianQuadratureRule.hpp"
#include "Node.hpp"
#include "Element.hpp"

class TestGaussianQuadratureRule : public CxxTest::TestSuite
{
public:

    /**
     * Check points and weights are in the right ranges.
     */
    void TestTheGaussianQuadratureRule()
    {
        // 0d
        GaussianQuadratureRule<0> quad_rule1(1);
        TS_ASSERT_EQUALS(quad_rule1.GetNumQuadPoints(),1U);

        GaussianQuadratureRule<0> quad_rule2(2);
        TS_ASSERT_EQUALS(quad_rule2.GetNumQuadPoints(),1U);

        TS_ASSERT_DELTA(quad_rule1.GetWeight(0),1,1e-12);

        // 1d
        for (int number_of_points=1; number_of_points<4; number_of_points++)
        {
            GaussianQuadratureRule<1> quad_rule(number_of_points);

            for (unsigned i=0; i<quad_rule.GetNumQuadPoints(); i++)
            {
                TS_ASSERT_LESS_THAN_EQUALS( quad_rule.GetWeight(i),1);
                TS_ASSERT_LESS_THAN_EQUALS(-quad_rule.GetWeight(i),0);

                TS_ASSERT_LESS_THAN( quad_rule.rGetQuadPoint(i)[0],1); // x<1
                TS_ASSERT_LESS_THAN(-quad_rule.rGetQuadPoint(i)[0],0); // x>0
            }
        }

        // 2d
        for (int number_of_points=1; number_of_points<4; number_of_points++)
        {
            GaussianQuadratureRule<2> quad_rule(number_of_points);

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
        for (int number_of_points=1; number_of_points<5; number_of_points++)
        {
            GaussianQuadratureRule<3> quad_rule(number_of_points);

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
        TS_ASSERT_THROWS_THIS(GaussianQuadratureRule<1> quad_rule(4),
                "Number of gauss points per dimension not supported.");
        TS_ASSERT_THROWS_THIS(GaussianQuadratureRule<2> quad_rule(4),
                "Number of gauss points per dimension not supported.");
        TS_ASSERT_THROWS_THIS(GaussianQuadratureRule<3> quad_rule(5),
                "Number of gauss points per dimension not supported.");
        TS_ASSERT_THROWS_THIS(GaussianQuadratureRule<4> quad_rule(1),
                "Gauss points not available for this dimension.");
    }

    /**
     * Test by integrating polynomials of the form x^n.
     * 1d case.
     */
    void TestGaussianQuadratureRuleIntegralOneD()
    {
        for (int num_quad_points=1; num_quad_points<4; num_quad_points++)
        {
            GaussianQuadratureRule<1> quad_rule(num_quad_points);

            for (int poly_degree=0; poly_degree<2*num_quad_points; poly_degree++)
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

                    double integrand_value = pow(transformed_quad_point[0],poly_degree);

                    integral+= integrand_value*jacobian_determinant
                               *quad_rule.GetWeight(quad_index);
                }

                TS_ASSERT_DELTA(integral, 1.0/(poly_degree+1.0)*(pow(3,poly_degree+1)-1), 1e-7);

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
        double expected[5][5] = { {1.0/2,  1.0/6,  1.0/12, 1.0/20, 1.0/30},
                                  {1.0/6,  1.0/24, 1.0/60, 1.0/120, 0},
                                  {1.0/12, 1.0/60, 1.0/180, 0, 0},
                                  {1.0/20, 1.0/120, 0,      0, 0},
                                  {1.0/30, 0,       0,      0, 0} };

        for (int num_quad_points=1; num_quad_points<4; num_quad_points++)
        {
            GaussianQuadratureRule<2> quad_rule(num_quad_points);

            for (int poly_degree_x=0; poly_degree_x<2*num_quad_points-1;
                 poly_degree_x++)
            {

                for (int poly_degree_y=0;
                     poly_degree_y<2*num_quad_points-1-poly_degree_x;
                     poly_degree_y++)
                {
                    double integral = 0.0;

                    for (unsigned quad_index=0;
                         quad_index<quad_rule.GetNumQuadPoints();
                         quad_index++)
                    {
                        const ChastePoint<2>& quad_point = quad_rule.rGetQuadPoint(quad_index);

                        integral += pow(quad_point[0], poly_degree_x)
                                   * pow(quad_point[1], poly_degree_y)
                                    *quad_rule.GetWeight(quad_index);
                    }

                    TS_ASSERT_DELTA(integral,
                                    expected[poly_degree_x][poly_degree_y],
                                    1e-7);
                }
            }
        }
    }

    /**
     * Test by integrating polynomials up to degree 2p-2, where p is the no. of
     * points in each dimension. This is the 3d case.
     *
     * We integrate things like x, y, z, x^2, y^2, z^2, xy, xz, yz...
     */
    void TestGaussianQuadratureRuleIntegralThreeD()
    {
        // Expected answers [degree_x][degree_y][degree_z]
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

        // Test 2 and 3 quadrature points per dimension in 3D
        for (int num_quad_points=2; num_quad_points<4; num_quad_points++)
        {
            GaussianQuadratureRule<3> quad_rule(num_quad_points);

            for (int poly_degree_x=0; poly_degree_x<2*num_quad_points-1;
                 poly_degree_x++)
            {

                for (int poly_degree_y=0;
                     poly_degree_y<2*num_quad_points-1/*-poly_degree_x*/;
                     poly_degree_y++)
                {

                    for (int poly_degree_z=0;
                         poly_degree_z<2*num_quad_points-1/*-poly_degree_y*/;
                         poly_degree_z++)
                    {
                        double integral = 0.0;

                        for (unsigned quad_index=0;
                             quad_index<quad_rule.GetNumQuadPoints();
                             quad_index++)
                        {
                            const ChastePoint<3>& quad_point = quad_rule.rGetQuadPoint(quad_index);

                            integral += pow(quad_point[0], poly_degree_x)
                                       * pow(quad_point[1], poly_degree_y)
                                       * pow(quad_point[2], poly_degree_z)
                                        *quad_rule.GetWeight(quad_index);
                        }

                        TS_ASSERT_DELTA(integral, expected[poly_degree_x][poly_degree_y][poly_degree_z],
                                        0.01);
                    }
                }
            }
        }
    }
};

#endif //_TESTGAUSSIANQUADRATURERULE_HPP_
