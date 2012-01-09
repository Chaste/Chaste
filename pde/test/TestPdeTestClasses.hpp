/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef TESTPDETESTCLASSES_HPP_
#define TESTPDETESTCLASSES_HPP_

#include <cxxtest/TestSuite.h>

#include "Node.hpp"

#include "HeatEquation.hpp"
#include "NonlinearEquationPde.hpp"
#include "SimplePoissonEquation.hpp"
#include "VaryingDiffusionAndSourceTermPde.hpp"

class TestPdeTestClasses : public CxxTest::TestSuite
{
public:
    void TestSimplePoissonEquation()
    {
        ChastePoint<1> zero1(0);
        ChastePoint<2> zero2(0,0);
        ChastePoint<3> zero3(0,0,0);

        SimplePoissonEquation<1,1> heat_equation1;
        SimplePoissonEquation<2,2> heat_equation2;
        SimplePoissonEquation<3,3> heat_equation3;

        // Diffusion matrices should be equal to identity
        c_matrix<double,1,1> diff1 = heat_equation1.ComputeDiffusionTerm(zero1);
        c_matrix<double,2,2> diff2 = heat_equation2.ComputeDiffusionTerm(zero2);
        c_matrix<double,3,3> diff3 = heat_equation3.ComputeDiffusionTerm(zero3);

        TS_ASSERT_DELTA(diff1(0,0),1,1e-12);

        TS_ASSERT_DELTA(diff2(0,0),1,1e-12);
        TS_ASSERT_DELTA(diff2(1,1),1,1e-12);
        TS_ASSERT_DELTA(diff2(0,1),0,1e-12);

        TS_ASSERT_DELTA(diff3(0,0),1,1e-12);
        TS_ASSERT_DELTA(diff3(1,1),1,1e-12);
        TS_ASSERT_DELTA(diff3(2,2),1,1e-12);
        TS_ASSERT_DELTA(diff3(0,1),0,1e-12);
        TS_ASSERT_DELTA(diff3(0,2),0,1e-12);
        TS_ASSERT_DELTA(diff3(1,2),0,1e-12);

        Node<2> zero(0);
        SimplePoissonEquation<2,2> heat_equation;
        TS_ASSERT_DELTA(heat_equation.ComputeConstantInUSourceTermAtNode(zero), 1.0, 1e-12);
    }

    void TestNonlinearEquationPde()
    {
        ChastePoint<1> zero1(0);
        ChastePoint<2> zero2(0,0);
        ChastePoint<3> zero3(0,0,0);
        double u = 2.0;

        NonlinearEquationPde<1> heat_equation1;
        NonlinearEquationPde<2> heat_equation2;
        NonlinearEquationPde<3> heat_equation3;

        TS_ASSERT_DELTA(heat_equation1.ComputeNonlinearSourceTerm(zero1,u),0.0,1e-12);
        TS_ASSERT_DELTA(heat_equation2.ComputeNonlinearSourceTerm(zero2,u),0.0,1e-12);
        TS_ASSERT_DELTA(heat_equation3.ComputeNonlinearSourceTerm(zero3,u),0.0,1e-12);

        // Diffusion matrices should be equal to identity * u;
        c_matrix<double, 1, 1> diff1 = heat_equation1.ComputeDiffusionTerm(zero1,u);
        c_matrix<double, 2, 2> diff2 = heat_equation2.ComputeDiffusionTerm(zero2,u);
        c_matrix<double, 3, 3> diff3 = heat_equation3.ComputeDiffusionTerm(zero3,u);

        TS_ASSERT_DELTA(diff1(0,0),u,1e-12);

        TS_ASSERT_DELTA(diff2(0,0),u,1e-12);
        TS_ASSERT_DELTA(diff2(1,1),u,1e-12);
        TS_ASSERT_DELTA(diff2(0,1),0,1e-12);

        TS_ASSERT_DELTA(diff3(0,0),u,1e-12);
        TS_ASSERT_DELTA(diff3(1,1),u,1e-12);
        TS_ASSERT_DELTA(diff3(2,2),u,1e-12);
        TS_ASSERT_DELTA(diff3(0,1),0,1e-12);
        TS_ASSERT_DELTA(diff3(0,2),0,1e-12);
        TS_ASSERT_DELTA(diff3(1,2),0,1e-12);
    }

    void TestHeatEquation()
    {
        ChastePoint<1> zero1(0);
        ChastePoint<2> zero2(0,0);
        ChastePoint<3> zero3(0,0,0);
        double u = 2.0;

        HeatEquation<1> pde1;
        HeatEquation<2> pde2;
        HeatEquation<3> pde3;

        TS_ASSERT_DELTA(pde1.ComputeSourceTerm(zero1,u), 0.0, 1e-12);
        TS_ASSERT_DELTA(pde2.ComputeSourceTerm(zero2,u), 0.0, 1e-12);
        TS_ASSERT_DELTA(pde3.ComputeSourceTerm(zero3,u), 0.0, 1e-12);

        TS_ASSERT_DELTA(pde1.ComputeDuDtCoefficientFunction(zero1), 1.0, 1e-12);
        TS_ASSERT_DELTA(pde2.ComputeDuDtCoefficientFunction(zero2), 1.0, 1e-12);
        TS_ASSERT_DELTA(pde3.ComputeDuDtCoefficientFunction(zero3), 1.0, 1e-12);

        // Diffusion matrices should be equal to identity
        c_matrix<double, 1, 1> diff1 = pde1.ComputeDiffusionTerm(zero1);
        c_matrix<double, 2, 2> diff2 = pde2.ComputeDiffusionTerm(zero2);
        c_matrix<double, 3, 3> diff3 = pde3.ComputeDiffusionTerm(zero3);

        TS_ASSERT_DELTA(diff1(0,0), 1, 1e-12);

        TS_ASSERT_DELTA(diff2(0,0), 1, 1e-12);
        TS_ASSERT_DELTA(diff2(1,1), 1, 1e-12);
        TS_ASSERT_DELTA(diff2(0,1), 0, 1e-12);

        TS_ASSERT_DELTA(diff3(0,0), 1, 1e-12);
        TS_ASSERT_DELTA(diff3(1,1), 1, 1e-12);
        TS_ASSERT_DELTA(diff3(2,2), 1, 1e-12);
        TS_ASSERT_DELTA(diff3(0,1), 0, 1e-12);
        TS_ASSERT_DELTA(diff3(0,2), 0, 1e-12);
        TS_ASSERT_DELTA(diff3(1,2), 0, 1e-12);

        Node<1> node(0, zero1);
        TS_ASSERT_DELTA(pde1.ComputeSourceTermAtNode(node,u), 0.0, 1e-12);
    }

    void TestVaryingPde1D()
    {
        VaryingDiffusionAndSourceTermPde<1> pde;
        ChastePoint<1> evaluation_point(2);
        TS_ASSERT_EQUALS(pde.ComputeConstantInUSourceTerm(evaluation_point, NULL),8.0);
        c_matrix<double, 1, 1> diffusion_term = pde.ComputeDiffusionTerm(evaluation_point);
        TS_ASSERT_EQUALS(diffusion_term(0,0), 4);
    }

    void TestVaryingPde2D()
    {
        VaryingDiffusionAndSourceTermPde<2> pde;
        ChastePoint<2> evaluation_point(3,4);
        TS_ASSERT_EQUALS(pde.ComputeConstantInUSourceTerm(evaluation_point, NULL),125.0);
        c_matrix<double, 2, 2> diffusion_term = pde.ComputeDiffusionTerm(evaluation_point);
        TS_ASSERT_EQUALS(diffusion_term(0,0), 25.0);
        TS_ASSERT_EQUALS(diffusion_term(0,1), 0.0);
        TS_ASSERT_EQUALS(diffusion_term(1,0), 0.0);
        TS_ASSERT_EQUALS(diffusion_term(1,1), 25.0);
    }
};

#endif // TESTPDETESTCLASSES_HPP_
