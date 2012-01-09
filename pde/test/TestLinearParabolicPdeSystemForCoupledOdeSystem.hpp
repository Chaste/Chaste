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

#ifndef TESTLINEARPARABOLICPDESYSTEMFORCOUPLEDODESYSTEM_HPP_
#define TESTLINEAPARABOLICRPDESYSTEMFORCOUPLEDODESYSTEM_HPP_

/**
 * TestSimpleLinearParabolicSystemForCoupledOdeSystem.hpp
 *
 * Test suite for the TestSimpleLinearParabolicSystemForCoupledOdeSystem class.
 */
#include <cxxtest/TestSuite.h>
#include "HeatEquationForCoupledOdeSystem.hpp"
#include "SchnackenbergCoupledPdeSystem.hpp"

class TestLinearParabolicPdeSystemForCoupledOdeSystem : public CxxTest::TestSuite
{
public:

    void TestHeatEquationForCoupledOdeSystem()
    {
        // Create PDE system object
        HeatEquationForCoupledOdeSystem<1> pde;

        ChastePoint<1> x(1.0);

        TS_ASSERT_DELTA(pde.ComputeDuDtCoefficientFunction(x,0), 1.0, 1e-6);

        c_vector<double,1> pde_solution;
        pde_solution(0) = 4.0;

        std::vector<double> ode_solution(1);
        ode_solution[0] = 5.0;

        TS_ASSERT_DELTA(pde.ComputeSourceTerm(x, pde_solution, ode_solution, 0), 0.0, 1e-6);

        Node<1> node(0);
        TS_ASSERT_DELTA(pde.ComputeSourceTermAtNode(node, pde_solution, ode_solution, 0), 0.0, 1e-6);

        c_matrix<double,1,1> diffusion_term = pde.ComputeDiffusionTerm(x, 0);
        TS_ASSERT_DELTA(diffusion_term(0,0), 1.0, 1e-6);
    }

    void TestSchnackenbergCoupledPdeSystem()
    {
        // Create PDE system object
        SchnackenbergCoupledPdeSystem<1> pde(1e-4, 1e-2, 0.1, 0.2, 0.3, 0.1);

        ChastePoint<1> x(1.0);

        TS_ASSERT_DELTA(pde.ComputeDuDtCoefficientFunction(x,0), 1.0, 1e-6);

        c_vector<double,2> pde_solution;
        pde_solution(0) = 2.0;
        pde_solution(1) = 0.75;

        std::vector<double> ode_solution(1);
        ode_solution[0] = 5.0;

        TS_ASSERT_DELTA(pde.ComputeSourceTerm(x, pde_solution, ode_solution, 0), 0.0, 1e-6);

        Node<1> node(0);
        TS_ASSERT_DELTA(pde.ComputeSourceTermAtNode(node, pde_solution, ode_solution, 0), 0.0, 1e-6);

        c_matrix<double,1,1> diffusion_term1 = pde.ComputeDiffusionTerm(x, 0);
        TS_ASSERT_DELTA(diffusion_term1(0,0), 1e-4, 1e-6);

        c_matrix<double,1,1> diffusion_term2 = pde.ComputeDiffusionTerm(x, 1);
        TS_ASSERT_DELTA(diffusion_term2(0,0), 1e-2, 1e-6);
    }
};

#endif /*TESTLINEARPARABOLICPDESYSTEMFORCOUPLEDODESYSTEM_HPP_*/
