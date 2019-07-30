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

#ifndef TESTLINEARPARABOLICPDESYSTEMFORCOUPLEDODESYSTEM_HPP_
#define TESTLINEARPARABOLICPDESYSTEMFORCOUPLEDODESYSTEM_HPP_

/**
 * TestSimpleLinearParabolicSystemForCoupledOdeSystem.hpp
 *
 * Test suite for the TestSimpleLinearParabolicSystemForCoupledOdeSystem class.
 */
#include <cxxtest/TestSuite.h>
#include "HeatEquationForCoupledOdeSystem.hpp"
#include "SchnackenbergCoupledPdeSystem.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

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
