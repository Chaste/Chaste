/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTSTOKESFLOWSOLVERLONG_HPP_
#define TESTSTOKESFLOWSOLVERLONG_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "StokesFlowAssembler.hpp"
#include "StokesFlowSolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraticMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "Warnings.hpp"
#include "NumericFileComparison.hpp"


class TestStokesFlowSolverLong : public CxxTest::TestSuite
{
public:
    // It is also slow because of #2083 - can be made non-nightly once #2083 is done
    void TestStokesWithLidDrivenCavity3d() throw(Exception)
    {
        unsigned num_elem = 5;
        QuadraticMesh<3> mesh(1.0/num_elem, 1.0, 1.0, 1.0);

        // Dynamic viscosity
        double mu = 1.0;

        // Boundary flow
        std::vector<unsigned> dirichlet_nodes;
        std::vector<c_vector<double,3> > dirichlet_flow;

        for ( TetrahedralMesh<3,3>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
              iter != mesh.GetBoundaryNodeIteratorEnd();
              ++iter)
        {
            double x = (*iter)->rGetLocation()[0];
            double y = (*iter)->rGetLocation()[1];
            double z = (*iter)->rGetLocation()[2];

            c_vector<double,3> flow = zero_vector<double>(3);

            if (fabs(z-1.0)<1e-6)
            {
                flow(0) = x*(1-x)*y*(1-y);
            }

            dirichlet_nodes.push_back((*iter)->GetIndex());
            dirichlet_flow.push_back(flow);
        }

        StokesFlowProblemDefinition<3> problem_defn(mesh);
        problem_defn.SetViscosity(mu);
        problem_defn.SetPrescribedFlowNodes(dirichlet_nodes, dirichlet_flow);

        StokesFlowSolver<3> solver(mesh, problem_defn, "LidDrivenCavityStokesFlow3d");

//        // Uncomment to make errors smaller
//        solver.SetKspAbsoluteTolerance(1e-10);

        solver.Solve();

        std::vector<c_vector<double,3> >& r_solution = solver.rGetVelocities();

        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            double z = mesh.GetNode(i)->rGetLocation()[2];
            if((fabs(x-0.6)<1e-6) && (fabs(y-0.6)<1e-6) && (fabs(z-0.6)<1e-6))
            {
            	TS_ASSERT_DELTA(r_solution[i](0), -8.5990e-03, 1e-4);
            	TS_ASSERT_DELTA(r_solution[i](1),  1.3326e-04, 1e-5);
            	TS_ASSERT_DELTA(r_solution[i](2), -5.1343e-03, 1e-4);
            }
        }
    }
};

#endif // TESTSTOKESFLOWSOLVER_HPP_
