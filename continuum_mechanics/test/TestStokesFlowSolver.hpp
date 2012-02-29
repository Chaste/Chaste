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

#ifndef TESTSTOKESFLOWSOLVER_HPP_
#define TESTSTOKESFLOWSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "StokesFlowAssembler.hpp"
#include "StokesFlowSolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraticMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "Warnings.hpp"
#include "NumericFileComparison.hpp"

// add 3d test, add test to TestStokesWithLidCavity
//
// Why not parallel
//
// Better preconditioner (mass matrix in C-block),


class TestStokesFlowSolver : public CxxTest::TestSuite
{
public:
    /*
     * Solution is u = [x, -y], p = const (= 1 as applying zero-Neumann on RHS).
     * Dirichlet BC applied on three sides, zero-Neumann on the other (so pressure is fully defined)
     * Just two elements.
     */
    void TestStokesWithDirichletVerySimple() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        for (unsigned run=0; run<2; run++)
        {
            // Set up a mesh on [0 1]x[0 1]
            unsigned num_elem = (run==0 ? 1 : 10);
            QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

            // Dynamic viscosity
            double mu = 1.0;

            // Boundary flow
            std::vector<unsigned> dirichlet_nodes;
            std::vector<c_vector<double,2> > dirichlet_flow;
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double x = mesh.GetNode(i)->rGetLocation()[0];
                double y = mesh.GetNode(i)->rGetLocation()[1];

                // Only apply on top, bottom and left boundaries
                if (x == 0.0 || y == 0.0 || y == 1.0)
                {
                    dirichlet_nodes.push_back(i);
                    c_vector<double,2> flow = zero_vector<double>(2);

                    flow(0) = x;
                    flow(1) = -y;
                    dirichlet_flow.push_back(flow);
                }
            }

            StokesFlowProblemDefinition<2> problem_defn(mesh);
            problem_defn.SetViscosity(mu);
            problem_defn.SetPrescribedFlowNodes(dirichlet_nodes, dirichlet_flow);

            StokesFlowSolver<2> solver(mesh, problem_defn, "SimpleStokesFlow");

            if(run==1)
            {
                // see comment above TS_ASSERTs, below
                solver.SetKspAbsoluteTolerance(1e-12);
            }

            solver.Solve();

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double x = mesh.GetNode(i)->rGetLocation()[0];
                double y = mesh.GetNode(i)->rGetLocation()[1];

                // solution is in finite element space, so FEM solution will be exact,
                // apart from linear solver errors
                TS_ASSERT_DELTA(solver.rGetVelocities()[i](0),  x, 1e-8);
                TS_ASSERT_DELTA(solver.rGetVelocities()[i](1), -y, 1e-8);
            }

            for (unsigned i=0; i<mesh.GetNumVertices(); i++)
            {
                // solution is in finite element space, so FEM solution will be exact,
                // apart from linear solver errors
                TS_ASSERT_DELTA(solver.rGetPressures()[i], 1.0, 1e-8);
            }

            // check the matrix is symmetric even after Dirichlet BCs have been applied
            TS_ASSERT(PetscMatTools::CheckSymmetry(solver.mSystemLhsMatrix));
        }
    }

    /*
     * Solution is u = [y(1-y), 0], p = 2(1-x).
     * Dirichlet BC applied on three sides, zero-Neumann on the other (so pressure is fully defined).
     */
    void TestStokesWithImposedPipeCondition() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Note: we could have num_elem=1 and test still pass, as FE solution is the same as the true
        // solution (analytic soln is in the FE space, ignoring linear solve errors. In fact, with
        // num_elem=1, the linear solve doesn't require the tolerance line below to be accurate
        // enough for the test).

        // set up a mesh on [0 1]x[0 1]
        unsigned num_elem = 10;
        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

        // Dynamic viscosity
        double mu = 1.0;

        // Boundary flow
        std::vector<unsigned> dirichlet_nodes;
        std::vector<c_vector<double,2> > dirichlet_flow;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            if (x == 0.0 || y == 0.0 || y == 1.0)
            {
                dirichlet_nodes.push_back(i);
                c_vector<double,2> flow = zero_vector<double>(2);

                flow(0) = y*(1-y);
                flow(1) = 0.0;
                dirichlet_flow.push_back(flow);
            }
        }

        assert(dirichlet_flow.size() == 6*num_elem +1);

        c_vector<double,2> body_force = zero_vector<double>(2);

        StokesFlowProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetViscosity(mu);
        problem_defn.SetPrescribedFlowNodes(dirichlet_nodes, dirichlet_flow);

        StokesFlowSolver<2> solver(mesh, problem_defn, "PipeStokesFlow");

        // see comment above TS_ASSERTs, below
        solver.SetKspAbsoluteTolerance(1e-12);

        solver.Solve();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double y = mesh.GetNode(i)->rGetLocation()[1];
            double exact_flow_x = y*(1-y);
            double exact_flow_y = 0.0;

            // solution is in finite element space, so FEM solution will be exact,
            // apart from linear solver errors
            TS_ASSERT_DELTA(solver.rGetVelocities()[i](0), exact_flow_x, 1e-8);
            TS_ASSERT_DELTA(solver.rGetVelocities()[i](1), exact_flow_y, 1e-8);
        }

        for (unsigned i=0; i<mesh.GetNumVertices(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double exact_pressure = 2*(1-x);

            // solution is in finite element space, so FEM solution will be exact,
            // apart from linear solver errors
            TS_ASSERT_DELTA( solver.rGetPressures()[i], exact_pressure, 1e-8);
        }

        // test output files
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "PipeStokesFlow";
        NumericFileComparison comp1(results_dir + "/flow_solution.nodes", "continuum_mechanics/test/data/PipeStokesFlow/flow_solution.nodes");
        TS_ASSERT(comp1.CompareFiles(1e-7));

        NumericFileComparison comp2(results_dir + "/pressure.txt", "continuum_mechanics/test/data/PipeStokesFlow/pressure.txt");
        TS_ASSERT(comp1.CompareFiles(1e-7));

        solver.WriteCurrentPressureSolution(10);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/pressure_10.txt " + results_dir + "/pressure.txt").c_str()), 0);
    }

    /*
     * Solution is u = [y(1-y), 0], p = 2(2-x).
     * Dirichlet BC applied on top and bottom, Neumann on the ends p(x=0)=3, p(x=1)=1.
     */
    void TestPoiseuilleFlow() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // set up a mesh on [0 1]x[0 1]
        unsigned num_elem = 10;
        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

        // Dynamic viscosity
        double mu = 1.0;

        // Boundary flow
        std::vector<unsigned> dirichlet_nodes;
        std::vector<c_vector<double,2> > dirichlet_flow;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double y=mesh.GetNode(i)->rGetLocation()[1];

            // Fix top and bottom
            if (y == 0.0 || y == 1.0)
            {
                dirichlet_nodes.push_back(i);
                c_vector<double,2> flow = zero_vector<double>(2);
                dirichlet_flow.push_back(flow);
            }
        }
        assert(dirichlet_flow.size()== 4*num_elem+2);

        // Normal stress boundary

        std::vector<BoundaryElement<1,2>*> boundary_elems;
        std::vector<c_vector<double,2> > normal_stresses;

        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
             iter != mesh.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            if (fabs((*iter)->CalculateCentroid()[0]) < 1e-4)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);

                c_vector<double,2> normal_stress = zero_vector<double>(2);
                normal_stress[0] = 3;

                normal_stresses.push_back(normal_stress);
            }
            else if (fabs((*iter)->CalculateCentroid()[0] - 1.0) < 1e-4)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);

                c_vector<double,2> normal_stress = zero_vector<double>(2);

                // This is negative because the outward pointing normal at this edge is opposite to the other edge
                normal_stress[0] = -1;

                normal_stresses.push_back(normal_stress);
            }
        }
        assert(boundary_elems.size() == 2*num_elem);

        c_vector<double,2> body_force = zero_vector<double>(2);

        StokesFlowProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetViscosity(mu);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, normal_stresses);
        problem_defn.SetPrescribedFlowNodes(dirichlet_nodes, dirichlet_flow);

        StokesFlowSolver<2> solver(mesh, problem_defn, "PoiseuilleFlow");

        solver.SetKspAbsoluteTolerance(1e-8);

        solver.Solve();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double y = mesh.GetNode(i)->rGetLocation()[1];

            double exact_flow_x = y*(1-y);
            double exact_flow_y = 0.0;

            // solution in FE space
            TS_ASSERT_DELTA(solver.rGetVelocities()[i](0), exact_flow_x, 1e-7);
            TS_ASSERT_DELTA(solver.rGetVelocities()[i](1), exact_flow_y, 1e-7);
        }

        for (unsigned i=0; i<mesh.GetNumVertices(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double exact_pressure = 2*(1-x) + 1;

            // solution is in FE space
            TS_ASSERT_DELTA( solver.rGetPressures()[i], exact_pressure, 1e-7);
        }

    }



///\todo #1959 This might not be a great test problem; discuss with Dave Kay.

    /*
     * Solution is u = [20xy^3, 5x^4-5y^4], p = 60x^2y-20y^3+const.
     * Dirichlet BC applied on all 4 sides so pressure is not fully defined.
     */
    void TestConvergenceWithAnalyticSolution() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        unsigned num_runs = 4;
        std::vector<double> L_inf_error_flow(num_runs, -1.0);
        std::vector<double> L_inf_error_p(num_runs, -1.0);
        std::vector<unsigned> num_elem(num_runs);

        // set up a mesh on [-1 1]x[-1 1]
        for(unsigned run=0; run<num_runs; run++)
        {
            num_elem[run] = unsigned(round(pow(2,run+1)));
            QuadraticMesh<2> mesh(2.0/num_elem[run], 2.0, 2.0);
            mesh.Translate(-1.0, -1.0);

            // Dynamic viscosity
            double mu = 1.0;

            // Boundary flow
            std::vector<unsigned> dirichlet_nodes;
            std::vector<c_vector<double,2> > dirichlet_flow;
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double x = mesh.GetNode(i)->rGetLocation()[0];
                double y = mesh.GetNode(i)->rGetLocation()[1];
                if (x == -1.0 || x == 1.0 || y == -1.0 || y == 1.0)
                {
                    dirichlet_nodes.push_back(i);
                    c_vector<double,2> flow = zero_vector<double>(2);

                    flow(0) = 20.0*x*y*y*y;
                    flow(1) = 5.0*x*x*x*x - 5.0*y*y*y*y;
                    dirichlet_flow.push_back(flow);
                }
            }

            assert(dirichlet_flow.size() == 8*num_elem[run]);

            c_vector<double,2> body_force = zero_vector<double>(2);

            StokesFlowProblemDefinition<2> problem_defn(mesh);
            problem_defn.SetViscosity(mu);
            problem_defn.SetPrescribedFlowNodes(dirichlet_nodes, dirichlet_flow);

            StokesFlowSolver<2> solver(mesh, problem_defn, "AnalyticalStokesFlow");

            solver.SetKspAbsoluteTolerance(1e-8);

            solver.Solve();


            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double x = mesh.GetNode(i)->rGetLocation()[0];
                double y = mesh.GetNode(i)->rGetLocation()[1];

                double exact_flow_x = 20.0*x*y*y*y;
                double exact_flow_y = 5.0*x*x*x*x - 5.0*y*y*y*y;

                //TS_ASSERT_DELTA(solver.rGetVelocities()[i](0), exact_flow_x, 5e-3);
                //TS_ASSERT_DELTA(solver.rGetVelocities()[i](1), exact_flow_y, 5e-3);

                double diff_x = fabs(solver.rGetVelocities()[i](0) - exact_flow_x);
                double diff_y = fabs(solver.rGetVelocities()[i](1) - exact_flow_y);
                double max_diff = std::max(diff_x,diff_y);
                L_inf_error_flow[run] = std::max(L_inf_error_flow[run], max_diff);
            }

            //Calculate the constant offset between the true solution and the numerical solution.
            double pressure_diff = 0.0;
            for (unsigned index = 0; index < mesh.GetNumVertices(); ++index)
            {
                double x = mesh.GetNode(index)->rGetLocation()[0];
                double y = mesh.GetNode(index)->rGetLocation()[1];
                double exact_pressure = 60.0*x*x*y -20.0*y*y*y;
                pressure_diff += solver.rGetPressures()[index] - exact_pressure;
            }
            pressure_diff /= mesh.GetNumVertices();

            for (unsigned i=0; i<mesh.GetNumVertices(); i++)
            {
                double x = mesh.GetNode(i)->rGetLocation()[0];
                double y = mesh.GetNode(i)->rGetLocation()[1];

                double exact_pressure = 60.0*x*x*y -20.0*y*y*y;

                //TS_ASSERT_DELTA( solver.rGetPressures()[i], exact_pressure + pressure_diff, 1e-0);

                L_inf_error_p[run] = std::max(L_inf_error_p[run], fabs(solver.rGetPressures()[i] - exact_pressure - pressure_diff));
            }
        }

        std::cout << "Num_elements Linf_error_flow Linf_error_p\n";
        for(unsigned i=0; i<num_runs; i++)
        {
            std::cout << num_elem[i] << " " << L_inf_error_flow[i] << " " << L_inf_error_p[i] << "\n";
        }

        /* results are:
         *
         * 2 1.15234 36.5909
         * 4 0.352193 9.51735
         * 8 0.0548267 1.86648
         * 16 0.00787528 0.504425
         * 32 0.00113687 0.117679
         * 64 0.000153278 0.0287536
         *
         * Everything is converging. Large errors in p down to this being an odd problem?
         */
        double res_flow[6] = { 1.15234, 0.352193, 0.0548267, 0.00787528, 0.00113687, 0.000153278 };
        double res_p[6] = { 36.5909, 9.51735, 1.86648, 0.504425, 0.117679, 0.0287536 };
        assert(num_runs <= 6);
        for(unsigned i=0; i<num_runs; i++)
        {
            TS_ASSERT_DELTA( L_inf_error_flow[i], res_flow[i], 1e-3);
            TS_ASSERT_DELTA( L_inf_error_p[i], res_p[i], 1e-3);
        }
    }

    /*
     * Simulation with regularised lid driven cavity u=1-x^4 on the top.
     */
    void TestStokesWithLidCavity() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Set up a mesh on [-1 1]x[-1 1]
        unsigned num_elem = 5;
        QuadraticMesh<2> mesh(2.0/num_elem, 2.0, 2.0);
        mesh.Translate(-1.0, -1.0);

        // Dynamic viscosity
        double mu = 1.0;

        // Boundary flow
        std::vector<unsigned> dirichlet_nodes;
        std::vector<c_vector<double,2> > dirichlet_flow;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            if (x == -1.0 || x == 1.0 || y == -1.0)
            {
                dirichlet_nodes.push_back(i);
                c_vector<double,2> flow = zero_vector<double>(2);
                dirichlet_flow.push_back(flow);
            }
            else if (y == 1.0) // this doesnt include corners
            {
                dirichlet_nodes.push_back(i);
                c_vector<double,2> flow = zero_vector<double>(2);

                flow(0) = 1-x*x*x*x;
                flow(1) = 0.0;
                dirichlet_flow.push_back(flow);
            }

        }

        assert(dirichlet_flow.size()== 8*num_elem);

        c_vector<double,2> body_force = zero_vector<double>(2);

        StokesFlowProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetViscosity(mu);
        problem_defn.SetPrescribedFlowNodes(dirichlet_nodes, dirichlet_flow);

        StokesFlowSolver<2> solver(mesh, problem_defn, "LidCavityStokesFlow");

        // Uncomment to make errors smaller
        //solver.SetKspAbsoluteTolerance(1e-12);

        solver.Solve();

        ///\todo Test something
    }
};

#endif // TESTSTOKESFLOWSOLVER_HPP_
