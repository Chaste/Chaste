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

            double x0 = mesh.GetNode(0)->rGetLocation()[0];
            double y0 = mesh.GetNode(0)->rGetLocation()[1];
            double exact_pressure_0 = 60.0*x0*x0*y0 -20.0*y0*y0*y0;
            double pressure_diff = solver.rGetPressures()[0] - exact_pressure_0;

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
         * 2 1.15234 73.1818
         * 4 0.352193 16.5878
         * 8 0.0548267 3.61898
         * 16 0.00787528 0.957452
         * 32 0.00113687 0.225023
         * 64 0.000153278, 0.0548126
         *
         * Everything is converging. Large errors in p down to this being an odd problem?
         */
        double res_flow[6] = { 1.15234, 0.352193, 0.0548267, 0.00787528, 0.00113687, 0.000153278 };
        double res_p[6] = { 73.1818, 16.5878, 3.61898, 0.957452, 0.225023, 0.0548126 };
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
