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
#include "FileComparison.hpp"

/*
 *  HOW_TO_TAG Continuum mechanics
 *  Solve Stokes' flow problems (this functionality is work-in-progress).
 *
 */
class TestStokesFlowSolver : public CxxTest::TestSuite
{
public:
    /*
     * Exact solution here is u = [x, -y], p = 20
     * Dirichlet BC applied on three sides, zero-stress on the other (so pressure is fully defined)
     * Just two elements.
     */
    void TestStokesExactSolutionSimple()
    {
        for (unsigned run=0; run<2; run++)
        {
            // Set up a mesh on [0 1]x[0 1]
            unsigned num_elem = (run==0 ? 1 : 10);
            QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

            // Dynamic viscosity
            double mu = 10.0;

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

            // coverage
            problem_defn.SetVerboseDuringSolve();

            StokesFlowSolver<2> solver(mesh, problem_defn, "SimpleStokesFlow");

            if (run==1)
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

            // test the pressures
            std::vector<double>& r_pressures = solver.rGetPressures();
            for (unsigned i=0; i<r_pressures.size(); i++)
            {
                // solution is in finite element space, so FEM solution will be exact,
                // apart from linear solver errors
                TS_ASSERT_DELTA(r_pressures[i], 20.0, 1e-6);
            }

            // check the matrix is symmetric even after Dirichlet BCs have been applied
            TS_ASSERT(PetscMatTools::CheckSymmetry(solver.mSystemLhsMatrix));
        }
    }

    // Exact solution for this problem is u = [y, -x], p = 0
    // For this flow: sigma = mu(grad u + (grad u)^T) - pI = mu*0 - pI = -pI
    // Again, exact solution is in FEM space so would work, up to
    // linear solve tolerance, with 1 element
    void TestStokesExactSolutionLessSimple()
    {
        // Set up a mesh on [0 1]x[0 1]
        unsigned num_elem = 3;
        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

        // Dynamic viscosity
        double mu = 10.0;

        // Boundary flow
        std::vector<unsigned> dirichlet_nodes;
        std::vector<c_vector<double,2> > dirichlet_flow;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];

            // Only apply on top and left boundaries (zero stress BCs on others)
            if (x == 0.0 || y == 0.0)
            {
                dirichlet_nodes.push_back(i);
                c_vector<double,2> flow = zero_vector<double>(2);

                flow(0) = y;
                flow(1) = -x;
                dirichlet_flow.push_back(flow);
            }
        }

        StokesFlowProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetViscosity(mu);
        problem_defn.SetPrescribedFlowNodes(dirichlet_nodes, dirichlet_flow);

        StokesFlowSolver<2> solver(mesh, problem_defn, "LessSimpleStokesFlow");

        solver.SetKspAbsoluteTolerance(1e-12);

        solver.Solve();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];

            // solution is in finite element space, so FEM solution will be exact,
            // apart from linear solver errors
            TS_ASSERT_DELTA(solver.rGetVelocities()[i](0),  y, 1e-8);
            TS_ASSERT_DELTA(solver.rGetVelocities()[i](1), -x, 1e-8);
        }

        // test the pressures
        std::vector<double>& r_pressures = solver.rGetPressures();
        for (unsigned i=0; i<r_pressures.size(); i++)
        {
            // solution is in finite element space, so FEM solution will be exact,
            // apart from linear solver errors
            TS_ASSERT_DELTA(r_pressures[i], 0.0, 1e-6);
        }
    }


    /*
     * Solution is u = [y(1-y), 0], p = 2(1-x) + const
     * Dirichlet BC applied on all four sides
     */
    void TestStokesWithImposedPipeCondition()
    {
        // Note: we could have num_elem=1 and test still pass, as FE solution is the same as the true
        // solution (analytic soln is in the FE space, ignoring linear solve errors.

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
            if (x == 0.0 || x==1.0 || y == 0.0 || y == 1.0) // should really be using a BoundaryNodeIterator
            {
                dirichlet_nodes.push_back(i);
                c_vector<double,2> flow = zero_vector<double>(2);

                flow(0) = y*(1-y);
                flow(1) = 0.0;
                dirichlet_flow.push_back(flow);
            }
        }

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


        std::vector<double>& r_pressures = solver.rGetPressures();

        // determine what the constant is
        double x = mesh.GetNode(0)->rGetLocation()[0];
        double constant = r_pressures[0] - 2*(1-x);
        // test the rest
        for (unsigned i=1; i<r_pressures.size(); i++)
        {
            x = mesh.GetNode(i)->rGetLocation()[0];
            double exact_pressure = 2*(1-x) + constant;
            // solution is in finite element space, so FEM solution will be exact,
            // apart from linear solver errors
            TS_ASSERT_DELTA(r_pressures[i], exact_pressure, 1e-8);
        }

        // test output files
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "PipeStokesFlow";
        NumericFileComparison comp1(results_dir + "/flow_solution.nodes", "continuum_mechanics/test/data/PipeStokesFlow/flow_solution.nodes");
        TS_ASSERT(comp1.CompareFiles(1e-7));

        NumericFileComparison comp2(results_dir + "/pressure.txt", "continuum_mechanics/test/data/PipeStokesFlow/pressure.txt");
        TS_ASSERT(comp2.CompareFiles(1e-7));

        solver.WriteCurrentPressureSolution(10);

        FileComparison comparer(results_dir + "/pressure_10.txt", results_dir + "/pressure.txt");
        TS_ASSERT(comparer.CompareFiles());

        NumericFileComparison comp3(results_dir + "/pressure_10.txt", results_dir + "/pressure.txt");
        TS_ASSERT(comp3.CompareFiles(1e-17));
    }


    // Exact solution for this problem is u = [y, -x], p = -3
    // For this flow: sigma = mu(grad u + (grad u)^T) - pI = mu*0 - pI = -pI
    // Again, exact solution is in FEM space so would work, up to
    // linear solve tolerance, with 1 element
    void TestStokesExactSolutionNonzeroNeumann()
    {
        // Set up a mesh on [0 1]x[0 1]
        unsigned num_elem = 3;
        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

        // Dynamic viscosity
        double mu = 10.0;

        // Boundary flow
        std::vector<unsigned> dirichlet_nodes;
        std::vector<c_vector<double,2> > dirichlet_flow;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];

            // Only apply on top and left boundaries
            if (x == 0.0 || y == 0.0)
            {
                dirichlet_nodes.push_back(i);
                c_vector<double,2> flow = zero_vector<double>(2);

                flow(0) = y;
                flow(1) = -x;
                dirichlet_flow.push_back(flow);
            }
        }

        // apply non-zero Neumann BCs on right and top sides
        std::vector<BoundaryElement<1,2>*> boundary_elems;
        std::vector<c_vector<double,2> > stresses;

        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
             iter != mesh.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            if (fabs((*iter)->CalculateCentroid()[0] - 1.0) < 1e-4)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);

                c_vector<double,2> stress = zero_vector<double>(2);
                stress(0) = 3.0; // stress = (3,0) = 3*normal
                stresses.push_back(stress);
            }
            else if (fabs((*iter)->CalculateCentroid()[1] - 1.0) < 1e-4)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);

                c_vector<double,2> stress = zero_vector<double>(2);
                stress(1) = 3.0; // stress = (0,3) = 3*normal
                stresses.push_back(stress);
            }
        }
        assert(boundary_elems.size() == 2*num_elem);

        StokesFlowProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetViscosity(mu);
        problem_defn.SetPrescribedFlowNodes(dirichlet_nodes, dirichlet_flow);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, stresses);

        StokesFlowSolver<2> solver(mesh, problem_defn, "StokesFlowNonZeroNeumann");

        solver.SetKspAbsoluteTolerance(1e-12);

        solver.Solve();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];

            // solution is in finite element space, so FEM solution will be exact,
            // apart from linear solver errors
            TS_ASSERT_DELTA(solver.rGetVelocities()[i](0),  y, 1e-8);
            TS_ASSERT_DELTA(solver.rGetVelocities()[i](1), -x, 1e-8);
        }

        // test the pressures
        std::vector<double>& r_pressures = solver.rGetPressures();
        for (unsigned i=0; i<r_pressures.size(); i++)
        {
            // solution is in finite element space, so FEM solution will be exact,
            // apart from linear solver errors
            TS_ASSERT_DELTA(r_pressures[i], -3.0, 1e-6);
        }
    }

    /*
     * Solution is u = [20xy^3, 5x^4-5y^4], p = 60x^2y-20y^3+const.
     * Dirichlet BC applied on all 4 sides so pressure is not fully defined.
     */
    void TestConvergenceWithAnalyticSolution()
    {
        unsigned num_runs = 4;
        std::vector<double> L_inf_error_flow(num_runs, -1.0);
        std::vector<double> L_inf_error_p(num_runs, -1.0);
        std::vector<unsigned> num_elem(num_runs);

        // set up a mesh on [-1 1]x[-1 1]
        for (unsigned run=0; run<num_runs; run++)
        {
            num_elem[run] = SmallPow(2u, run+1u);
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
                double max_difference = std::max(diff_x,diff_y);
                L_inf_error_flow[run] = std::max(L_inf_error_flow[run], max_difference);
            }

            //Calculate the constant offset between the true solution and the numerical solution.
            std::vector<double>& r_pressures = solver.rGetPressures();

            double pressure_difference = 0.0;
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double x = mesh.GetNode(i)->rGetLocation()[0];
                double y = mesh.GetNode(i)->rGetLocation()[1];
                double exact_pressure = 60.0*x*x*y -20.0*y*y*y;
                pressure_difference += r_pressures[i] - exact_pressure;
            }
            pressure_difference /= mesh.GetNumVertices();

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double x = mesh.GetNode(i)->rGetLocation()[0];
                double y = mesh.GetNode(i)->rGetLocation()[1];
                double exact_pressure = 60.0*x*x*y -20.0*y*y*y;
                L_inf_error_p[run] = std::max(L_inf_error_p[run], fabs(r_pressures[i] - exact_pressure - pressure_difference));
            }
        }

        std::cout << "Num_elements Linf_error_flow Linf_error_p\n";
        for (unsigned i=0; i<num_runs; i++)
        {
            std::cout << num_elem[i] << " " << L_inf_error_flow[i] << " " << L_inf_error_p[i] << "\n";
        }

        /* results are:
         *
         *
         * 2 1.51032 14.3527
         * 4 0.250461 5.12699
         * 8 0.0368664 1.95504
         * 16 0.00464591 0.576949
         * 32 0.000625753 0.155241
         * Everything is converging.
         */
        double res_flow[5] = { 1.51032, 0.250461, 0.0368664, 0.00464591, 0.000625753 };
        double res_p[5] = { 14.3527, 5.12699, 1.95504, 0.576949, 0.155241};
        assert(num_runs <= 5);
        for (unsigned i=0; i<num_runs; i++)
        {
            TS_ASSERT_DELTA( L_inf_error_flow[i], res_flow[i], 1e-3);
            TS_ASSERT_DELTA( L_inf_error_p[i], res_p[i], 1e-3);
        }
    }

    // Simulation with regularised lid driven cavity, u=x(1-x) on top
    // Using this top to compare with independently written code.
    // Alternatively, could use u=1-x^4 on the top (with geometry [-1,1]^2)
    // and compare with plots in Andy Wathen's fluids FEM book
    void TestStokesWithLidDrivenCavity()
    {
        unsigned num_elem = 5;
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

            if (y == 1.0)
            {
                dirichlet_nodes.push_back(i);
                c_vector<double,2> flow = zero_vector<double>(2);

                flow(0) = x*(1-x);
                flow(1) = 0.0;
                dirichlet_flow.push_back(flow);
            }
            else if (x == 0.0 || x == 1.0 || y == 0.0)
            {
                dirichlet_nodes.push_back(i);
                c_vector<double,2> flow = zero_vector<double>(2);
                flow(0) = 0.0;
                flow(1) = 0.0;
                dirichlet_flow.push_back(flow);
            }
        }

        StokesFlowProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetViscosity(mu);
        problem_defn.SetPrescribedFlowNodes(dirichlet_nodes, dirichlet_flow);

        StokesFlowSolver<2> solver(mesh, problem_defn, "LidDrivenCavityStokesFlow");

        // Uncomment to make errors smaller
        solver.SetKspAbsoluteTolerance(1e-10);

        solver.Solve();

        solver.CreateVtkOutput("Velocity");
#ifdef CHASTE_VTK
        //Check the VTK file exists
        FileFinder vtk_file("LidDrivenCavityStokesFlow/vtk/solution.vtu", RelativeTo::ChasteTestOutput);
        TS_ASSERT(vtk_file.Exists());
#endif


        double min_u = DBL_MAX;
        double max_u = -DBL_MAX;
        double min_v = DBL_MAX;
        double max_v = -DBL_MAX;
        double min_p = DBL_MAX;
        double max_p = -DBL_MAX;

        std::vector<c_vector<double,2> >& r_flow = solver.rGetVelocities();
        std::vector<double>& r_pressures = solver.rGetPressures();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            min_u = std::min(min_u, r_flow[i](0));
            max_u = std::max(max_u, r_flow[i](0));
            min_v = std::min(min_v, r_flow[i](1));
            max_v = std::max(max_v, r_flow[i](1));
            min_p = std::min(min_p, r_pressures[i]);
            max_p = std::max(max_p, r_pressures[i]);
        }

        // at the moment, one way of visualising this is to hack AbstractContinuumMechanicsSolver::WriteCurrentSpatialSolution()
        // to write node locations as well as the flow, then do for example, in matlab, s=load('LidDrivenCavityStokesFlow/flow_solution.nodes');
        // and quiver(s(:,1),s(:,2),s(:,3),s(:,4),0.2);

        // have visualised and compared to Raf's libmesh Stokes results, looks very similar

        // max u should be exactly 0.25 (dirichlet bc value)
        TS_ASSERT_DELTA(max_u, 0.25, 1e-8);

        // compare with Raf's libmesh Stokes, on same resolution mesh for which
        // min u = -0.04590, min v = -0.07382, max v = 0.07804
        TS_ASSERT_DELTA(min_u, -0.0459, 1e-3);
        TS_ASSERT_DELTA(min_v, -0.075, 5e-3);
        TS_ASSERT_DELTA(max_v,  0.075, 5e-3);

        // find a node in the interior for which u and v both not small
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            if (fabs(x-0.3)<1e-6 && fabs(y-0.6)<1e-6)
            {
                //at this node raf's solver returns (-0.033649 0.04652)
                TS_ASSERT_DELTA(r_flow[i](0), -0.035, 1e-3);
                TS_ASSERT_DELTA(r_flow[i](1),  0.047, 1e-3);
            }
        }

        TS_ASSERT_DELTA(min_p, -5.2119, 1e-2);
        TS_ASSERT_DELTA(max_p, 13.1247, 1e-2);
    }


    // solve problem for which solution is u=(x,y,-2z), p=const.
    void TestStokesSimple3d()
    {
        unsigned num_elem = 2;
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

            flow(0) = x;
            flow(1) = y;
            flow(2) = -2*z;

            dirichlet_nodes.push_back((*iter)->GetIndex());
            dirichlet_flow.push_back(flow);
        }

        StokesFlowProblemDefinition<3> problem_defn(mesh);
        problem_defn.SetViscosity(mu);
        problem_defn.SetPrescribedFlowNodes(dirichlet_nodes, dirichlet_flow);

        StokesFlowSolver<3> solver(mesh, problem_defn, "3dSimple");

        solver.SetKspAbsoluteTolerance(1e-12);
        solver.Solve();

        std::vector<c_vector<double,3> >& r_solution = solver.rGetVelocities();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            double z = mesh.GetNode(i)->rGetLocation()[2];

            TS_ASSERT_DELTA(r_solution[i](0), x,    1e-6);
            TS_ASSERT_DELTA(r_solution[i](1), y,    1e-6);
            TS_ASSERT_DELTA(r_solution[i](2), -2*z, 1e-6);
        }

        // test the pressures
        std::vector<double>& r_pressures = solver.rGetPressures();
        bool first = true;
        double value;
        for (unsigned i=0; i<r_pressures.size(); i++)
        {
            if (first)
            {
                value = r_pressures[i];
                first = false;
            }
            else
            {
                TS_ASSERT_DELTA(r_pressures[i], value, 5e-5);
            }
        }
    }

    void TestStokesWithLidDrivenCavity3d()
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
        solver.CreateVtkOutput("Velocity");
#ifdef CHASTE_VTK
        //Check the VTK file exists
        FileFinder vtk_file("LidDrivenCavityStokesFlow3d/vtk/solution.vtu", RelativeTo::ChasteTestOutput);
        TS_ASSERT(vtk_file.Exists());
#endif

        std::vector<c_vector<double,3> >& r_solution = solver.rGetVelocities();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            double z = mesh.GetNode(i)->rGetLocation()[2];
            if ((fabs(x-0.6)<1e-6) && (fabs(y-0.6)<1e-6) && (fabs(z-0.6)<1e-6))
            {
                TS_ASSERT_DELTA(r_solution[i](0), -8.5990e-03, 1e-4);
                TS_ASSERT_DELTA(r_solution[i](1),  1.1331e-04, 1e-5);
                TS_ASSERT_DELTA(r_solution[i](2), -5.1343e-03, 1e-4);
            }
        }
    }
};

#endif // TESTSTOKESFLOWSOLVER_HPP_
