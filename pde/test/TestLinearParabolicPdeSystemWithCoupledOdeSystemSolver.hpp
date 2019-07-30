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

#ifndef TESTLINEARPARABOLICPDESYSTEMWITHCOUPLEDODESYSTEMSOLVER_HPP_
#define TESTLINEARPARABOLICPDESYSTEMWITHCOUPLEDODESYSTEMSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include <vector>
#include <cmath>
#include "BoundaryConditionsContainerImplementation.hpp"
#include "AbstractBoundaryConditionsContainerImplementation.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OdeSystemForCoupledHeatEquation.hpp"
#include "OdeSystemForCoupledHeatEquationWithSource.hpp"
#include "HeatEquationForCoupledOdeSystem.hpp"
#include "HeatEquationWithSourceForCoupledOdeSystem.hpp"
#include "SchnackenbergCoupledPdeSystem.hpp"
#include "RandomNumberGenerator.hpp"
#include "LinearParabolicPdeSystemWithCoupledOdeSystemSolver.hpp"
#include "OutputFileHandler.hpp"
#include "NumericFileComparison.hpp"
#include "FileComparison.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscTools.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestLinearParabolicPdeSystemWithCoupledOdeSystemSolver : public CxxTest::TestSuite
{
public:

    void TestHeatEquationWithSourceWithCoupledOdeSystemIn1dWithZeroNeumann()
    {
        // Create mesh of domain [0,1]
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create PDE system object
        HeatEquationWithSourceForCoupledOdeSystem<1> pde;

        // Define zero Neumann boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
        iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

        // Create the correct number of ODE systems
        double a = 5.0;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> ode_systems;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            ode_systems.push_back(new OdeSystemForCoupledHeatEquationWithSource(a));
        }

        // Create PDE system solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<1,1,1> solver(&mesh, &pde, &bcc, ode_systems);

        // Test setting end time and timestep
        TS_ASSERT_THROWS_THIS(solver.SetTimes(1.0, 0.0), "Start time has to be less than end time");
        TS_ASSERT_THROWS_THIS(solver.SetTimeStep(0.0), "Time step has to be greater than zero");

        // Set end time and timestep
        double t_end = 0.1;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.001);

        // Set initial condition u(x,0) = 1 + cos(pi*x)
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            init_cond[i] = 1 + cos(M_PI*x);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE system and store result
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        /*
         * Test that solution is given by
         *
         * u(x,t) = 1 + (1 - exp(-a*t))/a + exp(-pi*pi*t)*cos(pi*x),
         * v(x,t) = exp(-a*t),
         *
         * with t = t_end.
         */
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];

            double u = 1 + (1 - exp(-a*t_end))/a + exp(-M_PI*M_PI*t_end)*cos(M_PI*x);
            TS_ASSERT_DELTA(result_repl[i], u, 0.1);

            double u_from_v = solver.GetOdeSystemAtNode(i)->rGetPdeSolution()[0];
            TS_ASSERT_DELTA(result_repl[i], u_from_v, 0.1);

            double v = exp(-a*t_end);
            TS_ASSERT_DELTA(ode_systems[i]->rGetStateVariables()[0], v, 0.1);
        }

        // Test the method GetOdeSystemAtNode()
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT(solver.GetOdeSystemAtNode(i) != NULL);
            TS_ASSERT_DELTA(static_cast<OdeSystemForCoupledHeatEquationWithSource*>(solver.GetOdeSystemAtNode(i))->GetA(), 5.0, 1e-6);
        }

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    void TestHeatEquationWithCoupledOdeSystemIn1dWithMixed()
    {
        // Create mesh of domain [0,1]
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create PDE system object
        HeatEquationForCoupledOdeSystem<1> pde;

        // Define non-zero Neumann boundary condition at x=0
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

        // Define zero Dirichlet boundary condition at x=1
        ConstBoundaryCondition<1>* p_boundary_condition2 = new ConstBoundaryCondition<1>(0.0);
        TetrahedralMesh<1,1>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorEnd();
        --node_iter;
        bcc.AddDirichletBoundaryCondition(*node_iter, p_boundary_condition2);

        // Create the correct number of ODE systems
        double a = 5.0;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> ode_systems;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            ode_systems.push_back(new OdeSystemForCoupledHeatEquation(a));
        }

        // Create PDE system solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<1,1,1> solver(&mesh, &pde, &bcc, ode_systems);

        // Set end time and timestep
        double t_end = 0.1;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.001);

        // Set initial condition u(x,0) = 1 - x
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            init_cond[i] = 1 - x;
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE system and store result
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        /*
         * Test that solution is given by
         *
         * u(x,t) = 1 - x,
         * v(x,t) = 1 + a*(1-x)*t,
         *
         * with t = t_end.
         */
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = 1 - x;
            TS_ASSERT_DELTA(result_repl[i], u, 0.1);

            double v = 1 + a*(1-x)*t_end;
            TS_ASSERT_DELTA(ode_systems[i]->rGetStateVariables()[0], v, 0.1);
        }

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    void TestHeatEquationWithCoupledOdeSystemIn2dWithZeroDirichlet()
    {
        // Create mesh of the domain [0,1]x[0,1]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create PDE system object
        HeatEquationForCoupledOdeSystem<2> pde;

        // Define zero Dirichlet boundary conditions on entire boundary
        BoundaryConditionsContainer<2,2,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

        // Create the correct number of ODE systems
        double a = 5.0;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> ode_systems;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            ode_systems.push_back(new OdeSystemForCoupledHeatEquation(a));
        }

        // Create PDE system solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<2,2,1> solver(&mesh, &pde, &bcc, ode_systems);

        // Set end time and timestep
        double t_end = 0.01;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.001);

        /*
         * Set initial condition
         *
         * u(x,y,0) = sin(pi*x)*sin(pi*y),
         *
         * which is an eigenfunction of the heat equation.
         */
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            init_cond[i] = sin(M_PI*x)*sin(M_PI*y);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE system and store result
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        /*
         * Test that solution is given by
         *
         * u(x,y,t) = e^{-2*pi*pi*t} sin(pi*x)*sin(pi*y),
         * v(x,y,t) = 1 + (1 - e^{-2*pi*pi*t})*sin(pi*x)*sin(pi*y)*a/(2*pi*pi),
         *
         * with t = t_end.
         */
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];

            double u = exp(-2*M_PI*M_PI*t_end)*sin(M_PI*x)*sin(M_PI*y);
            double v = 1.0 + (a/(2*M_PI*M_PI))*(1 - exp(-2*M_PI*M_PI*t_end))*sin(M_PI*x)*sin(M_PI*y);

            TS_ASSERT_DELTA(result_repl[i], u, 0.01);
            TS_ASSERT_DELTA(ode_systems[i]->rGetStateVariables()[0], v, 0.01);
        }

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    void TestSolveAndWriteResultsToFileMethod()
    {
        // Create mesh of the domain [0,1]x[0,1]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create PDE system object
        HeatEquationForCoupledOdeSystem<2> pde;

        // Define zero Dirichlet boundary conditions on entire boundary
        BoundaryConditionsContainer<2,2,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

        // Create the correct number of ODE systems
        double a = 5.0;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> ode_systems;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            ode_systems.push_back(new OdeSystemForCoupledHeatEquation(a));
        }

        // Create PDE system solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<2,2,1> solver(&mesh, &pde, &bcc, ode_systems);

        // Set end time and timestep (end time is not a multiple of timestep, for coverage)
        double t_end = 0.105;

        /*
         * Set initial condition
         *
         * u(x,y,0) = sin(pi*x)*sin(pi*y),
         *
         * which is an eigenfunction of the heat equation.
         */
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            init_cond[i] = sin(M_PI*x)*sin(M_PI*y);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);

        // Need an output folder
        TS_ASSERT_THROWS_THIS(solver.SolveAndWriteResultsToFile(),
                "SetOutputDirectory() must be called prior to SolveAndWriteResultsToFile()");
        solver.SetOutputDirectory("TestHeatEquationForCoupledOdeSystemIn2dWithZeroDirichletWithOutput");

        // Need a time interval
        TS_ASSERT_THROWS_THIS(solver.SolveAndWriteResultsToFile(),
                "SetTimes() must be called prior to SolveAndWriteResultsToFile()");
        solver.SetTimes(0, t_end);

        // Need a timestep
        TS_ASSERT_THROWS_THIS(solver.SolveAndWriteResultsToFile(),
                "SetTimeStep() must be called prior to SolveAndWriteResultsToFile()");
        solver.SetTimeStep(0.01);

        // Need sampling interval
        TS_ASSERT_THROWS_THIS(solver.SolveAndWriteResultsToFile(),
                "SetSamplingTimeStep() must be called prior to SolveAndWriteResultsToFile()");
        solver.SetSamplingTimeStep(0.1);

        // Need initial condition
        TS_ASSERT_THROWS_THIS(solver.SolveAndWriteResultsToFile(),
                "SetInitialCondition() must be called prior to SolveAndWriteResultsToFile()");
        solver.SetInitialCondition(initial_condition);

        solver.SolveAndWriteResultsToFile();
#ifdef CHASTE_VTK
        // Read back a sample of the output from VTK
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestHeatEquationForCoupledOdeSystemIn2dWithZeroDirichletWithOutput/";
        //Meta file
        FileComparison(results_dir + "/results.pvd",
                       "pde/test/data/results.pvd").CompareFiles();
        //Final time-step
        VtkMeshReader<2,2> vtk_reader(results_dir+"results_2.vtu");
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 128u);
        std::vector<double> pde_result, ode_result;
        vtk_reader.GetPointData("PDE variable 0", pde_result);
        vtk_reader.GetPointData("ODE variable 0", ode_result);
        TS_ASSERT_EQUALS(pde_result.size(), 81u );
        TS_ASSERT_EQUALS(ode_result.size(), 81u);
        // Centre
        TS_ASSERT_DELTA(pde_result[4], 0.1427, 1e-4);
        TS_ASSERT_DELTA(ode_result[4], 1.2545, 1e-4);
        // Top corner
        TS_ASSERT_DELTA(pde_result[2], 0.0, 1e-4);
        TS_ASSERT_DELTA(ode_result[2], 1.0, 1e-4);
#endif // CHASTE_VTK

        // Tidy up
        PetscTools::Destroy(initial_condition);
    }

    /**
     * This test provides an example of how to solve a coupled PDE system
     * where there is no coupled ODE system, and can be used as a template
     * for solving standard reaction-diffusion problems arising in the
     * study of pattern formation on fixed domains.
     */
    void TestSchnackenbergCoupledPdeSystemIn1dWithNonZeroDirichlet()
    {
        // Create mesh of domain [0,1]
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_1000_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create PDE system object
        SchnackenbergCoupledPdeSystem<1> pde(1e-4, 1e-2, 0.1, 0.2, 0.3, 0.1);

        // Create non-zero Dirichlet boundary conditions for each state variable
        BoundaryConditionsContainer<1,1,2> bcc;
        ConstBoundaryCondition<1>* p_bc_for_u = new ConstBoundaryCondition<1>(2.0);
        ConstBoundaryCondition<1>* p_bc_for_v = new ConstBoundaryCondition<1>(0.75);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_bc_for_u, 0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_bc_for_v, 1);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(mesh.GetNumNodes()-1), p_bc_for_u, 0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(mesh.GetNumNodes()-1), p_bc_for_v, 1);

        // Create PDE system solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<1,1,2> solver(&mesh, &pde, &bcc);

        // Set end time and time step
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-1);

        // Create initial conditions that are random perturbations of the uniform steady state
        std::vector<double> init_conds(2*mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            init_conds[2*i] = fabs(2.0 + RandomNumberGenerator::Instance()->ranf());
            init_conds[2*i + 1] = fabs(0.75 + RandomNumberGenerator::Instance()->ranf());
        }
        Vec initial_condition = PetscTools::CreateVec(init_conds);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE system and store result
        Vec solution = solver.Solve();
        ReplicatableVector solution_repl(solution);

        // Write results for visualization in gnuplot
        OutputFileHandler handler("TestSchnackenbergCoupledPdeSystemIn1dWithNonZeroDirichlet", false);
        out_stream results_file = handler.OpenOutputFile("schnackenberg.dat");
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double u = solution_repl[2*i];
            double v = solution_repl[2*i + 1];
            (*results_file) << x << "\t" << u << "\t" << v << "\n" << std::flush;
        }
        results_file->close();

        std::string results_filename = handler.GetOutputDirectoryFullPath() + "schnackenberg.dat";
        NumericFileComparison comp_results(results_filename, "pde/test/data/schnackenberg.dat");
        TS_ASSERT(comp_results.CompareFiles(1e-3));

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(solution);
    }
};

#endif /*TESTLINEARPARABOLICPDESYSTEMWITHCOUPLEDODESYSTEMSOLVER_HPP_*/
