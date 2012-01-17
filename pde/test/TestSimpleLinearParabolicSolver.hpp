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

#ifndef _TESTSIMPLELINEARPARABOLICSOLVER_HPP_
#define _TESTSIMPLELINEARPARABOLICSOLVER_HPP_

/**
 * TestSimpleLinearParabolicSolver.hpp
 *
 * Test suite for the SimpleLinearParabolicSolver class.
 *
 * Tests the class for the solution of parabolic pdes in 1D, 2D and 3D with and
 * without source terms with neumann and dirichlet booundary conditions.
 */

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include <petsc.h>
#include <vector>
#include <cmath>
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "ParallelColumnDataWriter.hpp"
#include "TrianglesMeshReader.hpp"
#include "FemlabMeshReader.hpp"
#include "HeatEquation.hpp"
#include "HeatEquationWithSourceTerm.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include "CompareHdf5ResultsFiles.hpp"

/*
 * Very simple toy time-adaptivity controller, for use in
 * Test1DProblemUsingTimeAdaptivityController. Returns
 * dt=0.001 for the first half of simulation, dt=0.01 for
 * the second half.
 */
class ToyController : public AbstractTimeAdaptivityController
{
    double ComputeTimeStep(double currentTime, Vec currentSolution)
    {
        if (currentTime < 0.05)
        {
            return 0.001;
        }
        else
        {
            return 0.01;
        }
    }

public:
    ToyController()
        : AbstractTimeAdaptivityController(0.001, 0.01)
    {
    }
};

class TestSimpleLinearParabolicSolver : public CxxTest::TestSuite
{
public:

    void TestSimpleLinearParabolicSolver1DZeroDirich()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<1> pde;

        // Create zero Dirichlet boundary conditions at first and last node
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode( mesh.GetNumNodes()-1 ), p_boundary_condition);

        // Create PDE solver
        SimpleLinearParabolicSolver<1,1> solver(&mesh,&pde,&bcc);

        // Set start and end times and timestep
        TS_ASSERT_THROWS_THIS(solver.SetTimes(1.0, 0.0),"Start time has to be less than end time");
        TS_ASSERT_THROWS_THIS(solver.SetTimeStep(0.0), "Time step has to be greater than zero");
        double t_end = 0.1;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.01);

        // Set initial condition u(0,x) = sin(x*pi)
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            init_cond[i] = sin(x*M_PI);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE and store solution
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be u = e^{-t*pi*pi} sin(x*pi), t=1
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = exp(-0.1*M_PI*M_PI)*sin(x*M_PI);
            TS_ASSERT_DELTA(result_repl[i], u, 0.1);
        }

        TS_ASSERT_EQUALS(solver.mMatrixIsAssembled, true);
        solver.SetMatrixIsNotAssembled();
        TS_ASSERT_EQUALS(solver.mMatrixIsAssembled, false);

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    void TestSimpleLinearParabolicSolver1DZeroDirichWithSourceTerm()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquationWithSourceTerm<1> pde;

        // Create zero Dirichlet boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        p_boundary_condition = new ConstBoundaryCondition<1>(-0.5);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode( mesh.GetNumNodes()-1 ), p_boundary_condition);

        // Create PDE solver
        SimpleLinearParabolicSolver<1,1> solver(&mesh,&pde,&bcc);

        // Set start and end times and timestep
        double t_end = 0.1;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.01);

        // Set initial condition u(0,x) = sin(x*pi)-0.5*x*x
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            init_cond[i] = sin(x*M_PI)-0.5*x*x;
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE and store solution
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be u = e^{-t*pi*pi} sin(x*pi) + 0.5*x^2, t=1
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = exp(-0.1*M_PI*M_PI)*sin(x*M_PI)-0.5*x*x;
            TS_ASSERT_DELTA(result_repl[i], u, 0.1);
        }

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    void TestSimpleLinearParabolicSolverNonzeroNeumannCondition()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<1> pde;

        // Create boundary conditions u(0)=0, u'(1)=1
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);

        ConstBoundaryCondition<1>* p_neumann_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, p_neumann_boundary_condition);

        // Create PDE solver
        SimpleLinearParabolicSolver<1,1> solver(&mesh,&pde,&bcc);

        // Set start and end times and timestep
        solver.SetTimes(0, 0.5);
        solver.SetTimeStep(0.01);

        // Set initial condition u(0,x) = x + sin(pi*x/2)
        const double PI_over_2 = M_PI/2.0;
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            init_cond[i] = x + sin(PI_over_2*x);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE and store solution
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Check result
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = x + exp(-0.5*PI_over_2*PI_over_2)*sin(x*PI_over_2);
            TS_ASSERT_DELTA(result_repl[i], u, 0.01);
        }

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    void TestSimpleLinearParabolicSolver2DZeroDirich()
    {
        // Read mesh on [0,1]x[0,1]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<2> pde;

        // Create zero Dirichlet boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

        // Create PDE solver
        SimpleLinearParabolicSolver<2,2> solver(&mesh,&pde,&bcc);

        // Set start and end times and timestep
        double t_end = 0.1;

        // Coverage of exception handling
        TS_ASSERT_THROWS_THIS(solver.Solve(),"SetTimes() has not been called");
        solver.SetTimes(0, t_end);
        TS_ASSERT_THROWS_THIS(solver.Solve(),"SetTimeStep() has not been called");
        solver.SetTimeStep(0.001);
        TS_ASSERT_THROWS_THIS(solver.Solve(),"SetInitialCondition() has not been called");

        /*
         * Set initial condition u(0,x,y) = sin(x*pi)*sin(y*pi) as this
         * is an eigenfunction of the heat equation.
         */
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            init_cond[i] = sin(x*M_PI)*sin(y*M_PI);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // Write output to HDF5, then to VTK
        solver.SetOutputToVtk(true);

        // Write output to HDF5, then to parallel VTK (.pvtu)
        solver.SetOutputToParallelVtk(true);

        // Write output to HDF5, then to a .txt format that is readable by Matlab
        solver.SetOutputToTxt(true);

        // Coverage of exception handling
        TS_ASSERT_THROWS_THIS(solver.Solve(),"Output directory or filename prefix has not been set");

        solver.SetOutputDirectoryAndPrefix("TestSimpleLinearParabolicSolver2DZeroDirich","results");

        // Solve PDE and store solution
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Check solution is u = e^{-2*t*pi*pi} sin(x*pi)*sin(y*pi), t=1
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = exp(-2*t_end*M_PI*M_PI)*sin(x*M_PI)*sin(y*M_PI);
            TS_ASSERT_DELTA(result_repl[i], u, 0.01);
        }
        // Test that there is an HDF5 file
        OutputFileHandler file_handler("TestSimpleLinearParabolicSolver2DZeroDirich", false);
        FileFinder h5_file = file_handler.FindFile("results.h5");
        TS_ASSERT(h5_file.Exists());

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestSimpleLinearParabolicSolver2DZeroDirich",
                                                "results", true,
                                                "pde/test/data", "results", false));

        // Test that there is a .vtu file
#ifdef CHASTE_VTK
        FileFinder vtk_file = file_handler.FindFile("vtk_output/results.vtu");
        TS_ASSERT(vtk_file.Exists());
#endif //CHASTE_VTK

        // Test that there are .txt files
        for (unsigned timestep=0; timestep<101; timestep++)
        {
            std::stringstream filename;
            filename << "txt_output/results_Variable_0_" << timestep << ".txt";
            FileFinder txt_file = file_handler.FindFile(filename.str());
            TS_ASSERT(txt_file.Exists());
        }

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    void TestSimpleLinearParabolicSolver2DNonzeroDirichWithSourceTerm()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquationWithSourceTerm<2> pde;

        // Create non-zero Dirichlet boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        while (iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            ConstBoundaryCondition<2>* p_dirichlet_boundary_condition =
                new ConstBoundaryCondition<2>(-0.25*(x*x+y*y));
            bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
            iter++;
        }

        // Create PDE solver
        SimpleLinearParabolicSolver<2,2> solver(&mesh,&pde,&bcc);

        // Set start and end times and timestep
        double t_end = 0.1;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.001);

        // Set initial condition u(0,x,y) = sin(x*pi)*sin(y*pi)-0.25*(x^2+y^2)
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            init_cond[i] = sin(x*M_PI)*sin(y*M_PI)-0.25*(x*x+y*y);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE and store solution
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Check solution is u = e^{-t*2*pi*pi} sin(x*pi) sin(y*pi) - 0.25(x^2+y^2), t=0.1
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = exp(-0.1*2*M_PI*M_PI)*sin(x*M_PI)*sin(y*M_PI)-0.25*(x*x+y*y);
            TS_ASSERT_DELTA(result_repl[i], u, 0.05);
        }

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    ///\todo This test fails with current tolerance
    void xTestSimpleLinearParabolicSolver2DNonzeroDirichletWithSourceTermOnFineMeshWithSmallDt()
    {
        // Create mesh from mesh reader
        FemlabMeshReader<2,2> mesh_reader("mesh/test/data/",
                                          "femlab_square_nodes.dat",
                                          "femlab_square_elements.dat",
                                          "femlab_square_edges.dat");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquationWithSourceTerm<2> pde;

        // Create non-zero Dirichlet boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorEnd();
        while (iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            ConstBoundaryCondition<2>* p_dirichlet_boundary_condition =
                new ConstBoundaryCondition<2>(-0.25*(x*x+y*y));
            bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
            iter++;
        }

        // Create PDE solver
        SimpleLinearParabolicSolver<2,2> solver(&mesh,&pde,&bcc);

        // Set start and end times and timestep
        double t_end = 0.1;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.001);

        // Set initial condition u(0,x,y) = sin(x*pi)*sin(y*pi)-0.25*(x^2+y^2)
        Vec initial_condition = PetscTools::CreateVec(mesh.GetNumNodes());
        double* p_initial_condition;
        VecGetArray(initial_condition, &p_initial_condition);
        int lo, hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNode(global_index)->GetPoint()[0];
            double y = mesh.GetNode(global_index)->GetPoint()[1];
            p_initial_condition[local_index] = sin(x*M_PI)*sin(y*M_PI)-0.25*(x*x+y*y);
        }
        VecRestoreArray(initial_condition, &p_initial_condition);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE and store solution
        Vec result = solver.Solve();

        // Check result
        double* p_result;
        VecGetArray(result, &p_result);

        // Check solution is u = e^{-t*2*pi*pi} sin(x*pi) sin(y*pi) - 0.25(x^2+y^2), t=0.1
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNode(global_index)->GetPoint()[0];
            double y = mesh.GetNode(global_index)->GetPoint()[1];
            double u = exp(-0.1*2*M_PI*M_PI)*sin(x*M_PI)*sin(y*M_PI)-0.25*(x*x+y*y);
            TS_ASSERT_DELTA(p_result[local_index], u, 0.001);
        }
        VecRestoreArray(result, &p_result);

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    void TestSimpleLinearParabolicSolver2DMixedOnCoarseMesh()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<2> pde;

        // Create mixed boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        while (iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];

            if ((fabs(y) < 0.01) || (fabs(y - 1.0) < 0.01) || (fabs(x) < 0.01))
            {
                ConstBoundaryCondition<2>* p_dirichlet_boundary_condition = new ConstBoundaryCondition<2>(x);
                bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
            }

            iter++;
        }

        TetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* p_neumann_boundary_condition =
            new ConstBoundaryCondition<2>(1.0);

        while (surf_iter < mesh.GetBoundaryElementIteratorEnd())
        {
            int node = (*surf_iter)->GetNodeGlobalIndex(0);
            double x = mesh.GetNode(node)->GetPoint()[0];

            if (fabs(x - 1.0) < 0.01)
            {
                bcc.AddNeumannBoundaryCondition(*surf_iter, p_neumann_boundary_condition);
            }
            surf_iter++;
        }

        // Create PDE solver
        SimpleLinearParabolicSolver<2,2> solver(&mesh,&pde,&bcc);

        // Set start and end times and timestep
        double t_end = 0.1;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.01);

        // Set initial condition u(0,x,y) = sin(0.5*M_PI*x)*sin(M_PI*y)+x
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            init_cond[i] = sin(0.5*M_PI*x)*sin(M_PI*y)+x;;
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE and store solution
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Check solution is u = e^{-5/4*M_PI*M_PI*t} sin(0.5*M_PI*x)*sin(M_PI*y)+x, t=0.1
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = exp((-5/4)*M_PI*M_PI*0.1) * sin(0.5*M_PI*x) * sin(M_PI*y) + x;
            TS_ASSERT_DELTA(result_repl[i], u, u*0.15);
        }

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    void TestSimpleLinearParabolicSolver2DMixed()
    {
        // Create mesh from mesh reader
        FemlabMeshReader<2,2> mesh_reader("mesh/test/data/",
                                          "femlab_square_nodes.dat",
                                          "femlab_square_elements.dat",
                                          "femlab_square_edges.dat");

        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<2> pde;

        // Create mixed boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();

        while (iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];

            if ((fabs(y) < 0.01) || (fabs(y - 1.0) < 0.01) || (fabs(x) < 0.01))
            {
                ConstBoundaryCondition<2>* p_dirichlet_boundary_condition =
                    new ConstBoundaryCondition<2>(x);
                bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
            }

            iter++;
        }

        TetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* p_neumann_boundary_condition = new ConstBoundaryCondition<2>(1.0);
        while (surf_iter != mesh.GetBoundaryElementIteratorEnd())
        {
            int node = (*surf_iter)->GetNodeGlobalIndex(0);
            double x = mesh.GetNode(node)->GetPoint()[0];

            if (fabs(x - 1.0) < 0.01)
            {
                bcc.AddNeumannBoundaryCondition(*surf_iter, p_neumann_boundary_condition);
            }

            surf_iter++;
        }

        // Create PDE solver
        SimpleLinearParabolicSolver<2,2> solver(&mesh,&pde,&bcc);

        // Set start and end times and timestep
        solver.SetTimes(0, 0.1);
        solver.SetTimeStep(0.01);

        // Set initial condition u(0,x,y) = sin(0.5*M_PI*x)*sin(M_PI*y)+x
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            init_cond[i] = sin(0.5*M_PI*x)*sin(M_PI*y)+x;
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE and store solution
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Check solution is u = e^{-5/4*M_PI*M_PI*t} sin(0.5*M_PI*x)*sin(M_PI*y)+x, t=0.1
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = exp((-5/4)*M_PI*M_PI*0.1) * sin(0.5*M_PI*x) * sin(M_PI*y) + x;
            TS_ASSERT_DELTA(result_repl[i], u, u*0.1);
        }

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    void TestHeatEquationSolutionDoesntDrift2D()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<2> pde;

        // Create non-zero constant Dirichlet boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        ConstBoundaryCondition<2>* dirichlet_bc = new ConstBoundaryCondition<2>(-84.5);
        while (iter < mesh.GetBoundaryNodeIteratorEnd())
        {
            bcc.AddDirichletBoundaryCondition(*iter, dirichlet_bc);
            iter++;
        }

        // Create PDE solver
        SimpleLinearParabolicSolver<2,2> solver(&mesh,&pde,&bcc);

        // Set start and end times and timestep
        double t_end = 1.0;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.01);

        // Set initial condition
        Vec initial_condition = PetscTools::CreateAndSetVec(mesh.GetNumNodes(), -84.5);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE and store solution
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Check solution is constant throughout the mesh
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            TS_ASSERT_DELTA(result_repl[i],-84.5, 0.0002);
        }

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    void TestHeatEquationSolutionDoesntDrift1D()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<1> pde;

        // Create non-zero constant Dirichlet boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        TetrahedralMesh<1,1>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        ConstBoundaryCondition<1>* dirichlet_bc = new ConstBoundaryCondition<1>(-84.5);
        while (iter < mesh.GetBoundaryNodeIteratorEnd())
        {
            bcc.AddDirichletBoundaryCondition(*iter, dirichlet_bc);
            iter++;
        }

        // Create PDE solver
        SimpleLinearParabolicSolver<1,1> solver(&mesh,&pde,&bcc);

        // Set start and end times and timestep
        double t_end = 1;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.01);

        solver.SetOutputDirectoryAndPrefix("HeatEquation","results");
        solver.SetOutputToTxt(true);
        solver.SetPrintingTimestepMultiple(20);

        // Set initial condition
        Vec initial_condition = PetscTools::CreateAndSetVec(mesh.GetNumNodes(), -84.5);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE and store solution
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Check solution is constant throughout the mesh
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            TS_ASSERT_DELTA(result_repl[i],-84.5, 0.0001);
        }

        FileFinder finder1("HeatEquation/txt_output/results_Variable_0_0.txt", RelativeTo::ChasteTestOutput);
        TS_ASSERT_EQUALS(finder1.Exists(), true);

        FileFinder finder2("HeatEquation/txt_output/results_Variable_0_5.txt", RelativeTo::ChasteTestOutput);
        TS_ASSERT_EQUALS(finder2.Exists(), true);

        FileFinder finder3("HeatEquation/txt_output/results_Variable_0_6.txt", RelativeTo::ChasteTestOutput);
        TS_ASSERT_EQUALS(finder3.Exists(), false);

        FileFinder finder4("HeatEquation/txt_output/results_Variable_0_100.txt", RelativeTo::ChasteTestOutput);
        TS_ASSERT_EQUALS(finder4.Exists(), false);

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    /*
     * Commented out heat equation with 2d mesh and initial condition non-zero at
     * centre, writing out data (doesn't test anything, wanted to see if we get a
     * circular diffusion pattern on such a small mesh, to compare with monodomain
     * with centre stimulus - result doesn't look like a circle).
     * !Need to change the diffusion coefficient to 0.001 if running this!
     */
    void DONOT_TestSimpleLinearParabolicSolver2DZeroNeumannNonZeroInCentre()
    {
        // Read mesh on [0,1]x[0,1]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<2> pde;

        // Create zero Neumann boundary condition
        BoundaryConditionsContainer<2,2,1> bcc;
        TetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* p_neumann_boundary_condition = new ConstBoundaryCondition<2>(0.0);
        while (surf_iter < mesh.GetBoundaryElementIteratorEnd())
        {
            bcc.AddNeumannBoundaryCondition(*surf_iter, p_neumann_boundary_condition);
            surf_iter++;
        }

        // Create PDE solver
        SimpleLinearParabolicSolver<2,2> solver(&mesh,&pde,&bcc);

        /*
         * Set  initial condition u(0,x,y) = sin(x*pi)*sin(y*pi) as this
         * is an eigenfunction of the heat equation.
         */
        Vec initial_condition = PetscTools::CreateVec(mesh.GetNumNodes());

        double* p_initial_condition;
        VecGetArray(initial_condition, &p_initial_condition);

        int lo, hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);

        // Stimulate
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            p_initial_condition[local_index] = 0;
            switch (global_index)
            {
                case 60:
                case 165:
                case 166:
                case 175:
                case 176:
                    p_initial_condition[local_index] = 100;
                    break;
            }
        }
        VecRestoreArray(initial_condition, &p_initial_condition);
        solver.SetInitialCondition(initial_condition);

        // Set start and end times and timestep
        double time = 0;
        double t_end = 0.1;
        double dt = 0.001;

        int time_var_id = 0;
        int heat_var_id = 0;

        ParallelColumnDataWriter* p_test_writer;
        p_test_writer = new ParallelColumnDataWriter("2DHeatEquation", "2DHeatEquation");

        p_test_writer->DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
        time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");

        heat_var_id = p_test_writer->DefineVariable("T","K");
        p_test_writer->EndDefineMode();

        p_test_writer->PutVariable(time_var_id, time);
        p_test_writer->PutVector(heat_var_id, initial_condition);
        p_test_writer->AdvanceAlongUnlimitedDimension();

        Vec result;

        solver.SetTimeStep(0.01);

        while (time < t_end)
        {
            time += dt;
            solver.SetTimes(time, time+dt);

            // Solve PDE and store solution
            result = solver.Solve();

            solver.SetInitialCondition(result);

            p_test_writer->PutVariable(time_var_id, time);
            p_test_writer->PutVector(heat_var_id, result);
            p_test_writer->AdvanceAlongUnlimitedDimension();
        }

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }

    /*
     * Identical problem to TestSimpleLinearParabolicSolver1DZeroDirich(), but uses
     * a time adaptivity controller to increase the timestep after a given time.
     */
    void Test1DProblemUsingTimeAdaptivityController()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<1> pde;

        // Create zero Dicihlet boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode( mesh.GetNumNodes()-1 ), p_boundary_condition);

        // Create PDE solver
        SimpleLinearParabolicSolver<1,1> solver(&mesh,&pde,&bcc);

        ToyController controller;

        // Set start and end times and timestep
        double t_end = 0.1;
        solver.SetTimes(0, t_end);

        // Set initial condition u(0,x) = sin(x*pi)
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            init_cond[i] = sin(x*M_PI);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // In this test we don't set a timestep, but instead give a timestep controller
        solver.SetTimeAdaptivityController(&controller);

        // Solve PDE and store solution
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Check solution is u = e^{-t*pi*pi} sin(x*pi), t=1
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = exp(-0.1*M_PI*M_PI)*sin(x*M_PI);
            TS_ASSERT_DELTA(result_repl[i], u, 0.1);
        }

        TS_ASSERT_EQUALS(solver.mMatrixIsAssembled, true);
        solver.SetMatrixIsNotAssembled();
        TS_ASSERT_EQUALS(solver.mMatrixIsAssembled, false);

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }
};

#endif //_TESTSIMPLELINEARPARABOLICSOLVER_HPP_
