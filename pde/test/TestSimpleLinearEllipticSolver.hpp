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

#ifndef _TESTSIMPLELINEARELLIPTICSOLVER_HPP_
#define _TESTSIMPLELINEARELLIPTICSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "SimplePoissonEquation.hpp"
#include "LinearPdeWithZeroSource.hpp"
#include "EllipticPdeWithLinearSource.hpp"
#include "EllipticPdeWithRadialLinearSource.hpp"
#include "SimpleLinearEllipticSolver.hpp"
#include <vector>
#include <cmath>
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "VaryingDiffusionAndSourceTermPde.hpp"
#include "TrianglesMeshReader.hpp"

/*
 * These are need for the nD problems in mD space (n!=m), as those
 * particular cases are not explicitly instantiated.
 */
#include "AbstractBoundaryConditionsContainerImplementation.hpp"
#include "BoundaryConditionsContainerImplementation.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestSimpleLinearEllipticSolver : public CxxTest::TestSuite
{
public:

    void TestWithPoissonsEquationAndMeshReader()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/trivial_1d_mesh");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        SimplePoissonEquation<1,1> pde;

        double value1 = pde.ComputeConstantInUSourceTermAtNode(*(mesh.GetNode(0)));
        double value2 = pde.ComputeConstantInUSourceTerm(mesh.GetNode(0)->GetPoint(), NULL);
        TS_ASSERT_DELTA(value1, value2, 1e-10);

        value1 = pde.ComputeLinearInUCoeffInSourceTermAtNode(*(mesh.GetNode(0)));
        value2 = pde.ComputeLinearInUCoeffInSourceTerm(mesh.GetNode(0)->GetPoint(), NULL);
        TS_ASSERT_DELTA(value1, value2, 1e-10);

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);

        // Solver
        SimpleLinearEllipticSolver<1,1> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be u = 0.5*x*(3-x)
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = 0.5*x*(3-x);
            TS_ASSERT_DELTA(result_repl[i], u, 0.001);
        }

        PetscTools::Destroy(result);
    }

    void TestWithHeatEquation2()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_mesh_5_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        SimplePoissonEquation<1,1> pde;

        // Boundary conditions u(-1)=1, u'(-3)=0
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);

        // Add Neumann condition to the left hand end
        ConstBoundaryCondition<1>* p_neumann_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, p_neumann_boundary_condition);

        // Solver
        SimpleLinearEllipticSolver<1,1> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = 1 - 0.5*(x+1)*(5+x);
            TS_ASSERT_DELTA(result_repl[i], u, 0.001);
        }

        PetscTools::Destroy(result);
    }

    void TestWithHeatEquationNonzeroNeumannCondition()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_mesh_5_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        SimplePoissonEquation<1,1> pde;

        // Boundary conditions u'(-3)=1, u(-1)=1
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], -1, 1e-12);

        // Note we pass -1 not 1; see comment for AddNeumannBoundaryCondition
        ConstBoundaryCondition<1>* p_neumann_boundary_condition = new ConstBoundaryCondition<1>(-1.0);

        // Add Neumann condition to the left hand end
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, p_neumann_boundary_condition);

        // Solver
        SimpleLinearEllipticSolver<1,1> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = -0.5*x*x - 2*x - 0.5;
            TS_ASSERT_DELTA(result_repl[i], u, 0.001);
        }

        PetscTools::Destroy(result);
    }

    void Test2dHeatEquationOnUnitSquare()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        DistributedTetrahedralMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        SimplePoissonEquation<2,2> pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        ConstBoundaryCondition<2>* p_boundary_condition = new ConstBoundaryCondition<2>(0.0);
        DistributedVectorFactory* p_factory = mesh.GetDistributedVectorFactory();
        for (unsigned i=0; i<4; i++)
        {
            if (p_factory->IsGlobalIndexLocal(i))
            {
                bcc.AddDirichletBoundaryCondition(mesh.GetNode(i), p_boundary_condition);
            }
        }

        // Solver
        SimpleLinearEllipticSolver<2,2> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);
        TS_ASSERT_DELTA(result_repl[4], 1.0/12.0, 0.001);

        PetscTools::Destroy(result);
    }

    void TestHeatEquationWithNeumannOnUnitDisc()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        SimplePoissonEquation<2,2> pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        // du/dn = -0.5 on r=1
        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* p_boundary_condition;
        p_boundary_condition = new ConstBoundaryCondition<2>(-0.5);
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
            iter++;
        }
        // u = 2 at some point on the boundary, say node 1
        p_boundary_condition = new ConstBoundaryCondition<2>(2.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(1), p_boundary_condition);

        // Solver
        SimpleLinearEllipticSolver<2,2> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            c_vector<double, 2> r;
            r(0) = mesh.GetNode(i)->GetPoint()[0];
            r(1) = mesh.GetNode(i)->GetPoint()[1];
            double u = -0.25 * inner_prod(r, r) + 2.25;
            TS_ASSERT_DELTA(result_repl[i], u, 0.01);
        }

        PetscTools::Destroy(result);
    }

    void TestVaryingPdeAndMeshReader1D()
    {
        /// Create mesh from mesh reader \todo set to correct mesh file?
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_mesh_1_to_3");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        VaryingDiffusionAndSourceTermPde<1> pde;

        // Boundary conditions u(1)=4
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_dirichlet_condition = new ConstBoundaryCondition<1>(4.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_dirichlet_condition);

        // Note we need to specify D * du/dx for the Neumann boundary condition
        ConstBoundaryCondition<1>* p_neumann_boundary_condition = new ConstBoundaryCondition<1>(7.0*9.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, p_neumann_boundary_condition);

        // Solver
        SimpleLinearEllipticSolver<1,1> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = -(x*x*x/12.0)-(333/(4*x))+4+1000.0/12.0;
            TS_ASSERT_DELTA(result_repl[i], u, 0.2);
        }

        PetscTools::Destroy(result);
    }

    /**
     * Test a simple PDE with nasty boundary conditions - du/dn has a discontinuity.
     * This test takes a little while to run, so is currently disabled.
     */
    void longTestKathrynHarrimanPage67EqFourPointOne()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        LinearPdeWithZeroSource<2> pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        // u = 0 on r<=1, z=0
        ConstBoundaryCondition<2>* p_boundary_dirichlet_condition = new ConstBoundaryCondition<2>(0.0);
        TetrahedralMesh<2,2>::BoundaryNodeIterator iter1 = mesh.GetBoundaryNodeIteratorBegin();
        while (iter1 != mesh.GetBoundaryNodeIteratorEnd())
        {
            if ((*iter1)->GetPoint()[0] <= 1.0 && fabs((*iter1)->GetPoint()[1]) < 0.0001)
            {
                bcc.AddDirichletBoundaryCondition(*iter1, p_boundary_dirichlet_condition);
            }
            iter1++;
        }
        // du/dn = 0 on r>1, z=0 and on r=0, z>=0
        TetrahedralMesh<2,2>::BoundaryElementIterator iter2 = mesh.GetBoundaryElementIteratorBegin();
        while (iter2 != mesh.GetBoundaryElementIteratorEnd())
        {
            // Condition is zero, so we don't actually have to do anything.
            iter2++;
        }
        // u=1 as r,z->infinity. We replace this by the exact solution on r=2 and on z=2
        iter1 = mesh.GetBoundaryNodeIteratorBegin();
        while (iter1 != mesh.GetBoundaryNodeIteratorEnd())
        {
            double r = (*iter1)->GetPoint()[0];
            double z = (*iter1)->GetPoint()[1];
            if (fabs(r - 2.0) <= 0.0001 || fabs(z - 2.0) < 0.0001)
            {
                double u = 1 - 2.0/M_PI*asin(2.0 / (sqrt(z*z + (1+r)*(1+r))
                                                    + sqrt(z*z + (1-r)*(1-r))));
                p_boundary_dirichlet_condition = new ConstBoundaryCondition<2>(u);
                bcc.AddDirichletBoundaryCondition(*iter1, p_boundary_dirichlet_condition);
            }
            iter1++;
        }

        // Solver
        SimpleLinearEllipticSolver<2,2> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();

        // Check result
        double* p_result;
        int lo, hi;
        VecGetOwnershipRange(result, &lo, &hi);
        VecGetArray(result, &p_result);
        for (unsigned global_index=0; global_index < mesh.GetNumNodes(); global_index++)
        {
            int local_index=global_index - lo;
            double r = mesh.GetNode(global_index)->GetPoint()[0];
            double z = mesh.GetNode(global_index)->GetPoint()[1];
            double u;
            if (z > 1e-12)
            {
                u = 1 - 2.0/M_PI*asin(2.0 / (sqrt(z*z + (1+r)*(1+r))
                                             + sqrt(z*z + (1-r)*(1-r))));
            }
            else if (r > 1.0)
            {
                u = 1 - 2.0/M_PI * asin(1.0/r);
            }
            else
            {
                u = 0;
            }
            TS_ASSERT_DELTA(p_result[local_index], u, 0.08);
        }

        VecRestoreArray(result, &p_result);
        PetscTools::Destroy(result);
    }

    // Test 3d data
    void Test3dEllipticEquationDirichletCondition()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        SimplePoissonEquation<3,3> pde;

        // Boundary conditions
        BoundaryConditionsContainer<3,3,1> bcc;
        TetrahedralMesh<3,3>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();

        while (iter < mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            double z = (*iter)->GetPoint()[2];

            ConstBoundaryCondition<3>* p_dirichlet_boundary_condition = new ConstBoundaryCondition<3>(-1.0/6*(x*x+y*y+z*z));
            bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
            iter++;
        }

        // Solver
        SimpleLinearEllipticSolver<3,3> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be -1/6*(x^2 + y^2 +z^2)
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double z = mesh.GetNode(i)->GetPoint()[2];
            double u = -1.0/6 * (x*x+y*y+z*z);
            TS_ASSERT_DELTA(result_repl[i], u, 0.01);
        }

        PetscTools::Destroy(result);
    }

    // Test 3d data
    void Test3dEllipticEquationNeumannCondition()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        SimplePoissonEquation<3,3> pde;

        // Boundary conditions
        BoundaryConditionsContainer<3,3,1> bcc;
        TetrahedralMesh<3,3>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();

        while (iter < mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            double z = (*iter)->GetPoint()[2];

            if (fabs(1-x)>=0.01)
            {
                //Dirichlet boundary condition
                ConstBoundaryCondition<3>* p_dirichlet_boundary_condition = new ConstBoundaryCondition<3>(-1.0/6*(x*x+y*y+z*z));
                bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
            }
            iter++;
        }

        TetrahedralMesh<3,3>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<3>* p_neumann_boundary_condition = new ConstBoundaryCondition<3>(-1.0/3);
        while (surf_iter < mesh.GetBoundaryElementIteratorEnd())
        {
            int node = (*surf_iter)->GetNodeGlobalIndex(0);
            double x = mesh.GetNode(node)->GetPoint()[0];
            // double y = mesh.GetNode(node)->GetPoint()[1];

            if (fabs(x - 1.0) < 0.01)
            {
                bcc.AddNeumannBoundaryCondition(*surf_iter, p_neumann_boundary_condition);
            }

            surf_iter++;
        }

        // Solver
        SimpleLinearEllipticSolver<3,3> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be -1/6*(x^2 + y^2 +z^2)
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double z = mesh.GetNode(i)->GetPoint()[2];
            double u = -1.0/6 * (x*x+y*y+z*z);
            TS_ASSERT_DELTA(result_repl[i], u, 0.1);
        }

        PetscTools::Destroy(result);
    }

    // Solve u_xx + 4*u = 0, u(0)=1, u(1)=2 => u = a sin(2x) + cos(2x), where a = (2-cos2)/sin2
    void TestWithLinearSourceTerm()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        EllipticPdeWithLinearSource<1> pde(4,0);

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);

        p_boundary_condition = new ConstBoundaryCondition<1>(2.0);
        unsigned last_node = (unsigned)(mesh.GetNumNodes()-1);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(last_node), p_boundary_condition);

        // Solver
        SimpleLinearEllipticSolver<1,1> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be u = a sin(x) + cos(x), where a = (2-cos1)/sin1
        double a = (2-cos(2))/sin(2);
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = a*sin(2*x) + cos(2*x);
            //std::cout << u << " " << result_repl[i] << "\n";
            TS_ASSERT_DELTA(result_repl[i], u, u*0.001);
        }

        PetscTools::Destroy(result);
    }

    /*
     * Picking the solution u=exp(xy), we solve the PDE u_xx + u_yy = (x^2+y^2) u,
     * with BCs u = exp(xy) on the boundary.
     */
    void TestWithLinearSourceTerm2d()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        EllipticPdeWithRadialLinearSource pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        for (TetrahedralMesh<2,2>::BoundaryNodeIterator iter =
              mesh.GetBoundaryNodeIteratorBegin();
            iter != mesh.GetBoundaryNodeIteratorEnd();
            iter++)
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            double val = exp(x*y);
            ConstBoundaryCondition<2>* p_boundary_condition = new ConstBoundaryCondition<2>(val);
            bcc.AddDirichletBoundaryCondition(*iter, p_boundary_condition);
        }

        // Solver
        SimpleLinearEllipticSolver<2,2> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be u = exp(xy)
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = exp(x*y);
            //std::cout << u << " " << result_repl[i] << "\n";
            TS_ASSERT_DELTA(result_repl[i], u, u*0.01);
        }

        PetscTools::Destroy(result);
    }

    /*
     * Test that the solver can read an ordering file and assign the correct
     * number of nodes to each processor.
     */
    void TestOrdering() throw(Exception)
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        if (PetscTools::GetNumProcs() == 2)
        {
            mesh.ReadNodesPerProcessorFile("mesh/test/data/nodes_per_processor_1.txt");

            // Instantiate PDE and BCC object, though not used
            SimplePoissonEquation<2,2> pde;
            BoundaryConditionsContainer<2,2,1> bcc;

            // Solver
            SimpleLinearEllipticSolver<2,2> solver(&mesh,&pde,&bcc);

            // test set up correctly
            PetscInt petsc_lo, petsc_hi;
            Vec vector = mesh.GetDistributedVectorFactory()->CreateVec();

            VecGetOwnershipRange(vector,&petsc_lo,&petsc_hi);

            if (PetscTools::GetMyRank() == 0)
            {
                TS_ASSERT_EQUALS(0, petsc_lo);
                TS_ASSERT_EQUALS(2, petsc_hi);
            }
            else
            {
                TS_ASSERT_EQUALS(2, petsc_lo);
                TS_ASSERT_EQUALS(5, petsc_hi);
            }
        }
    }

    // Solve a 1d problem in 2d space yet
    void TestWithPoissonsEquation1dMeshIn2dSpace()
    {
        const unsigned SPACE_DIM = 2;
        const unsigned ELEMENT_DIM = 1;

        // Create mesh from mesh reader
        TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader("mesh/test/data/trivial_1d_in_2d_mesh");
        TetrahedralMesh<ELEMENT_DIM,SPACE_DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        SimplePoissonEquation<ELEMENT_DIM,SPACE_DIM> pde;

        // Boundary conditions (u=0 on one end, u'=0 on other end)
        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1> bcc;
        ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition = new ConstBoundaryCondition<SPACE_DIM>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);

        // Solver
        SimpleLinearEllipticSolver<ELEMENT_DIM,SPACE_DIM> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be u = 0.5*x*(3-x)
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = 0.5*x*(3-x);
            TS_ASSERT_DELTA(result_repl[i], u, 0.001);
        }

        // Coverage
        TS_ASSERT(solver.GetLinearSystem()!=NULL);

        PetscTools::Destroy(result);
    }
};

#endif //_TESTSIMPLELINEARELLIPTICSOLVER_HPP_
