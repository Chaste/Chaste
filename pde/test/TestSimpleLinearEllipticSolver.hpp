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

#ifndef _TESTSIMPLELINEARELLIPTICSOLVER_HPP_
#define _TESTSIMPLELINEARELLIPTICSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "SimplePoissonEquation.hpp"
#include "LinearPdeWithZeroSource.hpp"
#include "EllipticPdeWithLinearSource.hpp"
#include "EllipticPdeWithLinearSourceCylindrical.hpp"
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


class SomePde : public AbstractLinearEllipticPde<2,2>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<2>& rX, Element<2,2>* )
    {
        double dist = sqrt( fabs(rX[0]-0.75)*fabs(rX[0]-0.75) + fabs(rX[1]-0.75)*fabs(rX[1]-0.75));
        return 0.1*exp(-10*dist*dist);
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>&, Element<2,2>* )
    {
        return 0.0;
    }

    c_matrix<double, 2, 2> ComputeDiffusionTerm(const ChastePoint<2>& )
    {
        return identity_matrix<double>(2);
    }
};

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
        double a = (2-cos(2.0))/sin(2.0);
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
    void TestOrdering()
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

    /**
     * #2033
     *
     * Periodic BCs are implemented but not being applied to the linear systems in any assembler.
     * This test will pass if the following two things are done:
     *   (1) The ApplyPeriodicBcsToLinearProblem() line is uncommented in AbstractAssemblerSolverHybrid.cpp
     *   (2) This test is run with '-ksp_type gmres' - periodic BCs are not being implemented at the
     *       moment in such a way that the matrix is kept symmetric, hence CG will diverge.
     *
     * \todo: Alter ApplyPeriodicBcsToLinearProblem() so that it keeps matrix symmetric, then
     *        replace ApplyDirichletToLinearProblem() calls in all assemblers to a call to a new method
     *        ApplyDirichletAndPeriodicBcsToLinearProblem(). Then this test should just pass.
     */
    void dont_Test2dHeatEquationWithPeriodicBcs()
    {
        TetrahedralMesh<2,2> mesh;
        double width = 1.0;
        mesh.ConstructRegularSlabMesh(0.1, width, width);

        // Instantiate PDE object
        SomePde pde;

        std::vector<std::pair<unsigned,unsigned> > identified_nodes; // for the test

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            if (fabs(x) < 1e-6)
            {
                for (unsigned j=0; j<mesh.GetNumNodes(); j++)
                {
                    double x2 = mesh.GetNode(j)->rGetLocation()[0];
                    double y2 = mesh.GetNode(j)->rGetLocation()[1];
                    if ((fabs(x2-width)<1e-6) && (fabs(y-y2)<1e-6))
                    {
                        std::pair<unsigned,unsigned> pair;
                        pair.first = i;
                        pair.second = j;
                        identified_nodes.push_back(pair);

                        bcc.AddPeriodicBoundaryCondition(mesh.GetNode(i), mesh.GetNode(j));
                    }
                }
            }

            if ((fabs(x-0.5)<1e-6) && (fabs(y)<1e-6))
            {
                ConstBoundaryCondition<2>* p_boundary_condition = new ConstBoundaryCondition<2>(0.0);
                bcc.AddDirichletBoundaryCondition(mesh.GetNode(i), p_boundary_condition);
            }
        }

        // Solver
        SimpleLinearEllipticSolver<2,2> solver(&mesh,&pde,&bcc);

        Vec result = solver.Solve();

        // write solution to file for visualisation
        ReplicatableVector res_repl(result);
        OutputFileHandler handler("PeriodicBcs");
        out_stream p_file = handler.OpenOutputFile("result.txt");
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            *p_file << mesh.GetNode(i)->rGetLocation()[0] << " " << mesh.GetNode(i)->rGetLocation()[1] << " " << res_repl[i] << "\n";
        }
        p_file->close();

        // test the periodicity has been enforced
        assert(identified_nodes.size()>0);
        for (unsigned i=0; i<identified_nodes.size(); i++)
        {
            TS_ASSERT_DELTA(res_repl[identified_nodes[i].first], res_repl[identified_nodes[i].second], 1e-8);
        }

        PetscTools::Destroy(result);
    }

    /*
     * Solve Del^2(u) = -4 in cylindrical polar co-ordinates
     * r in [0, 1], z in [0, Z], axisymmetric, so 2D.
     *
     *
     * BCs
     * u(1) = 0
     * du(0)/dR = 0
     * du(Z)/dz = 0
     * du(0)/dz = 1 - r^2
     *
     * u(R) => u(r,z) = 1-r^2
     *
     * We've done this here by simply hacking the PDE,
     * compare EllipticPdeWithLinearSource.hpp with
     * EllipticPdeWithLinearSourceCylindrical.hpp to see
     * how 'r' enters the terms that the finite element integrals see.
     */
    void TestCylinderAxisymmetricWithConstantSourceTerm()
    {
        std::cout << "Space step\tMax error\tTotal error\n";
        double max_error;
        double total_error;
        RunAnalyticProblemWithSpaceStep(0.025, max_error, total_error);
        TS_ASSERT_LESS_THAN(max_error, 1e-6);
        TS_ASSERT_LESS_THAN(total_error, 1e-4);
    }

    /*
     * Code to run a convergence analysis on the above problem.
     *
     * Not running every time as it takes a while...
     */
    void dontTestConvergenceOfCylinderAxisymmetricWithConstantSourceTerm()
    {
        std::cout << "Space step\tMax error\tTotal error\n";

        double max_error;
        double total_error; // Integrated up over small volume, so smaller than max error!

        for (unsigned i=0; i<7u; i++)
        {
            double node_spacing_in_mesh = 0.1/(double)(SmallPow(2u,i));

            RunAnalyticProblemWithSpaceStep(node_spacing_in_mesh, max_error, total_error);
        }
    }

private:
    /*
     * Method used by the tests:
     * TestCylinderAxisymmetricWithConstantSourceTerm
     * and
     * TestConvergenceOfCylinderAxisymmetricWithConstantSourceTerm
     */
    void RunAnalyticProblemWithSpaceStep(double node_spacing_in_mesh,
                                         double& rMaxError,
                                         double& rTotalError)
    {
        DistributedTetrahedralMesh<2,2> mesh;
        double R = 1.0;
        double Z = 1.0;
        mesh.ConstructRegularSlabMesh(node_spacing_in_mesh, R /*length*/, Z /*width*/);

        // Instantiate PDE object
        EllipticPdeWithLinearSourceCylindrical pde(0.0, 4.0);

        // Boundary conditions
        // (default is
        BoundaryConditionsContainer<2,2,1> bcc;
        ConstBoundaryCondition<2>* p_boundary_condition_1 = new ConstBoundaryCondition<2>(0.0);

        for (TetrahedralMesh<2u,2u>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorBegin();
             node_iter != mesh.GetBoundaryNodeIteratorEnd();
             ++node_iter)
        {
            double r = (*node_iter)->rGetLocation()[0];
            double z = (*node_iter)->rGetLocation()[1];

            if (fabs(r-R) < 1e-6)
            {
                //std::cout << "Applying R=0 BC" << std::endl;
                bcc.AddDirichletBoundaryCondition(*node_iter, p_boundary_condition_1);
            }

            if (fabs(z) < 1e-6)
            {
                //std::cout << "Applying 1-r^2 BC" << std::endl;
                ConstBoundaryCondition<2>* p_boundary_condition_2 = new ConstBoundaryCondition<2>(1 - r*r);
                bcc.AddDirichletBoundaryCondition(*node_iter, p_boundary_condition_2);
            }
        }

        // Solver
        SimpleLinearEllipticSolver<2,2> solver(&mesh,&pde,&bcc);

        /// \todo #2748
        //solver.SetOutputDirectoryAndPrefix("CylindricalAxisymmetricTestProblem","results");
        //solver.SetOutputToVtk(true);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be u(r,z) = 1 - r^2 (for all z)
        rMaxError = 0.0;
        rTotalError = 0.0;

        for (TetrahedralMesh<2u,2u>::ElementIterator elem_iter = mesh.GetElementIteratorBegin();
                elem_iter != mesh.GetElementIteratorEnd();
                ++elem_iter)
        {
            // Get analytic result at the centroid
            c_vector<double, 2u> location = elem_iter->CalculateCentroid();
            double analytic_result = 1.0 - (location[0]*location[0]);

            // Get simulation result interpolated to centroid too.
            unsigned node_0 = elem_iter->GetNodeGlobalIndex(0u);
            unsigned node_1 = elem_iter->GetNodeGlobalIndex(1u);
            unsigned node_2 = elem_iter->GetNodeGlobalIndex(2u);
            double sim_result = (result_repl[node_0] + result_repl[node_1] + result_repl[node_2])/3.0;

            // Error integrated over triangle (roughly)
            double error = 0.5*node_spacing_in_mesh*node_spacing_in_mesh*
                    sqrt((sim_result - analytic_result)*(sim_result - analytic_result));

            if (error > rMaxError)
            {
                rMaxError = error;
            }
            rTotalError += error;
        }

        PetscTools::Destroy(result);

//        // Output to file for visualizing somewhere
//        OutputFileHandler handler("CylindricalAxisymmetric", false);
//        out_stream results_file = handler.OpenOutputFile("results.dat");
//
//        for (TetrahedralMesh<2u,2u>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
//                     node_iter != mesh.GetNodeIteratorEnd();
//                     ++node_iter)
//        {
//            *results_file << node_iter->rGetLocation()[0] << "\t" << node_iter->rGetLocation()[1] << "\t" <<  result_repl[node_iter->GetIndex()] << std::endl;
//        }
//        results_file->close();

        std::cout << node_spacing_in_mesh << "\t" << rMaxError << "\t" << rTotalError << std::endl << std::flush;
    }
};

#endif //_TESTSIMPLELINEARELLIPTICSOLVER_HPP_
