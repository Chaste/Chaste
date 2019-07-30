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

#ifndef _TESTSIMPLENONLINEARELLIPTICSOLVER_HPP_
#define _TESTSIMPLENONLINEARELLIPTICSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include <petsc.h>
#include <vector>
#include <cmath>

#include "SimpleNonlinearEllipticSolver.hpp"
#include "SimplePetscNonlinearSolver.hpp"

#include "BoundaryConditionsContainer.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "ConstBoundaryCondition.hpp"

#include "NonlinearEquationPde.hpp"
#include "NonlinearEquation2Pde.hpp"
#include "NonlinearEquation3Pde.hpp"
#include "NonlinearEquation4Pde.hpp"
#include "NonlinearEquation5Pde.hpp"
#include "Example2DNonlinearEllipticPde.hpp"
#include "NonlinearLinearEquation.hpp"
#include "ExampleNasty2dNonlinearEllipticPde.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"

/**
 * For use in TestSimpleNonlinearEllipticSolver::Test2dOnUnitSquare.
 */
double bc_x1_func(const ChastePoint<2>& p)
{
    return 2*(2+p[1]*p[1]);
}
/**
 * For use in TestSimpleNonlinearEllipticSolver::Test2dOnUnitSquare.
 */
double bc_y1_func(const ChastePoint<2>& p)
{
    return 2*(2+p[0]*p[0]);
}

/**
 * For use in TestSimpleNonlinearEllipticSolver::TestNasty2dEquationOnUnitSquare.
 */
double bc_x1_func2(const ChastePoint<2>& p)
{
    return sin(2.0)*(sin(1.0)*sin(1.0)+1+p[1]*p[1]);
}
/**
 * For use in TestSimpleNonlinearEllipticSolver::TestNasty2dEquationOnUnitSquare.
 */
double bc_y1_func2(const ChastePoint<2>& p)
{
    return 2*(2+sin(p[0])*sin(p[0]));
}

/**
 * For use in TestSimpleNonlinearEllipticSolver::TestWithHeatEquation2DAndNeumannBCs
 */
double one_bc(const ChastePoint<2>& p)
{
    return p[1];
}

class TestSimpleNonlinearEllipticSolver : public CxxTest::TestSuite
{
public:
    void TestNumericalAgainstAnalyticJacobian()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquationPde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);

        // Solver
        SimpleNonlinearEllipticSolver<1,1> solver(&mesh, &pde, &bcc);

        TS_ASSERT( solver.VerifyJacobian(1e-3) );
    }

    void TestWithHeatEquation1D()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquationPde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);

        SimpleNonlinearEllipticSolver<1,1> solver(&mesh, &pde, &bcc);

        // Set up initial guess
        Vec initial_guess = PetscTools::CreateAndSetVec(mesh.GetNumNodes(),1.0);

        Vec answer = solver.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = sqrt(x*(1-x));
            TS_ASSERT_DELTA(answer_repl[i], u, 0.001);
        }

        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer);
    }

    void TestWithHeatEquation1DAndNeumannBCs()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquationPde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        // u(0) = 0
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        // u(1)*u'(1) = 1
        p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

        // Nonlinear solver to use
        SimpleNonlinearEllipticSolver<1,1> solver(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = PetscTools::CreateAndSetVec(mesh.GetNumNodes(),0.25);

        // Solve the PDE
        Vec answer = solver.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = sqrt(x*(4-x));
            TS_ASSERT_DELTA(answer_repl[i], u, 0.001);
        }

        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer);
    }

    void TestWithHeatEquation1D2()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquation2Pde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        p_boundary_condition = new ConstBoundaryCondition<1>(exp(1.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);

        SimpleNonlinearEllipticSolver<1,1> solver(&mesh, &pde, &bcc);

        // Set up initial Guess
        std::vector<double> init_guess(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            init_guess[i] = 1.0+0.01*i*i;
        }
        Vec initial_guess = PetscTools::CreateVec(init_guess);

        Vec answer = solver.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = exp(0.5*(3.0*x-x*x));
            TS_ASSERT_DELTA(answer_repl[i], u, 0.001);
        }

        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer);
    }

    void TestWithHeatEquation1D3()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquation3Pde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(sqrt(2.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);

        SimpleNonlinearEllipticSolver<1,1> solver(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = PetscTools::CreateVec(mesh.GetNumNodes());
        for (unsigned global_index=0; global_index<mesh.GetNumNodes(); global_index++)
        {
            PetscVecTools::SetElement(initial_guess, global_index, 1.5-0.15*global_index);
        }
        PetscVecTools::Finalise(initial_guess);

        Vec answer = solver.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            // note: 2.0*(exp(-x)-x*exp(-1.0) is always >= 0 for these values of x,
            // however when x=1, 2.0*(exp(-x)-x*exp(-1.0) should be
            // zero but can end up very slightly negative (-1e-17, for example), hence
            // the need for the fabs before the sqrt.
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = sqrt(fabs(2.0*(exp(-x)-x*exp(-1.0))));
            TS_ASSERT_DELTA(answer_repl[i], u, 0.001);
        }

        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer);
    }

    void TestWithHeatEquation1D4()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquation4Pde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        // u(1) = exp(1.0)
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(exp(-1.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);
        // u(0)^2*u'(0) = 0.0
        p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

        SimpleNonlinearEllipticSolver<1,1> solver(&mesh, &pde, &bcc);

        // Set up initial Guess
        std::vector<double> init_guess(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x1=0.1*(double)(i);
            init_guess[i] =  0.35*(1-x1*x1);
        }
        Vec initial_guess = PetscTools::CreateVec(init_guess);

        Vec answer = solver.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = x*exp(-x);
            TS_ASSERT_DELTA(answer_repl[i], u, 0.01);
        }

        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer);
    }

    void TestWithHeatEquation1D5()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquation5Pde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        // u(1) = exp(-1.0)
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(exp(-1.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);
        // u(0)^2*u'(0) = -1.0
        // Note that we specify 1 as the value, since figuring out which direction
        // the normal is in is hard in 1D.
        p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

        SimpleNonlinearEllipticSolver<1,1> solver(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = PetscTools::CreateVec(mesh.GetNumNodes());
        double x1;
        for (unsigned global_index=0; global_index<mesh.GetNumNodes(); global_index++)
        {
            x1=0.1*(double)(global_index);
            PetscVecTools::SetElement(initial_guess, global_index, 0.35*(1-x1*x1));
        }
        PetscVecTools::Finalise(initial_guess);

        Vec answer = solver.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = exp(-x);
            TS_ASSERT_DELTA(answer_repl[i], u, 0.01);
        }

        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer);
    }

    void TestWithHeatEquation1DAndNeumannBCs2()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquationPde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        // u(1) = sqrt(3.0)
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(sqrt(3.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);
        // u(0)*u'(0) = 2
        // Note that we specify -2 as the value, since figuring out which direction
        // the normal is in is hard in 1D.
        p_boundary_condition = new ConstBoundaryCondition<1>(-2.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

        SimpleNonlinearEllipticSolver<1,1> solver(&mesh, &pde, &bcc);

        // cover the bad size exception
        Vec badly_sized_init_guess = PetscTools::CreateAndSetVec(1,1.0); // size=1
        TS_ASSERT_THROWS_THIS( solver.Solve(badly_sized_init_guess, true),
                "Size of initial guess vector, 1, does not match size of problem, 11" );

        // Set up initial Guess
        Vec initial_guess = PetscTools::CreateAndSetVec(mesh.GetNumNodes(),1.0);

        // This problem seems unusally sensitive to the initial guess. Various other
        // choices failed to converge.
        Vec answer = solver.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = sqrt(x*(4-x));
            TS_ASSERT_DELTA(answer_repl[i], u, 0.002); //Error is likely to accumulate at the Neumann boundary (x=0)
        }

        PetscTools::Destroy(badly_sized_init_guess);
        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer);
    }

    void TestHeatEquationWithNeumannOnUnitDisc()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearLinearEquation<2> pde;

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

        SimpleNonlinearEllipticSolver<2,2> solver(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = PetscTools::CreateAndSetVec(mesh.GetNumNodes(),1.0);

        Vec answer = solver.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            c_vector<double, 2> r;
            r(0) = mesh.GetNode(i)->GetPoint()[0];
            r(1) = mesh.GetNode(i)->GetPoint()[1];
            double u = -0.25 * inner_prod(r, r) + 2.25;
            TS_ASSERT_DELTA(answer_repl[i], u, 0.01);
        }

        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer);
    }

    void TestWithHeatEquation2DAndNeumannBCs()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquationPde<2> pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        // u(y=0) = 0
        ConstBoundaryCondition<2>* zero_boundary_condition = new ConstBoundaryCondition<2>(0.0);
        TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorBegin();
        while (node_iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            if (fabs((*node_iter)->GetPoint()[1]) < 1e-12)
            {
                bcc.AddDirichletBoundaryCondition(*node_iter, zero_boundary_condition);
            }
            node_iter++;
        }

        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        FunctionalBoundaryCondition<2>* one_boundary_condition = new FunctionalBoundaryCondition<2>(&one_bc);
        AbstractBoundaryCondition<2>* p_boundary_condition;
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            double y = (*iter)->GetNodeLocation(0,1);
            if (fabs(y-1.0) < 1e-12)
            {
                // u(y=1)*u'(y=1) = 1
                p_boundary_condition = one_boundary_condition;
            }
            else
            {
                // No flux across left & right
                p_boundary_condition = zero_boundary_condition;
            }

            bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

            iter++;
        }

        SimpleNonlinearEllipticSolver<2,2> solver(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = PetscTools::CreateAndSetVec(mesh.GetNumNodes(),0.25);

        // Solve
        Vec answer = solver.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = sqrt(y*(4-y));
            TS_ASSERT_DELTA(answer_repl[i], u, 0.15);
        }

        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer);
    }

    void Test2dOnUnitSquare()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        Example2DNonlinearEllipticPde pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        ConstBoundaryCondition<2>* p_boundary_condition;
        TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorBegin();
        while (node_iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*node_iter)->GetPoint()[0];
            double y = (*node_iter)->GetPoint()[1];
            p_boundary_condition = NULL;
            if (fabs(x) < 1e-12)
            {
                // On x=0, u=1+y^2
                p_boundary_condition = new ConstBoundaryCondition<2>(1 + y*y);
            }
            else if (fabs(y) < 1e-12)
            {
                // On y=0, u=1+x^2
                p_boundary_condition = new ConstBoundaryCondition<2>(1 + x*x);
            }
            if (p_boundary_condition)
            {
                bcc.AddDirichletBoundaryCondition(*node_iter, p_boundary_condition);
            }

            node_iter++;
        }
        FunctionalBoundaryCondition<2>* p_functional_bc;
        TetrahedralMesh<2,2>::BoundaryElementIterator elt_iter = mesh.GetBoundaryElementIteratorBegin();
        while (elt_iter != mesh.GetBoundaryElementIteratorEnd())
        {
            double x = (*elt_iter)->GetNodeLocation(0,0);
            double y = (*elt_iter)->GetNodeLocation(0,1);
            p_functional_bc = NULL;
            if (fabs(y-1.0) < 1e-12)
            {
                // On y=1, Dgradu_dot_n = 2(2+x^2)
                p_functional_bc = new FunctionalBoundaryCondition<2>(&bc_y1_func);
            }
            else if (fabs(x-1.0) < 1e-12)
            {
                // On x=1, Dgradu_dot_n = 2(2+y^2)
                p_functional_bc = new FunctionalBoundaryCondition<2>(&bc_x1_func);
            }
            if (p_functional_bc)
            {
                bcc.AddNeumannBoundaryCondition(*elt_iter, p_functional_bc);
            }

            elt_iter++;
        }

        SimpleNonlinearEllipticSolver<2,2> solver(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = PetscTools::CreateAndSetVec(mesh.GetNumNodes(),4.0);

        // Numerical Jacobian
        Vec answer = solver.Solve(initial_guess, false);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = 1 + x*x + y*y;
            TS_ASSERT_DELTA(answer_repl[i], u, 0.01);
        }

        PetscTools::Destroy(answer);

        // Analytical Jacobian
        answer=solver.Solve(initial_guess, true);
        ReplicatableVector answer_repl2(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = 1 + x*x + y*y;
            TS_ASSERT_DELTA(answer_repl2[i], u, 0.01);
        }

        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer);
    }

    void TestNasty2dEquationOnUnitSquare()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        ExampleNasty2dNonlinearEllipticPde pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        ConstBoundaryCondition<2>* p_boundary_condition;
        TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorBegin();
        while (node_iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*node_iter)->GetPoint()[0];
            double y = (*node_iter)->GetPoint()[1];
            p_boundary_condition = NULL;
            if (fabs(x) < 1e-12)
            {
                // On x=0, u=1+y^2
                p_boundary_condition = new ConstBoundaryCondition<2>(1 + y*y);
            }
            else if (fabs(y) < 1e-12)
            {
                // On y=0, u=1+sin^2(x)
                p_boundary_condition = new ConstBoundaryCondition<2>(1 + sin(x)*sin(x));
            }
            if (p_boundary_condition)
            {
                bcc.AddDirichletBoundaryCondition(*node_iter, p_boundary_condition);
            }

            node_iter++;
        }
        FunctionalBoundaryCondition<2>* p_functional_bc;
        TetrahedralMesh<2,2>::BoundaryElementIterator elt_iter = mesh.GetBoundaryElementIteratorBegin();
        while (elt_iter != mesh.GetBoundaryElementIteratorEnd())
        {
            double x = (*elt_iter)->GetNodeLocation(0,0);
            double y = (*elt_iter)->GetNodeLocation(0,1);
            p_functional_bc = NULL;
            if (fabs(y-1.0) < 1e-12)
            {
                // On y=1, Dgradu_dot_n = 2(2+sin^2(x))
                p_functional_bc = new FunctionalBoundaryCondition<2>(&bc_y1_func2);
            }
            else if (fabs(x-1.0) < 1e-12)
            {
                // On x=1, Dgradu_dot_n = sin(2.0)(sin^2(1)+1+y^2)
                p_functional_bc = new FunctionalBoundaryCondition<2>(&bc_x1_func2);
            }
            if (p_functional_bc)
            {
                bcc.AddNeumannBoundaryCondition(*elt_iter, p_functional_bc);
            }

            elt_iter++;
        }

        SimpleNonlinearEllipticSolver<2,2> solver(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = PetscTools::CreateAndSetVec(mesh.GetNumNodes(),4.0);

        // Numerical Jacobian
        Vec answer = solver.Solve(initial_guess, false);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = 1 + sin(x)*sin(x) + y*y;
            TS_ASSERT_DELTA(answer_repl[i], u, 0.01);
        }

        PetscTools::Destroy(answer);

        // Analytical Jacobian
        answer=solver.Solve(initial_guess, true);
        ReplicatableVector answer_repl2(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = 1 + sin(x)*sin(x) + y*y;
            TS_ASSERT_DELTA(answer_repl2[i], u, 0.01);
        }

        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer);
    }
};

#endif //_TESTSIMPLENONLINEARELLIPTICSOLVER_HPP_
