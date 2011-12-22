/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTSOLVINGCOUPLEDNONLINEARPDES_HPP_
#define TESTSOLVINGCOUPLEDNONLINEARPDES_HPP_

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include <petsc.h>
#include <vector>
#include <cmath>
#include "PetscSetupAndFinalize.hpp"
#include "SimpleNonlinearEllipticSolver.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "TrianglesMeshReader.hpp"
#include "ReplicatableVector.hpp"
#include "NonlinearEquationPde.hpp"
#include "SimpleNewtonNonlinearSolver.hpp"

/**
 * A solver to solve the 'coupled' 2-unknown problem
 *    div.(u grad u) + 1  = 0,
 *    div.(v grav v) + lambda  = 0,
 * where lambda is taken in in the constructor.
 */
template <int DIM>
class MySimpleNonlinearCoupledSolver : public AbstractNonlinearAssemblerSolverHybrid<DIM,DIM,2>
{
private:
    double mLambda;
    virtual c_matrix<double,2*(DIM+1),2*(DIM+1) > ComputeMatrixTerm(c_vector<double, DIM+1>& rPhi,
            c_matrix<double, DIM, DIM+1>& rGradPhi,
            ChastePoint<DIM>& rX,
            c_vector<double,2>& rU,
            c_matrix<double,2,DIM>& rGradU,
            Element<DIM,DIM>* pElement)
    {
        c_matrix<double,2*(DIM+1),2*(DIM+1)> ret = zero_matrix<double>(2*(DIM+1), 2*(DIM+1));

        for (unsigned i=0; i<DIM+1; i++)
        {
            for (unsigned j=0; j<DIM+1; j++)
            {
                matrix_column<c_matrix<double,DIM,DIM+1> > grad_phi_i(rGradPhi,i);
                matrix_column<c_matrix<double,DIM,DIM+1> > grad_phi_j(rGradPhi,j);
                matrix_row<c_matrix<double,2,DIM> > gradU0(rGradU, 0);
                matrix_row<c_matrix<double,2,DIM> > gradU1(rGradU, 1);

                ret(2*i,  2*j)   = + rPhi(j)*inner_prod(gradU0,grad_phi_i) + rU(0)*inner_prod(grad_phi_j,grad_phi_i);
                ret(2*i+1,2*j+1) = + rPhi(j)*inner_prod(gradU1,grad_phi_i) + rU(1)*inner_prod(grad_phi_j,grad_phi_i);
            }
        }
        return ret;
    }


    virtual c_vector<double,2*(DIM+1)> ComputeVectorTerm(c_vector<double, DIM+1>& rPhi,
                                                         c_matrix<double, DIM, DIM+1>& rGradPhi,
                                                         ChastePoint<DIM>& rX,
                                                         c_vector<double,2>& rU,
                                                         c_matrix<double,2,DIM>& rGradU,
                                                         Element<DIM,DIM>* pElement)
    {
        c_vector<double,2*(DIM+1)> ret;

        for (unsigned i=0; i<DIM+1; i++)
        {
            matrix_column<c_matrix<double,DIM,DIM+1> > grad_phi_i(rGradPhi, i);
            matrix_row<c_matrix<double,2,DIM> > gradU0(rGradU, 0);
            matrix_row<c_matrix<double,2,DIM> > gradU1(rGradU, 1);

            ret(2*i)   = rU(0)*inner_prod(gradU0,grad_phi_i) - rPhi(i);
            ret(2*i+1) = rU(1)*inner_prod(gradU1,grad_phi_i) - mLambda*rPhi(i);
        }
        return ret;
    }

public:
    MySimpleNonlinearCoupledSolver(TetrahedralMesh<DIM,DIM>* pMesh,
                                      BoundaryConditionsContainer<DIM,DIM,2>* pBoundaryConditions,
                                      double lambda)
        : AbstractNonlinearAssemblerSolverHybrid<DIM,DIM,2>(pMesh,pBoundaryConditions)
    {
        mLambda = lambda;
    }
};

/**
 * A solver to solve the coupled 2-unknown problem
 *    div.(v gradu) = f(x,y),
 *    div.(u gradv) = g(x,y),
 * where f and g (and boundary conditions) are chosen such that the solution is
 *    u = x^2,  v = y.
 */
class AnotherCoupledNonlinearAssembler :  public AbstractNonlinearAssemblerSolverHybrid<2,2,2> // AnotherCoupledNonlinearAssembler>
{
private:
    double f(double x,double y)
    {
        return 2*y;
    }

    double g(double x,double y)
    {
        return 0;
    }

    virtual c_matrix<double,2*(2+1),2*(2+1)> ComputeMatrixTerm(c_vector<double, 2+1>& rPhi,
                                                               c_matrix<double, 2, 2+1>& rGradPhi,
                                                               ChastePoint<2>& rX,
                                                               c_vector<double,2>& rU,
                                                               c_matrix<double,2,2>& rGradU,
                                                               Element<2,2>* pElement)
    {
        c_matrix<double,2*(2+1),2*(2+1)> ret;

        for (unsigned i=0; i<2+1; i++)
        {
            for (unsigned j=0; j<2+1; j++)
            {
                matrix_column<c_matrix<double,2,2+1> > grad_phi_i(rGradPhi, i);
                matrix_column<c_matrix<double,2,2+1> > grad_phi_j(rGradPhi, j);
                matrix_row<c_matrix<double,2,2> > gradU0(rGradU, 0);
                matrix_row<c_matrix<double,2,2> > gradU1(rGradU, 1);

                ret(2*i,  2*j)   = rU(1)*inner_prod(grad_phi_j, grad_phi_i);
                ret(2*i,  2*j+1) = rPhi(j)*inner_prod(gradU0, grad_phi_i);
                ret(2*i+1,2*j)   = rPhi(j)*inner_prod(gradU1, grad_phi_i);
                ret(2*i+1,2*j+1) = rU(0)*inner_prod(grad_phi_j, grad_phi_i);
            }
        }
        return ret;
    }

    virtual c_vector<double,2*(2+1)> ComputeVectorTerm(c_vector<double, 2+1>& rPhi,
                                                       c_matrix<double, 2, 2+1>& rGradPhi,
                                                       ChastePoint<2>& rX,
                                                       c_vector<double,2>& rU,
                                                       c_matrix<double,2,2>& rGradU,
                                                       Element<2,2>* pElement)
    {
        c_vector<double,2*(2+1)> ret;

        for (unsigned i=0; i<2+1; i++)
        {
            matrix_column<c_matrix<double,2,2+1> > grad_phi_i(rGradPhi,i);
            matrix_row<c_matrix<double,2,2> > gradU0(rGradU, 0);
            matrix_row<c_matrix<double,2,2> > gradU1(rGradU, 1);

            ret(2*i)   = rU(1)*inner_prod(gradU0,grad_phi_i) + f(rX[0], rX[1])*rPhi(i);
            ret(2*i+1) = rU(0)*inner_prod(gradU1,grad_phi_i) + g(rX[0], rX[1])*rPhi(i);
        }
        return ret;
    }

public:
    AnotherCoupledNonlinearAssembler(TetrahedralMesh<2,2>* pMesh,
                                     BoundaryConditionsContainer<2,2,2>* pBoundaryConditions)
        : AbstractNonlinearAssemblerSolverHybrid<2,2,2>(pMesh,pBoundaryConditions)
    {
    }
};

//////////////////////////////////////////////////////////////////////////////
// Test class
//////////////////////////////////////////////////////////////////////////////
class TestSolvingCoupledNonlinearPdes : public CxxTest::TestSuite
{
private:

    template<int DIM>
    void runTestSimpleCoupledNonlinearPde()
    {
        std::string file;
        if (DIM==1)
        {
            file = "mesh/test/data/1D_0_to_1_10_elements";
        }
        else if (DIM==2)
        {
            file = "mesh/test/data/disk_522_elements";
        }
        else
        {
            file = "mesh/test/data/slab_395_elements";
        }
        TrianglesMeshReader<DIM,DIM> mesh_reader(file);


        TetrahedralMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ////////////////////////////////////////////////////////////////
        // Solve coupled system using solver defined above
        ////////////////////////////////////////////////////////////////

        // Boundary conditions for 2-unknown problem
        BoundaryConditionsContainer<DIM,DIM,2> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,0); // zero dirichlet for u
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,1); // zero dirichlet for v

        // For comparing residuals
        MySimpleNonlinearCoupledSolver<DIM> solver_lam_1(&mesh, &bcc, 1);

        // For comparing solutions
        MySimpleNonlinearCoupledSolver<DIM> solver(&mesh, &bcc, 4);

        ////////////////////////////////////////////
        // Store residual
        ////////////////////////////////////////////

        // Initialize 'solution' vector
        Vec guess = PetscTools::CreateAndSetVec(2*mesh.GetNumNodes(), 1.0);

        // Solve as well
        Vec result = solver.Solve(guess, true);
        ReplicatableVector result_repl(result);

        ///////////////////////////////////////////////////////////////////////
        // Now solve div.(u gradu) + 1 = 0 as an uncoupled 1-unknown problem
        ///////////////////////////////////////////////////////////////////////

        // Instantiate PDE object
        NonlinearEquationPde<DIM> pde;  //defined above

        // boundary conditions for 1-unknown problem
        BoundaryConditionsContainer<DIM,DIM,1> bcc_1unknown;
        bcc_1unknown.DefineZeroDirichletOnMeshBoundary(&mesh);

        // Assembler
        SimpleNonlinearEllipticSolver<DIM,DIM> solver_1unknown(&mesh,&pde,&bcc_1unknown);

        ////////////////////////////////////////////
        // Store residual
        ////////////////////////////////////////////
        Vec guess_1unknown = PetscTools::CreateAndSetVec(mesh.GetNumNodes(), 1.0);

        // Solve as well
        Vec result_1unknown = solver_1unknown.Solve(guess_1unknown, true);
        ReplicatableVector result_1unknown_repl(result_1unknown);

        /////////////////////////////////////////////////////////////////////////
        // check the residuals and solutions agree
        // (check the u solutions (result_repl[2*i]) is equal to the
        // solution of the 1-unknown problem and the v solutions
        // (result_repl[2*i+1]) equal to 2 times the 1-unknown
        // solution (as lambda=4))
        //////////////////////////////////////////////////////////////////////////
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(result_repl[2*i],   result_1unknown_repl[i], 1e-4);
            TS_ASSERT_DELTA(result_repl[2*i+1], 2*result_1unknown_repl[i], 1e-4);
        }
    }




public:

    /*
     * Solve:
     *     div.(u gradu) + 1  = 0,
     *     div.(v gradv) + 4  = 0,
     * with zero dirichlet on boundary.
     *
     * This is obviously really just two virtually identical uncoupled
     * problems
     */
    void TestSimpleCoupledNonlinearPde()   throw (Exception)
    {
        // run 1d version
        runTestSimpleCoupledNonlinearPde<1>();

        // run 2d version
        runTestSimpleCoupledNonlinearPde<2>();

        // run 3d version
        runTestSimpleCoupledNonlinearPde<3>();
    }

    /*
     * Solve:
     *     div.(u gradu) + 1  = 0,
     *     div.(v gradv) + 1  = 0,
     * with neumann boundary conditions (the same on both u and v)
     * on part of the boundary
     *
     * This is obviously two identical uncoupled problems
     */
    void TestSimpleCoupledNonlinearPdeWithNeumannBoundaryConditions() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ////////////////////////////////////////////////////////////////
        // Solve coupled system using solver defined above
        ////////////////////////////////////////////////////////////////

        // boundary conditions for 2-unknown problem
        BoundaryConditionsContainer<2,2,2> bcc;

        // du/dn = -0.5 on r=1
        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* p_boundary_condition = new ConstBoundaryCondition<2>(-0.5);
        ConstBoundaryCondition<2>* p_boundary_condition1 = new ConstBoundaryCondition<2>(-0.5);
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition,0);
            bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition1,1);
            iter++;
        }

        // u = 2 at some point on the boundary, say node 1
        p_boundary_condition = new ConstBoundaryCondition<2>(2.0);
        p_boundary_condition1 = new ConstBoundaryCondition<2>(2.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(1), p_boundary_condition,0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(1), p_boundary_condition1,1);

        // Use solver to solve (with lambda = 1)
        MySimpleNonlinearCoupledSolver<2> solver(&mesh,&bcc,1.0);

        Vec result = solver.Solve(PetscTools::CreateAndSetVec(2*mesh.GetNumNodes(),1.0),true);
        ReplicatableVector result_repl(result);

        ///////////////////////////////////////////////////////////////////////
        // Now solve div.(u gradu) + 1 = 0 as an uncoupled 1-unknown problem
        ///////////////////////////////////////////////////////////////////////

        // Instantiate PDE object
        NonlinearEquationPde<2> pde;

        // boundary conditions for 1-unknown problem
        BoundaryConditionsContainer<2,2,1> bcc_1unknown;

        iter = mesh.GetBoundaryElementIteratorBegin();
        p_boundary_condition = new ConstBoundaryCondition<2>(-0.5);
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc_1unknown.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
            iter++;
        }
        // u = 2 at some point on the boundary, say node 1
        p_boundary_condition = new ConstBoundaryCondition<2>(2.0);
        bcc_1unknown.AddDirichletBoundaryCondition(mesh.GetNode(1), p_boundary_condition);

        // Assembler
        SimpleNonlinearEllipticSolver<2,2> solver_1unknown(&mesh,&pde,&bcc_1unknown);

        Vec result_1unknown = solver_1unknown.Solve(PetscTools::CreateAndSetVec(mesh.GetNumNodes(),1.0), true);
        ReplicatableVector result_1unknown_repl(result_1unknown);

        // check the u solutions (result_repl[2*i]) is equal to the
        // solution of the 1-unknown problem and the v solutions
        // (result_repl[2*i+1]) are equal to the 1-unknown
        // solution
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(result_repl[2*i]  , result_1unknown_repl[i], 1e-6);
            TS_ASSERT_DELTA(result_repl[2*i+1], result_1unknown_repl[i], 1e-6);
        }
    }



    /* Solve the real coupled pde
     *
     *    v (u_xx+u_yy) = f
     *    u (v_xx+v_yy) = g
     *
     * where f,g and boundary conditions are chosen given the solution
     *
     *    u = x^2
     *    v = y
     */
    void TestRealCoupledPde() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");

// REPLACE ABOVE LINE WITH THIS AND THE ASSEMBLER SEEMS TO GET STUCK SOMEWHERE (applying
// boundary conditions maybe?)  ///\todo: find out why..
        //TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");

        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // boundary conditions for 2-unknown problem
        BoundaryConditionsContainer<2,2,2> bcc;

        TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        while (iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];

            // apply bc u=x^2
            ConstBoundaryCondition<2>* p_boundary_condition = new ConstBoundaryCondition<2>(x*x);
            bcc.AddDirichletBoundaryCondition(*iter, p_boundary_condition, 0);

            // apply bc v=x^2
            ConstBoundaryCondition<2>* p_boundary_condition1 = new ConstBoundaryCondition<2>(y);
            bcc.AddDirichletBoundaryCondition(*iter, p_boundary_condition1, 1);
            iter++;
        }


        // purpose-made solver for this problem:
        AnotherCoupledNonlinearAssembler solver(&mesh,&bcc);

        // use the newton solver (for coverage)
        SimpleNewtonNonlinearSolver newton_solver;
        newton_solver.SetTolerance(1e-10);
        newton_solver.SetWriteStats();
        solver.SetNonlinearSolver(&newton_solver);

        //// uncomment this to check whether ComputeMatrixTerm has been coded up correctly
        //// (by seeing whether the resulting analytic Jacobian matches the numerical one).
        //assert( solver.VerifyJacobian() );

        // IMPORTANT NOTE: both the petsc nonlinear solver and the Newton solver will FAIL
        // if an initial guess of zeroes is given.
        Vec result = solver.Solve( PetscTools::CreateAndSetVec(2*mesh.GetNumNodes(),1.0), false);

        int size;
        VecGetSize(result,&size);

        TS_ASSERT_EQUALS(size, (int)mesh.GetNumNodes()*2);

        ReplicatableVector result_repl(result);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];

            double u = result_repl[2*i];
            double v = result_repl[2*i+1];

            TS_ASSERT_DELTA(u, x*x, 1e-2);
            TS_ASSERT_DELTA(v, y, 1e-2);
        }

        VecDestroy(result);
    }
};
#endif /*TESTSOLVINGCOUPLEDNONLINEARPDES_HPP_*/

