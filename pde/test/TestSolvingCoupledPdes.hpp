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

#ifndef TESTSOLVINGCOUPLEDPDES_HPP_
#define TESTSOLVINGCOUPLEDPDES_HPP_

#include <cxxtest/TestSuite.h>

#include <vector>
#include <cmath>

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SimpleLinearEllipticSolver.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "AbstractAssemblerSolverHybrid.hpp"
#include "AbstractStaticLinearPdeSolver.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "ReplicatableVector.hpp"

#include <boost/mpl/if.hpp>
#include <boost/mpl/void.hpp>

//////////////////////////////////////////////////////////////////////////////
// A simple pde : u_xx + u_yy + x = 0
//////////////////////////////////////////////////////////////////////////////
class MySimplePde : public AbstractLinearEllipticPde<2,2>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<2>& rX, Element<2,2>* )
    {
        return rX[0];
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>&, Element<2,2>* )
    {
        return 0.0;
    }

    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& x)
    {
        return identity_matrix<double>(2);
    }
};

//////////////////////////////////////////////////////////////////////////////
// A solver to solve the 'coupled' 2-unknown problem
//    u_xx + u_yy +         x  = 0
//    v_xx + v_yy + \lambda x  = 0
//
//   \lambda is taken in in the constructor
//////////////////////////////////////////////////////////////////////////////
class MySimpleCoupledSolver
    : public AbstractAssemblerSolverHybrid<2,2,2,NORMAL>,
      public AbstractStaticLinearPdeSolver<2,2,2>
{
private:
    double mLambda;

    virtual c_matrix<double,2*(2+1),2*(2+1)> ComputeMatrixTerm(c_vector<double, 2+1>& rPhi,
                                                               c_matrix<double, 2, 2+1>& rGradPhi,
                                                               ChastePoint<2>& rX,
                                                               c_vector<double,2>& rU,
                                                               c_matrix<double,2,2>& rGradU,
                                                               Element<2,2>* pElement)
    {
        c_matrix<double,2*(2+1),2*(2+1)> ret = zero_matrix<double>(2*(2+1), 2*(2+1));

        /*
         * The following can be done more efficiently using matrix slices
         * and prods and so on (see BidomainAssembler) - efficiency not
         * needed for this test though.
         */
        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                for (unsigned k=0; k<2; k++)
                {
                    ret(2*i,  2*j)   += rGradPhi(k,i)*rGradPhi(k,j);
                    ret(2*i+1,2*j+1) += rGradPhi(k,i)*rGradPhi(k,j);
                }
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

        for (unsigned i=0; i<3; i++)
        {
            ret(2*i)   =         rX[0]*rPhi(i);
            ret(2*i+1) = mLambda*rX[0]*rPhi(i);
        }
        return ret;
    }

    void InitialiseForSolve(Vec initialSolution)
    {
        AbstractLinearPdeSolver<2,2,2>::InitialiseForSolve(initialSolution);
        assert(this->mpLinearSystem);
        this->mpLinearSystem->SetMatrixIsSymmetric(true);
        this->mpLinearSystem->SetKspType("cg");
    }

    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {
        SetupGivenLinearSystem(currentSolution, computeMatrix, this->mpLinearSystem);
    }

public:
    MySimpleCoupledSolver(TetrahedralMesh<2,2>* pMesh,
                          BoundaryConditionsContainer<2,2,2>* pBoundaryConditions,
                          double lambda)
        : AbstractAssemblerSolverHybrid<2,2,2,NORMAL>(pMesh,pBoundaryConditions),
          AbstractStaticLinearPdeSolver<2,2,2>(pMesh)
    {
        this->mpBoundaryConditions = pBoundaryConditions;
        mLambda = lambda;
    }
};

/*
 * A solver to solve the coupled 2-unknown problem
 *    u_xx + u_yy + v = f(x,y)
 *    v_xx + v_yy + u = g(x,y)
 *
 *  where f and g are chosen so that (with zero-dirichlet boundary conditions)
 *  the solution is
 *       u = sin(pi*x)sin(pi*x),   v = sin(2*pi*x)sin(2*pi*x)
 *
 *  The linear system that needs to be set up is, in block form
 *
 *   [ K   -M  ] [U]  =  [b1]
 *   [ -M   K  ] [V]     [b2]
 *
 *   where K is the stiffness matrix, M the mass matrix, U the vector of nodal values
 *   of u, V the vector of nodal values of v, b1_i = integral(f\phi_i dV) and
 *   b1_i = integral(g\phi_i dV), where the basis functions are phi_i
 *
 *   However, these Chaste solvers assume a STRIPED data format, ie that the unknown vector
 *   is [U_1 V_1 U_2 V_2 .. U_n V_n]  not  [ U_1 U_2 .. U_n V_1 V_2 .. V_n]
 *   Hence the matrix and vector coded below are the striped analogues of the
 *   matrix and vector written above
 */
class AnotherCoupledSolver
    : public AbstractAssemblerSolverHybrid<2,2,2,NORMAL>,
      public AbstractStaticLinearPdeSolver<2,2,2>
{
private:
    double f(double x,double y)
    {
        return -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y) + sin(2*M_PI*x)*sin(2*M_PI*y);
    }

    double g(double x,double y)
    {
        return -8*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y) + sin(M_PI*x)*sin(M_PI*y);
    }

    virtual c_matrix<double,2*(2+1),2*(2+1)> ComputeMatrixTerm(c_vector<double, 2+1>& rPhi,
                                                               c_matrix<double, 2, 2+1>& rGradPhi,
                                                               ChastePoint<2>& rX,
                                                               c_vector<double,2>& rU,
                                                               c_matrix<double,2,2>& rGradU,
                                                               Element<2,2>* pElement)
    {
        c_matrix<double,2*(2+1),2*(2+1)> ret = zero_matrix<double>(2*(2+1), 2*(2+1));

        /*
         * The following can be done more efficiently using matrix slices
         * and prods and so on (see BidomainAssembler) - efficiency not
         * needed for this test though.
         */
        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                for (unsigned k=0; k<2; k++)
                {
                    ret(2*i,  2*j)   += rGradPhi(k,i)*rGradPhi(k,j);
                    ret(2*i+1,2*j+1) += rGradPhi(k,i)*rGradPhi(k,j);
                }

                ret(2*i+1, 2*j)   = -rPhi(i)*rPhi(j);
                ret(2*i,   2*j+1) = -rPhi(i)*rPhi(j);
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

        for (unsigned i=0; i<3; i++)
        {
            ret(2*i)   = - f(rX[0],rX[1]) * rPhi(i);
            ret(2*i+1) = - g(rX[0],rX[1]) * rPhi(i);
        }
        return ret;
    }

    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {
        SetupGivenLinearSystem(currentSolution, computeMatrix, this->mpLinearSystem);
    }

public:
    AnotherCoupledSolver(TetrahedralMesh<2,2>* pMesh,
                         BoundaryConditionsContainer<2,2,2>* pBoundaryConditions)
        : AbstractAssemblerSolverHybrid<2,2,2,NORMAL>(pMesh,pBoundaryConditions),
          AbstractStaticLinearPdeSolver<2,2,2>(pMesh)
    {
    }
};

//////////////////////////////////////////////////////////////////////////////
// Test class
//////////////////////////////////////////////////////////////////////////////
class TestSolvingCoupledPdes : public CxxTest::TestSuite
{
public:

    /*
     * Solve:
     *     u_xx + u_yy + x  = 0
     *     v_xx + v_yy + 2x = 0
     *  with zero dirichlet on boundary
     *
     *  This is obviously really just two virtually identical uncoupled
     *  problems
     */
    void TestSimpleCoupledPde()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ////////////////////////////////////////////////////////////////
        // Solve coupled system using solver defined above
        ////////////////////////////////////////////////////////////////

        // Boundary conditions for 2-unknown problem
        BoundaryConditionsContainer<2,2,2> bcc_2unknowns;
        bcc_2unknowns.DefineZeroDirichletOnMeshBoundary(&mesh,0); // zero dirichlet for u
        bcc_2unknowns.DefineZeroDirichletOnMeshBoundary(&mesh,1); // zero dirichlet for v

        // Lambda in MySimpleCoupledSolver = 2
        MySimpleCoupledSolver solver_2unknowns(&mesh,&bcc_2unknowns,2.0);
        Vec result_2unknowns = solver_2unknowns.Solve();
        ReplicatableVector result_2unknowns_repl(result_2unknowns);

        ///////////////////////////////////////////////////////////////////
        // Now solve u_xx + u_yy + x = 0 as an uncoupled 1-unknown problem
        ///////////////////////////////////////////////////////////////////

        // Instantiate PDE object
        MySimplePde pde;  //defined above

        // Boundary conditions for 1-unknown problem
        BoundaryConditionsContainer<2,2,1> bcc_1unknown;
        bcc_1unknown.DefineZeroDirichletOnMeshBoundary(&mesh);

        // Assembler
        SimpleLinearEllipticSolver<2,2> solver_1unknown(&mesh,&pde,&bcc_1unknown);

        Vec result_1unknown = solver_1unknown.Solve();
        ReplicatableVector result_1unknown_repl(result_1unknown);

        /*
         * Check the u solutions (result_2unknowns_repl[2*i]) is equal
         * to the solution of the 1-unknown problem and the v solutions
         * (result_2unknowns_repl[2*i+1]) are equal to two times the
         * 1-unknown solution.
         */
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(result_2unknowns_repl[2*i]  ,   result_1unknown_repl[i], 1e-10);
            TS_ASSERT_DELTA(result_2unknowns_repl[2*i+1], 2*result_1unknown_repl[i], 1e-10);
            //std::cout << result_1unknown_repl[i] << " ";
        }

        PetscTools::Destroy(result_2unknowns);
        PetscTools::Destroy(result_1unknown);
    }

    /*  Solve:
     *     u_xx + u_yy + x = 0
     *     v_xx + v_yy + x = 0
     *  with neumann boundary conditions (the same on both u and v)
     *  on part of the boundary
     *
     *  This is obviously two identical uncoupled problems
     */
    void TestSimpleCoupledPdeWithNeumannBoundaryConditions()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ////////////////////////////////////////////////////////////////
        // Solve coupled system using solver defined above
        ////////////////////////////////////////////////////////////////

        // Boundary conditions for 2-unknown problem
        BoundaryConditionsContainer<2,2,2> bcc_2unknowns;

        // du/dn = -0.5 on r=1
        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* p_boundary_condition = new ConstBoundaryCondition<2>(-0.5);
        ConstBoundaryCondition<2>* p_boundary_condition1 = new ConstBoundaryCondition<2>(-0.5);
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc_2unknowns.AddNeumannBoundaryCondition(*iter, p_boundary_condition,0);
            bcc_2unknowns.AddNeumannBoundaryCondition(*iter, p_boundary_condition1,1);
            iter++;
        }
        // u = 2 at some point on the boundary, say node 1
        p_boundary_condition = new ConstBoundaryCondition<2>(2.0);
        p_boundary_condition1 = new ConstBoundaryCondition<2>(2.0);
        bcc_2unknowns.AddDirichletBoundaryCondition(mesh.GetNode(1), p_boundary_condition,0);
        bcc_2unknowns.AddDirichletBoundaryCondition(mesh.GetNode(1), p_boundary_condition1,1);

        // Use solver to solve (with lambda in MySimpleCoupledSolver = 1)

        MySimpleCoupledSolver solver_2unknowns(&mesh,&bcc_2unknowns,1.0);

        Vec result_2unknowns = solver_2unknowns.Solve();
        ReplicatableVector result_2unknowns_repl(result_2unknowns);

        ///////////////////////////////////////////////////////////////////
        // Now solve u_xx + u_yy + x = 0 as an uncoupled 1-unknown problem
        ///////////////////////////////////////////////////////////////////

        // Instantiate PDE object
        MySimplePde pde; // defined above

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
        SimpleLinearEllipticSolver<2,2> solver_1unknown(&mesh,&pde,&bcc_1unknown);

        Vec result_1unknown = solver_1unknown.Solve();
        ReplicatableVector result_1unknown_repl(result_1unknown);

        /*
         * Check the u solutions (result_2unknowns_repl[2*i]) is equal
         * to the solution of the 1-unknown problem and the v solutions
         * (result_2unknowns_repl[2*i+1]) are equal to the 1-unknown
         * solution.
         */
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(result_2unknowns_repl[2*i]  , result_1unknown_repl[i], 1e-6);
            TS_ASSERT_DELTA(result_2unknowns_repl[2*i+1], result_1unknown_repl[i], 1e-6);
            //std::cout << result_1unknown_repl[i] << " ";
        }

        PetscTools::Destroy(result_2unknowns);
        PetscTools::Destroy(result_1unknown);
    }

    /*
     * Solve a real coupled problem:
     *    u_xx + u_yy  + v = f(x,y)
     *    v_xx + v_yy  + u = g(x,y)
     *
     * where f and g are chosen so that (with zero-dirichlet boundary conditions)
     * the solution is
     *
     *    u = sin(pi*x)sin(pi*x),   v = sin(2*pi*x)sin(2*pi*x)
     */
    void TestRealCoupledPde()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Boundary conditions for 2-unknown problem
        BoundaryConditionsContainer<2,2,2> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,0); // zero dirichlet for u
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,1); // zero dirichlet for v

        // Purpose-made solver for this problem:
        AnotherCoupledSolver solver(&mesh,&bcc);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        //OutputFileHandler handler("Something");
        //out_stream p_file = handler.OpenOutputFile("file.txt");

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];

            double u = sin(M_PI*x)*sin(M_PI*y);
            double v = sin(2*M_PI*x)*sin(2*M_PI*y);

            //*p_file << x << " " << y << " " << result_repl[2*i] << " "
            //        <<  result_repl[2*i+1] << " " << u << " " << v << std::endl;

            TS_ASSERT_DELTA( result_repl[2*i]  , u, 0.002);
            TS_ASSERT_DELTA( result_repl[2*i+1], v, 0.007);
        }

        //p_file->close();
        PetscTools::Destroy(result);
    }
};

#endif /*TESTSOLVINGCOUPLEDPDES_HPP_*/
