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

/*
 *
 *
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 *
 *
 *
 */
#ifndef TESTWRITINGPDESOLVERSTUTORIAL_HPP_
#define TESTWRITINGPDESOLVERSTUTORIAL_HPP_

/* HOW_TO_TAG PDE
 * Write new PDE solvers (especially for linear coupled elliptic/parabolic systems)
 */


/*
 * == Introduction ==
 *
 * Chaste can be used to set up solvers for more general (coupled) PDEs. To do this the
 * user just needs to code up the integrands of any finite element (FE) matrices or vectors,
 * without having to deal with set-up, looping over elements, numerical quadrature, assembly
 * or solving the linear system. If you have a set of coupled linear PDEs for which it is
 * appropriate to use linear basis functions for each unknown (for example, a reaction-diffusion
 * system), then it is relatively straightforward to set up a solver that will be parallel and
 * reliable, since all the base components are heavily tested.
 *
 * Some solvers for general simple (uncoupled) linear PDEs are already provided in Chaste, such
 * as the `SimpleLinearEllipticSolver`. These are for PDEs that can be written in a generic
 * form (`SimpleLinearEllipticPde`, for example). However ''coupled'' PDEs can't
 * be easily written in generic form, so the user has to write their own solver. In this tutorial
 * we explain how to do this.
 *
 * For this tutorial the user needs to have read the solving-PDEs tutorials. It may also be
 * helpful to read the associated [wiki:ChasteGuides/NmoodLectureNotes lectures notes], in
 * particular the slides on solving equations using finite elements if you are not familiar
 * with this (lecture 2), the slides on the general design of the Chaste finite element solvers
 * (lecture 3), and the first part of lecture 4.
 *
 * Let us use the terminology "assembled in an FE manner" for any matrix or vector that is
 * defined via a volume/surface/line integral, and which is constructed by: looping over
 * elements (or surface elements, etc); computing the elemental contribution (i.e. a small
 * matrix/vector) using numerical quadrature; and adding to the full matrix/vector.
 *
 * We only consider linear problems here. In these problems the discretised FE problem leads
 * to a linear system, Ax=b, to be solved once in static problems and at each time step in
 * time-dependent problems. There are two cases to be distinguished. The first case is where
 * BOTH A and b are assembled in an FE manner, b possibly being composed of a volume integral
 * plus a surface integral. The other case is where this is not true, for example where b = Mc+d,
 * where the vector d and matrix M are assembled in an FE manner, but not the vector c.
 *
 * The Chaste PDE classes include ASSEMBLER classes, for setting up anything assembled in an
 * FE manner, and SOLVER classes, for setting up linear systems. In the general case, solvers
 * need to own assemblers for setting up each part of the linear system. However for the first
 * case described above (in which both A and b in Ax=b are assembled), we can use the design
 * where the solver IS AN assembler. We illustrate how to do this in the first tutorial.
 *
 * == Writing solvers ==
 *
 * Let us write a solver for the coupled 2-unknown problem
 * {{{
 *    Laplacian(u) + v = f(x,y)
 *    Laplacian(v) + u = g(x,y)
 * }}}
 * where `Laplacian(u)` of course represents u,,xx,,+u,,yy,, and where f and g are chosen so that,
 * with zero-dirichlet boundary conditions, the solution is given by u = sin(pi*x)sin(pi*x),
 * v = sin(2*pi*x)sin(2*pi*x).
 *
 *( As a brief aside, note that the solver we write will in fact work with general Dirichlet-Neumann
 * boundary conditions, though the test will only provide all-Dirichlet boundary conditions. We
 * save a discussion on general Dirichlet-Neumann boundary conditions for the second example.)
 *
 * Using linear basis functions, and a mesh with N nodes, the linear system that needs to be set up is
 * of size 2N by 2N, and in block form is:
 * {{{
 *   [ K   -M  ] [U]  =  [b1]
 *   [ -M   K  ] [V]     [b2]
 * }}}
 * where `K` is the stiffness matrix, `M` the mass matrix, `U` the vector of nodal values
 * of u, `V` the vector of nodal values of v, `b1` the vector with entries `integral(f\phi_i dV)` (i=1,..,N)
 * and `b2` has entries `integral(g\phi_i dV)` (here `phi_i` are the linear basis functions).
 *
 * This is the linear system which we now write a solver to set up. '''Note''', however, that
 * the main Chaste solvers assume a STRIPED data format, ie that the unknown vector
 * is `[U_1 V_1 U_2 V_2 .. U_n V_n]`, not `[ U_1 U_2 .. U_n V_1 V_2 .. V_n]`. ''We write down
 * equations in block form as it makes things clearer, but have to remember that the code
 * deals with STRIPED data structures.'' Striping is used in the code for parallelisation reasons.
 *
 * These are some basic includes as in the solving-PDEs tutorials
 */
#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "TrianglesMeshWriter.hpp"
/* We need to include the following two classes if we are going to use a combination of
 * (element_dim, space_dim, problem_dim) that isn't explicitly instantiated in
 * `BoundaryConditionsContainer.cpp` (see the bottom of that file) (without these includes this test will
 * fail to link).
 */
#include "BoundaryConditionsContainerImplementation.hpp"
#include "AbstractBoundaryConditionsContainerImplementation.hpp"
/* These two classes will be used in writing the solver */
#include "AbstractAssemblerSolverHybrid.hpp"
#include "AbstractStaticLinearPdeSolver.hpp"
/* We will solve a second problem, below, which will be time-dependent and will
 * require the following class */
#include "AbstractDynamicLinearPdeSolver.hpp"

/* The linear system is Ax=b where A and b are FE assembled, so we can use the solver-is-an-assembler
 * design. To construct our solver, we inherit from `AbstractAssemblerSolverHybrid` which links
 * the solver clases to the assembler classes, and `AbstractStaticLinearPdeSolver`. (For time-dependent
 * problems, the second parent would be `AbstractDynamicLinearPdeSolver`). Note the template
 * parameter `PROBLEM_DIM` below, in this case it is 2 as there are two unknowns.
 */
class MyTwoVariablePdeSolver
    : public AbstractAssemblerSolverHybrid<2/*elem_dim*/,2/*space_dim*/,2/*problem_dim*/,NORMAL/*amount of interpolation*/>,
      public AbstractStaticLinearPdeSolver<2,2,2>
{
private:
    /* The function f */
    double f(double x,double y)
    {
        return -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y) + sin(2*M_PI*x)*sin(2*M_PI*y);
    }
    /* The function g */
    double g(double x,double y)
    {
        return -8*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y) + sin(M_PI*x)*sin(M_PI*y);
    }

    /*
     *  The abstract assembler parent classes know how to assemble matrices and vectors, but the concrete
     *  class needs to provide the integrand of the elemental contribution to A and b. This first
     *  method returns the elemental contribution of the matrix A, given the provided bases
     *  (`rPhi`, `rGradPhi`). The '3's here represent the number of bases per element (ie the number
     *  of nodes as linear bases are being used). The returned matrix is 6 by 6 (problem_dim *
     *  num_bases_per_element = 2*3 = 6).
     */
    c_matrix<double,2*3,2*3> ComputeMatrixTerm(c_vector<double,3>& rPhi /* the three bases for the current element, evaluated at the current quad pt*/,
                                               c_matrix<double,2,3>& rGradPhi /* gradients of the three bases */,
                                               ChastePoint<2>& rX           /* physical coordinate of quad point */,
                                               c_vector<double,2>& rU       /* current solution (unused here as a linear static problem */,
                                               c_matrix<double,2,2>& rGradU /* current solution gradient (unused here as a linear static problem */,
                                               Element<2,2>* pElement)
    {
        /*
         * Set up the matrix, which corresponds to the elemental contribution for the matrix
         * written above, taking into account the striped nature of the matrices and vectors.
         * (Note: the following can be done more efficiently using matrix slices and products,
         * see `BidomainAssembler` for example).
         */
        c_matrix<double,2*3,2*3> ret = zero_matrix<double>(2*3, 2*3);

        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                for (unsigned k=0; k<2; k++)
                {
                    // stiffness matrix on diagonal 'blocks'
                    ret(2*i,  2*j)   += rGradPhi(k,i)*rGradPhi(k,j);
                    ret(2*i+1,2*j+1) += rGradPhi(k,i)*rGradPhi(k,j);
                }

                // (negative) mass matrix on off-diagonal 'blocks'
                ret(2*i+1, 2*j)   = -rPhi(i)*rPhi(j);
                ret(2*i,   2*j+1) = -rPhi(i)*rPhi(j);
            }
        }
        return ret;
    }

    /* Similarly compute the elemental contribution to the RHS vector */
    c_vector<double,2*3> ComputeVectorTerm(c_vector<double, 3>& rPhi,
                                           c_matrix<double, 2, 2+1>& rGradPhi,
                                           ChastePoint<2>& rX,
                                           c_vector<double,2>& rU,
                                           c_matrix<double,2,2>& rGradU,
                                           Element<2,2>* pElement)
    {
        c_vector<double,2*3> ret;

        for (unsigned i=0; i<3; i++)
        {
            ret(2*i)   = -f(rX[0],rX[1]) * rPhi(i);
            ret(2*i+1) = -g(rX[0],rX[1]) * rPhi(i);
        }
        return ret;
    }
    /* These classes which inherit from both assemblers and solvers must
     * provide the following method, which links the two. Just copy and paste
     * the following.
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {
        SetupGivenLinearSystem(currentSolution, computeMatrix, this->mpLinearSystem);
    }

public:
    /* The constructor takes in a mesh and boundary conditions container, and passes
     * them to the parent classes.
     */
    MyTwoVariablePdeSolver(TetrahedralMesh<2,2>* pMesh,
                           BoundaryConditionsContainer<2,2,2>* pBoundaryConditions)
        : AbstractAssemblerSolverHybrid<2,2,2,NORMAL>(pMesh,pBoundaryConditions),
          AbstractStaticLinearPdeSolver<2,2,2>(pMesh)
    {
    }
};
/*
 *  That is the solver written. The usage is the same as see the PDE solvers described in the
 *  previous tutorials - have a look at the first test below.
 *
 *  == A solver of 3 parabolic equations ==
 *
 *  Let us also write a solver for the following problem, which is composed of 3 parabolic PDEs
 *  {{{
 *    u_t = Laplacian(u) + v
 *    v_t = Laplacian(v) + u + 2w
 *    w_t = Laplacian(w) + g(t,x,y)
 *  }}}
 *  where g(t,x,y) = t if x>0.5 and 0 otherwise. This time we assume general
 *  Dirichlet-Neumann boundary conditions will be specified.
 *
 *  The `AbstractAssemblerSolverHybrid` deals with the Dirichlet and Neumann boundary parts of the implementation,
 *  so, we, the writer of the solver, don't have to worry about this. The user has to realise though that they are
 *  specifying NATURAL Neumann BCs, which are whatever appears naturally in the weak form of the problem. In this case, natural
 *  Neumann BCs are specification of: `du/dn = s1, dv/dn = s2, dw/dn = s3`, which coincide with usual Neumann BCs. However,
 *  suppose the last equation was `w_t = Laplacian(w) + Div(D grad(u))`, then the natural BCs would be:
 *  `du/dn = s1, dv/dn = s2, dw/dn + (Dgradu).n = s3`.
 *
 *  We need to choose a time-discretisation. Let us choose an implicit discretisation, ie
 *  {{{
 *  (u^{n+1} - u^{n})/dt = Laplacian(u^{n+1}) + v^{n+1}
 *  (v^{n+1} - v^{n})/dt = Laplacian(v^{n+1}) + u^{n+1} + 2w^{n+1}
 *  (w^{n+1} - w^{n})/dt = Laplacian(w^{n+1}) + g(t^{n+1},x)
 *  }}}
 *
 *  Using linear basis functions, and a mesh with N nodes, the linear system that needs to be set up is
 *  of size 3N by 3N, and in block form is:
 *  {{{
 *    [ M/dt+K     -M       0    ] [U^{n+1}]  =  [b1]  +  [c1]
 *    [   -M     M/dt+K    -2M   ] [V^{n+1}]     [b2]  +  [c2]
 *    [    0        0     M/dt+K ] [W^{n+1}]     [b3]  +  [c3]
 *  }}}
 * where `K` is the stiffness matrix, `M` the mass matrix, `U^n` the vector of nodal values
 * of u at time t_n, etc, `b1` has entries `integral( (u^n/dt)phi_i dV )`, and similarly for
 * `b2` and `b3`. Writing the Neumann boundary conditions for
 *  u as `du/dn = s1(x,y)` on Gamma, a subset of the boundary, then `c1` has entries
 * `integral_over_Gamma (s1*phi_i dS)`, and similarly for `c2` and `c3`.
 *
 * Let us create a solver for this linear system, which will be written in a way in which the RHS
 * vector is assembled in an FE manner, so that the solver-is-an-assembler design can be used.
 * Note that this solver inherits from `AbstractDynamicLinearPdeSolver` and PROBLEM_DIM is now equal
 * to 3. We don't have to worry about setting up [c1 c2 c3] (we just need to take in a `BoundaryConditionsContainer`
 * and the parent will use it in assembling this vector). We do however have to tell it
 * how to assemble the volume integral part of the RHS vector, and the LHS matrix.
 */
class ThreeParabolicPdesSolver
    : public AbstractAssemblerSolverHybrid<2,2,3,NORMAL>,
      public AbstractDynamicLinearPdeSolver<2,2,3>
{
private:
    /* Define the function g(t,x,y) */
    double g(double t, ChastePoint<2>& rX)
    {
        return t*(rX[0]>0.5);
    }

    /* Provide the (elemental contribution to the) LHS matrix. The matrix is 9 by 9, where
     * 9 = 3*3 = PROBLEM_DIM * NUM_NODES_PER_ELEMENT */
    c_matrix<double,3*3,3*3> ComputeMatrixTerm(c_vector<double,3>& rPhi,
                                               c_matrix<double,2,3>& rGradPhi,
                                               ChastePoint<2>& rX,
                                               c_vector<double,3>& rU,
                                               c_matrix<double,3,2>& rGradU,
                                               Element<2,2>* pElement)
    {
        c_matrix<double,9,9> ret = zero_matrix<double>(9,9);

        // this is how to get the current timestep
        double dt = PdeSimulationTime::GetPdeTimeStep();

        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                // mass matrix on the diagonal blocks
                ret(3*i,  3*j)   =  rPhi(i)*rPhi(j)/dt;
                ret(3*i+1,3*j+1) =  rPhi(i)*rPhi(j)/dt;
                ret(3*i+2,3*j+2) =  rPhi(i)*rPhi(j)/dt;

                // mass matrix on some off-diagonal blocks
                ret(3*i,  3*j+1) =  -rPhi(i)*rPhi(j);
                ret(3*i+1,3*j) =  -rPhi(i)*rPhi(j);
                ret(3*i+1,3*j+2) =  -2*rPhi(i)*rPhi(j);

                // stiffness matrix on the diagonal blocks
                for (unsigned dim=0; dim<2; dim++)
                {
                    ret(3*i,  3*j)   += rGradPhi(dim,i)*rGradPhi(dim,j);
                    ret(3*i+1,3*j+1) += rGradPhi(dim,i)*rGradPhi(dim,j);
                    ret(3*i+2,3*j+2) += rGradPhi(dim,i)*rGradPhi(dim,j);
                }
            }
        }
        return ret;
    }

    /* Provide the volume elemental contribution to the RHS vector, ie the vector `[b1 b2 b3]` above */
    c_vector<double,3*3> ComputeVectorTerm(c_vector<double, 3>& rPhi,
                                           c_matrix<double, 2, 3>& rGradPhi,
                                           ChastePoint<2>& rX,
                                           c_vector<double,3>& rU,
                                           c_matrix<double,3,2>& rGradU,
                                           Element<2,2>* pElement)
    {
        c_vector<double,3*3> ret;

        // get u,v,w out of the provided parameters
        double u = rU(0);
        double v = rU(1);
        double w = rU(2);

        // this is how to get the current time, timestep and inverse of the timestep
        double t = PdeSimulationTime::GetTime();
        double dt = PdeSimulationTime::GetPdeTimeStep();
        double inverse_dt = PdeSimulationTime::GetPdeTimeStepInverse();

        for (unsigned i=0; i<3; i++)
        {
            ret(3*i)   =  u* inverse_dt * rPhi(i);
            ret(3*i+1) =  v* inverse_dt * rPhi(i);
            ret(3*i+2) =  (w * inverse_dt + g(t+dt,rX)) * rPhi(i);
        }
        return ret;
    }

    /* Define this method as before */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {
        SetupGivenLinearSystem(currentSolution, computeMatrix, this->mpLinearSystem);
    }

public:
    /* The constructor is similar to before. However: '''important''' - by default the dynamic solvers
     * will reassemble the matrix each timestep. In this (and most other) problems the matrix is constant
     * and only needs to be assembled once. Make sure we tell the solver this, otherwise performance
     * will be destroyed.
     */
    ThreeParabolicPdesSolver(TetrahedralMesh<2,2>* pMesh,
                             BoundaryConditionsContainer<2,2,3>* pBoundaryConditions)
        : AbstractAssemblerSolverHybrid<2,2,3,NORMAL>(pMesh,pBoundaryConditions),
          AbstractDynamicLinearPdeSolver<2,2,3>(pMesh)
    {
        this->mMatrixIsConstant = true;
    }
};

/* Now the tests using the two solvers */
class TestWritingPdeSolversTutorial : public CxxTest::TestSuite
{
public:
   /* Use the first solver to solve the static PDE. We apply zero Dirichlet boundary conditions
    * on the whole of the boundary for both variables.
    */
    void TestMyTwoVariablePdeSolver()
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.01 /*h*/, 1.0 /*width*/, 1.0 /*height*/);

        // Boundary conditions for 2-unknown problem
        BoundaryConditionsContainer<2,2,2> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,0); // zero dirichlet for u
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,1); // zero dirichlet for v

        // Use our purpose-made solver for this problem:
        MyTwoVariablePdeSolver solver(&mesh,&bcc);

        /* The `AbstractStaticLinearPdeSolver` class from which our solver
         * inherits, provides a `Solve` method.
         */
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        /* Compare against the exact solution. */
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];

            double u = result_repl[2*i];
            double v = result_repl[2*i+1];

            double u_exact = sin(M_PI*x)*sin(M_PI*y);
            double v_exact = sin(2*M_PI*x)*sin(2*M_PI*y);

            TS_ASSERT_DELTA(u, u_exact, 0.002);
            TS_ASSERT_DELTA(v, v_exact, 0.007);
        }
        PetscTools::Destroy(result);
    }

    /* Now run a test solving the parabolic-parabolic-parabolic PDE system. */
    void TestMyParaEllipticSetOfPdesSolver()
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.05 /*h*/, 1.0 /*width*/, 1.0 /*height*/);

        /* Set up the boundary conditions. v and w are zero on the entire boundary,
         * and du/dn=1 on the LHS and 0 otherwise.
         */
        BoundaryConditionsContainer<2,2,3> bcc;

        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,1 /*index of unknown, ie v*/);
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,2 /*index of unknown, ie w*/);

        ConstBoundaryCondition<2>* p_neumann_bc = new ConstBoundaryCondition<2>(1.0);
        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        while (iter < mesh.GetBoundaryElementIteratorEnd())
        {
            if (fabs((*iter)->CalculateCentroid()[0])<1e-6)
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_neumann_bc, 0 /*index of unknown, ie u*/);
            }
            iter++;
        }

        /* Use our solver */
        ThreeParabolicPdesSolver solver(&mesh,&bcc);

        /* The interface is exactly the same as the `SimpleLinearParabolicSolver`. */
        solver.SetTimeStep(0.01);
        solver.SetTimes(0.0, 2.0);

        Vec initial_condition = PetscTools::CreateAndSetVec(3*mesh.GetNumNodes(), 0.0);
        solver.SetInitialCondition(initial_condition);

        /* For this test we show how to output results to file for multiple sampling times. We start by
         * specifying an output directory and filename prefix for our results file:
         */
        solver.SetOutputDirectoryAndPrefix("ThreeVarCoupledProblem","results");

        /* When an output directory has been specified, the solver writes output in HDF5 format. To
         * convert this to another output format, we call the relevant method. Here, we convert
         * the output to plain txt files. We also say how often to write the data, telling the
         * solver to output results to file every tenth timestep. The solver will create
         * one file for each variable and for each time, so for example, the file
         * results_Variable_0_10 is the results for u, over all nodes, at the 11th printed time.
         * Have a look in the output directory after running the test. (Note: the HDF5 data can also be
         * converted to VTK or cmgui formats).
         */
        solver.SetOutputToTxt(true);
        solver.SetPrintingTimestepMultiple(10);

        /* We are now ready to solve the system. */
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        /* The plain txt output can be loaded into matlab for easy visualisation. For this we
         * also need the mesh data - at the very least the nodal locations - so we also write out
         * the mesh
         */
        TrianglesMeshWriter<2,2> mesh_writer("ThreeVarCoupledProblem", "mesh", false /*don't clean (ie delete everything in) directory!*/);
        mesh_writer.WriteFilesUsingMesh(mesh);

        /* Note that we need to destroy the initial condition vector as well as the solution. */
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
        /*
         *
         * '''Visualisation:''' To visualise in matlab/octave, you can load the node file,
         * and then the data files. However, the node file needs to be edited to remove any
         * comment lines (lines starting with '#') and the header line (which says how
         * many nodes there are). (The little matlab script 'anim/matlab/read_chaste_node_file.m'
         * can be used to do this, though it is not claimed to be very robust). As an
         * example matlab visualisation script, the following, if run from the output
         * directory, plots w:
         * {{{
         * pos = read_chaste_node_file('mesh.node');
         * for i=0:200
         *   file = ['txt_output/results_Variable_2_',num2str(i),'.txt'];
         *   w = load(file);
         *   plot3(pos(:,1),pos(:,2),w,'*');
         *   pause;
         * end;
         * }}}
         *
         */
    }
};

#endif // TESTWRITINGPDESOLVERSTUTORIAL_HPP_
