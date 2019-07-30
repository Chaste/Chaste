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
#ifndef TESTWRITINGPDESOLVERSTWOTUTORIAL_HPP_
#define TESTWRITINGPDESOLVERSTWOTUTORIAL_HPP_


/*
 *  == Introduction ==
 *
 *  In the previous tutorial we showed how a PDE solver could be written for the
 *  'simple' case in which the FEM discretisation leads to a linear system Ax=b where
 *  both A and b are 'assembled'. In this tutorial, we consider the more general case,
 *  and show to write assembler classes which assemble one particular matrix or vector,
 *  and how to write solver classes which ''use'' assemblers to create and solve the FEM
 *  linear system.
 *
 *  We will take as the test problem the heat equation, `u_t = u_{xx}`, with Dirichlet
 *  BCs `u = u*` on `Gamma1` and `du/dn = g` on `Gamma2`.
 *
 *  We write a solver which uses an '''explicit''' time-discretisation (as opposed to the implicit
 *  discretisations used throughout the rest of the code). The FEM linear system that needs to be set up is
 *  {{{
 *  M U^{n+1} = (M + dt K) U^{n}  +  c
 *  }}}
 *  where `M` is the mass matrix, `K` the stiffness matrix, and `U^{n}` the vector of nodal
 *  values of u at timestep n. c is the surface integral term coming from the Neumann BCs,
 *  ie `c_i = integral_over_Gamma2 (g * phi_i dS)`. (This can be compared with an
 *  implicit time-discretisation, for which we solve `(M - dt K) U^{n+1} = M U^{n} + c`).
 *
 *  Let us call `M + dt*K` the 'RHS matrix'. We will write a solver, inheriting from
 *  `AbstractDynamicLinearPdeSolver`, which is going to ''use'' three assemblers: (i) an assembler of
 *  the mass matrix (already written); (ii) an assembler of the RHS matrix (we have to write this ourselves);
 *  and (iii) an assembler of surface term, c (already written).
 *
 *  Firstly, include `AbstractFeVolumeIntegralAssembler` which the assembler we write will inherit from,
 *  and `AbstractDynamicLinearPdeSolver`, which the solver we write will inherit from.
 */
#include <cxxtest/TestSuite.h>
#include "AbstractFeVolumeIntegralAssembler.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
/* Some standard includes */
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
/* The two assemblers that we can use */
#include "MassMatrixAssembler.hpp"
#include "NaturalNeumannSurfaceTermAssembler.hpp"
/* Ignore these for the time being */
//#include "HeatEquation.hpp"
//#include "SimpleLinearParabolicSolver.hpp"

/* == Writing assemblers ==
 *
 * We need to write an assembler for setting up the matrix `M + dt K`.
 *
 * Any new assembler should inherit from `AbstractFeVolumeIntegralAssembler`, which deals with looping over
 * elements, looping over quadrature points, etc. Concrete classes need to provide the integrand for the matrix
 * or vector being assembled (exactly as in the previous tutorials). However, in general, the assembler
 * class can be used to assemble a matrix OR a vector OR both. The class we write here needs to assemble
 * a matrix but not a vector. Note that the parent class `AbstractFeVolumeIntegralAssembler` has two booleans
 * in the template list (as well as the dimension template parameters as normal) - these booleans say
 * whether this class will be assembling a vector or a matrix (or both).
 */
template<unsigned DIM>
class RhsMatrixAssembler
    : public AbstractFeVolumeIntegralAssembler<DIM,DIM,1/*problem dim*/,false /*doesn't assemble vectors*/,true/*assembles a matrix*/,NORMAL /*amount of interpolation*/>
{
private:
    /* Even when a class isn't being written for a very general dimensions sometimes it is a good idea
     * to define the following, and then use `ELEMENT_DIM` etc in the below, as it can make the code a
     * bit easier to understand.
     */
    static const unsigned ELEMENT_DIM = DIM;
    static const unsigned SPACE_DIM = DIM;
    static const unsigned PROBLEM_DIM = 1;

    /* We are assembling a matrix, we means we need to provide a `ComputeMatrixTerm()` method, to return the
     * elemental contribution to the RHS matrix. Note that `ELEMENT_DIM+1` is the number of
     * nodes in the element (=number of basis functions).
     */
    c_matrix<double,PROBLEM_DIM*(ELEMENT_DIM+1),PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeMatrixTerm(
                                                                                c_vector<double, ELEMENT_DIM+1> &rPhi,
                                                                                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                                                                                ChastePoint<SPACE_DIM> &rX,
                                                                                c_vector<double,PROBLEM_DIM> &rU,
                                                                                c_matrix<double, PROBLEM_DIM, SPACE_DIM> &rGradU /* not used */,
                                                                                Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> ret = zero_matrix<double>(ELEMENT_DIM+1,ELEMENT_DIM+1);

        double dt = PdeSimulationTime::GetPdeTimeStep();

        for (unsigned i=0; i<ELEMENT_DIM+1; i++) // essentially a loop over the basis functions
        {
            for (unsigned j=0; j<ELEMENT_DIM+1; j++) // essentially a loop over the basis functions
            {
                // mass matrix
                ret(i,j) = rPhi(i)*rPhi(j);
                // -dt * stiffness matrix
                for (unsigned dim=0; dim<SPACE_DIM; dim++)
                {
                    ret(i,j) -= dt * rGradPhi(dim,i)*rGradPhi(dim,j);
                }
            }
        }
        return ret;
        // this could been done more efficiently and succinctly
        // using outer_prod(rPhi, rPhi) and prod(trans(rGradPhi), rGradPhi);
    }

    /* (If we were (also) assembling a vector, we would also have to provide a `ComputeVectorTerm()` method, which is
     * very similar).
     *
     * Now write the constructor. */
public:
    RhsMatrixAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
        : AbstractFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,1,false,true,NORMAL>(pMesh)
    {
    }
};
/* That's the assembler written. The following solver class will show how to use it.
 *
 * == Writing the solver class ==
 *
 * The parent class here is `AbstractDynamicLinearPdeSolver`, which contains a linear system
 * (`this->mpLinearSystem`), and will deal with allocating memory and solving the linear system.
 * The concrete class needs to implement a `SetupLinearSystem()` method which completely sets
 * up the linear system. In this case, it needs to set the LHS matrix in the linear system to
 * be M, and set the RHS vector to be `rhs_matrix * current_soln`.
 */
template<unsigned DIM>
class ExplicitHeatEquationSolver : public AbstractDynamicLinearPdeSolver<DIM,DIM,1>
{
private:
    /* The constuctor will take in a mesh and a BCC, the latter will be stored as a member variable */
    BoundaryConditionsContainer<DIM,DIM,1>* mpBoundaryConditions;
    /* Declare a matrix for the RHS matrix */
    Mat mRhsMatrix;

    /* This is the main method which needs to be implemented. It takes in the current solution, and a
     * boolean saying whether the matrix (ie A in Ax=b) is being computed or not.
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {
        /* This is how to use assemblers to set up matrices. We declare a mass matrix assembler,
         * pass it the LHS matrix of the linear system, and tell it to assemble. We also declare
         * one of our purpose-built `RhsMatrixAssemblers`, pass it the matrix `mRhsMatrix`, and
         * tell it to assemble.
         *
         * '''Important note''': if any of the assemblers will require the current solution (ie solution
         * at the current timestep), this needs to be passed to the assembler, as in the commented
         * line below.
         */
        if (computeMatrix)
        {
            MassMatrixAssembler<DIM,DIM> mass_matrix_assembler(this->mpMesh);
            RhsMatrixAssembler<DIM> rhs_matrix_assembler(this->mpMesh);

            mass_matrix_assembler.SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
            mass_matrix_assembler.AssembleMatrix();

            rhs_matrix_assembler.SetMatrixToAssemble(mRhsMatrix);
            //rhs_matrix_assembler.SetCurrentSolution(currentSolution);
            rhs_matrix_assembler.AssembleMatrix();


            this->mpLinearSystem->FinaliseLhsMatrix(); // (Petsc communication)
            PetscMatTools::Finalise(mRhsMatrix);       // (Petsc communication)
        }

        /* Use the RHS matrix to set up the RHS vector, ie set `b=(M+dtK)U^n` */
        MatMult(mRhsMatrix, currentSolution, this->mpLinearSystem->rGetRhsVector());

        /* The third assembler we use is the `NaturalNeumannSurfaceTermAssembler`, which assembles
         * the vector `c` defined above, using the Neumann BCs stored in the `BoundaryConditionsContainer`
         * which is passed in in the constructor
         */
        NaturalNeumannSurfaceTermAssembler<DIM,DIM,1> surface_integral_assembler(this->mpMesh, mpBoundaryConditions);
        surface_integral_assembler.SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false /*don't zero vector before assembling!*/);
        surface_integral_assembler.Assemble();

        /* Some necessary PETSc communication before applying Dirichet BCs */
        this->mpLinearSystem->FinaliseRhsVector();         // (Petsc communication)
        this->mpLinearSystem->SwitchWriteModeLhsMatrix();  // (Petsc communication - needs to called when going from adding entries to inserting entries)

        /* Apply the dirichlet BCs from the BCC to the linear system */
        mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);

        /* Some necessary PETSc communication to finish */
        this->mpLinearSystem->FinaliseRhsVector();
        this->mpLinearSystem->FinaliseLhsMatrix();
    }
public:
    /* The constructor needs to call the parent constructor, save the BCC, ''say that the (LHS) matrix is constant
     * in time'' (so it is only computed once), and allocate memory for the RHS matrix.
     */
    ExplicitHeatEquationSolver(TetrahedralMesh<DIM,DIM>* pMesh,
                               BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions)
         : AbstractDynamicLinearPdeSolver<DIM,DIM,1>(pMesh),
           mpBoundaryConditions(pBoundaryConditions)
    {
        this->mMatrixIsConstant = true;
        PetscTools::SetupMat(mRhsMatrix, this->mpMesh->GetNumNodes(), this->mpMesh->GetNumNodes(), 9);
    }

    /* Destructor */
    ~ExplicitHeatEquationSolver()
    {
        PetscTools::Destroy(mRhsMatrix);
    }
};
/* That's all that needs to be written to write your own solver using the solver hierarchy
 *
 * = A test using the solver =
 *
 * The following test uses the new solver. Since the interface is exactly the same as the
 * other solvers, except for not taking in a PDE (the fact that it solves a parameterless
 * heat equation is hardcoded into the solver), all of the below should be recognisable.
 * Note however the tiny timestep - this is needed for stability as this is an explicit scheme.
 * Also, to compare with the implicit solver, comment out the appropriate lines below. Note that
 * the implicit solver may seem quite slow in comparison - this is because the linear system is
 * much harder to solve (linear system is Ax=b, for explicit A=M, for implicit A=M-dt*K), but
 * remember that the implicit solver can use much larger timesteps.
 */

class TestWritingPdeSolversTwoTutorial : public CxxTest::TestSuite
{
public:
    void TestExplicitSolver()
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.05 /*h*/, 1.0 /*width*/, 1.0 /*height*/);

        // Set up BCs u=0 on entire boundary
        BoundaryConditionsContainer<2,2,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

        ExplicitHeatEquationSolver<2> solver(&mesh,&bcc);
        //// To use the old solver instead, comment out the above line
        //// and use these instead (also uncomment the appropriate includes).
        //HeatEquation<2> pde;
        //SimpleLinearParabolicSolver<2,2> solver(&mesh,&pde,&bcc);

        /* The interface is exactly the same as the `SimpleLinearParabolicSolver`. */
        solver.SetTimeStep(0.0001);
        solver.SetTimes(0.0, 0.2);

        std::vector<double> init_cond(mesh.GetNumNodes(), 0.0);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            double distance_from_centre = sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) );
            if (distance_from_centre < 1.0/3.0)
            {
                init_cond[i] = 1.0;
            }
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        solver.SetOutputDirectoryAndPrefix("ExplicitHeatEquationSolver","results");

        solver.SetOutputToTxt(true);

        solver.SetPrintingTimestepMultiple(100);
        /* We are now ready to solve the system. */
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Check nothing has changed in this tutorial
        TS_ASSERT_DELTA(result_repl[220], 0.019512, 1e-4);

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }
};

#endif // TESTWRITINGPDESOLVERSTWOTUTORIAL_HPP_
