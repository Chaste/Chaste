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

#ifndef TESTINCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_
#define TESTINCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "IncompressibleNonlinearElasticitySolver.hpp"
#include "ExponentialMaterialLaw.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "NonlinearElasticityTools.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "NumericFileComparison.hpp"
#include "VtkNonlinearElasticitySolutionWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

double MATERIAL_PARAM = 1.0;
double ALPHA = 0.2;

// Body force corresponding to the deformation
// x = X+0.5*alpha*X^2, y=Y/(1+alpha*X), with p=2c
//
//   TO TEST THE TIME-DEPENDENCE, THIS REQUIRES THE
//   CURRENT TIME TO BE SET TO 1.0
//
c_vector<double,2> MyBodyForce(c_vector<double,2>& rX, double t)
{
    assert(rX(0)>=0 && rX(0)<=1 && rX(1)>=0 && rX(1)<=1);

    c_vector<double,2> body_force;
    double lam = 1+ALPHA*rX(0);
    body_force(0) = -2*MATERIAL_PARAM * ALPHA;
    body_force(1) = -2*MATERIAL_PARAM * 2*ALPHA*ALPHA*rX(1)/(lam*lam*lam);

    // Make sure the time has been passed through to here correctly.
    // This function requires t=1 for the test to pass
    body_force(0) += (t-1)*5723485;
    return body_force;
}

// Surface traction on three sides of a cube, corresponding to
// x = X+0.5*alpha*X^2, y=Y/(1+alpha*X), with p=2c
//
//   TO TEST THE TIME-DEPENDENCE, THIS REQUIRES THE
//   CURRENT TIME TO BE SET TO 1.0
//
c_vector<double,2> MyTraction(c_vector<double,2>& location, double t)
{
    c_vector<double,2> traction = zero_vector<double>(2);

    double lam = 1+ALPHA*location(0);
    if (fabs(location(0)-1.0) <= DBL_EPSILON) //Right edge
    {
        traction(0) =  2*MATERIAL_PARAM * (lam - 1.0/lam);
        traction(1) = -2*MATERIAL_PARAM * location(1)*ALPHA/(lam*lam);
    }
    else if (fabs(location(1))  <= DBL_EPSILON) //Bottom edge
    {
        traction(0) =  2*MATERIAL_PARAM * location(1)*ALPHA/(lam*lam);
        traction(1) = -2*MATERIAL_PARAM * (-lam + 1.0/lam);
    }
    else if (fabs(location(1) - 1.0) <= DBL_EPSILON)//Top edge
    {
        traction(0) = -2*MATERIAL_PARAM * location(1)*ALPHA/(lam*lam);
        traction(1) =  2*MATERIAL_PARAM * (-lam + 1.0/lam);
    }
    else
    {
        NEVER_REACHED;
    }

    // Make sure the time has been passed through to here correctly.
    // This function requires t=1 for the test to pass
    traction(0) += (t-1.0)*543548;

    return traction;
}

// see TestSolveWithPressureBcsOnDeformedSurface()
double MyPressureFunction(double t)
{
    // the hardcoded value 1.32317 comes from outputting
    // double pressure = (2*c1*(pow(lambda,-1) - lambda*lambda*lambda))/lambda;
    return 1.32317 + 100*(t-1.0);
}

class TestIncompressibleNonlinearElasticitySolver : public CxxTest::TestSuite
{
public:
    // This is purely for coverage of assembling a 3D system...
    void TestAssembleSystem3D()
    {
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader1("mesh/test/data/3D_Single_tetrahedron_element_quadratic",2,1,false);

        mesh.ConstructFromMeshReader(mesh_reader1);

        ExponentialMaterialLaw<3> law(2.0, 3.0);
        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);

        SolidMechanicsProblemDefinition<3> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);

        IncompressibleNonlinearElasticitySolver<3> solver(mesh,
                                                          problem_defn,
                                                          "");

        solver.AssembleSystem(true, true);

        // cover exception
        CompressibleMooneyRivlinMaterialLaw<3> comp_law(2.0,1.0);
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&comp_law);
        TS_ASSERT_THROWS_THIS(IncompressibleNonlinearElasticitySolver<3> another_solver(mesh,problem_defn,""), "SolidMechanicsProblemDefinition object contains compressible material laws");
    }

    void TestAssembleSystem()
    {
        QuadraticMesh<2> mesh(0.5, 1.0, 1.0);
        ExponentialMaterialLaw<2> law(2.0, 3.0);
        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);

        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                          problem_defn,
                                                          "");

        ///////////////////////////////////////////////////////////////////
        // test whether residual vector is currently zero (as
        // current solution should have been initialised to u=0, p=p0
        ///////////////////////////////////////////////////////////////////
        solver.AssembleSystem(true, true);
        ReplicatableVector rhs_vec0(solver.mResidualVector);

        TS_ASSERT_EQUALS( rhs_vec0.GetSize(), 3U*25U);

        for (unsigned i=0; i<rhs_vec0.GetSize(); i++)
        {
            TS_ASSERT_DELTA(rhs_vec0[i], 0.0, 1e-12);
        }

        ///////////////////////////////////////////////////////////////////
        // Include pressure-on-deformed-body boundary conditions, which
        // add contributions to residual and jacobian
        ///////////////////////////////////////////////////////////////////
        std::vector<BoundaryElement<1,2>*> boundary_elems;
        double pressure = 3.0;

        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            BoundaryElement<1,2>* p_element = *iter;
            boundary_elems.push_back(p_element);
        }

        problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, pressure);

        solver.mLastDampingValue = 1.0;
        solver.AssembleSystem(true, true); // See comments in AbstractNonlinearElasticitySolver::ShouldAssembleMatrixTermForPressureOnDeformedBc
        ReplicatableVector rhs_vec(solver.mResidualVector);


        ///////////////////////////////////////////////////////////////////
        // compute numerical Jacobian and compare with analytic jacobian
        // (about u=0, p=p0)
        ///////////////////////////////////////////////////////////////////
        unsigned num_dofs = rhs_vec.GetSize();
        double h = 1e-6;

        int lo, hi;
        MatGetOwnershipRange(solver.mrJacobianMatrix, &lo, &hi);

        for (unsigned j=0; j<num_dofs; j++)
        {
            solver.rGetCurrentSolution().clear();
            solver.FormInitialGuess();
            solver.rGetCurrentSolution()[j] += h;

            solver.AssembleSystem(true, false);

            ReplicatableVector perturbed_rhs( solver.mResidualVector );

            for (unsigned i=0; i<num_dofs; i++)
            {
                if ((lo<=(int)i) && ((int)i<hi))
                {
                    double analytic_matrix_val = PetscMatTools::GetElement(solver.mrJacobianMatrix,i,j);
                    double numerical_matrix_val = (perturbed_rhs[i] - rhs_vec[i])/h;
                    if ((fabs(analytic_matrix_val)>1e-6) && (fabs(numerical_matrix_val)>1e-6))
                    {
                        // relative error
                        TS_ASSERT_DELTA( (analytic_matrix_val-numerical_matrix_val)/analytic_matrix_val, 0.0, 1e-2);
                    }
                    else
                    {
                        // absolute error
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-2);
                    }
                }
            }
        }
        PetscTools::Barrier();

        //////////////////////////////////////////////////////////
        // compare numerical and analytic jacobians again, this
        // time using a non-zero displacement, u=lambda x, v = mu y
        // (lambda not equal to 1/nu), p = x*y
        //////////////////////////////////////////////////////////
        double lambda = 1.2;
        double mu = 1.0/1.3;

        solver.rGetCurrentSolution().clear();
        solver.FormInitialGuess();
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            solver.rGetCurrentSolution()[3*i]   = (lambda-1)*mesh.GetNode(i)->rGetLocation()[0];
            solver.rGetCurrentSolution()[3*i+1] = (mu-1)*mesh.GetNode(i)->rGetLocation()[1];
            // should really be only setting this on vertices..
            solver.rGetCurrentSolution()[3*i+2] = mesh.GetNode(i)->rGetLocation()[0]*mesh.GetNode(i)->rGetLocation()[1];
        }

        solver.AssembleSystem(true, true);
        ReplicatableVector rhs_vec2(solver.mResidualVector);

        h=1e-8; // needs to be smaller for this one

        // check identity block
        for (unsigned i=mesh.GetNumVertices(); i<mesh.GetNumNodes(); i++)
        {
            assert(mesh.GetNode(i)->IsInternal());
            unsigned row = 3*i+2;
            if ((lo<=(int)row) && ((int)row<hi))
            {
                for (unsigned col=0; col<num_dofs; col++)
                {
                    double val = PetscMatTools::GetElement(solver.mrJacobianMatrix,row,col);
                    TS_ASSERT_DELTA(val, (double)(row==col), 1e-12);
                }
            }
        }

        // compare with numerical jacobian
        for (unsigned j=0; j<num_dofs; j++)
        {
            solver.rGetCurrentSolution()[j] += h;
            solver.AssembleSystem(true, false);

            ReplicatableVector perturbed_rhs( solver.mResidualVector );


            for (unsigned i=0; i<num_dofs; i++)
            {
                if ((lo<=(int)i) && ((int)i<hi))
                {
                    double analytic_matrix_val = PetscMatTools::GetElement(solver.mrJacobianMatrix,i,j);
                    double numerical_matrix_val = (perturbed_rhs[i] - rhs_vec2[i])/h;

                    if ((fabs(analytic_matrix_val)>1e-6) && (fabs(numerical_matrix_val)>1e-6))
                    {
                        // relative error
                        TS_ASSERT_DELTA( (analytic_matrix_val-numerical_matrix_val)/analytic_matrix_val, 0.0, 1e-2);
                    }
                    else
                    {
                        // absolute error
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-2);
                    }
                }

            }
            solver.rGetCurrentSolution()[j] -= h;
        }
    }

    void TestComputeResidualAndGetNormWithBadDeformation()
    {
        // There are 10 nodes (including internals) which means that the (non)linear system
        // can't be distributed over more than 10 processes.  This may cause deadlock in
        // VecNorm when calculating the norm of the residual at the end of this test.
        if (PetscTools::GetNumProcs() > 10u)
        {
            TS_TRACE("This test is designed for fewer processes!");
            return;
        }
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader1("mesh/test/data/3D_Single_tetrahedron_element_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader1);
        NashHunterPoleZeroLaw<3> law;

        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);

        SolidMechanicsProblemDefinition<3> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        IncompressibleNonlinearElasticitySolver<3> solver(mesh,
                                                          problem_defn,
                                                          "");

        // compute the residual norm - should be zero as no force or tractions
        TS_ASSERT_DELTA( solver.ComputeResidualAndGetNorm(false), 0.0, 1e-7);

        // the change the current solution (=displacement) to correspond to a small stretch
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                solver.rGetCurrentSolution()[3*i+j] = 0.01*mesh.GetNode(i)->rGetLocation()[j];
            }
        }

        // compute the residual norm - check computes fine
        TS_ASSERT_LESS_THAN( 0.0, solver.ComputeResidualAndGetNorm(false));

        // the change the current solution (=displacement) to correspond to a large stretch
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                solver.rGetCurrentSolution()[3*i+j] = mesh.GetNode(i)->rGetLocation()[j];
            }
        }

        // compute the residual norm - material law should complain - exception thrown...
        TS_ASSERT_THROWS_CONTAINS( solver.ComputeResidualAndGetNorm(false), "strain unacceptably large");
        // ..unless we set the allowException parameter to be true, in which case we should get
        // infinity returned.
        TS_ASSERT_EQUALS( solver.ComputeResidualAndGetNorm(true), DBL_MAX);
    }


    // A test where the solution should be zero displacement
    // It mainly tests that the initial guess was set up correctly to
    // the final correct solution, ie u=0, p=zero_strain_pressure (!=0)
    void TestWithZeroDisplacement()
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        double c1 = 3.0;
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(c1);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&mooney_rivlin_law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);


        IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                          problem_defn,
                                                          "");

        // for coverage
        TS_ASSERT_THROWS_THIS(solver.SetWriteOutput(true),
                "Can\'t write output if no output directory was given in constructor");
        solver.SetWriteOutput(false);

        solver.Solve();
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 0u);

        TS_ASSERT_THROWS_CONTAINS(solver.CreateCmguiOutput(), "No output directory was given");
        TS_ASSERT_THROWS_CONTAINS(solver.CreateVtkOutput("Displacement"), "No output directory was given");

        // get deformed position
        std::vector<c_vector<double,2> >& r_deformed_position
            = solver.rGetDeformedPosition();

        for (unsigned i=0; i<r_deformed_position.size(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[0], r_deformed_position[i](0), 1e-8);
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[1], r_deformed_position[i](1), 1e-8);
        }

        // check the final pressure
        std::vector<double>& r_pressures = solver.rGetPressures();
        TS_ASSERT_EQUALS(r_pressures.size(), mesh.GetNumNodes());
        for (unsigned i=0; i<r_pressures.size(); i++)
        {
            TS_ASSERT_DELTA(r_pressures[i], 2*c1, 1e-6);
        }
    }

    void TestSettingUpHeterogeneousProblem()
    {
        // two element quad mesh on the square
        QuadraticMesh<2> mesh(1.0, 1.0, 1.0);

        MooneyRivlinMaterialLaw<2> law_1(5.0);
        MooneyRivlinMaterialLaw<2> law_2(1.0);
        std::vector<AbstractMaterialLaw<2>*> laws;
        laws.push_back(&law_1);
        laws.push_back(&law_2);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,laws);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                          problem_defn,
                                                          "");

        // Recall variables in the solution vector are ordered
        // [U0 V0 P0 U1 V1 P1 .. Un Vn Pn]
        // with P_i a dummy pressure variable if node i is not a vertex
        // pressure for node 0 at (0,0) in element 1 (c1=1.0 above)
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(solver.rGetCurrentSolution()[2], 2.0, 1e-6); //0*3+2
        // pressure for node 3 at (1,1) in element 0 (c1=5.0 above)
        TS_ASSERT_DELTA(mesh.GetNode(3)->rGetLocation()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(3)->rGetLocation()[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(solver.rGetCurrentSolution()[11], 10.0, 1e-6); //3*3+2
    }

    void TestSolve()
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        MooneyRivlinMaterialLaw<2> law(1.0);
        c_vector<double,2> body_force;
        body_force(0) = 3.0;
        body_force(1) = 0.0;

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetBodyForce(body_force);

        IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                          problem_defn,
                                                          "simple_nonlin_elas");

        solver.Solve();
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 4u); // 'hardcoded' answer, protects against Jacobian getting messed up

        std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

        double xend = 1.17199;
        double yend = 0.01001;

        ///////////////////////////////////////////////////////////
        // compare the solution at the corners with the values
        // obtained using the dealii finite elasticity solver
        //
        // Results have been visually checked to see they agree
        // (they do, virtually or completely overlapping.
        ///////////////////////////////////////////////////////////

        // bottom lhs corner should still be at (0,0)
        assert( fabs(mesh.GetNode(0)->rGetLocation()[0] - 0) < 1e-9 );
        assert( fabs(mesh.GetNode(0)->rGetLocation()[1] - 0) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[0](0), 0.0, 1e-6 );
        TS_ASSERT_DELTA( r_solution[0](1), 0.0, 1e-6 );

        // top lhs corner should still be at (0,1)
        assert( fabs(mesh.GetNode(3)->rGetLocation()[0] - 0) < 1e-9 );
        assert( fabs(mesh.GetNode(3)->rGetLocation()[1] - 1) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[3](0), 0.0, 1e-6 );
        TS_ASSERT_DELTA( r_solution[3](1), 1.0, 1e-6 );

        // DEALII value for bottom rhs corner is (1.17199,0.01001)
        assert( fabs(mesh.GetNode(1)->rGetLocation()[0] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(1)->rGetLocation()[1] - 0) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[1](0), xend, 1e-3 );
        TS_ASSERT_DELTA( r_solution[1](1), yend, 1e-3 );

        // DEALII value for top rhs corner is (1.17199,0.98999)
        assert( fabs(mesh.GetNode(2)->rGetLocation()[0] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(2)->rGetLocation()[1] - 1) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[2](0), xend,   1e-3 );
        TS_ASSERT_DELTA( r_solution[2](1), 1-yend, 1e-3 );

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }

    /**
     *  Solve a problem with non-zero dirichlet boundary conditions
     *  and non-zero tractions. THIS TEST COMPARES AGAINST AN EXACT SOLUTION.
     *
     *  Choosing the deformation x=X/lambda, y=lambda*Y, with a
     *  Mooney-Rivlin material, then
     *   F = [1/lam 0; 0 lam], T = [2*c1-p*lam^2, 0; 0, 2*c1-p/lam^2],
     *   sigma = [2*c1/lam^2-p, 0; 0, 2*c1*lam^2-p].
     *  Choosing p=2*c1*lam^2, then sigma = [2*c1/lam^2-p 0; 0 0].
     *  The surface tractions are then
     *   TOP and BOTTOM SURFACE: 0
     *   RHS: s = SN = J*invF*sigma*N = [lam 0; 0 1/lam]*sigma*[1,0]
     *          = [2*c1(1/lam-lam^3), 0]
     *
     *  So, we have to specify displacement boundary conditions (y=lam*Y) on
     *  the LHS (X=0), and traction bcs (s=the above) on the RHS (X=1), and can
     *  compare the computed displacement and pressure against the true solution.
     *
     */
    void TestSolveWithNonZeroBoundaryConditions()
    {
        double lambda = 0.85;
        double c1 = 1.0;
        unsigned num_elem = 5;

        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

//        TrianglesMeshReader<2,2> reader("mesh/test/data/square_128_elements_quadratic_reordered",2,1,false);
//        mesh.ConstructFromMeshReader(reader);


        MooneyRivlinMaterialLaw<2> law(c1);

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,2> > locations;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if (fabs(mesh.GetNode(i)->rGetLocation()[0]) < 1e-6)
            {
                fixed_nodes.push_back(i);
                c_vector<double,2> new_position;
                new_position(0) = 0;
                new_position(1) = lambda*mesh.GetNode(i)->rGetLocation()[1];
                locations.push_back(new_position);
            }
        }

        std::vector<BoundaryElement<1,2>*> boundary_elems;
        std::vector<c_vector<double,2> > tractions;
        c_vector<double,2> traction;
        traction(0) = 2*c1*(pow(lambda,-1) - lambda*lambda*lambda);
        traction(1) = 0;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            if (fabs((*iter)->CalculateCentroid()[0] - 1.0)<1e-4)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
                tractions.push_back(traction);
            }
        }
        assert(boundary_elems.size()==num_elem);


        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetFixedNodes(fixed_nodes, locations);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, tractions);

        IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                          problem_defn,
                                                          "nonlin_elas_non_zero_bcs");

        /////////////////////////////////////////////////////////////////
        // Provide the exact solution as the initial guess and check
        // residual is (nearly) exactly zero (as the solution is in
        // the FEM space, the FEM solution, assuming the nonlinear
        // system was solved exactly, is the exact solution
        /////////////////////////////////////////////////////////////////

        std::vector<double> old_current_soln = solver.rGetCurrentSolution();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double exact_x = (1.0/lambda)*mesh.GetNode(i)->rGetLocation()[0];
            double exact_y = lambda*mesh.GetNode(i)->rGetLocation()[1];

            solver.rGetCurrentSolution()[3*i] = exact_x - mesh.GetNode(i)->rGetLocation()[0];
            solver.rGetCurrentSolution()[3*i+1] = exact_y - mesh.GetNode(i)->rGetLocation()[1];

            if (mesh.GetNode(i)->IsInternal())
            {
                solver.rGetCurrentSolution()[3*i+2] =  0.0;
            }
            else
            {
                solver.rGetCurrentSolution()[3*i+2] =  2*c1*lambda*lambda;
            }
        }

        /* HOW_TO_TAG Continuum mechanics
         * Get or output stresses during a solve
         */

        // get the solver to save the stresses on each element (averaged over quad point stresses)
        solver.SetComputeAverageStressPerElementDuringSolve();


        solver.Solve();

        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 0u); // initial guess was solution

        // test stresses. The 1st PK stress should satisfy S = [s(0) 0 ; 0 0], where s is the
        // applied traction. This has to be multiplied by F^{-T} to get the 2nd PK stress.
        assert(solver.mAverageStressesPerElement.size()==mesh.GetNumElements()); //Will fail when we move to DistributedQuadraticMesh
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            if (mesh.CalculateDesignatedOwnershipOfElement(i))
            {
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,0), lambda*traction(0), 1e-8);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(1,0), 0.0, 1e-8);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,1), 0.0, 1e-8);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(1,1), 0.0, 1e-8);
            }
        }

        ///////////////////////////////////////////////////////////////////////////
        // Now solve properly
        ///////////////////////////////////////////////////////////////////////////

        solver.rGetCurrentSolution() = old_current_soln;
        // coverage
        solver.SetKspAbsoluteTolerance(1e-10);

        solver.Solve();

        // write the stresses
        solver.WriteCurrentAverageElementStresses("solution");


        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 3u); // 'hardcoded' answer, protects against Jacobian getting messed up

        std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

        for (unsigned i=0; i<fixed_nodes.size(); i++)
        {
            unsigned index = fixed_nodes[i];
            TS_ASSERT_DELTA(r_solution[index](0), locations[i](0), 1e-8);
            TS_ASSERT_DELTA(r_solution[index](1), locations[i](1), 1e-8);
        }

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double exact_x = (1.0/lambda)*mesh.GetNode(i)->rGetLocation()[0];
            double exact_y = lambda*mesh.GetNode(i)->rGetLocation()[1];

            TS_ASSERT_DELTA( r_solution[i](0), exact_x, 1e-5 );
            TS_ASSERT_DELTA( r_solution[i](1), exact_y, 1e-5 );
        }

        std::vector<double>& r_pressures = solver.rGetPressures();
        TS_ASSERT_EQUALS(r_pressures.size(), mesh.GetNumNodes());
        for (unsigned i=0; i<r_pressures.size(); i++)
        {
            TS_ASSERT_DELTA(r_pressures[i], 2*c1*lambda*lambda, 1e-5);
        }

        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            if (mesh.CalculateDesignatedOwnershipOfElement(i))
            {
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,0), lambda*traction(0), 1e-3);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(1,0), 0.0, 1e-3);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,1), 0.0, 1e-3);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(1,1), 0.0, 1e-3);
            }
        }

        // check the written stresses
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        NumericFileComparison comparison(test_output_directory + "/nonlin_elas_non_zero_bcs/solution.stress", "continuum_mechanics/test/data/exact.stress");
        TS_ASSERT(comparison.CompareFiles(2e-4));

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }

    /**
     *  Test with functional (rather than constant) body force and surface traction, against a known
     *  solution. Since a non-zero body force is used here and a known solution, this is the MOST
     *  IMPORTANT TEST.
     *
     *  See equations and finite element implementations document for more info.
     *
     *  Choose x=X+0.5*alpha*X*X, y=Y/(1+alpha*X), p=2c, then F has determinant 1, and S can be shown to
     *  be, where lam = 1 +alpha X (ie dx/dX)
     *    S = 2c[lam-1/lam,   -Y*alpha*lam^{-2}; -Y*alpha*lam^{-2}, 1/lam - lam]
     *  in which case the required body force and surface traction can be computed to be
     *
     *  b = 2c/density [ -alpha, -2*Y*alpha^2 * lam^{-3} ]
     *  s = 2c[lam-1/lam, -Y*alpha/(lam^2)] on X=1
     *  s = 2c[0, lam - 1/lam]              on Y=0
     *  s = 2c[Y*alpha/lam^2, 1/lam - lam]  on Y=1
     *
     */
    void TestAgainstExactSolution()
    {
        for (unsigned run=0; run<2; run++)
        {
            MechanicsEventHandler::Reset();

            unsigned num_elem = 5;
            QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

            MooneyRivlinMaterialLaw<2> law(MATERIAL_PARAM);

            std::vector<unsigned> fixed_nodes
              = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);


            std::vector<BoundaryElement<1,2>*> boundary_elems;
            for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
                  = mesh.GetBoundaryElementIteratorBegin();
                iter != mesh.GetBoundaryElementIteratorEnd();
                ++iter)
            {
                // get all boundary elems except those on X=0
                if (fabs((*iter)->CalculateCentroid()[0])>1e-6)
                {
                    BoundaryElement<1,2>* p_element = *iter;
                    boundary_elems.push_back(p_element);
                }
            }
            assert(boundary_elems.size()==3*num_elem);

            SolidMechanicsProblemDefinition<2> problem_defn(mesh);


            problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
            problem_defn.SetZeroDisplacementNodes(fixed_nodes);
            problem_defn.SetBodyForce(MyBodyForce);
            problem_defn.SetTractionBoundaryConditions(boundary_elems, MyTraction);

            if (run==1)
            {
                problem_defn.SetVerboseDuringSolve(); // coverage
                problem_defn.SetSolveUsingSnes();
            }

            IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                              problem_defn,
                                                              "nonlin_elas_functional_data");

            // this test requires the time to be set to t=1 to pass (see comment
            // in and MyBodyForce() and MyTraction()
            solver.SetCurrentTime(1.0);

            // cover the option of writing output for each iteration
            solver.SetWriteOutputEachNewtonIteration();

            solver.Solve();

            if (run==0)
            {
                // matrix might have (small) errors introduced if this fails
                TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 3u); // 'hardcoded' answer, protects against Jacobian getting messed up
            }

            // check CreateCmguiOutput() - call and check output files were written.
            solver.CreateCmguiOutput();

            // Create .vtu file
            VtkNonlinearElasticitySolutionWriter<2> vtk_writer(solver);
            vtk_writer.Write();

            std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double X = mesh.GetNode(i)->rGetLocation()[0];
                double Y = mesh.GetNode(i)->rGetLocation()[1];

                double exact_x = X + 0.5*ALPHA*X*X;
                double exact_y = Y/(1+ALPHA*X);

                TS_ASSERT_DELTA(r_solution[i](0), exact_x, 1e-4);
                TS_ASSERT_DELTA(r_solution[i](1), exact_y, 1e-4);
            }

            std::vector<double>& r_pressures = solver.rGetPressures();
            for (unsigned i=0; i<r_pressures.size(); i++)
            {
                TS_ASSERT_DELTA( r_pressures[i]/(2*MATERIAL_PARAM), 1.0, 1e-3);
            }

            MechanicsEventHandler::Headings();
            MechanicsEventHandler::Report();

            if (run==0)
            {
                // Check output files were created: the standard output files initial.nodes and solution.nodes, the extra newton iteration
                // output files created as SetWriteOutputEachNewtonIteration() was called above, and the cmgui files created as
                // CreateCmguiOutput() was called above.
                FileFinder init_file("nonlin_elas_functional_data/initial.nodes", RelativeTo::ChasteTestOutput);
                TS_ASSERT(init_file.Exists());
                FileFinder nodes_file("nonlin_elas_functional_data/solution.nodes", RelativeTo::ChasteTestOutput);
                TS_ASSERT(nodes_file.Exists());
                FileFinder newton_file("nonlin_elas_functional_data/newton_iteration_3.nodes", RelativeTo::ChasteTestOutput);
                TS_ASSERT(newton_file.Exists());
                FileFinder exelem_file("nonlin_elas_functional_data/cmgui/solution_0.exelem", RelativeTo::ChasteTestOutput);
                TS_ASSERT(exelem_file.Exists());
                FileFinder exnode0_file("nonlin_elas_functional_data/cmgui/solution_0.exnode", RelativeTo::ChasteTestOutput);
                TS_ASSERT(exnode0_file.Exists());
                FileFinder exnode1_file("nonlin_elas_functional_data/cmgui/solution_1.exnode", RelativeTo::ChasteTestOutput);
                TS_ASSERT(exnode1_file.Exists());

#ifdef CHASTE_VTK
                //Check the VTK file exists
                FileFinder vtk_file("nonlin_elas_functional_data/vtk/solution.vtu", RelativeTo::ChasteTestOutput);
                TS_ASSERT(vtk_file.Exists());
#endif

                solver.rGetCurrentSolution().clear();
                solver.rGetCurrentSolution().resize(solver.mNumDofs, 0.0);
                solver.SetKspAbsoluteTolerance(1); // way too high
                TS_ASSERT_THROWS_CONTAINS(solver.Solve(), "KSP Absolute tolerance was too high");
            }
        }
    }

    /*
     *  Test the functionality for specifying that a pressure should act in the normal direction on the
     *  DEFORMED SURFACE.
     *
     *  The deformation is based on that in TestSolveWithNonZeroBoundaryConditions (x=X/lambda,
     *  y=Y*lam; see comments for this test), but the exact solution here is this rotated by
     *  45 degrees anticlockwise. We choose dirichlet boundary conditions on the X=0 surface to
     *  match this, and a pressure on the opposite surface (similar to the traction provided in
     *  TestSolveWithNonZeroBoundaryConditions, except we don't provide the direction, the code
     *  needs to work this out), and it is scaled by 1.0/lambda as it acts on a smaller surface
     *  than would on the undeformed surface.
     */
    void TestSolveWithPressureBcsOnDeformedSurface()
    {
        std::vector<double> soln_first_run;

        for (unsigned run=0; run<2; run++)
        {
            double lambda = 0.85;
            double c1 = 1.0;
            unsigned num_elem = 10;

            QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);
            MooneyRivlinMaterialLaw<2> law(c1);

            std::vector<unsigned> fixed_nodes;
            std::vector<c_vector<double,2> > locations;
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                if (fabs(mesh.GetNode(i)->rGetLocation()[0]) < 1e-6)
                {
                    fixed_nodes.push_back(i);
                    c_vector<double,2> new_position;
                    new_position(0) = (lambda/sqrt(2.0)) * mesh.GetNode(i)->rGetLocation()[1];
                    new_position(1) = (lambda/sqrt(2.0)) * mesh.GetNode(i)->rGetLocation()[1];
                    locations.push_back(new_position);
                }
            }

            std::vector<BoundaryElement<1,2>*> boundary_elems;
            double pressure = (2*c1*(pow(lambda,-1) - lambda*lambda*lambda))/lambda;

            for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
                  = mesh.GetBoundaryElementIteratorBegin();
                iter != mesh.GetBoundaryElementIteratorEnd();
                ++iter)
            {
                if (fabs((*iter)->CalculateCentroid()[0] - 1.0)<1e-4)
                {
                    BoundaryElement<1,2>* p_element = *iter;
                    boundary_elems.push_back(p_element);
                }
            }
            assert(boundary_elems.size()==num_elem);

            SolidMechanicsProblemDefinition<2> problem_defn(mesh);

            problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
            problem_defn.SetFixedNodes(fixed_nodes, locations);

            // the two variants of SetApplyNormalPressureOnDeformedSurface(). MyPressureFunction
            // will return the same value as that in the variable pressure IF the current time
            // is set to 1.0.
            if (run==0)
            {
                problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, pressure);
            }
            else
            {
                problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, MyPressureFunction);
            }

            IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                              problem_defn,
                                                              "nonlin_elas_pressure_on_deformed");

            if (run==1)
            {
                solver.SetCurrentTime(1.0);
                // To speed up this test, we provide the correct answer as the initial guess for
                // the second run. Note: the second run may still take an iteration or two, as the
                // final newton solve tolerance may be different (eg if relative tolerance)
                solver.rGetCurrentSolution() = soln_first_run;
            }

            solver.Solve();

            if (run==0)
            {
                soln_first_run = solver.rGetCurrentSolution();
            }

            std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double exact_x_before_rotation = (1.0/lambda)*mesh.GetNode(i)->rGetLocation()[0];
                double exact_y_before_rotation = lambda*mesh.GetNode(i)->rGetLocation()[1];

                double exact_x = (1.0/sqrt(2.0))*( exact_x_before_rotation + exact_y_before_rotation);
                double exact_y = (1.0/sqrt(2.0))*(-exact_x_before_rotation + exact_y_before_rotation);

                TS_ASSERT_DELTA( r_solution[i](0), exact_x, 1e-3 );
                TS_ASSERT_DELTA( r_solution[i](1), exact_y, 1e-3 );
            }

            // check the final pressure
            std::vector<double>& r_pressures = solver.rGetPressures();
            TS_ASSERT_EQUALS(r_pressures.size(), mesh.GetNumNodes());
            for (unsigned i=0; i<r_pressures.size(); i++)
            {
                TS_ASSERT_DELTA(r_pressures[i], 2*c1*lambda*lambda, 5e-2 );
            }
        }
    }

    /*
     *  Repeat of previous test but with SNES solver from PETSc.
     *
     *  Test the functionality for specifying that a pressure should act in the normal direction on the
     *  DEFORMED SURFACE.
     *
     *  The deformation is based on that in TestSolveWithNonZeroBoundaryConditions (x=X/lambda,
     *  y=Y*lam; see comments for this test), but the exact solution here is this rotated by
     *  45 degrees anticlockwise. We choose dirichlet boundary conditions on the X=0 surface to
     *  match this, and a pressure on the opposite surface (similar to the traction provided in
     *  TestSolveWithNonZeroBoundaryConditions, except we don't provide the direction, the code
     *  needs to work this out), and it is scaled by 1.0/lambda as it acts on a smaller surface
     *  than would on the undeformed surface.
     */
    void TestSolveWithPressureBcsOnDeformedSurfaceSnes()
    {
        std::vector<double> soln_first_run;

        for (unsigned run=0; run<2; run++)
        {
            double lambda = 0.85;
            double c1 = 1.0;
            unsigned num_elem = 10;

            QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);
            MooneyRivlinMaterialLaw<2> law(c1);

            std::vector<unsigned> fixed_nodes;
            std::vector<c_vector<double,2> > locations;
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                if (fabs(mesh.GetNode(i)->rGetLocation()[0]) < 1e-6)
                {
                    fixed_nodes.push_back(i);
                    c_vector<double,2> new_position;
                    new_position(0) = (lambda/sqrt(2.0)) * mesh.GetNode(i)->rGetLocation()[1];
                    new_position(1) = (lambda/sqrt(2.0)) * mesh.GetNode(i)->rGetLocation()[1];
                    locations.push_back(new_position);
                }
            }

            std::vector<BoundaryElement<1,2>*> boundary_elems;
            double pressure = (2*c1*(pow(lambda,-1) - lambda*lambda*lambda))/lambda;

            for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
                  = mesh.GetBoundaryElementIteratorBegin();
                iter != mesh.GetBoundaryElementIteratorEnd();
                ++iter)
            {
                if (fabs((*iter)->CalculateCentroid()[0] - 1.0)<1e-4)
                {
                    BoundaryElement<1,2>* p_element = *iter;
                    boundary_elems.push_back(p_element);
                }
            }
            assert(boundary_elems.size()==num_elem);

            SolidMechanicsProblemDefinition<2> problem_defn(mesh);

            problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
            problem_defn.SetFixedNodes(fixed_nodes, locations);

            /*** Switch on SNES this time ***/
            problem_defn.SetSolveUsingSnes();
            /*** Switch on SNES this time ***/

            // the two variants of SetApplyNormalPressureOnDeformedSurface(). MyPressureFunction
            // will return the same value as that in the variable pressure IF the current time
            // is set to 1.0.
            if (run==0)
            {
                problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, pressure);
            }
            else
            {
                problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, MyPressureFunction);
            }

            IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                              problem_defn,
                                                              "nonlin_elas_pressure_on_deformed");

            if (run==1)
            {
                solver.SetCurrentTime(1.0);
                // To speed up this test, we provide the correct answer as the initial guess for
                // the second run. Note: the second run may still take an iteration or two, as the
                // final newton solve tolerance may be different (eg if relative tolerance)
                solver.rGetCurrentSolution() = soln_first_run;
            }

            solver.Solve();

            if (run==0)
            {
                soln_first_run = solver.rGetCurrentSolution();
            }

            std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double exact_x_before_rotation = (1.0/lambda)*mesh.GetNode(i)->rGetLocation()[0];
                double exact_y_before_rotation = lambda*mesh.GetNode(i)->rGetLocation()[1];

                double exact_x = (1.0/sqrt(2.0))*( exact_x_before_rotation + exact_y_before_rotation);
                double exact_y = (1.0/sqrt(2.0))*(-exact_x_before_rotation + exact_y_before_rotation);

                TS_ASSERT_DELTA( r_solution[i](0), exact_x, 1e-3 );
                TS_ASSERT_DELTA( r_solution[i](1), exact_y, 1e-3 );
            }

            // check the final pressure
            std::vector<double>& r_pressures = solver.rGetPressures();
            TS_ASSERT_EQUALS(r_pressures.size(), mesh.GetNumNodes());
            for (unsigned i=0; i<r_pressures.size(); i++)
            {
                TS_ASSERT_DELTA(r_pressures[i], 2*c1*lambda*lambda, 5e-2 );
            }
        }
    }

    /**
     *  Same as TestSolveWithNonZeroBoundaryConditions()
     *  but using sliding boundary conditions, in order to test
     *  that sliding boundary conditions are implemented correctly
     *
     *  Choosing the deformation x=X/lambda, y=lambda*Y, with a
     *  Mooney-Rivlin material, then
     *   F = [1/lam 0; 0 lam], T = [2*c1-p*lam^2, 0; 0, 2*c1-p/lam^2],
     *   sigma = [2*c1/lam^2-p, 0; 0, 2*c1*lam^2-p].
     *  Choosing p=2*c1*lam^2, then sigma = [2*c1/lam^2-p 0; 0 0].
     *  The surface tractions are then
     *   TOP and BOTTOM SURFACE: 0
     *   RHS: s = SN = J*invF*sigma*N = [lam 0; 0 1/lam]*sigma*[1,0]
     *          = [2*c1(1/lam-lam^3), 0]
     *
     *  So, we have to specify partial displacement boundary conditions
     *  (x=0, y=FREE) on the LHS (X=0), and traction bcs (s=the above) on
     *  the RHS (X=1), and can
     *  compare the computed displacement and pressure against the true solution.
     *
     */
    void TestSlidingBoundaryConditions()
    {
        double lambda = 0.85;
        double c1 = 1.0;
        unsigned num_elem = 5;

        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);
        MooneyRivlinMaterialLaw<2> law(c1);

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,2> > locations;

        // fix node 0 (located at the origin) fully
        fixed_nodes.push_back(0);
        locations.push_back(zero_vector<double>(2));

        // for the rest of the nodes, if X=0, set x=0, leave y free.
        for (unsigned i=1; i<mesh.GetNumNodes(); i++)
        {
            if (fabs(mesh.GetNode(i)->rGetLocation()[0]) < 1e-6)
            {
                fixed_nodes.push_back(i);
                c_vector<double,2> new_position;
                new_position(0) = 0;
                new_position(1) = SolidMechanicsProblemDefinition<2>::FREE;
                locations.push_back(new_position);
            }
        }


        std::vector<BoundaryElement<1,2>*> boundary_elems;
        std::vector<c_vector<double,2> > tractions;
        c_vector<double,2> traction;
        traction(0) = 2*c1*(pow(lambda,-1) - lambda*lambda*lambda);
        traction(1) = 0;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            if (fabs((*iter)->CalculateCentroid()[0] - 1.0)<1e-4)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
                tractions.push_back(traction);
            }
        }
        assert(boundary_elems.size()==num_elem);


        SolidMechanicsProblemDefinition<2> problem_defn(mesh);

        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetFixedNodes(fixed_nodes, locations);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, tractions);

        // coverage
        problem_defn.SetVerboseDuringSolve();

        IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                          problem_defn,
                                                          "TestSlidingBoundaryConditions");


        // coverage
        solver.SetKspAbsoluteTolerance(1e-8);

        solver.Solve();
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 3u); // 'hardcoded' answer, protects against Jacobian getting messed up

        std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

        for (unsigned i=0; i<fixed_nodes.size(); i++)
        {
            unsigned index = fixed_nodes[i];
            TS_ASSERT_DELTA(r_solution[index](0), 0.0, 1e-8);

            double exact_y = lambda*mesh.GetNode(index)->rGetLocation()[1];
            TS_ASSERT_DELTA(r_solution[index](1), exact_y, 1e-5);
        }

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double exact_x = (1.0/lambda)*mesh.GetNode(i)->rGetLocation()[0];
            double exact_y = lambda*mesh.GetNode(i)->rGetLocation()[1];

            TS_ASSERT_DELTA( r_solution[i](0), exact_x, 1e-5 );
            TS_ASSERT_DELTA( r_solution[i](1), exact_y, 1e-5 );
        }


        std::vector<double>& r_pressures = solver.rGetPressures();
        for (unsigned i=0; i<r_pressures.size(); i++)
        {
            TS_ASSERT_DELTA(r_pressures[i], 2*c1*lambda*lambda, 1e-4 );
        }

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }


    // same set up as TestSolveWithNonZeroBoundaryConditions, here we test
    // that the simulation works fine if the internal nodes are not assumed to
    // have indices greater than vertex nodes
    void TestWithReordering()
    {
        double lambda = 0.85;
        double c1 = 1.0;

        std::vector<double> soln_normal;
        std::vector<double> soln_reordered;

        double end_residual_norm_normal =0.0;
        double end_residual_norm_reordered = 0.0;

        for (unsigned run=0; run<2; run++)
        {
            std::string mesh_file = (run==0 ? "mesh/test/data/square_128_elements_quadratic" : "mesh/test/data/square_128_elements_quadratic_reordered");

            QuadraticMesh<2> mesh;
            TrianglesMeshReader<2,2> reader(mesh_file,2,1,false);
            mesh.ConstructFromMeshReader(reader);

            // same BCs as in test mentioned above
            MooneyRivlinMaterialLaw<2> law(c1);
            std::vector<unsigned> fixed_nodes;
            std::vector<c_vector<double,2> > locations;
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                if (fabs(mesh.GetNode(i)->rGetLocation()[0]) < 1e-6)
                {
                    fixed_nodes.push_back(i);
                    c_vector<double,2> new_position;
                    new_position(0) = 0;
                    new_position(1) = lambda*mesh.GetNode(i)->rGetLocation()[1];
                    locations.push_back(new_position);
                }
            }

            std::vector<BoundaryElement<1,2>*> boundary_elems;
            std::vector<c_vector<double,2> > tractions;
            c_vector<double,2> traction;
            traction(0) = 2*c1*(pow(lambda,-1) - lambda*lambda*lambda);
            traction(1) = 0;
            for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
                  = mesh.GetBoundaryElementIteratorBegin();
                iter != mesh.GetBoundaryElementIteratorEnd();
                ++iter)
            {
                if (fabs((*iter)->CalculateCentroid()[0] - 1.0)<1e-4)
                {
                    BoundaryElement<1,2>* p_element = *iter;
                    boundary_elems.push_back(p_element);
                    tractions.push_back(traction);
                }
            }

            SolidMechanicsProblemDefinition<2> problem_defn(mesh);
            problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
            problem_defn.SetFixedNodes(fixed_nodes, locations);
            problem_defn.SetTractionBoundaryConditions(boundary_elems, tractions);

            IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                              problem_defn,
                                                              "");

            solver.Solve();

            if (run==0)
            {
                soln_normal = solver.rGetCurrentSolution();
                end_residual_norm_normal = solver.ComputeResidualAndGetNorm(false);
            }
            else
            {
                soln_reordered = solver.rGetCurrentSolution();
                end_residual_norm_reordered = solver.ComputeResidualAndGetNorm(false);
            }
        }

        // check the end residual norms were the same
        // Need to do relative error as they will both be very small;
        double rel_error = fabs(end_residual_norm_normal-end_residual_norm_reordered)/end_residual_norm_reordered;
        TS_ASSERT_LESS_THAN(rel_error, 0.011);


        // The two meshes are the same except the nodes 4 and 81 have been swapped around.
        // Hence, the solutions should be the same except for the unknowns at these nodes
        for (unsigned i=0; i<289 /*num total nodes*/; i++)
        {
            if (i!=4 && i!=81)
            {
                for (unsigned j=0; j<3; j++) // [u,v,p] unknowns
                {
                    TS_ASSERT_DELTA(soln_normal[3*i+j], soln_reordered[3*i+j], 7e-7);
                }
            }
        }

        // Check the solution at nodes 4 and 81
        for (unsigned j=0; j<3; j++) // [u,v,p] unknowns
        {
            TS_ASSERT_DELTA(soln_normal[3*4+j],  soln_reordered[3*81+j], 4e-7);
            TS_ASSERT_DELTA(soln_normal[3*81+j], soln_reordered[3*4+j],  5e-7);
        }
    }

    // same set up as SolveWithNonZeroBoundaryConditions, here we test
    // that the simulation works fine with a DistributedQuadraticMesh
    void TestSolveDistributedQuadraticMeshWithNonZeroBoundaryConditions()
    {
        double lambda = 0.85;
        double c1 = 1.0;

        DistributedQuadraticMesh<2> mesh;

        TrianglesMeshReader<2,2> reader("mesh/test/data/square_128_elements_quadratic_reordered",2,1,false);
        mesh.ConstructFromMeshReader(reader);


        MooneyRivlinMaterialLaw<2> law(c1);

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,2> > locations;
        for (AbstractTetrahedralMesh<2,2>::NodeIterator iter = mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            if (fabs(iter->rGetLocation()[0]) < 1e-6)
            {
                fixed_nodes.push_back(iter->GetIndex());
                c_vector<double,2> new_position;
                new_position(0) = 0;
                new_position(1) = lambda*iter->rGetLocation()[1];
                locations.push_back(new_position);
            }
        }

        std::vector<BoundaryElement<1,2>*> boundary_elems;
        std::vector<c_vector<double,2> > tractions;
        c_vector<double,2> traction;
        traction(0) = 2*c1*(pow(lambda,-1) - lambda*lambda*lambda);
        traction(1) = 0;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            if (fabs((*iter)->CalculateCentroid()[0] - 1.0)<1e-4)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
                tractions.push_back(traction);
            }
        }

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetFixedNodes(fixed_nodes, locations);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, tractions);

        IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                          problem_defn,
                                                          "nonlin_elas_non_zero_bcs");

        /////////////////////////////////////////////////////////////////
        // Provide the exact solution as the initial guess and check
        // residual is (nearly) exactly zero (as the solution is in
        // the FEM space, the FEM solution, assuming the nonlinear
        // system was solved exactly, is the exact solution
        /////////////////////////////////////////////////////////////////

        std::vector<double> old_current_soln = solver.rGetCurrentSolution();

        for (AbstractTetrahedralMesh<2,2>::NodeIterator iter = mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            double exact_x = (1.0/lambda)*iter->rGetLocation()[0];
            double exact_y = lambda*iter->rGetLocation()[1];

            solver.rGetCurrentSolution()[3*iter->GetIndex()] = exact_x - iter->rGetLocation()[0];
            solver.rGetCurrentSolution()[3*iter->GetIndex()+1] = exact_y - iter->rGetLocation()[1];

            if (iter->IsInternal())
            {
                solver.rGetCurrentSolution()[3*iter->GetIndex()+2] =  0.0;
            }
            else
            {
                solver.rGetCurrentSolution()[3*iter->GetIndex()+2] =  2*c1*lambda*lambda;
            }
        }

        // get the solver to save the stresses on each element (averaged over quad point stresses)
        solver.SetComputeAverageStressPerElementDuringSolve();

///\todo #2223 Trips exception  No Dirichlet boundary conditions` in `ContinuumMechanicsProblemDefinition`
//        solver.Solve();
//
//        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 0u); // initial guess was solution
//
//        // test stresses. The 1st PK stress should satisfy S = [s(0) 0 ; 0 0], where s is the
//        // applied traction. This has to be multiplied by F^{-T} to get the 2nd PK stress.
//        assert(solver.mAverageStressesPerElement.size()==mesh.GetNumElements()); //Will fail when we move to DistributedQuadraticMesh
//        for (unsigned i=0; i<mesh.GetNumElements(); i++)
//        {
//            if (mesh.CalculateDesignatedOwnershipOfElement(i))
//            {
//                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,0), lambda*traction(0), 1e-8);
//                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(1,0), 0.0, 1e-8);
//                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,1), 0.0, 1e-8);
//                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(1,1), 0.0, 1e-8);
//            }
//        }
//
//
//
//        ///////////////////////////////////////////////////////////////////////////
//        // Now solve properly
//        ///////////////////////////////////////////////////////////////////////////
//
//        solver.rGetCurrentSolution() = old_current_soln;
//        // coverage
//        solver.SetKspAbsoluteTolerance(1e-10);
//
//        solver.Solve();
//
//        // write the stresses
//        solver.WriteCurrentAverageElementStresses("solution");
//
//
//        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 3u); // 'hardcoded' answer, protects against Jacobian getting messed up
//
//        std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();
//
//        for (unsigned i=0; i<fixed_nodes.size(); i++)
//        {
//            unsigned index = fixed_nodes[i];
//            TS_ASSERT_DELTA(r_solution[index](0), locations[i](0), 1e-8);
//            TS_ASSERT_DELTA(r_solution[index](1), locations[i](1), 1e-8);
//        }
//
//        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
//        {
//            double exact_x = (1.0/lambda)*mesh.GetNode(i)->rGetLocation()[0];
//            double exact_y = lambda*mesh.GetNode(i)->rGetLocation()[1];
//
//            TS_ASSERT_DELTA( r_solution[i](0), exact_x, 1e-5 );
//            TS_ASSERT_DELTA( r_solution[i](1), exact_y, 1e-5 );
//        }
//
//        std::vector<double>& r_pressures = solver.rGetPressures();
//        TS_ASSERT_EQUALS(r_pressures.size(), mesh.GetNumNodes());
//        for (unsigned i=0; i<r_pressures.size(); i++)
//        {
//            TS_ASSERT_DELTA(r_pressures[i], 2*c1*lambda*lambda, 1e-5);
//        }
//
//        for (unsigned i=0; i<mesh.GetNumElements(); i++)
//        {
//            if (mesh.CalculateDesignatedOwnershipOfElement(i))
//            {
//                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,0), lambda*traction(0), 1e-3);
//                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(1,0), 0.0, 1e-3);
//                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,1), 0.0, 1e-3);
//                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(1,1), 0.0, 1e-3);
//            }
//        }
//
//        // check the written stresses
//        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
//        NumericFileComparison comparison(test_output_directory + "/nonlin_elas_non_zero_bcs/solution.stress", "continuum_mechanics/test/data/exact.stress");
//        TS_ASSERT(comparison.CompareFiles(2e-4));
//
//        MechanicsEventHandler::Headings();
//        MechanicsEventHandler::Report();
    }


    // This is purely to ensure coverage of vtk output in 3d
    void TestVtkCoverage3d()
    {
       QuadraticMesh<3> mesh;
       TrianglesMeshReader<3,3> mesh_reader1("mesh/test/data/3D_Single_tetrahedron_element_quadratic",2,1,false);
       mesh.ConstructFromMeshReader(mesh_reader1);
       NashHunterPoleZeroLaw<3> law;

       std::vector<unsigned> fixed_nodes;
       fixed_nodes.push_back(0);

       SolidMechanicsProblemDefinition<3> problem_defn(mesh);
       problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
       problem_defn.SetZeroDisplacementNodes(fixed_nodes);

       IncompressibleNonlinearElasticitySolver<3> solver(mesh,
                                                         problem_defn,
                                                         "TestIncompressibleNonlinearElasticitySolver");
       solver.Solve();

       // Create .vtu file
       VtkNonlinearElasticitySolutionWriter<3> vtk_writer(solver);
       vtk_writer.Write();

#ifdef CHASTE_VTK
       //Check the VTK file exists
       FileFinder vtk_file("nonlin_elas_functional_data/vtk/solution.vtu", RelativeTo::ChasteTestOutput);
       TS_ASSERT(vtk_file.Exists());
#endif
    }


    // The incompressible elasticity (and fluids) solvers remove the dummy pressure solutions (p=0 at internal nodes)
    // after the solve by linearly interpolating from vertices. Here we test this directly - in particular in 3d as the
    // 3d tests use constant pressure solutions so it is important to check with non-constant pressure
    // at the vertices.
    void TestRemoveDummyPressure()
    {
        // 2d version
        {
            QuadraticMesh<2> mesh(1.0/5, 1.0, 1.0);

            std::vector<unsigned> fixed_nodes;
            fixed_nodes.push_back(0);
            MooneyRivlinMaterialLaw<2> law(1);
            SolidMechanicsProblemDefinition<2> problem_defn(mesh);
            problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
            problem_defn.SetZeroDisplacementNodes(fixed_nodes);

            IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                              problem_defn,
                                                              "");


            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                if (! mesh.GetNode(i)->IsInternal())
                {
                    double x = mesh.GetNode(i)->rGetLocation()[0];
                    double y = mesh.GetNode(i)->rGetLocation()[1];
                    // p defined as a linear function of (x,y) at vertices
                    solver.rGetCurrentSolution()[3*i+2] = 5.6234534 + 2.2432*x + 7.3432*y;
                }
            }

            solver.RemovePressureDummyValuesThroughLinearInterpolation();

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double x = mesh.GetNode(i)->rGetLocation()[0];
                double y = mesh.GetNode(i)->rGetLocation()[1];

                TS_ASSERT_DELTA(solver.rGetCurrentSolution()[3*i  ], 0.0, 1e-8);
                TS_ASSERT_DELTA(solver.rGetCurrentSolution()[3*i+1], 0.0, 1e-8);

                // p should be correct at all nodes, since linearly interpolated a linear function
                TS_ASSERT_DELTA(solver.rGetCurrentSolution()[3*i+2], 5.6234534 + 2.2432*x + 7.3432*y, 1e-8);
            }
        }

        // 3d version
        {
            QuadraticMesh<3> mesh(1.0/5, 1.0, 1.0, 2.0);

            std::vector<unsigned> fixed_nodes;
            fixed_nodes.push_back(0);
            MooneyRivlinMaterialLaw<3> law(1,1);
            SolidMechanicsProblemDefinition<3> problem_defn(mesh);
            problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
            problem_defn.SetZeroDisplacementNodes(fixed_nodes);

            IncompressibleNonlinearElasticitySolver<3> solver(mesh,
                                                              problem_defn,
                                                              "");

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                if (!mesh.GetNode(i)->IsInternal())
                {
                    double x = mesh.GetNode(i)->rGetLocation()[0];
                    double y = mesh.GetNode(i)->rGetLocation()[1];
                    double z = mesh.GetNode(i)->rGetLocation()[2];
                    solver.rGetCurrentSolution()[4*i+3] = 5.6234534 + 2.2432*x + 7.3432*y + 6.24523*z;
                }
            }

            solver.RemovePressureDummyValuesThroughLinearInterpolation();

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double x = mesh.GetNode(i)->rGetLocation()[0];
                double y = mesh.GetNode(i)->rGetLocation()[1];
                double z = mesh.GetNode(i)->rGetLocation()[2];

                TS_ASSERT_DELTA(solver.rGetCurrentSolution()[4*i  ], 0.0, 1e-8);
                TS_ASSERT_DELTA(solver.rGetCurrentSolution()[4*i+1], 0.0, 1e-8);
                TS_ASSERT_DELTA(solver.rGetCurrentSolution()[4*i+2], 0.0, 1e-8);
                TS_ASSERT_DELTA(solver.rGetCurrentSolution()[4*i+3], 5.6234534 + 2.2432*x + 7.3432*y + 6.24523*z, 1e-8);
            }
        }
    }
};

#endif /*TESTINCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_*/
