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

#ifndef TESTCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_
#define TESTCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "CompressibleNonlinearElasticitySolver.hpp"
#include "IncompressibleNonlinearElasticitySolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ToyCompressibleMaterialLaw.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "NonlinearElasticityTools.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "CompressibleExponentialLaw.hpp"
#include "NumericFileComparison.hpp"
#include "FileComparison.hpp"

/*
 * All these are for the MyBodyForce and MySurfaceTraction functions below.
 * See TestAgainstExactNonlinearSolution().
 */
static const double C_PARAM = 1.0;
static const double D_PARAM = 0.1;
static const double A_PARAM = 0.1;
static const double Q_PARAM = 0.9;
static const double m = -0.5; // -1/DIM
static const double w1 = C_PARAM*pow(Q_PARAM,2*m);

double ComputeLambda(double X)
{
    return 1+A_PARAM*X;
}

double ComputeI1(double X, double Y)
{
    double lam = ComputeLambda(X);
    return Q_PARAM*Q_PARAM*lam*lam + A_PARAM*A_PARAM*Y*Y/(lam*lam*lam*lam) + 1.0/(lam*lam);
}

double ComputeW3(double X, double Y)
{
    return m*C_PARAM*ComputeI1(X,Y)*pow(Q_PARAM,2*m-2) + D_PARAM*(1-1.0/Q_PARAM);
}

double ComputeDW3dX(double X, double Y)
{
    double lam = ComputeLambda(X);
    return m*C_PARAM*pow(Q_PARAM,2*m-2)*(2*Q_PARAM*Q_PARAM*lam - 4*A_PARAM*A_PARAM*Y*Y*pow(lam,-5) - 2/(lam*lam*lam))*A_PARAM;
}

double ComputeDW3dY(double X, double Y)
{
    double lam = ComputeLambda(X);
    return m*C_PARAM*pow(Q_PARAM,2*m-2)*A_PARAM*A_PARAM*2*Y/(lam*lam*lam*lam);
}

c_matrix<double,2,2> Compute1stPkStress(double X, double Y)
{
    c_matrix<double,2,2> S;

    double w3 = ComputeW3(X,Y);
    double lam = ComputeLambda(X);

    S(0,0) = 2*( w1*Q_PARAM*lam + w3*Q_PARAM/(lam) );
    S(0,1) = -2* w1*Y*A_PARAM/(lam*lam);
    S(1,0) = 2* w3*Y*A_PARAM*Q_PARAM/(lam*lam);
    S(1,1) = 2*( w1/lam + w3*lam*Q_PARAM*Q_PARAM );

    return S;
}

c_vector<double,2> MyBodyForce(c_vector<double,2>& rX, double t)
{
    assert(rX(0)>=0 && rX(0)<=1 && rX(1)>=0 && rX(1)<=1);

    double lam = ComputeLambda(rX(0));
    double a = A_PARAM;
    double q = Q_PARAM;
    double w3 = ComputeW3(rX(0),rX(1));
    double dw3dX = ComputeDW3dX(rX(0),rX(1));
    double dw3dY = ComputeDW3dY(rX(0),rX(1));

    double dS00dX = 2*(w1*q*a - w3*q*a/(lam*lam) + dw3dX*q/(lam));
    double dS01dX = 2*(2*w1*a*a*rX(1)/(lam*lam*lam));
    double dS10dY = 2*(w3*a*q/(lam*lam) + a*rX(1)*q*dw3dY/(lam*lam));
    double dS11dY = 2*lam*dw3dY*q*q;

    c_vector<double,2> body_force;
    body_force(0) = -dS00dX-dS10dY;
    body_force(1) = -dS01dX-dS11dY;
    return body_force;
}

c_vector<double,2> MyTraction(c_vector<double,2>& rX, double t)
{
    c_matrix<double,2,2> S = Compute1stPkStress(rX(0), rX(1));

    c_vector<double,2> traction = zero_vector<double>(2);

    if (fabs(rX(0)-1.0) <= 1e-12) //Right edge
    {
        traction(0) = S(0,0);
        traction(1) = S(0,1);
    }
    else if (fabs(rX(1))  <= 1e-12) //Bottom edge
    {
        traction(0) = -S(1,0);
        traction(1) = -S(1,1);
    }
    else if (fabs(rX(1) - 1.0) <= 1e-12)//Top edge
    {
        traction(0) = S(1,0);
        traction(1) = S(1,1);
    }
    else
    {
        NEVER_REACHED;
    }
    return traction;
}

class TestCompressibleNonlinearElasticitySolver : public CxxTest::TestSuite
{
public:

    /*
     * This is purely for coverage of assembling a 3D system (and also uses
     * alternative, heterogeneous constructor, also for coverage).
     */
    void TestAssembleSystem3D()
    {
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader1("mesh/test/data/3D_Single_tetrahedron_element_quadratic",2,1,false);

        mesh.ConstructFromMeshReader(mesh_reader1);

        ToyCompressibleMaterialLaw<3> law(1.0, 0.0, -1.0);
        std::vector<AbstractMaterialLaw<3>*> laws;
        laws.push_back(&law);

        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);

        SolidMechanicsProblemDefinition<3> problem_defn(mesh);
        problem_defn.SetMaterialLaw(COMPRESSIBLE,laws);

        CompressibleNonlinearElasticitySolver<3> solver(mesh,
                                                        problem_defn,
                                                        "");
        solver.AssembleSystem(true, true);

        TS_ASSERT(PetscMatTools::CheckSymmetry(solver.mrJacobianMatrix));

        // cover exception
        MooneyRivlinMaterialLaw<3> incomp_law(2.0,2.0);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&incomp_law);
        TS_ASSERT_THROWS_THIS( CompressibleNonlinearElasticitySolver<3> another_solver(mesh,problem_defn,""), "SolidMechanicsProblemDefinition object contains incompressible material laws");
    }

    // compare computed Jacobian against a numerically computed
    // Jacobian
    void TestAssembleSystem()
    {
        QuadraticMesh<2> mesh(1.0/2.0, 1.0, 1.0);
        CompressibleExponentialLaw<2> law;

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);

        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        "");
        solver.AssembleSystem(true, true);

        ///////////////////////////////////////////////////////////////////
        // test whether residual vector is currently zero (as
        // current solution should have been initialised to u=0
        ///////////////////////////////////////////////////////////////////
        ReplicatableVector rhs_vec(solver.mResidualVector);
        TS_ASSERT_EQUALS( rhs_vec.GetSize(), 2U*25U );
        for (unsigned i=0; i<rhs_vec.GetSize(); i++)
        {
            TS_ASSERT_DELTA(rhs_vec[i], 0.0, 1e-12);
        }

        ///////////////////////////////////////////////////////////////////
        // compute numerical Jacobian and compare with analytic jacobian
        // (about u=0)
        ///////////////////////////////////////////////////////////////////
        unsigned num_dofs = rhs_vec.GetSize();
        double h = 1e-6;

        int lo, hi;
        MatGetOwnershipRange(solver.mrJacobianMatrix, &lo, &hi);

        for (unsigned j=0; j<num_dofs; j++)
        {
            solver.rGetCurrentSolution().clear();
            solver.rGetCurrentSolution().resize(num_dofs, 0.0);
            solver.rGetCurrentSolution()[j] += h;

            solver.AssembleSystem(true, false);

            ReplicatableVector perturbed_rhs( solver.mResidualVector );

            for (unsigned i=0; i<num_dofs; i++)
            {
                if ((lo<=(int)i) && ((int)i<hi))
                {
                    double analytic_matrix_val = PetscMatTools::GetElement(solver.mrJacobianMatrix,i,j);
                    double numerical_matrix_val = (perturbed_rhs[i] - rhs_vec[i])/h;

                    double average_val = (fabs(analytic_matrix_val)+fabs(numerical_matrix_val))/2.0;

                    if (average_val > 1.0)
                    {
                        // relative error
                        TS_ASSERT_DELTA( (analytic_matrix_val-numerical_matrix_val)/average_val, 0.0, 1e-3);
                    }
                    else
                    {
                        // absolute error - tolerance chosen for finer grid: QuadraticMesh<2> mesh(1.0/4.0, 1.0, 1.0);
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 6e-4);
                    }

//                    double difference =  analytic_matrix_val - numerical_matrix_val;
//                    if (fabs(difference)<1e-6)
//                    {
//                        difference = 0.0;
//                    }
//                    std::cout << difference << " ";
                }
            }
//            std::cout << "\n";
        }
        PetscTools::Barrier();


        //////////////////////////////////////////////////////////
        // compare numerical and analytic jacobians again, this
        // time using a non-zero displacement, u=lambda x, v = mu y
        // (lambda not equal to 1/nu)
        //////////////////////////////////////////////////////////
        double lambda = 1.2;
        double mu = 1.0/1.3;

        solver.rGetCurrentSolution().clear();
        solver.rGetCurrentSolution().resize(num_dofs, 0.0);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            solver.rGetCurrentSolution()[2*i]   = (lambda-1)*mesh.GetNode(i)->rGetLocation()[0];
            solver.rGetCurrentSolution()[2*i+1] = (mu-1)*mesh.GetNode(i)->rGetLocation()[1];
        }

        solver.AssembleSystem(true, true);
        ReplicatableVector rhs_vec2(solver.mResidualVector);

        h=1e-8; // needs to be smaller for this one [COMMENT COPIED FROM INCOMPRESSIBLE VERSION OF THIS TEST]

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
                        TS_ASSERT_DELTA( (analytic_matrix_val-numerical_matrix_val)/analytic_matrix_val, 0.0, 1e-3);
                    }
                    else
                    {
                        // absolute error
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-4);
                    }

//                    double difference =  analytic_matrix_val - numerical_matrix_val;
//
//                    if (fabs(difference)<1e-5)
//                    {
//                        difference = 0.0;
//                    }
//                    std::cout << difference << " ";
                }
            }
//            std::cout << "\n";

            solver.rGetCurrentSolution()[j] -= h;
        }


    }


    // It just tests that nothing happens if zero force and tractions are given
    void TestWithZeroDisplacement()
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        ToyCompressibleMaterialLaw<2> law(1.0, 0.0, -1.0);

        std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        "");

        solver.Solve();
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 0u);

        // Get deformed position
        std::vector<c_vector<double,2> >& r_deformed_position
            = solver.rGetDeformedPosition();

        for (unsigned i=0; i<r_deformed_position.size(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[0], r_deformed_position[i](0), 1e-8);
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[1], r_deformed_position[i](1), 1e-8);
        }
    }

    /**
     * Test against an exact solution.
     *
     * Suppose the deformation is given by  x = (alpha X, beta Y), with a nonlinear Mooney-Rivlin material law
     * W(I1,I2,I3) = c(I1*I3^{-1/2} -3) - d(I3^{1/2} - 1)^2
     *
     * On the unit square we specify displacement boundaries on the X=0 which match this deformation, we assume
     * zero body force, zero traction boundary conditions on the top/bottom surfaces, and fixed traction value, s, on
     * the X=1 surface. Using the above deformation and material law we can compute S by hand, in terms of alpha and
     * beta, and then S11 defines s, and S22=0 gives a relationship between alpha and beta.
     *
     * In this case (writing a for alpha, etc), 0.5*S22 = c/a + a^2 b(-0.5c*(a^2+b^2)/(ab)^3 + d(1 - 1/(ab))
     *
     * which gives a cubic equation to determine beta given alpha, c and d. Let D=d/c, then:
     *  (2Da^3)b^3 + (1-2Da^2)b^2  - a^2 = 0
     *
     * For a given alpha we can use matlab to get the solution (choosing the positive real root):
     *
     *  >> a=0.9; D=0.5;
     *  >> roots([2*D*a*a*a, 1-2*D*a*a, 0.0, -a*a])
     *  ans =
     *    -0.608190204001744 + 0.890314286611269i
     *    -0.608190204001744 - 0.890314286611269i
     *     0.955749406631746
     */
    void TestSolveForSimpleDeformationWithCompMooneyRivlin()
    {
        double c = 2.2;
        double d = 1.1;
        double alpha = 0.9;
        double beta = 0.955749406631746;

        double w1 = c/(alpha*beta); // dW_dI1
        double w3 = -0.5*c*(alpha*alpha+beta*beta)*pow(alpha*beta,-3) + d*(1.0 - 1.0/(alpha*beta)); // dW_dI3

        double traction_value = 2*w1*alpha + 2*w3*alpha*beta*beta;

        unsigned num_elem = 5;

        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);
        CompressibleMooneyRivlinMaterialLaw<2> law(c, d);

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,2> > locations;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if (fabs(mesh.GetNode(i)->rGetLocation()[0]) < 1e-6)
            {
                fixed_nodes.push_back(i);
                c_vector<double,2> new_position;
                new_position(0) = 0;
                new_position(1) = beta*mesh.GetNode(i)->rGetLocation()[1];
                locations.push_back(new_position);
            }
        }

        std::vector<BoundaryElement<1,2>*> boundary_elems;
        std::vector<c_vector<double,2> > tractions;
        c_vector<double,2> traction;
        traction(0) = traction_value;
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
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
        problem_defn.SetFixedNodes(fixed_nodes, locations);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, tractions);


        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        "comp_nonlin_compMR_simple");

        // Coverage
        solver.SetKspAbsoluteTolerance(1e-10);


        solver.SetComputeAverageStressPerElementDuringSolve();


        /////////////////////////////////////////////////////////////////
        // Provide the exact solution as the initial guess and check
        // residual is (nearly) exactly zero (as the solution is in
        // the FEM space, the FEM solution, assuming the nonlinear
        // system was solved exactly, is the exact solution
        /////////////////////////////////////////////////////////////////

        std::vector<double> old_current_soln = solver.rGetCurrentSolution();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double exact_x = alpha*mesh.GetNode(i)->rGetLocation()[0];
            double exact_y = beta*mesh.GetNode(i)->rGetLocation()[1];

            solver.rGetCurrentSolution()[2*i]   = exact_x - mesh.GetNode(i)->rGetLocation()[0];
            solver.rGetCurrentSolution()[2*i+1] = exact_y - mesh.GetNode(i)->rGetLocation()[1];
        }

        solver.Solve();

        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 0u); // initial guess was solution


        // test stresses. The 1st PK stress should satisfy S = [s(0) 0 ; 0 0], where s is the
        // applied traction. This has to be multiplied by F^{-T} to get the 2nd PK stress.
        assert(solver.mAverageStressesPerElement.size()==mesh.GetNumElements());
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            if (mesh.CalculateDesignatedOwnershipOfElement(i))
            {
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,0), traction(0)/alpha, 1e-8);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(1,0), 0.0, 1e-8);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,1), 0.0, 1e-8);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(1,1), 0.0, 1e-8);
            }
        }

        ///////////////////////////////////////////////////////////////////////////
        // Now solve properly
        ///////////////////////////////////////////////////////////////////////////
        solver.rGetCurrentSolution() = old_current_soln;

        solver.Solve();

        std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

        for (unsigned i=0; i<fixed_nodes.size(); i++)
        {
            unsigned index = fixed_nodes[i];
            TS_ASSERT_DELTA(r_solution[index](0), locations[i](0), 1e-8);
            TS_ASSERT_DELTA(r_solution[index](1), locations[i](1), 1e-8);
        }

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double exact_x = alpha*mesh.GetNode(i)->rGetLocation()[0];
            double exact_y = beta*mesh.GetNode(i)->rGetLocation()[1];

            TS_ASSERT_DELTA( r_solution[i](0), exact_x, 1e-5 );
            TS_ASSERT_DELTA( r_solution[i](1), exact_y, 1e-5 );
        }


        // check the stresses (averaged over each quad point). The alpha below is for converting
        // from 1st PK stress (for which we have SN=s => S(0,0) = traction_value) to 2nd PK stress,
        // using T = SF^{-T}
        assert(solver.mAverageStressesPerElement.size()==mesh.GetNumElements());
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            //if (mesh.CalculateDesignatedOwnershipOfElement(i)) ///\todo #2223 This is correct for a distributed mesh
            if (mesh.GetElement(i)->GetOwnership()) // For a shared mesh (not distributed) this is a larger set than the commented line above
            {
                TS_ASSERT_DELTA((solver.GetAverageStressPerElement(i)(0,0)*alpha - traction_value)/traction_value, 0.0, 5e-4);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(1,0), 0.0, 5e-4);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,1), 0.0, 5e-4);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,1), 0.0, 5e-4);
            }
            else
            {
                //The element was ignored by the assembler and thus the stress is unchanged from its initial value
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,0), 0.0, 5e-4);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(1,0), 0.0, 5e-4);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,1), 0.0, 5e-4);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(i)(0,1), 0.0, 5e-4);
            }
        }


        // write the stresses
        solver.WriteCurrentAverageElementStresses("solution");

        // check the written stresses
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        NumericFileComparison comparison(test_output_directory + "/comp_nonlin_compMR_simple/solution.stress", "continuum_mechanics/test/data/nonlin_comp.stress");
        TS_ASSERT(comparison.CompareFiles(2e-4));

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }

    /**
     * Test using a nonlinear material law and for a nonlinear deformation
     *
     * We take the unit square as before, fix one side, use the compressible Mooney-Rivlin
     * law,
     *  W(I1,I2,I3) = c(I1*I3^{-1/2} -3) - d(I3^{1/2} - 1)^2,
     * and want to prescribe tractions on the remaining sides and a body force such that
     *
     *  x = [ q(X + aX^2/2), Y/(1+aX) ]
     *
     * Note that when q=1, this is an incompressible nonlinear deformation (see similar test of
     * incompressible solver). q adds some compressibility
     *
     * Then after a page of algebra, we can derive what the 1st PK stress is, which allows us to
     * determine the required traction and body force.
     *
     * The calculation is written out fully in the FiniteElementImplementations document.
     */
    void TestAgainstExactNonlinearSolution()
    {
        for (unsigned run = 0; run < 2; run++)
        {
            unsigned num_elem = 10;
            QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

            CompressibleMooneyRivlinMaterialLaw<2> law(C_PARAM,D_PARAM);

            std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

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
            problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
            problem_defn.SetZeroDisplacementNodes(fixed_nodes);
            problem_defn.SetBodyForce(MyBodyForce);
            problem_defn.SetTractionBoundaryConditions(boundary_elems, MyTraction);

            if (run==1)
            {
                problem_defn.SetSolveUsingSnes();
            }

            CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                            problem_defn,
                                                            "comp_nonlin_elas_exact_soln");


            solver.Solve();

            std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double X = mesh.GetNode(i)->rGetLocation()[0];
                double Y = mesh.GetNode(i)->rGetLocation()[1];

                double exact_x = Q_PARAM*(X + 0.5*A_PARAM*X*X);
                double exact_y = Y/(1+A_PARAM*X);

                TS_ASSERT_DELTA(r_solution[i](0), exact_x, 1e-4);
                TS_ASSERT_DELTA(r_solution[i](1), exact_y, 1e-4);
            }

            if (run==0)
            {
                // Check that the last matrix was symmetric
                TS_ASSERT(PetscMatTools::CheckSymmetry(solver.mrJacobianMatrix));

                ////////////////////////////////////////////////////////////////////
                // Completely separately, we now test the AssembleOnBoundaryElement
                // method for the situation where a normal pressure on the deformed surface
                // is chosen. This is tested in the incompressible case by solving
                // a full problem and testing against an exact solution, here we just
                // check that AssembleOnBoundaryElement in the compressible solver matches
                // AssembleOnBoundaryElement in the incompressible solve
                ////////////////////////////////////////////////////////////////////
                double pressure = 12.32423;
                problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, pressure);

                c_matrix<double,6,6> a_elem;
                c_vector<double,6> b_elem;
                solver.AssembleOnBoundaryElement(*(boundary_elems[0]), a_elem, b_elem, true, false, 0);

                MooneyRivlinMaterialLaw<2> mooney_rivlin_incompressible(1.0);
                problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&mooney_rivlin_incompressible);

                IncompressibleNonlinearElasticitySolver<2> incompressible_solver(mesh,
                                                                                 problem_defn,
                                                                                 "");

                c_matrix<double,6,6> a_elem_incompressible;
                c_vector<double,6> b_elem_incompressible;

                for (unsigned i=0; i<mesh.GetNumNodes(); i++)
                {
                    // spatial variables
                    for (unsigned j=0; j<2; j++)
                    {
                        incompressible_solver.mCurrentSolution[3*i+j] = solver.mCurrentSolution[2*i+j];
                    }
                    // pressure variable
                    incompressible_solver.mCurrentSolution[3*i+2] = 0.0;
                }

                incompressible_solver.AssembleOnBoundaryElement(*(boundary_elems[0]), a_elem_incompressible, b_elem_incompressible, true, false, 0);

                for (unsigned i=0; i<6; i++)
                {
                    TS_ASSERT_DELTA( b_elem_incompressible(i), b_elem(i), 1e-12 );
                }
            }
        }
    }


    void TestCheckPositiveDefinitenessOfJacobianMatrix()
    {
        unsigned num_elem = 10;

        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);
        CompressibleExponentialLaw<2> law;

        std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        "");

        solver.AssembleSystem(true,true);
        unsigned N = solver.mNumDofs;

        // create some vectors, but make sure they have the same ownership range as
        // the Jacobian matrix
        Vec test_vec;
        VecDuplicate(solver.mResidualVector, &test_vec);
        Vec product_vec;
        VecDuplicate(solver.mResidualVector, &product_vec);
        PetscVecTools::Zero(test_vec);
        PetscVecTools::Zero(product_vec);

        for (unsigned i=0; i<N; i++)
        {
            PetscVecTools::SetElement(test_vec, i, 1.0);

            MatMult(solver.mrJacobianMatrix,test_vec,product_vec);
            double vT_J_v = 0.0;
            VecDot(product_vec, test_vec, &vT_J_v);
            // vT_J_v is the i-th term on the diagonal of the Jacobian matrix (and should be positive)
            TS_ASSERT_LESS_THAN(0.0, vT_J_v);

            PetscVecTools::SetElement(test_vec, i, 0.0);
        }

        PetscTools::Destroy(test_vec);
        PetscTools::Destroy(product_vec);
    }

    // Solve using an exponential material law. Doesn't test against an exact solution, just that check that the
    // solver converges. Doesn't seem very robust.
    void TestSolveForSimpleDeformationWithExponentialLaw()
    {
        unsigned num_elem = 5;

        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);
        CompressibleExponentialLaw<2> law;

        std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        // works ok with g1 = -1,1,2,3. Doesn't newton converge for g1=-2, and gets worse as num_elem increases..
        c_vector<double,2> gravity;
        gravity(0) = 2.0;
        gravity(1) = 0.0;
        problem_defn.SetBodyForce(gravity);

        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        "CompressibleExponentialLawSolve");

        // In parallel, this fails with the current solver / preconditioner combination (cg + bjacobi)
        // although checked that matrix is positive definite. Fails because one of the blocks is
        // not positive def? Get round this by using just jacobi preconditioning.
        PetscTools::SetOption("-pc_type","jacobi");

        solver.Solve();

        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 4u);

        // check node 5 is (1,0)
        assert( fabs(mesh.GetNode(5)->rGetLocation()[0] - 1.0) < 1e-8 );
        assert( fabs(mesh.GetNode(5)->rGetLocation()[1] - 0.0) < 1e-8 );


        std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

        TS_ASSERT_DELTA(r_solution[5](0), 1.0360, 1e-4);
        TS_ASSERT_DELTA(r_solution[5](1), 0.0021, 1e-4);
    }

    /**
     * Same as TestSolveForSimpleDeformationWithCompMooneyRivlin (see comments for this),
     * except the y position of the fixed nodes is left free, i.e. sliding boundary conditions
     * are given
     */
    void TestSolveUsingSlidingBoundaryConditions2d()
    {
        double c = 2.2;
        double d = 1.1;
        double alpha = 0.9;
        double beta = 0.955749406631746;

        double w1 = c/(alpha*beta); // dW_dI1
        double w3 = -0.5*c*(alpha*alpha+beta*beta)*pow(alpha*beta,-3) + d*(1.0 - 1.0/(alpha*beta)); // dW_dI3

        double traction_value = 2*w1*alpha + 2*w3*alpha*beta*beta;

        unsigned num_elem = 5;

        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);
        CompressibleMooneyRivlinMaterialLaw<2> law(c, d);

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
        traction(0) = traction_value;
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
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
        problem_defn.SetFixedNodes(fixed_nodes, locations);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, tractions);


        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        "CompressibleMechanicsSlidingBcs2d");

        // Coverage
        solver.SetKspAbsoluteTolerance(1e-10);

        solver.Solve();

        std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

        for (unsigned i=0; i<fixed_nodes.size(); i++)
        {
            unsigned index = fixed_nodes[i];
            TS_ASSERT_DELTA(r_solution[index](0), 0.0, 1e-8);

            double exact_y = beta*mesh.GetNode(index)->rGetLocation()[1];
            TS_ASSERT_DELTA(r_solution[index](1), exact_y, 1e-5);
        }

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double exact_x = alpha*mesh.GetNode(i)->rGetLocation()[0];
            double exact_y = beta*mesh.GetNode(i)->rGetLocation()[1];

            TS_ASSERT_DELTA( r_solution[i](0), exact_x, 1e-5 );
            TS_ASSERT_DELTA( r_solution[i](1), exact_y, 1e-5 );
        }

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }

    // 3d sliding boundary conditions test
    void TestSolveUsingSlidingBoundaryConditions3d()
    {
        unsigned num_elem = 2;

        QuadraticMesh<3> mesh(1.0/num_elem, 1.0, 1.0, 1.0);
        CompressibleMooneyRivlinMaterialLaw<3> law(1.0,1.0);

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,3> > locations;

        // fix node 0 (located at the origin) fully
        fixed_nodes.push_back(0);
        locations.push_back(zero_vector<double>(3));

        // for the rest of the nodes, if Y=0, set y=0, leave x,z free.
        for (unsigned i=1; i<mesh.GetNumNodes(); i++)
        {
            if (fabs(mesh.GetNode(i)->rGetLocation()[1]) < 1e-6)
            {
                fixed_nodes.push_back(i);
                c_vector<double,3> new_position;
                new_position(0) = SolidMechanicsProblemDefinition<3>::FREE;
                new_position(1) = 0;
                new_position(2) = SolidMechanicsProblemDefinition<3>::FREE;
                locations.push_back(new_position);
            }
        }

        SolidMechanicsProblemDefinition<3> problem_defn(mesh);
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
        problem_defn.SetFixedNodes(fixed_nodes,locations);

        // gravity pushing down against the fixed plane
        c_vector<double,3> gravity;
        gravity(0) = 0.0;
        gravity(1) = -1.0;
        gravity(2) = 0.0;
        problem_defn.SetBodyForce(gravity);

        CompressibleNonlinearElasticitySolver<3> solver(mesh,
                                                        problem_defn,
                                                        "CompressibleMechanicsSlidingBcs3d");


        solver.Solve();

        std::vector<c_vector<double,3> >& r_solution = solver.rGetDeformedPosition();

        // just check the Y=0 nodes still have y=0 but have moved in X and Z (except for node at origin)
        for (unsigned i=1; i<mesh.GetNumNodes(); i++)
        {
            if (fabs(mesh.GetNode(i)->rGetLocation()[1]) < 1e-6)
            {
                TS_ASSERT_DELTA(r_solution[i](1), 0.0, 1e-8);
                TS_ASSERT_DIFFERS(r_solution[i](0), mesh.GetNode(i)->rGetLocation()[0]);
                TS_ASSERT_DIFFERS(r_solution[i](2), mesh.GetNode(i)->rGetLocation()[2]);
            }
        }
    }

    /* HOW_TO_TAG Continuum mechanics
     * Write strain after solve
     */
    void TestWritingStrain()
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1.0);

        CompressibleMooneyRivlinMaterialLaw<2> law(1.0, 1.0);

        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        "TestWritingStrain");

        c_matrix<double,2,2> F;
        solver.GetElementCentroidStrain(DEFORMATION_GRADIENT_F,*(mesh.GetElement(0)), F);

        TS_ASSERT_DELTA(F(0,0), 1.0, 1e-8);
        TS_ASSERT_DELTA(F(0,1), 0.0, 1e-8);
        TS_ASSERT_DELTA(F(1,0), 0.0, 1e-8);
        TS_ASSERT_DELTA(F(1,1), 1.0, 1e-8);

        double alpha = 0.324;

        // Apply a shear to the mesh:
        // move the nodes on Y=1, nodes, 2, 8, 3, across to the right
        solver.mCurrentSolution[4]  = alpha;
        solver.mCurrentSolution[16] = alpha;
        solver.mCurrentSolution[6]  = alpha;
        // move the nodes on Y=0.5, nodes, 5, 6, 7, across to the right by half as much
        solver.mCurrentSolution[10]  = alpha/2.0;
        solver.mCurrentSolution[12] = alpha/2.0;
        solver.mCurrentSolution[14] = alpha/2.0;

        solver.GetElementCentroidStrain(DEFORMATION_GRADIENT_F,*(mesh.GetElement(0)), F);

        TS_ASSERT_DELTA(F(0,0), 1.0, 1e-8);
        TS_ASSERT_DELTA(F(0,1), alpha, 1e-8);
        TS_ASSERT_DELTA(F(1,0), 0.0, 1e-8);
        TS_ASSERT_DELTA(F(1,1), 1.0, 1e-8);

        solver.GetElementCentroidStrain(DEFORMATION_GRADIENT_F,*(mesh.GetElement(1)), F);

        TS_ASSERT_DELTA(F(0,0), 1.0, 1e-8);
        TS_ASSERT_DELTA(F(0,1), alpha, 1e-8);
        TS_ASSERT_DELTA(F(1,0), 0.0, 1e-8);
        TS_ASSERT_DELTA(F(1,1), 1.0, 1e-8);

        c_matrix<double,2,2> C;
        solver.GetElementCentroidStrain(DEFORMATION_TENSOR_C,*(mesh.GetElement(1)), C);
        TS_ASSERT_DELTA(C(0,0), 1.0, 1e-8);
        TS_ASSERT_DELTA(C(0,1), alpha, 1e-8);
        TS_ASSERT_DELTA(C(1,0), alpha, 1e-8);
        TS_ASSERT_DELTA(C(1,1), 1.0+alpha*alpha, 1e-8);

        c_matrix<double,2,2> E;
        solver.GetElementCentroidStrain(LAGRANGE_STRAIN_E,*(mesh.GetElement(1)), E);
        TS_ASSERT_DELTA(E(0,0), 0.0, 1e-8);
        TS_ASSERT_DELTA(E(0,1), alpha/2, 1e-8);
        TS_ASSERT_DELTA(E(1,0), alpha/2, 1e-8);
        TS_ASSERT_DELTA(E(1,1), alpha*alpha/2, 1e-8);

        solver.WriteCurrentStrains(DEFORMATION_GRADIENT_F,"shear_2d",0);

        FileFinder generated_file("TestWritingStrain/shear_2d_0.strain", RelativeTo::ChasteTestOutput);
        FileFinder reference_file("continuum_mechanics/test/data/shear_2d_0.strain", RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated_file,reference_file);
        TS_ASSERT(comparer.CompareFiles());
    }

    void TestWritingStrain3d()
    {
        QuadraticMesh<3> mesh(1.0, 1.0, 1.0, 1.0);

        CompressibleMooneyRivlinMaterialLaw<3> law(1.0, 1.0);

        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);
        SolidMechanicsProblemDefinition<3> problem_defn(mesh);
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        CompressibleNonlinearElasticitySolver<3> solver(mesh,
                                                        problem_defn,
                                                        "");
        double alpha = 1.1;
        double beta = 0.15;
        double gamma = 0.53;

        // Apply a deformation
        //
        // (x,y,z) = (aX, bY-cX, cZ+aY)
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double X = mesh.GetNode(i)->rGetLocation()[0];
            double Y = mesh.GetNode(i)->rGetLocation()[1];
            double Z = mesh.GetNode(i)->rGetLocation()[2];


            solver.mCurrentSolution[3*i]   = alpha*X - X;
            solver.mCurrentSolution[3*i+1] = beta*Y - gamma*X - Y;
            solver.mCurrentSolution[3*i+2] = gamma*Z + alpha*Y - Z;
        }

        // so deformation gradient should be
        // F = [a  0  0]
        //     [-c b  0]
        //     [0  a  c]
        c_matrix<double,3,3> F;

        for (unsigned i=0; i<1 /*mesh.GetNumElements()*/; i++)
        {
            solver.GetElementCentroidStrain(DEFORMATION_GRADIENT_F,*(mesh.GetElement(i)), F);

            TS_ASSERT_DELTA(F(0,0), alpha, 1e-8);
            TS_ASSERT_DELTA(F(0,1), 0.0, 1e-8);
            TS_ASSERT_DELTA(F(0,2), 0.0, 1e-8);
            TS_ASSERT_DELTA(F(1,0), -gamma, 1e-8);
            TS_ASSERT_DELTA(F(1,1), beta, 1e-8);
            TS_ASSERT_DELTA(F(1,2), 0.0, 1e-8);
            TS_ASSERT_DELTA(F(2,0), 0.0, 1e-8);
            TS_ASSERT_DELTA(F(2,1), alpha, 1e-8);
            TS_ASSERT_DELTA(F(2,2), gamma, 1e-8);
        }

        // no output directory given, cover a return statement in WriteCurrentStrains()
        solver.WriteCurrentStrains(DEFORMATION_GRADIENT_F,"wont_be_written",0);
    }


    void TestWritingStress3dAndExceptions()
    {
        QuadraticMesh<3> mesh(1.0, 1.0, 1.0, 1.0);

        CompressibleMooneyRivlinMaterialLaw<3> law(1.0, 1.0);

        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);

        SolidMechanicsProblemDefinition<3> problem_defn(mesh);
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        // cover a return statement where nothing happens if no output directory set
        CompressibleNonlinearElasticitySolver<3> solver_0(mesh,
                                                          problem_defn,
                                                          "");

        solver_0.SetComputeAverageStressPerElementDuringSolve();
        solver_0.WriteCurrentAverageElementStresses("wont_be_written",0);

        CompressibleNonlinearElasticitySolver<3> solver(mesh,
                                                        problem_defn,
                                                        "TestWritingStress");

        // cover exceptions
        TS_ASSERT_THROWS_CONTAINS(solver.GetAverageStressPerElement(0), "Call SetComputeAverageStressPerElementDuringSolve() before solve if calling GetAverageStressesPerElement()");
        TS_ASSERT_THROWS_CONTAINS(solver.WriteCurrentAverageElementStresses("wont_be_written",0), "Call SetComputeAverageStressPerElementDuringSolve() before solve if calling WriteCurrentAverageElementStresses()");

        // test method, but not through calling Solve. Tests where SetComputeAverageStressPerElementDuringSolve() is
        // called and then Solve() is called are TestSolveForSimpleDeformationWithCompMooneyRivlin above (homogeneous
        // stress/strain) and the corresponding homogeneous stress/strain incompressible test.
        //
        solver.SetComputeAverageStressPerElementDuringSolve();

        assert(solver.mAverageStressesPerElement.size()==6u);
        for (unsigned i=0; i<6; i++)
        {
            solver.mAverageStressesPerElement[i] = zero_vector<double>(6);
        }

        c_matrix<double,3,3> T;
        T(0,0) = 0.043;
        T(1,1) = 0.564;
        T(2,2) = 0.243;
        T(0,1) = 0.43;
        T(1,0) = 0.43;
        T(0,2) = 0.03;
        T(2,0) = 0.03;
        T(1,2) = 1.03;
        T(2,1) = 1.03;

        solver.AddStressToAverageStressPerElement(T,0);

        solver.AddStressToAverageStressPerElement(T,1);
        solver.AddStressToAverageStressPerElement(T,1);

        T(2,2) = 10.243;

        solver.AddStressToAverageStressPerElement(T,2);

        for (unsigned i=0; i<6; i++)
        {
            solver.mAverageStressesPerElement[3](i) = i;
        }

        T(2,2) = 0.243;
        double corrrect_T3[3][3] = { {0, 1, 2}, {1, 3, 4}, {2, 4, 5} };


        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(0)(i,j), T(i,j), 1e-12);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(1)(i,j), 2*T(i,j), 1e-12);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(3)(i,j), corrrect_T3[i][j], 1e-12);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(4)(i,j), 0.0, 1e-12);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(5)(i,j), 0.0, 1e-12);
            }
        }

        T(2,2) = 10.243;

        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(2)(i,j), T(i,j), 1e-12);
            }
        }

        solver.WriteCurrentAverageElementStresses("some_stresses_3d",10);

        FileFinder generated_file("TestWritingStress/some_stresses_3d_10.stress", RelativeTo::ChasteTestOutput);
        FileFinder reference_file("continuum_mechanics/test/data/some_stresses_3d_10.stress", RelativeTo::ChasteSourceRoot);
        NumericFileComparison comparison(generated_file, reference_file);
        TS_ASSERT(comparison.CompareFiles());

    }

    // quick 2d test that complements above test
    void TestWritingStress2d()
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1.0);
        CompressibleMooneyRivlinMaterialLaw<2> law(1.0, 1.0);

        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        "");

        solver.SetComputeAverageStressPerElementDuringSolve();

        assert(solver.mAverageStressesPerElement.size()==2u);
        for (unsigned i=0; i<2; i++)
        {
            solver.mAverageStressesPerElement[i] = zero_vector<double>(3);
        }

        c_matrix<double,2,2> T;
        T(0,0) = 0.043;
        T(1,1) = 0.564;
        T(0,1) = 0.43;
        T(1,0) = 0.43;
        solver.AddStressToAverageStressPerElement(T,0);

        for (unsigned i=0; i<2; i++)
        {
            for (unsigned j=0; j<2; j++)
            {
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(0)(i,j), T(i,j), 1e-12);
                TS_ASSERT_DELTA(solver.GetAverageStressPerElement(1)(i,j), 0.0, 1e-12);
            }
        }
    }

    //Covers exceptions that are thrown in AbstractContinuumMechanicsSolver if used with a non-quadratic mesh
    void TestAbstractContinuumMechanicsSolverMeshType()
    {
        TetrahedralMesh<2,2> mesh;
        SolidMechanicsProblemDefinition<2> problem_defn(mesh);

        TS_ASSERT_THROWS_CONTAINS(CompressibleNonlinearElasticitySolver<2> solver(mesh, problem_defn, ""),
                                  "Continuum mechanics solvers require a quadratic mesh");

        TetrahedralMesh<3,3> mesh3d;
        SolidMechanicsProblemDefinition<3> problem_defn3d(mesh3d);

        TS_ASSERT_THROWS_CONTAINS(CompressibleNonlinearElasticitySolver<3> solver3d(mesh3d, problem_defn3d, ""),
                                  "Continuum mechanics solvers require a quadratic mesh");

    }
};

#endif /* TESTCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_ */
