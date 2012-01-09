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


#ifndef TESTIMPLICITCARDIACMECHANICSSOLVER_HPP_
#define TESTIMPLICITCARDIACMECHANICSSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "ImplicitCardiacMechanicsSolver.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraturePointsGroup.hpp"
#include "NonlinearElasticityTools.hpp"
#include "ReplicatableVector.hpp"


// useful typedef
typedef ImplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<2>,2> IncompressibleImplicitSolver2d;



// helper function - the frobenius norm of a matrix (though any norm would have done).
double MatrixNorm(c_matrix<double,2,2> mat)
{
    return sqrt(mat(0,0)*mat(0,0)+mat(0,1)*mat(0,1)+mat(1,0)*mat(1,0)+mat(1,1)*mat(1,1));
}


class TestImplicitCardiacMechanicsSolver : public CxxTest::TestSuite
{
public:
    void TestCompareJacobians() throw(Exception)
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1.0);
        MooneyRivlinMaterialLaw<2> law(0.02);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        IncompressibleImplicitSolver2d solver(NHS,mesh,problem_defn,"");

        std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints(), 0.0);
        std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);

        for(unsigned i=0; i<calcium_conc.size(); i++)
        {
            calcium_conc[i] = 0.05;
        }

        solver.SetCalciumAndVoltage(calcium_conc, voltages);

        // NOTE: calling CompareJacobians below bypasses calling Solve(t0,t1,dt), hence the
        // time info will not be set. We therefore must explicitly set them here.
        solver.mCurrentTime = 0.0;
        solver.mNextTime = 0.01;
        solver.mOdeTimestep = 0.01;


        ///////////////////////////////////////////////////////////////////
        // compute numerical jacobian and compare with analytic jacobian
        // (about u=0, p=p0)
        ///////////////////////////////////////////////////////////////////
        solver.AssembleSystem(true, true);
        ReplicatableVector rhs_vec(solver.mResidualVector);
        unsigned num_dofs = rhs_vec.GetSize();
        double h = 1e-6;
        int lo, hi;
        MatGetOwnershipRange(solver.mJacobianMatrix, &lo, &hi);

        for(unsigned j=0; j<num_dofs; j++)
        {
            solver.mCurrentSolution.clear();
            solver.FormInitialGuess();
            solver.mCurrentSolution[j] += h;

            solver.AssembleSystem(true, false);

            ReplicatableVector perturbed_rhs( solver.mResidualVector );

            for(unsigned i=0; i<num_dofs; i++)
            {
                if((lo<=(int)i) && ((int)i<hi))
                {
                    double analytic_matrix_val = PetscMatTools::GetElement(solver.mJacobianMatrix,i,j);
                    double numerical_matrix_val = (perturbed_rhs[i] - rhs_vec[i])/h;
                    if((fabs(analytic_matrix_val)>1e-6) && (fabs(numerical_matrix_val)>1e-6))
                    {
                        // relative error
                        TS_ASSERT_DELTA( (analytic_matrix_val-numerical_matrix_val)/analytic_matrix_val, 0.0, 1e-2);
                    }
                    else
                    {
                        // absolute error
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-4);
                    }
                }
            }
        }
        PetscTools::Barrier("TestCompareJacobians");

//// default no longer allowed
//        // coverage - test default material law works ok
//        IncompressibleImplicitSolver2d another_solver(NHS, mesh,"",fixed_nodes);
//        c_matrix<double,2,2> F = zero_matrix<double>(2,2);
//        F(0,0)=F(1,1)=1.1;
//        double pressure = 1;
//        c_matrix<double,2,2> S;
//        another_solver.mMaterialLaws[0]->Compute1stPiolaKirchoffStress(F,pressure,S);
//        TS_ASSERT_DELTA(S(0,0), 1.5805, 1e-3);
    }

    // A test where we specify the 'resting' intracellular calcium concentration
    // for which the active tension should be zero, so should solve in 0 newton
    // iterations
    void TestWithZeroActiveTension() throw(Exception)
    {
        QuadraticMesh<2> mesh(0.125, 1.0, 1.0);
        MooneyRivlinMaterialLaw<2> law(0.02);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        IncompressibleImplicitSolver2d solver(NHS,mesh,problem_defn,"ImplicitCardiacMech/ZeroActiveTension");

        TS_ASSERT_EQUALS(solver.GetTotalNumQuadPoints(), mesh.GetNumElements()*9u);

        // 0.0002 is the initial Ca conc in Lr91
        std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints(), 0.0002);
        std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);
        solver.SetCalciumAndVoltage(calcium_conc, voltages);

        solver.Solve(0,0.1,0.01);

        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 0u);
    }



    // Specifies a non-constant active tension and checks the lambda behaves
    // as it should do. Also has hardcoded tests
    void TestSpecifiedCalciumCompression() throw(Exception)
    {
        // NOTE: test hardcoded for num_elem = 4
        QuadraticMesh<2> mesh(0.25, 1.0, 1.0);
        MooneyRivlinMaterialLaw<2> law(0.02);

        // need to leave the mesh as unfixed as possible
        std::vector<unsigned> fixed_nodes(2);
        fixed_nodes[0] = 0;
        fixed_nodes[1] = 5;

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        IncompressibleImplicitSolver2d solver(NHS,mesh,problem_defn,"ImplicitCardiacMech/SpecifiedCaCompression");
        QuadraturePointsGroup<2> quad_points(mesh, *(solver.GetQuadratureRule()));

        std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints());
        for(unsigned i=0; i<calcium_conc.size(); i++)
        {
            double Y = quad_points.Get(i)(1);
            // 0.0002 is the initial Ca conc in Lr91, 0.001 is the greatest Ca conc
            // value in one of the Lr91 TestIonicModel tests
            calcium_conc[i] = 0.0002 + 0.001*Y;
        }

        std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);
        solver.SetCalciumAndVoltage(calcium_conc, voltages);

        // solve for quite a long time to get some deformation
        solver.Solve(0,10,1);

        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 7u); // hardcoded 7, this check is to make sure the jac is correctly computed

        // have visually checked the answer and seen that it looks ok, so have
        // a hardcoded test here. Node that 24 is the top-right corner node,
        TS_ASSERT_DELTA( solver.rGetDeformedPosition()[24](0), 0.9480, 1e-2); //different results in 3dp with different preconditioners
        TS_ASSERT_DELTA( solver.rGetDeformedPosition()[24](1), 1.0516, 1e-2); //different results in 3dp with different preconditioners

        std::vector<double>& lambda = solver.rGetFibreStretches();

//// removing this test, its a pain to maintain as it requires refitting the cubic
//        // the lambdas should be less than 1 (ie compression), and also
//        // should be near the same for any particular value of Y, ie the
//        // same along any fibre. Lambda should decrease nonlinearly.
//        // Uncomment trace and view in matlab (plot y against lambda)
//        // to observe this - SEE FIGURE ATTACHED TO TICKET #757.
//        // The lambda are constant for given Y if Y>0.1.5 (ie not near fixed nodes)
//        // and a cubic polynomial can be fitted with matlab
//        for(unsigned i=0; i<lambda.size(); i++)
//        {
//            TS_ASSERT_LESS_THAN(lambda[i], 1.0);
//
//            // Get the value of Y for the point
//            double Y = quad_points.Get(i)(1);
//            // Lambda should be near a value obtained by fitting a
//            // cubic polynomial of lambda against Y
//            // Matlab code:
//            //  %% load data from file to variable data, data=(Y,lambda)
//            //  i=find(data(:,1)>0.15)
//            //  data = data = data(i,:)
//            //  c = polyfit(data(:,1),data(:,2),3)
//            //  %% To plot
//            //  yy = 0:0.01:1
//            //  fit = c(1)*yy.^3 + c(2)*yy.^2 + c(3)*yy + c(4);
//            //  plot(data(:,1),data(:,2),'*')
//            //  hold on
//            //  plot(yy,fit,'r')
//            double lam_fit = -0.026920*Y*Y*Y + 0.066128*Y*Y - 0.056929*Y + 0.978174;
//
//            if(Y>0.6)
//            {
//                double error = 0.0005;
//                TS_ASSERT_DELTA(lambda[i], lam_fit, error);
//            }
//            else if(Y>0.15)
//            {
//                double error = 0.0030;
//                TS_ASSERT_DELTA(lambda[i], lam_fit, error);
//            }
//
//            //// don't delete:
//            //std::cout << quad_points.Get(i)(0) << " " << quad_points.Get(i)(1) << " " << lambda[i] << "\n";
//        }

        // hardcoded test
        TS_ASSERT_DELTA(lambda[34], 0.9737, 2e-3);
    }

    // Same as above test but has fibres in Y-direction (and bottom surface fixed - so results are the same),
    // through either setting a constant fibre direction or using a file of different (though in this case
    // all equal) fibre directions for each element.
    // The stretches should be identical to in the above test, the deformation rotated.
    void TestDifferentFibreDirections() throw(Exception)
    {
        for (unsigned run=1; run<=2; run++)
        {
            // NOTE: test hardcoded for num_elem = 4
            QuadraticMesh<2> mesh(0.25, 1.0, 1.0);
            MooneyRivlinMaterialLaw<2> law(0.02);

            // need to leave the mesh as unfixed as possible
            std::vector<unsigned> fixed_nodes(2);
            fixed_nodes[0] = 0;
            fixed_nodes[1] = 1; // was 5 in the above test, {0,1}=small part of X=0 surface, {0,5}=small part of Y=0 surface

            SolidMechanicsProblemDefinition<2> problem_defn(mesh);
            problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
            problem_defn.SetZeroDisplacementNodes(fixed_nodes);

            IncompressibleImplicitSolver2d solver(NHS,mesh,problem_defn,"ImplicitCardiacMech/FibresInYDirection");

            if(run==1)
            {
                c_matrix<double,2,2> non_orthogonal_mat = zero_matrix<double>(2,2);
                non_orthogonal_mat(0,0) = 1.0;
                TS_ASSERT_THROWS_CONTAINS(solver.SetConstantFibreSheetDirections(non_orthogonal_mat), "not orthogonal");

                // ortho matrix = [0 1; 1 0], ie fibres in Y direction
                c_matrix<double,2,2> ortho_matrix = zero_matrix<double>(2,2);
                ortho_matrix(0,1) = 1.0;
                ortho_matrix(1,0) = 1.0;
                solver.SetConstantFibreSheetDirections(ortho_matrix);
            }
            else
            {
                TS_ASSERT_THROWS_CONTAINS(solver.SetVariableFibreSheetDirections("heart/test/data/fibre_tests/badheader_4by4mesh_fibres.ortho", false), "found 32342, expected 32");
                solver.SetVariableFibreSheetDirections("heart/test/data/fibre_tests/4by4mesh_fibres.ortho", false);
            }

            QuadraturePointsGroup<2> quad_points(mesh, *(solver.GetQuadratureRule()));

            std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints());
            for(unsigned i=0; i<calcium_conc.size(); i++)
            {
                double X = quad_points.Get(i)(0);
                // 0.0002 is the initial Ca conc in Lr91, 0.001 is the greatest Ca conc
                // value in one of the Lr91 TestIonicModel tests
                calcium_conc[i] = 0.0002 + 0.001*X;
            }

            std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);
            solver.SetCalciumAndVoltage(calcium_conc, voltages);

            // solve for quite a long time to get some deformation
            solver.Solve(0,10,1);

            TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 7u); // hardcoded 7, this check is to make sure the jac is correctly computed

            // have visually checked the answer and seen that it looks ok, so have
            // a hardcoded test here. Node that 24 is the top-right corner node,
            TS_ASSERT_DELTA( solver.rGetDeformedPosition()[24](1), 0.9429, 1e-2);
            TS_ASSERT_DELTA( solver.rGetDeformedPosition()[24](0), 1.0565, 1e-2);

            std::vector<double>& lambda = solver.rGetFibreStretches();

            //for(unsigned i=0; i<lambda.size(); i++)
            //{
            //    std::cout << quad_points.Get(i)(0) << " " << quad_points.Get(i)(1) << " " << lambda[i] << "\n";
            //}

            // hardcoded test
            TS_ASSERT_DELTA(lambda[34], 0.9693, 1e-3);  // ** different value to previous test - attributing the difference in results to the fact mesh isn't rotation-invariant
        }
    }



    // cover all other contraction model options which are allowed but not been used in a test
    // so far (or in TestExplicitCardiacMechanicsSolver)
    void TestCoverage() throw(Exception)
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1.0);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        IncompressibleImplicitSolver2d impl_solver1(KERCHOFFS2003,mesh,problem_defn,"");
        IncompressibleImplicitSolver2d impl_solver2(NONPHYSIOL3,mesh,problem_defn,"");

        // call with TS_ASSERT_THROWS_CONTAINS with any disallowed contraction models here:
    }


    void TestComputeDeformationGradientAndStretchesEachElement() throw(Exception)
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1.0);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        IncompressibleImplicitSolver2d solver(KERCHOFFS2003,mesh,problem_defn,"");

        // compute the stretches, they should be 1
        std::vector<double> stretches(mesh.GetNumElements());
        std::vector<c_matrix<double,2,2> > deformation_gradients(mesh.GetNumElements());

        // initialise to junk
        for(unsigned i=0; i<stretches.size(); i++)
        {
            stretches[i] = 13482.534578;
            deformation_gradients[i](0,0)
             = deformation_gradients[i](0,1)
               = deformation_gradients[i](1,0)
                 = deformation_gradients[i](1,1) = 7777.777727777777777;
        }


        solver.ComputeDeformationGradientAndStretchInEachElement(deformation_gradients, stretches);
        for(unsigned i=0; i<stretches.size(); i++)
        {
            TS_ASSERT_DELTA(stretches[i], 1.0, 1e-6);
            double err = MatrixNorm(deformation_gradients[i]-identity_matrix<double>(2));
            TS_ASSERT_DELTA(err, 0.0, 1e-10);
        }


        // get the current solution (displacement), and contract in the non-fibre direction
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned j=1;
            solver.mCurrentSolution[2*i+j] = -0.1*mesh.GetNode(i)->rGetLocation()[j];
        }

        // stretches should still be 1, F should be equal to [1,0;0,0.9]
        solver.ComputeDeformationGradientAndStretchInEachElement(deformation_gradients, stretches);

        c_matrix<double,2,2> correct_F = identity_matrix<double>(2);
        correct_F(1,1) = 0.9;
        for(unsigned i=0; i<stretches.size(); i++)
        {
            TS_ASSERT_DELTA(stretches[i], 1.0, 1e-6);
            double err = MatrixNorm(deformation_gradients[i]-correct_F);
            TS_ASSERT_DELTA(err, 0.0, 1e-10);
        }

        // contract in the fibre direction
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned j=0;
            solver.mCurrentSolution[2*i+j] = -0.2*mesh.GetNode(i)->rGetLocation()[j];
        }

        // stretches should now be 0.8, F should be equal to [0.8,0;0,0.9]
        solver.ComputeDeformationGradientAndStretchInEachElement(deformation_gradients, stretches);
        correct_F(0,0) = 0.8;
        for(unsigned i=0; i<stretches.size(); i++)
        {
            TS_ASSERT_DELTA(stretches[i], 0.8, 1e-3);
            double err = MatrixNorm(deformation_gradients[i]-correct_F);
            TS_ASSERT_DELTA(err, 0.0, 1e-10);
        }
    }

    // this test defines a fibre direction for each quadrature point - but the fibres directions are all constant
    // so the results can be checked against the other versions.. See CardiacElectroMechanicsProblem tests for
    // varying fibre directions per quad point.
    void TestDefineFibresPerQuadraturePoint()
    {
        QuadraticMesh<2> mesh(0.25, 1.0, 1.0);
        MooneyRivlinMaterialLaw<2> law(0.02);

        // need to leave the mesh as unfixed as possible
        std::vector<unsigned> fixed_nodes(2);
        fixed_nodes[0] = 0;
        fixed_nodes[1] = 1;

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        IncompressibleImplicitSolver2d solver(NHS, mesh, problem_defn, "ImplicitCardiacMech/FibresInYDirectionDefinePerQuadPoint");

        TS_ASSERT_THROWS_CONTAINS( solver.SetVariableFibreSheetDirections("heart/test/data/fibre_tests/badheader_4by4mesh_fibres.orthoquad", true), "found 45430, expected 288");
        solver.SetVariableFibreSheetDirections("heart/test/data/fibre_tests/4by4mesh_fibres.orthoquad", true);

        QuadraturePointsGroup<2> quad_points(mesh, *(solver.GetQuadratureRule()));
        std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints());
        for(unsigned i=0; i<calcium_conc.size(); i++)
        {
            double X = quad_points.Get(i)(0);
            // 0.0002 is the initial Ca conc in Lr91, 0.001 is the greatest Ca conc
            // value in one of the Lr91 TestIonicModel tests
            calcium_conc[i] = 0.0002 + 0.001*X;
        }

        std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);
        solver.SetCalciumAndVoltage(calcium_conc, voltages);

        // solve for quite a long time to get some deformation
        solver.Solve(0,10,1);

        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 7u); // hardcoded 7, this check is to make sure the jac is correctly computed

        // have visually checked the answer and seen that it looks ok, so have
        // a hardcoded test here. Node that 24 is the top-right corner node,
        TS_ASSERT_DELTA( solver.rGetDeformedPosition()[24](1), 0.9429, 1e-2);
        TS_ASSERT_DELTA( solver.rGetDeformedPosition()[24](0), 1.0565, 1e-2);

        std::vector<double>& lambda = solver.rGetFibreStretches();
        TS_ASSERT_DELTA(lambda[34], 0.9693, 1e-3);
    }



//    //
//    // This test compares the implicit solver with the old dead explicit solver in
//    // dealii. It was written before the explicit chaste solver. It fails in parallel and
//    // can be deleted really, but I'm keeping it here (without ever running it) precisely
//    // because of the strange behaviour in parallel.
//    //
//    // This passes in sequential and fails in parallel (newton failing to converge), and
//    // appears to be because of the first linear system not being solved correctly in parallel.
//    // In sequential 75 iters are reported on the linear system, in parallel 10000 iters (the max)
//    // are reported (but no divergence reported?!) (using PCBJACOBI). So something bad goes
//    // wrong in parallel. In fact, with PCJACOBI, 10000 iters are reported for sequential.
//    // It seems to be related to the fixed nodes: changing the fixed nodes from the 5 set below
//    // to the entire X=0 edge makes it runs in sequential (~60 iters PCBJACOBI) and parallel
//    // (~400 iters PCBJACOBI).
//    //
//    void strangefailingbehaviourinparallelTestCompareWithDeadExplicitSolver() throw(Exception)
//    {
//        // note 8 elements is assumed in the fixed nodes
//        QuadraticMesh<2> mesh(0.125, 1.0, 1.0);
//        MooneyRivlinMaterialLaw<2> law(0.02);
//
//        // fix all nodes on elements containing the origin (as was done in
//        // dealii test)
//        /* THESE CAUSE PROBLEMS.... */
//        std::vector<unsigned> fixed_nodes;
//        fixed_nodes.push_back(0);
//        fixed_nodes.push_back(82);
//        fixed_nodes.push_back(9);
//        fixed_nodes.push_back(1);
//        fixed_nodes.push_back(81);
//
////        /* ....WHEREAS DOING WOULD BE OK */
////        std::vector<unsigned> fixed_nodes
////          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);
//
//        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
//        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
//
//        IncompressibleImplicitSolver2d solver(NHS, mesh, problem_defn, "ImplicitCardiacMech/CompareWithExplicit",&law);
//
//        std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints(), 1); // unrealistically large Ca (but note random material law used)
//        std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);
//        solver.SetCalciumAndVoltage(calcium_conc, voltages);
//        solver.Solve(0,0.01,0.01);
//
//
//        // Visually compared results, they are identical to the dealii results
//        // Hardcoded value for (1,1) node
//        TS_ASSERT_DELTA(solver.rGetDeformedPosition()[80](0),  0.98822 /*dealii*/, 5e-4);
//        TS_ASSERT_DELTA(solver.rGetDeformedPosition()[80](1),  1.01177 /*dealii*/, 3e-4);
//        // Hardcoded value for (0,1) node
//        TS_ASSERT_DELTA(solver.rGetDeformedPosition()[72](0), -0.00465 /*dealii*/, 4e-4);
//        TS_ASSERT_DELTA(solver.rGetDeformedPosition()[72](1),  1.00666 /*dealii*/, 1e-4);
//    }
};

#endif /*TESTIMPLICITCARDIACMECHANICSSOLVER_HPP_*/
