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
#include "TetrahedralMesh.hpp"

// useful typedef
typedef ImplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<2>,2> IncompressibleImplicitSolver2d;
typedef ImplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<3>,3> IncompressibleImplicitSolver3d;

// helper function - the frobenius norm of a matrix (though any norm would have done).
double MatrixNorm(c_matrix<double,2,2> mat)
{
    return sqrt(mat(0,0)*mat(0,0)+mat(0,1)*mat(0,1)+mat(1,0)*mat(1,0)+mat(1,1)*mat(1,1));
}


class TestImplicitCardiacMechanicsSolver : public CxxTest::TestSuite
{
public:
    void TestCompareJacobians()
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1.0);
        MooneyRivlinMaterialLaw<2> law(0.02);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        for (unsigned run=0; run<2; run++)
        {
            // Two runs - one with cross fibre tension applied, one without.

            ElectroMechanicsProblemDefinition<2> problem_defn(mesh);
            problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
            problem_defn.SetZeroDisplacementNodes(fixed_nodes);
            problem_defn.SetContractionModel(NHS,0.01);
            problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass
            if (run==1)
            {
                bool apply_cross_fibre_tension = true;
                double cross_fibre_tension_fraction = 0.25;
                problem_defn.SetApplyIsotropicCrossFibreTension(apply_cross_fibre_tension, cross_fibre_tension_fraction);
            }

            IncompressibleImplicitSolver2d solver(mesh,problem_defn,"");

            //The following lines are not relevant to this test but need to be there
            TetrahedralMesh<2,2>* p_fine_mesh = new TetrahedralMesh<2,2>();//unused in this test
            p_fine_mesh->ConstructRegularSlabMesh(1.0,1.0,1.0);
            TetrahedralMesh<2,2>* p_coarse_mesh = new TetrahedralMesh<2,2>();//unused in this test
            p_coarse_mesh->ConstructRegularSlabMesh(1.0,1.0,1.0);
            FineCoarseMeshPair<2>* p_pair = new FineCoarseMeshPair<2>(*p_fine_mesh, *p_coarse_mesh);//also unused in this test
            p_pair->SetUpBoxesOnFineMesh();
            p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(solver.GetQuadratureRule()), false);
            p_pair->DeleteFineBoxCollection();
            solver.SetFineCoarseMeshPair(p_pair);
            ///////////////////////////////////////////////////////////////////////////

            solver.Initialise();
            std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints(), 0.0);
            std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);

            for (unsigned i=0; i<calcium_conc.size(); i++)
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
            MatGetOwnershipRange(solver.mrJacobianMatrix, &lo, &hi);

            for (unsigned j=0; j<num_dofs; j++)
            {
                solver.mCurrentSolution.clear();
                solver.FormInitialGuess();
                solver.mCurrentSolution[j] += h;

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
                            TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-4);
                        }
                    }
                }
            }

            //in need of deletion even if all these 3 have no influence at all on this test
            delete p_fine_mesh;
            delete p_coarse_mesh;
            delete p_pair;
            PetscTools::Barrier("TestCompareJacobians");
        }
    }

    // A test where we specify the 'resting' intracellular calcium concentration
    // for which the active tension should be zero, so should solve in 0 newton
    // iterations
    void TestWithZeroActiveTension()
    {
        QuadraticMesh<2> mesh(0.125, 1.0, 1.0);
        MooneyRivlinMaterialLaw<2> law(0.02);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetContractionModel(NHS,0.01);
        problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass

        IncompressibleImplicitSolver2d solver(mesh,problem_defn,"ImplicitCardiacMech/ZeroActiveTension");

        //The following lines are not relevant to this test but need to be there
        TetrahedralMesh<2,2>* p_fine_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_fine_mesh->ConstructRegularSlabMesh(0.125, 1.0, 1.0);
        TetrahedralMesh<2,2>* p_coarse_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_coarse_mesh->ConstructRegularSlabMesh(0.125, 1.0, 1.0);
        FineCoarseMeshPair<2>* p_pair = new FineCoarseMeshPair<2>(*p_fine_mesh, *p_coarse_mesh);//also unused in this test
        p_pair->SetUpBoxesOnFineMesh();
        p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(solver.GetQuadratureRule()), false);
        p_pair->DeleteFineBoxCollection();
        solver.SetFineCoarseMeshPair(p_pair);
        ///////////////////////////////////////////////////////////////////////////
        solver.Initialise();

        TS_ASSERT_EQUALS(solver.GetTotalNumQuadPoints(), mesh.GetNumElements()*6u);

        // 0.0002 is the initial Ca conc in Lr91
        std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints(), 0.0002);
        std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);
        solver.SetCalciumAndVoltage(calcium_conc, voltages);

        solver.Solve(0,0.1,0.01);

        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 0u);

        //in need of deletion even if all these 3 have no influence at all on this test
        delete p_fine_mesh;
        delete p_coarse_mesh;
        delete p_pair;
    }

    // Specifies a non-constant active tension and checks the lambda behaves
    // as it should do. Also has hardcoded tests
    void TestSpecifiedCalciumCompression()
    {
        // NOTE: test hardcoded for num_elem = 4
        QuadraticMesh<2> mesh(0.25, 1.0, 1.0);
        MooneyRivlinMaterialLaw<2> law(0.02);

        // need to leave the mesh as unfixed as possible
        std::vector<unsigned> fixed_nodes(2);
        fixed_nodes[0] = 0;
        fixed_nodes[1] = 5;

        ElectroMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetContractionModel(NHS,0.01);
        problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass

        IncompressibleImplicitSolver2d solver(mesh,problem_defn,"ImplicitCardiacMech/SpecifiedCaCompression");
        QuadraturePointsGroup<2> quad_points(mesh, *(solver.GetQuadratureRule()));

        //The following lines are not relevant to this test but need to be there
        TetrahedralMesh<2,2>* p_fine_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_fine_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0);
        TetrahedralMesh<2,2>* p_coarse_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_coarse_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0);
        FineCoarseMeshPair<2>* p_pair = new FineCoarseMeshPair<2>(*p_fine_mesh, *p_coarse_mesh);//also unused in this test
        p_pair->SetUpBoxesOnFineMesh();
        p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(solver.GetQuadratureRule()), false);
        p_pair->DeleteFineBoxCollection();
        solver.SetFineCoarseMeshPair(p_pair);
        ///////////////////////////////////////////////////////////////////////////

        solver.Initialise();
        std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints());
        for (unsigned i=0; i<calcium_conc.size(); i++)
        {
            double Y = quad_points.rGet(i)(1);
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

//// removing this test, its a pain to maintain as it requires refitting the cubic
//        // the lambdas should be less than 1 (ie compression), and also
//        // should be near the same for any particular value of Y, ie the
//        // same along any fibre. Lambda should decrease nonlinearly.
//        // Uncomment trace and view in matlab (plot y against lambda)
//        // to observe this - SEE FIGURE ATTACHED TO TICKET #757.
//        // The lambda are constant for given Y if Y>0.1.5 (ie not near fixed nodes)
//        // and a cubic polynomial can be fitted with matlab
//        for (unsigned i=0; i<lambda.size(); i++)
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
//            if (Y>0.6)
//            {
//                double error = 0.0005;
//                TS_ASSERT_DELTA(lambda[i], lam_fit, error);
//            }
//            else if (Y>0.15)
//            {
//                double error = 0.0030;
//                TS_ASSERT_DELTA(lambda[i], lam_fit, error);
//            }
//
//            //// don't delete:
//            std::cout << quad_points.Get(i)(0) << " " << quad_points.Get(i)(1) << " " << lambda[i] << "\n";
//        }

        // hardcoded test
        // Was quad point 34 = 3*9 + 7 (quad 7 in element 3) when there were 9 quads per element
        // Investigate quad point 21 = 3*6 + 3 (quad 3 in element 3)

        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(quad_points.rGet(21)), 3u);

        std::map<unsigned,DataAtQuadraturePoint>::iterator iter = solver.rGetQuadPointToDataAtQuadPointMap().find(19);
        if (iter != solver.rGetQuadPointToDataAtQuadPointMap().end()) //ie because some processes won't own this in parallel
        {
            TS_ASSERT_DELTA(iter->second.Stretch, 0.9737, 2e-3);
        }

        //in need of deletion even if all these 3 have no influence at all on this test
        delete p_fine_mesh;
        delete p_coarse_mesh;
        delete p_pair;
    }

    // Same as above test but has fibres in Y-direction (and bottom surface fixed - so results are the same),
    // through either setting a constant fibre direction or using a file of different (though in this case
    // all equal) fibre directions for each element.
    // The stretches should be identical to in the above test, the deformation rotated.
    void TestDifferentFibreDirections()
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

            ElectroMechanicsProblemDefinition<2> problem_defn(mesh);
            problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
            problem_defn.SetZeroDisplacementNodes(fixed_nodes);
            problem_defn.SetContractionModel(NHS,0.01);
            problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass

            IncompressibleImplicitSolver2d solver(mesh,problem_defn,"ImplicitCardiacMech/FibresInYDirection");

            //The following lines are not relevant to this test but need to be there
            TetrahedralMesh<2,2>* p_fine_mesh = new TetrahedralMesh<2,2>();//unused in this test
            p_fine_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0);
            TetrahedralMesh<2,2>* p_coarse_mesh = new TetrahedralMesh<2,2>();//unused in this test
            p_coarse_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0);
            FineCoarseMeshPair<2>* p_pair = new FineCoarseMeshPair<2>(*p_fine_mesh, *p_coarse_mesh);//also unused in this test
            p_pair->SetUpBoxesOnFineMesh();
            p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(solver.GetQuadratureRule()), false);
            p_pair->DeleteFineBoxCollection();
            solver.SetFineCoarseMeshPair(p_pair);
            ///////////////////////////////////////////////////////////////////////////
            solver.Initialise();

            if (run==1)
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
                FileFinder bad_file("heart/test/data/fibre_tests/badheader_4by4mesh_fibres.ortho", RelativeTo::ChasteSourceRoot);
                TS_ASSERT_THROWS_CONTAINS(solver.SetVariableFibreSheetDirections(bad_file, false),
                                          "found 32342, expected 32");
                FileFinder good_file("heart/test/data/fibre_tests/4by4mesh_fibres.ortho", RelativeTo::ChasteSourceRoot);
                solver.SetVariableFibreSheetDirections(good_file, false);
            }

            QuadraturePointsGroup<2> quad_points(mesh, *(solver.GetQuadratureRule()));

            std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints());
            for (unsigned i=0; i<calcium_conc.size(); i++)
            {
                double X = quad_points.rGet(i)(0);
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


            // Was quad point 34 = 3*9 + 7 (quad 7 in element 3) when there were 9 quads per element
            // Investigate quad point 19 = 3*6 + 1 (quad 3 in element 3)
            std::map<unsigned,DataAtQuadraturePoint>::iterator iter = solver.rGetQuadPointToDataAtQuadPointMap().find(19);
            if (iter != solver.rGetQuadPointToDataAtQuadPointMap().end()) //ie because some processes won't own this in parallel
            {
                TS_ASSERT_DELTA(iter->second.Stretch, 0.9682, 1e-3);  // ** different value to previous test - attributing the difference in results to the fact mesh isn't rotation-invariant
            }

            //in need of deletion even if all these 3 have no influence at all on this test
            delete p_fine_mesh;
            delete p_coarse_mesh;
            delete p_pair;
        }
    }


    // cover all other contraction model options which are allowed but not been used in a test
    // so far (or in TestExplicitCardiacMechanicsSolver)
    void TestCoverage()
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1.0);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetContractionModel(KERCHOFFS2003,0.01);
        problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass

        //The following lines are not relevant to this test but need to be there
        TetrahedralMesh<2,2>* p_fine_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_fine_mesh->ConstructRegularSlabMesh(1.0, 1.0, 1.0);
        TetrahedralMesh<2,2>* p_coarse_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_coarse_mesh->ConstructRegularSlabMesh(1.0, 1.0, 1.0);
        FineCoarseMeshPair<2>* p_pair = new FineCoarseMeshPair<2>(*p_fine_mesh, *p_coarse_mesh);//also unused in this test
        p_pair->SetUpBoxesOnFineMesh();

        //////////////////////////////////////////////////////////////////
        IncompressibleImplicitSolver2d impl_solver1(mesh,problem_defn,"");
        p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(impl_solver1.GetQuadratureRule()), false);
        p_pair->DeleteFineBoxCollection();
        impl_solver1.SetFineCoarseMeshPair(p_pair);
        impl_solver1.Initialise();

        // Prevent an assertion being thrown about setting the cell factory more than once.
        // (just for testing)
        delete problem_defn.mpContractionCellFactory;
        problem_defn.mpContractionCellFactory=NULL;

        problem_defn.SetContractionModel(NONPHYSIOL3,0.01);
        IncompressibleImplicitSolver2d impl_solver2(mesh,problem_defn,"");
        impl_solver2.SetFineCoarseMeshPair(p_pair);
        impl_solver2.Initialise();

        // call with TS_ASSERT_THROWS_CONTAINS with any disallowed contraction models here:

        //in need of deletion even if all these 3 have no influence at all on this test
        delete p_fine_mesh;
        delete p_coarse_mesh;
        delete p_pair;
    }


    void TestComputeDeformationGradientAndStretchesEachElement()
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1.0);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetContractionModel(KERCHOFFS2003,0.01);
        problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass

        IncompressibleImplicitSolver2d solver(mesh,problem_defn,"");

        //The following lines are not relevant to this test but need to be there
        TetrahedralMesh<2,2>* p_fine_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_fine_mesh->ConstructRegularSlabMesh(1.0, 1.0, 1.0);
        TetrahedralMesh<2,2>* p_coarse_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_coarse_mesh->ConstructRegularSlabMesh(1.0, 1.0, 1.0);
        FineCoarseMeshPair<2>* p_pair = new FineCoarseMeshPair<2>(*p_fine_mesh, *p_coarse_mesh);//also unused in this test
        p_pair->SetUpBoxesOnFineMesh();
        p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(solver.GetQuadratureRule()), false);
        p_pair->DeleteFineBoxCollection();
        solver.SetFineCoarseMeshPair(p_pair);
        ///////////////////////////////////////////////////////////////////////////
        solver.Initialise();

        // compute the stretches, they should be 1
        std::vector<double> stretches(mesh.GetNumElements());
        std::vector<c_matrix<double,2,2> > deformation_gradients(mesh.GetNumElements());

        // initialise to junk
        for (unsigned i=0; i<stretches.size(); i++)
        {
            stretches[i] = 13482.534578;
            deformation_gradients[i](0,0)
             = deformation_gradients[i](0,1)
               = deformation_gradients[i](1,0)
                 = deformation_gradients[i](1,1) = 7777.777727777777777;
        }


        solver.ComputeDeformationGradientAndStretchInEachElement(deformation_gradients, stretches);
        for (unsigned i=0; i<stretches.size(); i++)
        {
            TS_ASSERT_DELTA(stretches[i], 1.0, 1e-6);
            double err = MatrixNorm(deformation_gradients[i]-identity_matrix<double>(2));
            TS_ASSERT_DELTA(err, 0.0, 1e-10);
        }


        // get the current solution (displacement), and contract in the non-fibre direction
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned j=1;
            // note: the 3 here is DIM+1, the problem dim for incompressible problems
            solver.mCurrentSolution[3*i+j] = -0.1*mesh.GetNode(i)->rGetLocation()[j];
        }

        // stretches should still be 1, F should be equal to [1,0;0,0.9]
        solver.ComputeDeformationGradientAndStretchInEachElement(deformation_gradients, stretches);

        c_matrix<double,2,2> correct_F = identity_matrix<double>(2);
        correct_F(1,1) = 0.9;        //in need of deletion even if all these 3 have no influence at all on this test

        for (unsigned i=0; i<stretches.size(); i++)
        {
            TS_ASSERT_DELTA(stretches[i], 1.0, 1e-6);
            double err = MatrixNorm(deformation_gradients[i]-correct_F);
            TS_ASSERT_DELTA(err, 0.0, 1e-10);
        }

        // contract in the fibre direction
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned j=0;
            // note: the 3 here is DIM+1, the problem dim for incompressible problems
            solver.mCurrentSolution[3*i+j] = -0.2*mesh.GetNode(i)->rGetLocation()[j];
        }

        // stretches should now be 0.8, F should be equal to [0.8,0;0,0.9]
        solver.ComputeDeformationGradientAndStretchInEachElement(deformation_gradients, stretches);
        correct_F(0,0) = 0.8;
        for (unsigned i=0; i<stretches.size(); i++)
        {
            TS_ASSERT_DELTA(stretches[i], 0.8, 1e-3);
            double err = MatrixNorm(deformation_gradients[i]-correct_F);
            TS_ASSERT_DELTA(err, 0.0, 1e-10);
        }

        //in need of deletion even if all these 3 have no influence at all on this test
        delete p_fine_mesh;
        delete p_coarse_mesh;
        delete p_pair;
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

        ElectroMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetContractionModel(NHS,0.01);
        problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass

        IncompressibleImplicitSolver2d solver(mesh, problem_defn, "ImplicitCardiacMech/FibresInYDirectionDefinePerQuadPoint");

        //The following lines are not relevant to this test but need to be there
        TetrahedralMesh<2,2>* p_fine_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_fine_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0);
        TetrahedralMesh<2,2>* p_coarse_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_coarse_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0);
        FineCoarseMeshPair<2>* p_pair = new FineCoarseMeshPair<2>(*p_fine_mesh, *p_coarse_mesh);//also unused in this test
        p_pair->SetUpBoxesOnFineMesh();
        p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(solver.GetQuadratureRule()), false);
        p_pair->DeleteFineBoxCollection();
        solver.SetFineCoarseMeshPair(p_pair);
        ///////////////////////////////////////////////////////////////////////////
        solver.Initialise();

        FileFinder bad_file("heart/test/data/fibre_tests/badheader_4by4mesh_fibres.orthoquad", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_CONTAINS(solver.SetVariableFibreSheetDirections(bad_file, true),
                                  "found 45430, expected 192");

        FileFinder good_file("heart/test/data/fibre_tests/4by4mesh_fibres.orthoquad", RelativeTo::ChasteSourceRoot);
        solver.SetVariableFibreSheetDirections(good_file, true);

        QuadraturePointsGroup<2> quad_points(mesh, *(solver.GetQuadratureRule()));
        std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints());
        for (unsigned i=0; i<calcium_conc.size(); i++)
        {
            double X = quad_points.rGet(i)(0);
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

        std::map<unsigned,DataAtQuadraturePoint>::iterator iter = solver.rGetQuadPointToDataAtQuadPointMap().find(19);
        if (iter != solver.rGetQuadPointToDataAtQuadPointMap().end()) //ie because some processes won't own this in parallel
        {
            TS_ASSERT_DELTA(iter->second.Stretch, 0.9682, 1e-3);
        }

        //in need of deletion even if all these 3 have no influence at all on this test
        delete p_fine_mesh;
        delete p_coarse_mesh;
        delete p_pair;
    }

    void TestCrossFibreTensionWithSimpleContractionModel()
    {
        QuadraticMesh<2> mesh(0.25, 1.0, 1.0);
        MooneyRivlinMaterialLaw<2> law(1);

        std::vector<unsigned> fixed_nodes
        = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetContractionModel(NHS,0.01);
        problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass


        //Cross fibre tension fractions to be tested
        c_vector<double, 4> tension_fractions;
        tension_fractions[0] = 0.0;
        tension_fractions[1] = 0.1;
        tension_fractions[2] = 0.5;
        tension_fractions[3] = 1.0;

        //Expected resulting deformed location of Node 4.
        c_vector<double, 4> x;
        c_vector<double, 4> y;
        x[0] = 0.9984;
        x[1] = 0.9985;
        x[2] = 0.9990;
        x[3] = 1.0; //Note, for cross_tension == 1.0 there should be no deformation of the tissue square (tensions balance)
        y[0] = -0.0013;
        y[1] = -0.0013;
        y[2] = -0.0008;
        y[3] = 0.0;

        for (unsigned i=0; i < tension_fractions.size();i++)
        {
            problem_defn.SetApplyIsotropicCrossFibreTension(true,tension_fractions[i]);

            TS_ASSERT_EQUALS(problem_defn.GetApplyCrossFibreTension(), true);
            TS_ASSERT_DELTA(problem_defn.GetSheetTensionFraction(),tension_fractions[i], 1e-6);
            TS_ASSERT_DELTA(problem_defn.GetSheetNormalTensionFraction(),tension_fractions[i], 1e-6);

            // NONPHYSIOL1 => NonphysiologicalContractionModel 1
            IncompressibleImplicitSolver2d solver(mesh,problem_defn,"TestImplicitCardiacMech");

            // The following lines are not relevant to this test but need to be there
            // as the solver is expecting an electrics node to be paired up with each mechanics node.
            TetrahedralMesh<2,2>* p_fine_mesh = new TetrahedralMesh<2,2>();//electrics ignored in this test
            p_fine_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0);
            FineCoarseMeshPair<2>* p_pair = new FineCoarseMeshPair<2>(*p_fine_mesh, mesh);
            p_pair->SetUpBoxesOnFineMesh();
            p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(solver.GetQuadratureRule()), false);
            p_pair->DeleteFineBoxCollection();
            solver.SetFineCoarseMeshPair(p_pair);
            ///////////////////////////////////////////////////////////////////////////
            solver.Initialise();

            // Set up a voltage and calcium level at each quadrature point.
            QuadraturePointsGroup<2> quad_points(mesh, *(solver.GetQuadratureRule()));
            std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints());
            for (unsigned j=0; j<calcium_conc.size(); j++)
            {
                double X = quad_points.rGet(j)(0);
                // 0.0002 is the initial Ca conc in Lr91, 0.001 is the greatest Ca conc
                // value in one of the Lr91 TestIonicModel tests
                calcium_conc[j] = 0.0005 + 0.001*X;
            }
            std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);
            solver.SetCalciumAndVoltage(calcium_conc, voltages);

            // solve UP TO t=0. So Ta(lam_n,t_{n+1})=5*sin(0)=0, ie no deformation
            solver.Solve(-0.01,0.0,0.01);
            TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(),0u);

            //solver.Solve(0.24,0.25,0.01);
            solver.Solve(0.0,2,0.01);
//// These fail due to changes from #2185 but no point fixing as these numbers will change later anyway (#2180)
            TS_ASSERT_DELTA(solver.rGetDeformedPosition()[4](0),  x[i], 1e-3);
            TS_ASSERT_DELTA(solver.rGetDeformedPosition()[4](1), y[i], 1e-3);

            // clean up memory
            delete p_fine_mesh;
            delete p_pair;
        }
    }


    /**
     * Here we apply the same tension in the fibre direction,
     * and in both sheet and sheet-normal (cross fibre) directions.
     *
     * Therefore nothing should happen!
     */
    void TestIsotropicCrossFibreTensions()
    {
        /*
         * Expected resulting deformed location of Nodes 4, 24, 104, 124:
         * 4: 1, 0, 0
         * 24: 1, 1, 0
         * 104: 1, 0, 1
         * 124: 1, 1, 1
         */
        c_vector<double, 4> x;
        c_vector<double, 4> y;
        c_vector<double, 4> z;
        x[0] = 1;
        x[1] = 1;
        x[2] = 1;
        x[3] = 1;
        y[0] = 0;
        y[1] = 1;
        y[2] = 0;
        y[3] = 1;
        z[0] = 0;
        z[1] = 0;
        z[2] = 1;
        z[3] = 1;

        QuadraticMesh<3> mesh(0.25, 1.0, 1.0, 1.0);
        MooneyRivlinMaterialLaw<3> law(1,1);

        std::vector<unsigned> fixed_nodes
        = NonlinearElasticityTools<3>::GetNodesByComponentValue(mesh,0,0.0);

        ElectroMechanicsProblemDefinition<3> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetContractionModel(NONPHYSIOL1,0.01);
        problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass
        double tension_fraction=1;
        problem_defn.SetApplyIsotropicCrossFibreTension(true,tension_fraction);

        // NONPHYSIOL1 => NonphysiologicalContractionModel 1
        IncompressibleImplicitSolver3d solver(mesh,problem_defn,"TestIsotropicCrossFibreImplicit");

        // The following lines are not relevant to this test but need to be there
        // as the solver is expecting an electrics node to be paired up with each mechanics node.
        TetrahedralMesh<3,3>* p_fine_mesh = new TetrahedralMesh<3,3>();//electrics ignored in this test
        p_fine_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0, 1.0);
        FineCoarseMeshPair<3>* p_pair = new FineCoarseMeshPair<3>(*p_fine_mesh, mesh);
        p_pair->SetUpBoxesOnFineMesh();
        p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(solver.GetQuadratureRule()), false);
        p_pair->DeleteFineBoxCollection();
        solver.SetFineCoarseMeshPair(p_pair);
        ///////////////////////////////////////////////////////////////////////////
        solver.Initialise();

        // coverage
        QuadraturePointsGroup<3> quad_points(mesh, *(solver.GetQuadratureRule()));

        std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints(), 0.0);
        std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);

        solver.SetCalciumAndVoltage(calcium_conc, voltages);

        // solve UP TO t=0. So Ta(lam_n,t_{n+1})=5*sin(0)=0, ie no deformation
        solver.Solve(-0.01,0.0,0.01);
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(),0u);

        solver.Solve(0,0.01,0.01);

        std::vector<unsigned> nodes;
        nodes.push_back(4);
        nodes.push_back(24);
        nodes.push_back(104);
        nodes.push_back(124);

        for (unsigned node=0; node<4; node++)
        {
            std::cout << "Node: " << nodes[node] << "\n";
            TS_ASSERT_DELTA(solver.rGetDeformedPosition()[nodes[node]](0), x[node], 1e-4);
            TS_ASSERT_DELTA(solver.rGetDeformedPosition()[nodes[node]](1), y[node], 1e-4);
            TS_ASSERT_DELTA(solver.rGetDeformedPosition()[nodes[node]](2), z[node], 1e-4);
        }
        //in need of deletion even if all these 3 have no influence at all on this test
        delete p_fine_mesh;
        delete p_pair;
    }

    /**
     * This time we will make x (fibres) and z (sheet-normal) contract,
     * and y will not contract (so nodes will expand out for y, shrink in for x,z).
     */
    void TestAnisotropicCrossFibreTensions()
    {
        /*
         * Expected resulting deformed location of Nodes 4, 24, 104, 124:
         * 4: 1, 0, 0
         * 24: 1, 1, 0
         * 104: 1, 0, 1
         * 124: 1, 1, 1
         */
        c_vector<double, 4> x;
        c_vector<double, 4> y;
        c_vector<double, 4> z;
        x[0] = 0.9905;
        x[1] = 0.9908;
        x[2] = 0.9904;
        x[3] = 0.9907;
        y[0] = -0.0113;
        y[1] = 1.0105;
        y[2] = -0.0113;
        y[3] = 1.0105;
        z[0] = 0.0056;
        z[1] = 0.0056;
        z[2] = 0.9946;
        z[3] = 0.9946;

        QuadraticMesh<3> mesh(0.25, 1.0, 1.0, 1.0);
        MooneyRivlinMaterialLaw<3> law(1,1);

        std::vector<unsigned> fixed_nodes
        = NonlinearElasticityTools<3>::GetNodesByComponentValue(mesh,0,0.0);

        ElectroMechanicsProblemDefinition<3> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetContractionModel(NONPHYSIOL1,0.01);
        problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass
        double sheet_tension_fraction=0;
        double sheet_normal_tension_fraction=1;
        problem_defn.SetApplyAnisotropicCrossFibreTension(true,sheet_tension_fraction, sheet_normal_tension_fraction);

        // NONPHYSIOL1 => NonphysiologicalContractionModel 1
        IncompressibleImplicitSolver3d solver(mesh,problem_defn,"TestAnisotropicCrossFibreImplicit");

        // The following lines are not relevant to this test but need to be there
        // as the solver is expecting an electrics node to be paired up with each mechanics node.
        TetrahedralMesh<3,3>* p_fine_mesh = new TetrahedralMesh<3,3>();//electrics ignored in this test
        p_fine_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0, 1.0);
        FineCoarseMeshPair<3>* p_pair = new FineCoarseMeshPair<3>(*p_fine_mesh, mesh);
        p_pair->SetUpBoxesOnFineMesh();
        p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(solver.GetQuadratureRule()), false);
        p_pair->DeleteFineBoxCollection();
        solver.SetFineCoarseMeshPair(p_pair);
        ///////////////////////////////////////////////////////////////////////////
        solver.Initialise();

        // coverage
        QuadraturePointsGroup<3> quad_points(mesh, *(solver.GetQuadratureRule()));

        std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints(), 0.0);
        std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);

        solver.SetCalciumAndVoltage(calcium_conc, voltages);

        // solve UP TO t=0. So Ta(lam_n,t_{n+1})=5*sin(0)=0, ie no deformation
        solver.Solve(-0.01,0.0,0.01);
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(),0u);

        solver.Solve(0.24,0.25,0.01);

        std::vector<unsigned> nodes;
        nodes.push_back(4);
        nodes.push_back(24);
        nodes.push_back(104);
        nodes.push_back(124);

        for (unsigned node=0; node<4; node++)
        {
            std::cout << "Node: " << nodes[node] << "\n";
            TS_ASSERT_DELTA(solver.rGetDeformedPosition()[nodes[node]](0), x[node], 1e-3);
            TS_ASSERT_DELTA(solver.rGetDeformedPosition()[nodes[node]](1), y[node], 1e-3);
            TS_ASSERT_DELTA(solver.rGetDeformedPosition()[nodes[node]](2), z[node], 1e-3);
        }

        // tidy up memory
        delete p_fine_mesh;
        delete p_pair;
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
//    void strangefailingbehaviourinparallelTestCompareWithDeadExplicitSolver()
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
//          //The following lines are not relevant to this test but need to be there
//          TetrahedralMesh<2,2>* p_fine_mesh = new TetrahedralMesh<2,2>();//unused in this test
//          p_fine_mesh->ConstructRegularSlabMesh(0.125, 1.0, 1.0);
//          TetrahedralMesh<2,2>* p_coarse_mesh = new TetrahedralMesh<2,2>();//unused in this test
//          p_coarse_mesh->ConstructRegularSlabMesh(0.125, 1.0, 1.0);
//          FineCoarseMeshPair<2>* p_pair = new FineCoarseMeshPair<2>(*p_fine_mesh, *p_coarse_mesh);//also unused in this test
//          p_pair->SetUpBoxesOnFineMesh();
//          p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(solver.GetQuadratureRule()), false);
//          p_pair->DeleteFineBoxCollection();
//          solver.SetFineCoarseMeshPair(p_pair);
//          ///////////////////////////////////////////////////////////////////////////
//          solver.Initialise();

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
//
//        //in need of deletion even if all these 3 have no influence at all on this test
//        delete p_fine_mesh;
//        delete p_coarse_mesh;
//        delete p_pair;
//    }
};

#endif /*TESTIMPLICITCARDIACMECHANICSSOLVER_HPP_*/
