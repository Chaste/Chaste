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


#ifndef TESTEXPLICITCARDIACMECHANICSSOLVER_HPP_
#define TESTEXPLICITCARDIACMECHANICSSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "ExplicitCardiacMechanicsSolver.hpp"
#include "ImplicitCardiacMechanicsSolver.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraturePointsGroup.hpp"
#include "NonlinearElasticityTools.hpp"
#include "ReplicatableVector.hpp"

// some useful typedefs
typedef ExplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<2>,2> IncompressibleExplicitSolver2d;
typedef ExplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<3>,3> IncompressibleExplicitSolver3d;
typedef ImplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<2>,2> IncompressibleImplicitSolver2d;

class TestExplicitCardiacMechanicsSolver : public CxxTest::TestSuite
{
public:
    void TestWithSimpleContractionModel()
    {
        QuadraticMesh<2> mesh(0.25, 1.0, 1.0);
        MooneyRivlinMaterialLaw<2> law(1);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetContractionModel(NONPHYSIOL1,0.01);
        problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass

        // NONPHYSIOL1 => NonphysiologicalContractionModel 1
        IncompressibleExplicitSolver2d solver(mesh,problem_defn,"TestExplicitCardiacMech");

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

        // coverage
        QuadraturePointsGroup<2> quad_points(mesh, *(solver.GetQuadratureRule()));

        std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints(), 0.0);
        std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);

        solver.SetCalciumAndVoltage(calcium_conc, voltages);

        // solve UP TO t=0. So Ta(lam_n,t_{n+1})=5*sin(0)=0, ie no deformation
        solver.Solve(-0.01,0.0,0.01);
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(),0u);

        solver.Solve(0.24,0.25,0.01);

        TS_ASSERT_DELTA(solver.rGetDeformedPosition()[4](0),  0.9732, 1e-2);
        TS_ASSERT_DELTA(solver.rGetDeformedPosition()[4](1), -0.0156, 1e-2);

        //in need of deletion even if all these 3 have no influence at all on this test
        delete p_fine_mesh;
        delete p_coarse_mesh;
        delete p_pair;
    }

    // with stretch (and stretch-rate) independent contraction models the implicit and explicit schemes
    // are identical
    void TestCompareImplicitAndExplicitWithStretchIndependentContractionModel()
    {
        QuadraticMesh<2> mesh(0.25, 1.0, 1.0);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetContractionModel(NONPHYSIOL1,0.01);
        problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass

        //The following lines are not relevant to this test but need to be there
        TetrahedralMesh<2,2>* p_fine_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_fine_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0);
        TetrahedralMesh<2,2>* p_coarse_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_coarse_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0);
        FineCoarseMeshPair<2>* p_pair = new FineCoarseMeshPair<2>(*p_fine_mesh, *p_coarse_mesh);//also unused in this test
        p_pair->SetUpBoxesOnFineMesh();
        /////////////////////////////////////////////////////////////////////

        // NONPHYSIOL 1 - contraction model is of the form sin(t)
        IncompressibleExplicitSolver2d expl_solver(mesh,problem_defn,""/*"TestCompareExplAndImplCardiacSolvers_Exp"*/);
        p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(expl_solver.GetQuadratureRule()), false);
        p_pair->DeleteFineBoxCollection();
        expl_solver.SetFineCoarseMeshPair(p_pair);
        expl_solver.Initialise();

        IncompressibleImplicitSolver2d impl_solver(mesh,problem_defn,""/*"TestCompareExplAndImplCardiacSolvers_Imp"*/);
        impl_solver.SetFineCoarseMeshPair(p_pair);
        impl_solver.Initialise();

        double dt = 0.25;
        for (double t=0; t<3; t+=dt)
        {
            expl_solver.Solve(t,t+dt,dt);
            impl_solver.Solve(t,t+dt,dt);

            // computations should be identical
            TS_ASSERT_EQUALS(expl_solver.GetNumNewtonIterations(), impl_solver.GetNumNewtonIterations());
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(expl_solver.rGetDeformedPosition()[i](0),  impl_solver.rGetDeformedPosition()[i](0), 1e-9);
                TS_ASSERT_DELTA(expl_solver.rGetDeformedPosition()[i](1),  impl_solver.rGetDeformedPosition()[i](1), 1e-9);
            }
        }

        //in need of deletion even if all these 3 have no influence at all on this test
        delete p_fine_mesh;
        delete p_coarse_mesh;
        delete p_pair;
    }


    // with stretch-dependent contraction models the implicit and explicit schemes can be similar
    void TestCompareImplicitAndExplicitWithStretchDependentContractionModel()
    {
        QuadraticMesh<2> mesh(0.25, 1.0, 1.0);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetContractionModel(NONPHYSIOL2,0.01);
        problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass

        //The following lines are not relevant to this test but need to be there
        TetrahedralMesh<2,2>* p_fine_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_fine_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0);
        TetrahedralMesh<2,2>* p_coarse_mesh = new TetrahedralMesh<2,2>();//unused in this test
        p_coarse_mesh->ConstructRegularSlabMesh(0.25, 1.0, 1.0);
        FineCoarseMeshPair<2>* p_pair = new FineCoarseMeshPair<2>(*p_fine_mesh, *p_coarse_mesh);//also unused in this test
        p_pair->SetUpBoxesOnFineMesh();
        /////////////////////////////////////////////////////////////////////

        // NONPHYSIOL 2 - contraction model is of the form lam*sin(t)
        IncompressibleExplicitSolver2d expl_solver(mesh,problem_defn,"TestCompareExplAndImplCardiacSolversStretch_Exp");
        p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(expl_solver.GetQuadratureRule()), false);
        p_pair->DeleteFineBoxCollection();
        expl_solver.SetFineCoarseMeshPair(p_pair);
        expl_solver.Initialise();

        IncompressibleImplicitSolver2d impl_solver(mesh,problem_defn,"TestCompareExplAndImplCardiacSolversStretch_Imp");
        impl_solver.SetFineCoarseMeshPair(p_pair);
        impl_solver.Initialise();

        expl_solver.WriteCurrentSpatialSolution("solution","nodes",0);
        impl_solver.WriteCurrentSpatialSolution("solution","nodes",0);

        unsigned counter = 1;

        double t0 = 0.0;
        double t1 = 0.25; // to be quite quick (min stretch ~=0.88), make this 5/4 (?) say for min stretch < 0.7
        double dt = 0.025;

        for (double t=t0; t<t1; t+=dt)
        {
            //std::cout << "\n **** t = " << t << " ****\n" << std::flush;

            expl_solver.SetWriteOutput(false);
            expl_solver.Solve(t,t+dt,dt);
            expl_solver.SetWriteOutput();
            expl_solver.WriteCurrentSpatialSolution("solution","nodes",counter);

            impl_solver.SetWriteOutput(false);
            impl_solver.Solve(t,t+dt,dt);
            impl_solver.SetWriteOutput();
            impl_solver.WriteCurrentSpatialSolution("solution","nodes",counter);

            // the solutions turn out to be very close to each other
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(expl_solver.rGetDeformedPosition()[i](0),  impl_solver.rGetDeformedPosition()[i](0), 2e-3);
                TS_ASSERT_DELTA(expl_solver.rGetDeformedPosition()[i](1),  impl_solver.rGetDeformedPosition()[i](1), 2e-3);
            }

            // visualisation in matlab or octave
            // run from CHASTETESTOUTPUT
            //
            // close all; figure; hold on
            // for i=0:10, x1 = load(['TestCompareExplAndImplCardiacSolversStretch_Exp/solution_',num2str(i),'.nodes']); plot(i,x1(5,1),'b*'); end
            // for i=0:10, x2 = load(['TestCompareExplAndImplCardiacSolversStretch_Imp/solution_',num2str(i),'.nodes']); plot(i,x2(5,1),'r*'); end
            //
            // close all; figure; hold on
            // for i=0:10, x2 = load(['TestCompareExplAndImplCardiacSolversStretch_Imp/solution_',num2str(i),'.nodes']); plot(x2(:,1),x2(:,2),'r*'); end

            counter++;
        }

        // check there was significant deformation - node 4 is (1,0)
        TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_LESS_THAN(impl_solver.rGetDeformedPosition()[4](0), 0.9);

        //in need of deletion even if all these 3 have no influence at all on this test
        delete p_fine_mesh;
        delete p_coarse_mesh;
        delete p_pair;
    }

    // cover all other contraction model options which are allowed but not been used in a test so far
    void TestCoverage()
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1.0);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetContractionModel(NONPHYSIOL3,0.01);
        problem_defn.SetMechanicsSolveTimestep(0.01); //This is only set to make ElectroMechanicsProblemDefinition::Validate pass

        TS_ASSERT_THROWS_THIS(problem_defn.SetApplyAnisotropicCrossFibreTension(true,1.0,1.0),
                              "You can only apply anisotropic cross fibre tensions in a 3D simulation.");

        TetrahedralMesh<2,2>* p_fine_mesh = new TetrahedralMesh<2,2>();
        p_fine_mesh->ConstructRegularSlabMesh(1.0, 1.0, 1.0);
        TetrahedralMesh<2,2>* p_coarse_mesh = new TetrahedralMesh<2,2>();
        p_coarse_mesh->ConstructRegularSlabMesh(1.0, 1.0, 1.0);
        FineCoarseMeshPair<2>* p_pair = new FineCoarseMeshPair<2>(*p_fine_mesh, *p_coarse_mesh);
        p_pair->SetUpBoxesOnFineMesh();

        TetrahedralMesh<2,2>* p_coarse_mesh_big = new TetrahedralMesh<2,2>();
        p_coarse_mesh_big->ConstructRegularSlabMesh(1.0, 3.0, 3.0);
        FineCoarseMeshPair<2>* p_pair_wrong = new FineCoarseMeshPair<2>(*p_fine_mesh, *p_coarse_mesh_big);

        // Some strange memory thing happens if we don't have the test with the
        // exceptions below in a separate block. I have no idea why.
        {
            IncompressibleExplicitSolver2d expl_solver(mesh,problem_defn,"");
            p_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(expl_solver.GetQuadratureRule()), false);
            p_pair->DeleteFineBoxCollection();
            expl_solver.SetFineCoarseMeshPair(p_pair);
            expl_solver.Initialise();

            // Prevent an assertion being thrown about setting the cell factory more than once.
            delete problem_defn.mpContractionCellFactory;
            problem_defn.mpContractionCellFactory = NULL;

            problem_defn.SetContractionModel(NASH2004,0.01);
            IncompressibleExplicitSolver2d expl_solver_with_nash(mesh,problem_defn,"");
            expl_solver_with_nash.SetFineCoarseMeshPair(p_pair);
            expl_solver_with_nash.Initialise();

            // Prevent an assertion being thrown about setting the cell factory more than once.
            delete problem_defn.mpContractionCellFactory;
            problem_defn.mpContractionCellFactory = NULL;

            problem_defn.SetContractionModel(KERCHOFFS2003,0.01);
            IncompressibleExplicitSolver2d expl_solver_with_kerchoffs(mesh,problem_defn,"");
            expl_solver_with_kerchoffs.SetFineCoarseMeshPair(p_pair);
            expl_solver_with_kerchoffs.Initialise();

        }

        {
            // bad contraction model
            ElectroMechanicsProblemDefinition<2> problem_defn2(mesh);
            problem_defn2.SetMaterialLaw(INCOMPRESSIBLE,&law);
            problem_defn2.SetContractionModel(NHS,0.01);
            IncompressibleExplicitSolver2d solver(mesh,problem_defn2,"");
            TS_ASSERT_THROWS_THIS(solver.SetFineCoarseMeshPair(p_pair_wrong),
                    "When setting a mesh pair into the solver, the coarse mesh of the mesh pair must be the same as the quadratic mesh");
            solver.SetFineCoarseMeshPair(p_pair);

            TS_ASSERT_THROWS_THIS(solver.Initialise(),
                    "stretch-rate-dependent contraction model requires an IMPLICIT cardiac mechanics solver.");
        }

        delete p_fine_mesh;
        delete p_coarse_mesh;
        delete p_pair;
        delete p_coarse_mesh_big;
        delete p_pair_wrong;
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
        problem_defn.SetContractionModel(NONPHYSIOL1,0.01);
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
        x[0] = 0.9732;
        x[1] = 0.9759;
        x[2] = 0.9865;
        x[3] = 1.0; //Note, for cross_tension == 1.0 there should be no deformation of the tissue square (tensions balance)
        y[0] = -0.0156;
        y[1] = -0.0140;
        y[2] = -0.0077;
        y[3] = 0.0;

        for (unsigned i=0; i < tension_fractions.size();i++)
        {
            problem_defn.SetApplyIsotropicCrossFibreTension(true,tension_fractions[i]);

            TS_ASSERT_EQUALS(problem_defn.GetApplyCrossFibreTension(), true);
            TS_ASSERT_DELTA(problem_defn.GetSheetTensionFraction(),tension_fractions[i], 1e-6);
            TS_ASSERT_DELTA(problem_defn.GetSheetNormalTensionFraction(),tension_fractions[i], 1e-6);

            // NONPHYSIOL1 => NonphysiologicalContractionModel 1
            IncompressibleExplicitSolver2d solver(mesh,problem_defn,"TestExplicitCardiacMech");

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

            // coverage
            QuadraturePointsGroup<2> quad_points(mesh, *(solver.GetQuadratureRule()));

            std::vector<double> calcium_conc(solver.GetTotalNumQuadPoints(), 0.0);
            std::vector<double> voltages(solver.GetTotalNumQuadPoints(), 0.0);

            solver.SetCalciumAndVoltage(calcium_conc, voltages);

            // solve UP TO t=0. So Ta(lam_n,t_{n+1})=5*sin(0)=0, ie no deformation
            solver.Solve(-0.01,0.0,0.01);
            TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(),0u);

            solver.Solve(0.24,0.25,0.01);

//// These fail due to changes from #2185 but no point fixing as these numbers will change later anyway (#2180)
            TS_ASSERT_DELTA(solver.rGetDeformedPosition()[4](0), x[i], 1e-3);
            TS_ASSERT_DELTA(solver.rGetDeformedPosition()[4](1), y[i], 1e-3);

            //in need of deletion even if all these 3 have no influence at all on this test
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
        IncompressibleExplicitSolver3d solver(mesh,problem_defn,"TestIsotropicCrossFibreExplicit");

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
        IncompressibleExplicitSolver3d solver(mesh,problem_defn,"TestAnisotropicCrossFibreExplicit");

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

        // Tidy up memory
        delete p_fine_mesh;
        delete p_pair;
    }
};

#endif /*TESTEXPLICITCARDIACMECHANICSSOLVER_HPP_*/
