/*

Copyright (c) 2005-2012, University of Oxford.
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


#ifndef TESTCARDIACELECTROMECHANICSPROBLEM_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEM_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "LuoRudy1991.hpp"
#include "NonlinearElasticityTools.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "NumericFileComparison.hpp"
#include "Hdf5DataReader.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"

// cell factory which stimulates everything at once
class EntirelyStimulatedTissueCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    EntirelyStimulatedTissueCellFactory()
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-100000.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
    }
};



class TestCardiacElectroMechanicsProblem : public CxxTest::TestSuite
{
public:
    // In this test all nodes are stimulated at the same time, hence
    // exactly the same active tension is generated at all nodes,
    // so the internal force in entirely homogeneous.
    // We fix the x-coordinate of the X=0, but leave the y-coordinate
    // free - ie sliding boundary conditions. With a homogeneous
    // force this means the solution should be a perfect rectangle,
    // ie x=alpha*X, y=beta*Y for some alpha, beta.
    void TestWithHomogeneousEverythingCompressible() throw(Exception)
    {
        EntirelyStimulatedTissueCellFactory cell_factory;

        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01, 0.05, 0.05);

        QuadraticMesh<2> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.025, 0.05, 0.05);

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,2> > fixed_node_locations;

        // fix the node at the origin so that the solution is well-defined (ie unique)
        fixed_nodes.push_back(0);
        fixed_node_locations.push_back(zero_vector<double>(2));

        // for the rest of the nodes, if they lie on X=0, fix x=0 but leave y free.
        for(unsigned i=1 /*not 0*/; i<mechanics_mesh.GetNumNodes(); i++)
        {
            if(fabs(mechanics_mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                c_vector<double,2> new_position;
                new_position(0) = 0.0;
                new_position(1) = SolidMechanicsProblemDefinition<2>::FREE;
                fixed_nodes.push_back(i);
                fixed_node_locations.push_back(new_position);
            }
        }

        ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(COMPRESSIBLE);
        problem_defn.SetFixedNodes(fixed_nodes, fixed_node_locations);
        problem_defn.SetMechanicsSolveTimestep(1.0);

        // the following is just for coverage - applying a zero pressure so has no effect on deformation
        std::vector<BoundaryElement<1,2>*> boundary_elems;
        boundary_elems.push_back(* (mechanics_mesh.GetBoundaryElementIteratorBegin()));
        problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, 0.0);

        HeartConfig::Instance()->SetSimulationDuration(10.0);

        CardiacElectroMechanicsProblem<2> problem(COMPRESSIBLE,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  &cell_factory,
                                                  &problem_defn,
                                                  "TestCardiacEmHomogeneousEverythingCompressible");
        problem.Solve();
        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();

        // not sure how easy is would be determine what the deformation should be
        // exactly, but it certainly should be constant squash in X direction, constant
        // stretch in Y.

        // first, check node 8 starts is the far corner
        assert(fabs(mechanics_mesh.GetNode(8)->rGetLocation()[0] - 0.05)<1e-8);
        assert(fabs(mechanics_mesh.GetNode(8)->rGetLocation()[1] - 0.05)<1e-8);

        double X_scale_factor = r_deformed_position[8](0)/0.05;
        double Y_scale_factor = r_deformed_position[8](1)/0.05;

        std::cout << "Scale_factors = " << X_scale_factor << " " << Y_scale_factor << ", product = " << X_scale_factor*Y_scale_factor<<"\n";

        for(unsigned i=0; i<mechanics_mesh.GetNumNodes(); i++)
        {
            double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
            double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];

            TS_ASSERT_DELTA( r_deformed_position[i](0), X * X_scale_factor, 1e-6);
            TS_ASSERT_DELTA( r_deformed_position[i](1), Y * Y_scale_factor, 1e-6);
        }
    }

    // This test is virtually identical to the above test, except it uses
    // incompressible solid mechanics. Since the solution should be
    // x=alpha*X, y=beta*Y for some alpha, beta (see comments for above test),
    // we also test that alpha*beta = 1.0
    void TestWithHomogeneousEverythingIncompressible() throw(Exception)
    {
        EntirelyStimulatedTissueCellFactory cell_factory;

        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01, 0.05, 0.05);

        QuadraticMesh<2> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.025, 0.05, 0.05);

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,2> > fixed_node_locations;
        fixed_nodes.push_back(0);
        fixed_node_locations.push_back(zero_vector<double>(2));
        for(unsigned i=1 /*not 0*/; i<mechanics_mesh.GetNumNodes(); i++)
        {
            if(fabs(mechanics_mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                c_vector<double,2> new_position;
                new_position(0) = 0.0;
                new_position(1) = SolidMechanicsProblemDefinition<2>::FREE;
                fixed_nodes.push_back(i);
                fixed_node_locations.push_back(new_position);
            }
        }

        ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        problem_defn.SetFixedNodes(fixed_nodes, fixed_node_locations);
        problem_defn.SetMechanicsSolveTimestep(1.0);

        // coverage, this file is just X-direction fibres
        problem_defn.SetVariableFibreSheetDirectionsFile("heart/test/data/fibre_tests/2by2mesh_fibres.ortho", false);

        HeartConfig::Instance()->SetSimulationDuration(10.0);

        CardiacElectroMechanicsProblem<2> problem(INCOMPRESSIBLE,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  &cell_factory,
                                                  &problem_defn,
                                                  "TestCardiacEmHomogeneousEverythingIncompressible");

        problem.Solve();


        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();

        // not sure how easy is would be determine what the deformation should be
        // exactly, but it certainly should be constant squash in X direction, constant
        // stretch in Y.

        // first, check node 8 starts is the far corner
        assert(fabs(mechanics_mesh.GetNode(8)->rGetLocation()[0] - 0.05)<1e-8);
        assert(fabs(mechanics_mesh.GetNode(8)->rGetLocation()[1] - 0.05)<1e-8);

        double X_scale_factor = r_deformed_position[8](0)/0.05;
        double Y_scale_factor = r_deformed_position[8](1)/0.05;

        std::cout << "Scale_factors = " << X_scale_factor << " " << Y_scale_factor << ", product = " << X_scale_factor*Y_scale_factor<<"\n";

        for(unsigned i=0; i<mechanics_mesh.GetNumNodes(); i++)
        {
            double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
            double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];

            TS_ASSERT_DELTA( r_deformed_position[i](0), X * X_scale_factor, 1e-6);
            TS_ASSERT_DELTA( r_deformed_position[i](1), Y * Y_scale_factor, 1e-6);
        }

        // here, we also test the deformation is incompressible
        TS_ASSERT_DELTA(X_scale_factor * Y_scale_factor, 1.0, 1e-6);
    }

    // These tests are older than the above tests..
    void TestImplicitNhs2dOneMechanicsElement() throw(Exception)
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

        HeartConfig::Instance()->SetSimulationDuration(10.0);

        CardiacElectroMechProbRegularGeom<2> problem(INCOMPRESSIBLE,
                                                     0.05, /* width (cm) */
                                                     1,    /* mech mesh size*/
                                                     5,    /* elec elem each dir */
                                                     &cell_factory,
                                                     NHS,
                                                     1.0,  /* mechanics solve timestep */
                                                     0.01, /* contraction model ode timestep */
                                                     "TestCardiacElectroMechOneElement");
        c_vector<double,2> pos;
        pos(0) = 0.05;
        pos(1) = 0.0;

        problem.SetWatchedPosition(pos);

        if(PetscTools::IsSequential()) // #2084
        {
            TS_ASSERT_THROWS_CONTAINS(problem.SetOutputDeformationGradientsAndStress(3.4),"not a multiple");
        	problem.SetOutputDeformationGradientsAndStress(3.0);
        }
        problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();
        TS_ASSERT_DELTA(r_deformed_position[1](0), 0.0497, 1e-4);

        OutputFileHandler handler("TestCardiacElectroMechOneElement",false);

        NumericFileComparison comparer(handler.GetOutputDirectoryFullPath() + "watched.txt","heart/test/data/good_watched.txt");
        TS_ASSERT(comparer.CompareFiles(1e-2));

        // check electrics output was written
        std::string command = "ls " + handler.GetOutputDirectoryFullPath() + "/electrics";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        if(PetscTools::IsSequential())
        {
			command = "ls " + handler.GetOutputDirectoryFullPath() + "/deformation/deformation_gradient_0.strain";
			TS_ASSERT_EQUALS(system(command.c_str()), 0);

			command = "ls " + handler.GetOutputDirectoryFullPath() + "/deformation/deformation_gradient_3.strain";
			TS_ASSERT_EQUALS(system(command.c_str()), 0);

			command = "ls " + handler.GetOutputDirectoryFullPath() + "/deformation/deformation_gradient_6.strain";
			TS_ASSERT_EQUALS(system(command.c_str()), 0);

			command = "ls " + handler.GetOutputDirectoryFullPath() + "/deformation/second_PK_0.stress";
			TS_ASSERT_EQUALS(system(command.c_str()), 0);

			command = "ls " + handler.GetOutputDirectoryFullPath() + "/deformation/second_PK_3.stress";
			TS_ASSERT_EQUALS(system(command.c_str()), 0);

			command = "ls " + handler.GetOutputDirectoryFullPath() + "/deformation/second_PK_6.stress";
			TS_ASSERT_EQUALS(system(command.c_str()), 0);
        }


        // coverage

        HeartConfig::Instance()->SetSimulationDuration(10.0); // has to be reset after a solve, it seems..
        CardiacElectroMechProbRegularGeom<2> prob_with_bad_model(INCOMPRESSIBLE,0.05,1,5,&cell_factory,NONPHYSIOL1,1,0.01,"");
        TS_ASSERT_THROWS_CONTAINS(prob_with_bad_model.Solve(),"Invalid contraction model");


        CardiacElectroMechProbRegularGeom<2> prob_with_bad_model_comp(COMPRESSIBLE,0.05,1,5,&cell_factory,NONPHYSIOL1,1,0.01,"");
        TS_ASSERT_THROWS_CONTAINS(prob_with_bad_model_comp.Solve(),"Invalid contraction model");

        CardiacElectroMechProbRegularGeom<2> prob_with_bad_timesteps(INCOMPRESSIBLE,0.05,1,5,&cell_factory,NHS,0.025,0.01,"");
        TS_ASSERT_THROWS_CONTAINS(prob_with_bad_timesteps.Initialise(),"does not divide");


        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }

    void TestWithKerchoffs() throw(Exception)
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

        HeartConfig::Instance()->SetSimulationDuration(20.0);

        CardiacElectroMechProbRegularGeom<2> problem(INCOMPRESSIBLE,
                                                     0.05, /* width (cm) */
                                                     1,    /* mech mesh size*/
                                                     5,    /* elec elem each dir */
                                                     &cell_factory,
                                                     KERCHOFFS2003,
                                                     1.0,   /* mechanics solve timestep */
                                                     0.01,  /* Kerchoffs ode timestep */
                                                     "TestCardiacEmWithKerchoffs");

        c_vector<double,2> pos;
        pos(0) = 0.05;
        pos(1) = 0.0;

        problem.SetWatchedPosition(pos);
        problem.SetNoElectricsOutput();
        problem.Initialise();

        problem.Solve();

        //visualise to verify

        // hardcoded result
        TS_ASSERT_EQUALS(problem.mWatchedMechanicsNodeIndex, 1u);
        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[1](0), 0.0479, 0.0002);
        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[1](1),-0.0003, 0.0002);
    }

    //
    //  BAD test - fails with HYPRE (for some reason HYPRE can't solve the one of the linear systems, and
    //  the search direction in the end doesn't decrease the residual), and also with ILU if you increase
    //  the number of elements (whether LR91 or N98 is used). Probably the active tension is too high.
    //
    //
    void TestExplicitSolverWithNash2004() throw(Exception)
    {
#ifdef MECH_USE_HYPRE
        TS_FAIL("This test is known to fail with HYPRE - see comments in test");
        return;
#endif

        PlaneStimulusCellFactory<CML_noble_varghese_kohl_noble_1998_basic_with_sac, 2> cell_factory(-1000*1000);

        HeartConfig::Instance()->SetSimulationDuration(20.0);

        CardiacElectroMechProbRegularGeom<2> problem(INCOMPRESSIBLE,
                                                     0.05, /* width (cm) */
                                                     1,    /* mech mesh size*/
                                                     5,    /* elec elem each dir */
                                                     &cell_factory,
                                                     NASH2004,
                                                     1.0,   /* mechanics solve timestep */
                                                     0.01,  /* nash ode timestep */
                                                     "TestExplicitWithNash");

        c_vector<double,2> pos;
        pos(0) = 0.05;
        pos(1) = 0.0;

        problem.SetWatchedPosition(pos);
        problem.SetNoElectricsOutput();
        problem.Initialise();

        problem.Solve();

        //visualise to verify

        // hardcoded result
        TS_ASSERT_EQUALS(problem.mWatchedMechanicsNodeIndex, 1u);
        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[1](0), 0.0419, 0.0002);
    }

    void TestWithCompressibleApproach() throw(Exception)
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

        HeartConfig::Instance()->SetSimulationDuration(20.0);

        CardiacElectroMechProbRegularGeom<2> problem(COMPRESSIBLE,
                                                     0.05, /* width (cm) */
                                                     1,    /* mech mesh size*/
                                                     5,    /* elec elem each dir */
                                                     &cell_factory,
                                                     KERCHOFFS2003,
                                                     1.0,   /* mechanics solve timestep */
                                                     0.01,  /* Kerchoffs ode timestep */
                                                     "TestCompressibleWithKerchoffs");

        problem.Solve();

        // Mainly just testing no errors when Solve was called.
        // The results of this test can be visually compared with the results of the
        // equivalent incompressible simulation in TestWithKerchoffs.

        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[1](0), 0.0465, 0.0002);
        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[1](1),-0.0012, 0.0002);

        // create and initialise an incompressible NASH2004 problem, just for coverage..
        HeartConfig::Instance()->SetSimulationDuration(20.0);
        CardiacElectroMechProbRegularGeom<2> problem2(COMPRESSIBLE,0.05,1,5,&cell_factory,
                                                      NASH2004,
                                                      1.0, 0.01,"");
        problem2.Initialise();
    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEM_HPP_*/
