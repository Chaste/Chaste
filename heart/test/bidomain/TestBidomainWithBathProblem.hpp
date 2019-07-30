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

#ifndef TESTBIDOMAINWITHBATH_HPP_
#define TESTBIDOMAINWITHBATH_HPP_

#include "CardiacSimulationArchiver.hpp"

#include <cxxtest/TestSuite.h>
#include <vector>

#include "BidomainWithBathProblem.hpp"

#include "LuoRudy1991.hpp"
#include "TenTusscher2006EpiBackwardEuler.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "ConstBoundaryCondition.hpp"
#include "HeartEventHandler.hpp"
#include "HeartRegionCodes.hpp"
#include "Timer.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SimpleBathProblemSetup.hpp"
#include "HeartConfig.hpp"
#include "NumericFileComparison.hpp"

class TestBidomainWithBathProblem : public CxxTest::TestSuite
{
public:
    void tearDown()
    {
        HeartConfig::Reset();
    }

    void TestLabellingNodes()
    {
        HeartConfig::Instance()->SetSimulationDuration(0.01);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_10_elements_with_two_attributes");
        HeartConfig::Instance()->SetOutputDirectory("bidomain_bath");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> bidomain_cell_factory;
        BidomainWithBathProblem<1> bidomain_problem( &bidomain_cell_factory );
        bidomain_problem.Initialise();

        AbstractTetrahedralMesh<1,1>* p_mesh = &(bidomain_problem.rGetMesh());

        // the middle 4 elements are 'heart' elements (ie region=0),
        // so the middle 5 nodes should be heart nodes
        char expected_node_regions[11]={ 'B', 'B', 'B',
                       'T', 'T', 'T', 'T', 'T',
                       'B','B','B'};

        for (unsigned i=0; i<11; i++)
        {
            if (p_mesh->GetDistributedVectorFactory()->IsGlobalIndexLocal(i))
            {
                if (expected_node_regions[i] == 'B')
                {
                    TS_ASSERT(HeartRegionCode::IsRegionBath( p_mesh->GetNode(i)->GetRegion()) );
                }
                else
                {
                    TS_ASSERT_EQUALS(expected_node_regions[i], 'T' );
                    TS_ASSERT(HeartRegionCode::IsRegionTissue( p_mesh->GetNode(i)->GetRegion()) );
                }
            }
        }
        // we need to call solve as otherwise an EventHandler exception is thrown
        bidomain_problem.Solve();
    }


    void TestFailsIfNoBathElements()
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("bidomain_bath");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> bidomain_cell_factory;
        BidomainWithBathProblem<1> bidomain_problem( &bidomain_cell_factory );
        // Fails because no bath
        TS_ASSERT_THROWS_THIS(bidomain_problem.Initialise(), "No bath element found");

        // Prevent an EventHandling exception in later tests
        HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
    }

    void TestCheckForBathElementsNoDeadlock()
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("bidomain_bath");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> bidomain_cell_factory;
        BidomainWithBathProblem<1> bidomain_problem( &bidomain_cell_factory );

        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        try
        {
            mesh.GetElement(0)->SetAttribute(HeartRegionCode::GetValidBathId());
        }
        catch(Exception&)
        {
            // I don't own element 0
        }

        bidomain_problem.SetMesh(&mesh);

        // Fails because no bath
        TS_ASSERT_THROWS_NOTHING(bidomain_problem.Initialise());

        // Prevent an EventHandling exception in later tests
        HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
    }


    void TestBathIntracellularStimulation()
    {
        HeartConfig::Instance()->SetSimulationDuration(10.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_1d");

        c_vector<double,1> centre;
        centre(0) = 0.5;
        BathCellFactory<1> cell_factory(-1e6, centre); // stimulates x=0.5 node

        BidomainWithBathProblem<1> bidomain_problem( &cell_factory );

        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        // set the x<0.25 and x>0.75 regions as the bath region
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            double x = mesh.GetElement(i)->CalculateCentroid()[0];
            if ((x<0.25) || (x>0.75))
            {
                mesh.GetElement(i)->SetAttribute(HeartRegionCode::GetValidBathId());
            }
        }

        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();

        bidomain_problem.Solve();

        Vec sol = bidomain_problem.GetSolution();
        ReplicatableVector sol_repl(sol);

        // test V = 0 for all bath nodes
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if (HeartRegionCode::IsRegionBath( mesh.GetNode(i)->GetRegion() )) // bath
            {
                TS_ASSERT_DELTA(sol_repl[2*i], 0.0, 1e-12);
            }
        }

        // test symmetry of V and phi_e
        for (unsigned i=0; i<=(mesh.GetNumNodes()-1)/2; i++)
        {
            unsigned opposite = mesh.GetNumNodes()-i-1;
            TS_ASSERT_DELTA(sol_repl[2*i], sol_repl[2*opposite], 2e-3);      // V
            TS_ASSERT_DELTA(sol_repl[2*i+1], sol_repl[2*opposite+1], 2e-3);  // phi_e
        }

        // a couple of hardcoded values
        TS_ASSERT_DELTA(sol_repl[2*50], 3.7684, 1e-3);
        TS_ASSERT_DELTA(sol_repl[2*70], 5.1777, 1e-3);
    }


    // In this test we have no cardiac tissue, so that the equations are just sigma * phi_e''=0
    // throughout the domain (with a Neumann boundary condition on x=1 and a dirichlet boundary
    // condition (ie grounding) on x=0), so the exact solution can be calculated and compared
    // against.
    void Test1dProblemOnlyBathGroundedOneSide()
    {
        HeartConfig::Instance()->SetSimulationDuration(0.5);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBathOnlyBath");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath");

        c_vector<double,1> centre;
        centre(0) = 0.5;
        BathCellFactory<1> cell_factory(-1e6, centre);

        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            mesh.GetElement(i)->SetAttribute(HeartRegionCode::GetValidBathId());
        }

        // create boundary conditions container
        double boundary_val = 1.0;
        boost::shared_ptr<BoundaryConditionsContainer<1,1,2> > p_bcc(new BoundaryConditionsContainer<1,1,2>);
        ConstBoundaryCondition<1>* p_bc_stim = new ConstBoundaryCondition<1>(boundary_val);
        ConstBoundaryCondition<1>* p_zero_stim = new ConstBoundaryCondition<1>(0.0);

        // loop over boundary elements and set (sigma\gradphi).n = 1.0 on RHS edge
        for (TetrahedralMesh<1,1>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
           iter != mesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            if (((*iter)->GetNodeLocation(0))[0]==1.0)
            {
                /// \todo: I think you need to provide a boundary condition for unknown#1 if you are gonig to provide one for unknown#2?
                p_bcc->AddNeumannBoundaryCondition(*iter, p_zero_stim, 0);
                p_bcc->AddNeumannBoundaryCondition(*iter, p_bc_stim,   1);
            }
        }

        BidomainWithBathProblem<1> bidomain_problem( &cell_factory );

        bidomain_problem.SetBoundaryConditionsContainer(p_bcc);
        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();

        // fix phi=0 on LHS edge
        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);
        bidomain_problem.SetFixedExtracellularPotentialNodes(fixed_nodes);

        bidomain_problem.Solve();

        Vec sol = bidomain_problem.GetSolution();
        ReplicatableVector sol_repl(sol);

        // test phi = x*boundary_val/sigma (solution of phi''=0, phi(0)=0, sigma*phi'(1)=boundary_val
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double bath_cond = HeartConfig::Instance()->GetBathConductivity();
            double x = mesh.GetNode(i)->rGetLocation()[0];
            TS_ASSERT_DELTA(sol_repl[2*i],   0.0,   1e-12);               // V
            TS_ASSERT_DELTA(sol_repl[2*i+1], x*boundary_val/bath_cond, 1e-4);   // phi_e
        }
    }

    void Test2dBathIntracellularStimulation()
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath2d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_2d");

        c_vector<double,2> centre;
        centre(0) = 0.05;
        centre(1) = 0.05;
        BathCellFactory<2> cell_factory(-5e6, centre); // stimulates x=0.05 node

        BidomainWithBathProblem<2> bidomain_problem( &cell_factory );

        DistributedTetrahedralMesh<2,2>* p_mesh = Load2dMeshAndSetCircularTissue<DistributedTetrahedralMesh<2,2> >(
            "mesh/test/data/2D_0_to_1mm_400_elements", 0.05, 0.05, 0.04);

        bidomain_problem.SetMesh(p_mesh);
        bidomain_problem.Initialise();

        bidomain_problem.Solve();

        Vec sol = bidomain_problem.GetSolution();
        ReplicatableVector sol_repl(sol);

        // test V = 0 for all bath nodes
        for (AbstractTetrahedralMesh<2,2>::NodeIterator iter=p_mesh->GetNodeIteratorBegin();
             iter != p_mesh->GetNodeIteratorEnd(); ++iter)
        {
            if (HeartRegionCode::IsRegionBath( (*iter).GetRegion() )) // bath
            {
                unsigned index=(*iter).GetIndex();

                TS_ASSERT_DELTA(sol_repl[2*index], 0.0, 1e-12);
            }
        }

        // A couple of hardcoded values. We would normally have to look up the index in the
        // permutation vector when using a DistributedTetrahedralMesh, however the call to
        // Load2dMeshAndSetCircularTissue above imposes a DUMB partition, so no need.
        TS_ASSERT_DELTA(sol_repl[2*50], 28.3912, 1e-3); // node 50
        TS_ASSERT_DELTA(sol_repl[2*70], 28.3912, 1e-3); // node 70

        delete p_mesh;
    }

    void Test2dBathInputFluxEqualsOutputFlux()
    {
        HeartConfig::Instance()->SetSimulationDuration(3.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath2dFluxCompare");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_2d_fluxes");

        // Coverage of Hdf5ToCmguiConverter with bath
        HeartConfig::Instance()->SetVisualizeWithCmgui(true);

        HeartConfig::Instance()->SetOdeTimeStep(0.001);  //ms

        // need to create a cell factory but don't want any intra stim, so magnitude
        // of stim is zero.
        c_vector<double,2> centre;
        centre(0) = 0.05; // cm
        centre(1) = 0.05; // cm
        BathCellFactory<2> cell_factory( 0.0, centre);

        BidomainWithBathProblem<2> bidomain_problem( &cell_factory );

        // Coverage
        TS_ASSERT(bidomain_problem.GetHasBath());

        TetrahedralMesh<2,2>* p_mesh = Load2dMeshAndSetCircularTissue<TetrahedralMesh<2,2> >(
            "mesh/test/data/2D_0_to_1mm_400_elements", 0.05, 0.05, 0.02);

        //boundary flux for Phi_e. -10e3 is under threshold, -14e3 crashes the cell model
        double boundary_flux = -11.0e3;
        double start_time = 0.5;
        double duration = 1.9; // of the stimulus, in ms

        HeartConfig::Instance()->SetElectrodeParameters(false,0, boundary_flux, start_time, duration);


        bidomain_problem.SetMesh(p_mesh);
        bidomain_problem.Initialise();

        bidomain_problem.Solve();

        Vec sol = bidomain_problem.GetSolution();
        ReplicatableVector sol_repl(sol);

        bool ap_triggered = false;
        /*
         * We are checking the last time step. This test will only make sure that an upstroke is triggered.
         * We ran longer simulation for 350 ms and a nice AP was observed.
         */

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            // test V = 0 for all bath nodes and that an AP is triggered in the tissue
            if (HeartRegionCode::IsRegionBath( p_mesh->GetNode(i)->GetRegion() )) // bath
            {
                TS_ASSERT_DELTA(sol_repl[2*i], 0.0, 1e-12);
            }
            else if (sol_repl[2*i] > 0.0)//at the last time step
            {
                ap_triggered = true;
            }
        }

        TS_ASSERT_EQUALS(bidomain_problem.mpElectrodes->mAreActive, false); // should be switched off by now..
        TS_ASSERT(ap_triggered);

        delete p_mesh;
    }

    void Test2dBathMultipleBathConductivities()
    {
        HeartConfig::Instance()->SetSimulationDuration(2.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath2dMultipleBathConductivities");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_2d");

        HeartConfig::Instance()->SetOdeTimeStep(0.001);  //ms ???

        std::set<unsigned> tissue_ids;
        tissue_ids.insert(0); // Same as default value defined in HeartConfig

        std::set<unsigned> bath_ids;
        bath_ids.insert(2);
        bath_ids.insert(3);
        bath_ids.insert(4);
        HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids);

        // need to create a cell factory but don't want any intra stim, so magnitude
        // of stim is zero.
        c_vector<double,2> centre;
        centre(0) = 0.05; // cm
        centre(1) = 0.05; // cm
        BathCellFactory<2> cell_factory( 0.0, centre);

        BidomainWithBathProblem<2> bidomain_problem( &cell_factory );

        DistributedTetrahedralMesh<2,2> mesh;

        mesh.ConstructRegularSlabMesh(0.05, 0.9, 0.9);

        // set the x<0.25 and x>0.75 regions as the bath region
        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            double x = iter->CalculateCentroid()[0];
            double y = iter->CalculateCentroid()[1];
            if ((x>0.3) && (x<0.6) && (y>0.3) && (y<0.6))
            {
                iter->SetAttribute(0);
            }
            else
            {
                if (y < 0.2)
                {
                    iter->SetAttribute(2);
                }
                else if (y < 0.7)
                {
                    iter->SetAttribute(3);
                }
                else if (y < 0.9)
                {
                    iter->SetAttribute(4);
                }
            }
        }

        std::map<unsigned, double> multiple_bath_conductivities;
        multiple_bath_conductivities[2] = 7.0;
        multiple_bath_conductivities[3] = 1.0;
        multiple_bath_conductivities[4] = 0.001;

        HeartConfig::Instance()->SetBathMultipleConductivities(multiple_bath_conductivities);

        double boundary_flux = -3.0e3;
        double start_time = 0.0;
        double duration = 1.0; // of the stimulus, in ms

        HeartConfig::Instance()->SetElectrodeParameters(false, 0, boundary_flux, start_time, duration);

        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();

        bidomain_problem.Solve();

        DistributedVector distributed_solution = bidomain_problem.GetSolutionDistributedVector();
        DistributedVector::Stripe voltage(distributed_solution, 0);

        /*
         * We are checking the last time step. This test will only make sure that an AP is triggered.
         */
        bool ap_triggered = false;

        for (DistributedVector::Iterator index = distributed_solution.Begin();
             index!= distributed_solution.End();
             ++index)
        {
            // test V = 0 for all bath nodes and that an AP is triggered in the tissue
            if (HeartRegionCode::IsRegionBath( mesh.GetNode(index.Global)->GetRegion() )) // bath
            {
                TS_ASSERT_DELTA(voltage[index], 0.0, 1e-12);
            }
            else if (voltage[index] > 0.0)//at the last time step
            {
                ap_triggered = true;
            }
        }

        TS_ASSERT(PetscTools::ReplicateBool(ap_triggered));
    }

    void TestProblemChecksUsingBathWithMultipleBathConductivities()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_400_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        std::set<unsigned> tissue_ids;
        tissue_ids.insert(0);

        std::set<unsigned> bath_ids;
        bath_ids.insert(1);
        bath_ids.insert(2); // non-default identifier!

        HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids);

        BathCellFactory<2> cell_factory( 0.0, Create_c_vector(0.0, 0.0) );
        BidomainProblem<2> bidomain_problem( &cell_factory ); // non-bath problem, despite specifying bath stuff above!
        bidomain_problem.SetMesh( &mesh );
        TS_ASSERT_THROWS_THIS( bidomain_problem.Initialise() , "User has set bath identifiers, but the BidomainProblem isn't expecting a bath. Did you mean to use BidomainProblem(..., true)? Or alternatively, BidomainWithBathProblem(...)?");
    }


    void Test2dBathGroundedElectrodeStimulusSwitchesOnOff()
    {
        // Total execution time is 5 ms. Electrodes are on in [1.0, 3.0]
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath2dGroundedOnOff");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_2d_grounded_on_off");
        HeartConfig::Instance()->SetOdeTimeStep(0.001);  //ms

        // need to create a cell factory but don't want any intra stim, so magnitude
        // of stim is zero.
        c_vector<double,2> centre;
        centre(0) = 0.05; // cm
        centre(1) = 0.05; // cm
        BathCellFactory<2> cell_factory( 0.0, centre);

        BidomainWithBathProblem<2> bidomain_problem( &cell_factory );

        TetrahedralMesh<2,2>* p_mesh = Load2dMeshAndSetCircularTissue<TetrahedralMesh<2,2> >(
            "mesh/test/data/2D_0_to_1mm_400_elements", 0.05, 0.05, 0.02);

        //boundary flux for Phi_e. -10e3 is under threshold, -14e3 crashes the cell model
        double boundary_flux = -11.0e3;
        double start_time = 1.0;
        double duration = 2.0; // of the stimulus, in ms

        HeartConfig::Instance()->SetElectrodeParameters( true, 0, boundary_flux,start_time, duration );


        bidomain_problem.SetMesh(p_mesh);
        bidomain_problem.Initialise();

        /*
         *  While t in [0.0, 1.0) electrodes are off
         */
        {
            HeartConfig::Instance()->SetSimulationDuration(0.5);  //ms
            bidomain_problem.Solve();

            /// \todo: we don't need a ReplicatableVector here. Every processor can check locally
            Vec sol = bidomain_problem.GetSolution();
            ReplicatableVector sol_repl(sol);

            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                // test phi_e close to 0 for all bath nodes since electrodes are off
                if (HeartRegionCode::IsRegionBath( p_mesh->GetNode(i)->GetRegion() )) // bath
                {
                    TS_ASSERT_DELTA(sol_repl[2*i+1], 0.0, 0.5);
                }
            }

            TS_ASSERT_EQUALS(bidomain_problem.mpElectrodes->mAreActive, false); // should be switched off by now..
        }


        /*
         *  At the end of the simulation AP has been triggered
         */
        {
            HeartConfig::Instance()->SetSimulationDuration(5.0);  //ms
            bidomain_problem.Solve();

            Vec sol = bidomain_problem.GetSolution();
            ReplicatableVector sol_repl(sol);

            bool ap_triggered = false;
            /*
             * We are checking the last time step. This test will only make sure that an upstroke is triggered.
             * We ran longer simulation for 350 ms and a nice AP was observed.
             */

            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                // test V = 0 for all bath nodes and that an AP is triggered in the tissue
                if (HeartRegionCode::IsRegionBath( p_mesh->GetNode(i)->GetRegion() )) // bath
                {
                    TS_ASSERT_DELTA(sol_repl[2*i], 0.0, 1e-12);
                }
                else if (sol_repl[2*i] > 0.0)//at the last time step
                {
                    ap_triggered = true;
                }
            }

            // Check that grounded electrode has been successfully removed and therefore phi_e !=0.
            // Nodes defining grounded electrode are 10, 21, 32, 43, 54, ... , 120
            TS_ASSERT_DELTA(sol_repl[21], -80.2794, 1e-3);
            TS_ASSERT_DELTA(sol_repl[43], -80.2794, 1e-3);
            TS_ASSERT_DELTA(sol_repl[241], -80.2794, 1e-3);

            TS_ASSERT_EQUALS(bidomain_problem.mpElectrodes->mAreActive, false); // should be switched off by now..
            TS_ASSERT(ap_triggered);
        }

        delete p_mesh;
    }


    void TestMatrixBasedAssembledBath(void)
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms

        // need to create a cell factory but don't want any intra stim, so magnitude
        // of stim is zero.
        c_vector<double,2> centre;
        centre(0) = 0.05;
        centre(1) = 0.05;
        BathCellFactory<2> cell_factory( 0.0, centre);

        //boundary flux for Phi_e
        double boundary_flux = -4e2;
        double duration = 0.2; //ms

        DistributedTetrahedralMesh<2,2>* p_mesh = Load2dMeshAndSetCircularTissue<DistributedTetrahedralMesh<2,2> >(
            "mesh/test/data/2D_0_to_1mm_400_elements", 0.05, 0.05, 0.02);

        ///////////////////////////////////////////////////////////////////
        // matrix based
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("BidomainBathMatrixBased");
        HeartConfig::Instance()->SetOutputFilenamePrefix("matrix_based");

        BidomainWithBathProblem<2> matrix_based_bido( &cell_factory );

        HeartConfig::Instance()->SetElectrodeParameters(true,0,boundary_flux, 0.0, duration);

        {
            Timer::Reset();

            matrix_based_bido.SetMesh(p_mesh);
            matrix_based_bido.Initialise();
            matrix_based_bido.Solve();

            Timer::Print("2D Matrix based");
        }

        ///////////////////////////////////////////////////////////////////
        // non matrix based
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("BidomainBathNonMatrixBased");
        HeartConfig::Instance()->SetOutputFilenamePrefix("non_matrix_based");

        BidomainWithBathProblem<2> non_matrix_based_bido( &cell_factory);

        {
            Timer::Reset();

            non_matrix_based_bido.SetMesh(p_mesh);
            non_matrix_based_bido.Initialise();
            non_matrix_based_bido.Solve();

            Timer::Print("2D non matrix based");
        }

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector matrix_based_solution = matrix_based_bido.GetSolutionDistributedVector();
        DistributedVector non_matrix_based_solution = non_matrix_based_bido.GetSolutionDistributedVector();

        DistributedVector::Stripe matrix_based_voltage(matrix_based_solution, 0);
        DistributedVector::Stripe non_matrix_based_voltage(non_matrix_based_solution, 0);

        DistributedVector::Stripe matrix_based_ex_pot(matrix_based_solution, 1);
        DistributedVector::Stripe non_matrix_based_ex_pot(non_matrix_based_solution, 1);

        for (DistributedVector::Iterator index = matrix_based_solution.Begin();
             index != matrix_based_solution.End();
             ++index)
        {
            TS_ASSERT_DELTA(matrix_based_voltage[index], non_matrix_based_voltage[index], 1e-7);
            //TS_ASSERT_DELTA(matrix_based_ex_pot[index], non_matrix_based_ex_pot[index], 1e-7);
            //std::cout << matrix_based_voltage[index] << std::endl;
        }

        delete p_mesh;
    }

    void TestArchivingBidomainProblemWithElectrodes(void)
    {
        std::string archive_dir = "BidomainWithElectrodesArchiving";

        // Create the mesh outside the save scope, so we can compare with the loaded version.
        TetrahedralMesh<2,2>* p_mesh = Load2dMeshAndSetCircularTissue<TetrahedralMesh<2,2> >(
            "mesh/test/data/2D_0_to_1mm_400_elements", 0.05, 0.05, 0.02);

        { // save
            HeartConfig::Instance()->SetSimulationDuration(3.0);  // ms
            HeartConfig::Instance()->SetOutputDirectory(archive_dir + "Output");
            HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_2d_fluxes");
            HeartConfig::Instance()->SetOdeTimeStep(0.001);  // ms

            // need to create a cell factory but don't want any intra stim.
            ZeroStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;

            BidomainWithBathProblem<2> bidomain_problem( &cell_factory );

            //boundary flux for Phi_e. -10e3 is under threshold, -14e3 crashes the cell model
            double boundary_flux = -11.0e3;
            double duration = 1.9; // of the stimulus, in ms

            HeartConfig::Instance()->SetElectrodeParameters(false,0,boundary_flux, 0.0, duration);

            bidomain_problem.SetMesh(p_mesh);
            bidomain_problem.Initialise();

            // Save using helper class
            CardiacSimulationArchiver<BidomainWithBathProblem<2> >::Save(bidomain_problem, archive_dir, false);
        }

        { // load
            AbstractCardiacProblem<2,2,2>* p_abstract_problem = CardiacSimulationArchiver<BidomainWithBathProblem<2> >::Load(archive_dir);

            // get the new mesh
            AbstractTetrahedralMesh<2,2>& r_mesh = p_abstract_problem->rGetMesh();
            TS_ASSERT_EQUALS(p_mesh->GetNumElements(), r_mesh.GetNumElements());

            // Check that the bath is in the right place
            for (unsigned i=0; i<r_mesh.GetNumElements(); i++)
            {
                double x = r_mesh.GetElement(i)->CalculateCentroid()[0];
                double y = r_mesh.GetElement(i)->CalculateCentroid()[1];
                if (sqrt((x-0.05)*(x-0.05) + (y-0.05)*(y-0.05)) > 0.02)
                {
                    TS_ASSERT(HeartRegionCode::IsRegionBath(r_mesh.GetElement(i)->GetUnsignedAttribute()));
                }
                else
                {
                    TS_ASSERT(HeartRegionCode::IsRegionTissue(r_mesh.GetElement(i)->GetUnsignedAttribute()));
                }
                // Compare mesh before & after
                TS_ASSERT_EQUALS(r_mesh.GetElement(i)->GetUnsignedAttribute(), p_mesh->GetElement(i)->GetUnsignedAttribute());
                TS_ASSERT_DELTA(x, p_mesh->GetElement(i)->CalculateCentroid()[0], 1e-12);
                TS_ASSERT_DELTA(y, p_mesh->GetElement(i)->CalculateCentroid()[1], 1e-12);
            }

            // Check that there's an exact correspondence between bath nodes and fake cells
            for (unsigned i=r_mesh.GetDistributedVectorFactory()->GetLow(); i<r_mesh.GetDistributedVectorFactory()->GetHigh(); i++)
            {
                TS_ASSERT_EQUALS(r_mesh.GetNode(i)->GetRegion(), p_mesh->GetNode(i)->GetRegion());
                FakeBathCell* p_fake = dynamic_cast<FakeBathCell*>(p_abstract_problem->GetTissue()->GetCardiacCell(i));
                if (HeartRegionCode::IsRegionBath(r_mesh.GetNode(i)->GetRegion()))
                {
                    TS_ASSERT(p_fake != NULL);
                }
                else
                {
                    TS_ASSERT(p_fake == NULL);
                }
            }

            // This should only generate action potential if the electrodes were correctly saved and restored.
            p_abstract_problem->Solve();

            Vec sol = p_abstract_problem->GetSolution();
            ReplicatableVector sol_repl(sol);

            bool ap_triggered = false;
            /*
             * We are checking the last time step. This test will only make sure that an upstroke is triggered.
             * We ran longer simulation for 350 ms and a nice AP was observed.
             */
            for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
            {
                // test V = 0 for all bath nodes and that an AP is triggered in the tissue
                if (HeartRegionCode::IsRegionBath(r_mesh.GetNode(i)->GetRegion()))
                {
                    TS_ASSERT_DELTA(sol_repl[2*i], 0.0, 1e-12);
                }
                else if (sol_repl[2*i] > 0.0) //at the last time step
                {
                    ap_triggered = true;
                }
            }

            // We can get away with the following line only because this is a friend class and test.
            boost::shared_ptr<Electrodes<2> > p_electrodes = static_cast<BidomainWithBathProblem<2>* >(p_abstract_problem)->mpElectrodes;

            TS_ASSERT_EQUALS(p_electrodes->mAreActive, false); // should be switched off by now..
            TS_ASSERT_EQUALS(ap_triggered, true);

            delete p_abstract_problem;
        }

        delete p_mesh;
    }

    void TestSettingElectrodesOnResumedSimulation(void)
    {
        std::string archive_dir = "BidomainWithElectrodesArchiving";

        // Create the mesh outside the save scope, so we can compare with the loaded version.
        TetrahedralMesh<2,2>* p_mesh = Load2dMeshAndSetCircularTissue<TetrahedralMesh<2,2> >(
            "mesh/test/data/2D_0_to_1mm_400_elements", 0.05, 0.05, 0.02);

        { // save
            HeartConfig::Instance()->SetSimulationDuration(3.0);  // ms
            HeartConfig::Instance()->SetOutputDirectory(archive_dir + "Output");
            HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_2d_fluxes");
            HeartConfig::Instance()->SetOdeTimeStep(0.001);  // ms

            // need to create a cell factory but don't want any intra stim.
            ZeroStimulusCellFactory<CellTenTusscher2006EpiFromCellMLBackwardEuler, 2> cell_factory;

            BidomainWithBathProblem<2> bidomain_problem( &cell_factory );

            bidomain_problem.SetMesh(p_mesh);
            bidomain_problem.Initialise();

            // Save using helper class
            CardiacSimulationArchiver<BidomainWithBathProblem<2> >::Save(bidomain_problem, archive_dir, false);
        }

        { // load
            AbstractCardiacProblem<2,2,2>* p_abstract_problem = CardiacSimulationArchiver<BidomainWithBathProblem<2> >::Load(archive_dir);

            // get the new mesh
            AbstractTetrahedralMesh<2,2>& r_mesh = p_abstract_problem->rGetMesh();
            TS_ASSERT_EQUALS(p_mesh->GetNumElements(), r_mesh.GetNumElements());

            //boundary flux for Phi_e. -10e3 is under threshold, -14e3 crashes the cell model
            double boundary_flux = -11.0e3;
            double duration = 1.9; // of the stimulus, in ms
            double start_time = 0.5; // of the stimulus, in ms

            HeartConfig::Instance()->SetElectrodeParameters(false,0,boundary_flux, start_time, duration);
            p_abstract_problem->SetElectrodes();

            // This should only generate action potential if the electrodes were correctly saved and restored.
            p_abstract_problem->Solve();

            Vec sol = p_abstract_problem->GetSolution();
            ReplicatableVector sol_repl(sol);

            bool ap_triggered = false;
            /*
             * We are checking the last time step. This test will only make sure that an upstroke is triggered.
             * We ran longer simulation for 350 ms and a nice AP was observed.
             */
            for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
            {
                // test V = 0 for all bath nodes and that an AP is triggered in the tissue
                if (HeartRegionCode::IsRegionBath(r_mesh.GetNode(i)->GetRegion()))
                {
                    TS_ASSERT_DELTA(sol_repl[2*i], 0.0, 1e-12);
                }
                else if (sol_repl[2*i] > 0.0) //at the last time step
                {
                    ap_triggered = true;
                }
            }

            // We can get away with the following line only because this is a friend class and test.
            boost::shared_ptr<Electrodes<2> > p_electrodes = static_cast<BidomainWithBathProblem<2>* >(p_abstract_problem)->mpElectrodes;

            TS_ASSERT_EQUALS(p_electrodes->mAreActive, false); // should be switched off by now..
            TS_ASSERT_EQUALS(ap_triggered, true);

            delete p_abstract_problem;
        }

        delete p_mesh;
    }

    void TestArchivingMeshFileWithAttributes()
    {
        std::string archive_dir = "TestArchivingMeshFileWithAttributes";

        { // save...
            HeartConfig::Instance()->SetSimulationDuration(0.01);  //ms
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_10_elements_with_three_attributes");
            HeartConfig::Instance()->SetOutputDirectory(archive_dir + "Output");
            HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

            std::set<unsigned> tissue_ids;
            tissue_ids.insert(0u); // (The default)
            std::set<unsigned> bath_ids;
            bath_ids.insert(1u); // (The default)
            bath_ids.insert(2u); // Some other type of bath
            HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids,bath_ids);

            std::map<unsigned, double> multiple_bath_conductivities;
            multiple_bath_conductivities[1] = 3.14;
            multiple_bath_conductivities[2] = 2.72;
            HeartConfig::Instance()->SetBathMultipleConductivities(multiple_bath_conductivities);

            ZeroStimulusCellFactory<CellLuoRudy1991FromCellML, 1> bidomain_cell_factory;
            BidomainWithBathProblem<1> bidomain_problem( &bidomain_cell_factory );
            bidomain_problem.Initialise();

            // Save using helper class
            CardiacSimulationArchiver<BidomainWithBathProblem<1> >::Save(bidomain_problem, archive_dir, false);
        }

        HeartConfig::Instance()->Reset(); // Forget these IDs/conductivities, they should be loaded from the archive.

        { // load...
            AbstractCardiacProblem<1,1,2>* p_abstract_problem = CardiacSimulationArchiver<BidomainWithBathProblem<1> >::Load(archive_dir);

            // Check the identifiers have made it
            std::set<unsigned> tissue_ids = HeartConfig::Instance()->rGetTissueIdentifiers();
            TS_ASSERT( tissue_ids.size() == 1u );
            TS_ASSERT( *(tissue_ids.begin()) == 0u );
            std::set<unsigned> bath_ids = HeartConfig::Instance()->rGetBathIdentifiers();
            TS_ASSERT( bath_ids.size() == 2u );
            TS_ASSERT( *bath_ids.begin() == 1u );
            TS_ASSERT( *(++bath_ids.begin()) == 2u );

            AbstractTetrahedralMesh<1,1>* p_mesh = &(p_abstract_problem->rGetMesh());

            /* Check the element ids have been loaded properly. The middle 4 elements are 'heart' elements
             * (region=0). The edges are different "flavours" of bath... */
            unsigned expected_element_regions[10]={ 2,1,1,0,0,0,0,1,1,2 };
            for (AbstractTetrahedralMesh<1,1>::ElementIterator iter = p_mesh->GetElementIteratorBegin();
                 iter != p_mesh->GetElementIteratorEnd();
                 ++iter)
            {
                unsigned element_index = iter->GetIndex();
                unsigned element_attribute = iter->GetUnsignedAttribute();
                switch ( expected_element_regions[element_index] )
                {
                case 0:
                    TS_ASSERT(HeartRegionCode::IsRegionTissue( element_attribute ));
                    break;
                case 1:
                    TS_ASSERT_EQUALS( expected_element_regions[element_index], 1u);
                    TS_ASSERT(HeartRegionCode::IsRegionBath( element_attribute ));
                    TS_ASSERT_DELTA(HeartConfig::Instance()->GetBathConductivity( element_attribute ), 3.14, 1e-9);
                    break;
                case 2:
                    TS_ASSERT_EQUALS( expected_element_regions[element_index], 2u);
                    TS_ASSERT(HeartRegionCode::IsRegionBath( element_attribute ));
                    TS_ASSERT_DELTA(HeartConfig::Instance()->GetBathConductivity( element_attribute ), 2.72, 1e-9);
                    break;
                default:
                    NEVER_REACHED;
                    break;
                };

            }
            /* ...so the middle 5 nodes should be heart nodes */
            char expected_node_regions[11]={ 'B', 'B', 'B',
                       'T', 'T', 'T', 'T', 'T',
                       'B','B','B'};
            for (unsigned i=0; i<10; i++)
            {
                if (p_mesh->GetDistributedVectorFactory()->IsGlobalIndexLocal(i))
                {
                    if (expected_node_regions[i] == 'T')
                    {
                         TS_ASSERT(HeartRegionCode::IsRegionTissue( p_mesh->GetNode(i)->GetRegion() ));
                    }
                    else
                    {
                         TS_ASSERT_EQUALS( expected_node_regions[i], 'B')
                         TS_ASSERT(HeartRegionCode::IsRegionBath( p_mesh->GetNode(i)->GetRegion() ));
                    }
                }
            }

            delete p_abstract_problem;
        }
    }

    void TestSwitchesOffAtCorrectTime()
    {
        // zero stim cell factory
        c_vector<double,2> centre;
        centre(0) = 0.05; // cm
        centre(1) = 0.05; // cm
        BathCellFactory<2> cell_factory( 0.0, centre);

        // boundary flux for Phi_e. -10e3 is under threshold, -14e3 crashes the cell model
        //
        // Will use printing dt = 0.01 in second run below, so choose start and end times of the
        // electrode which don't coincide with printing times
        double boundary_flux = -11.0e3;
        double start_time = 0.5;
        double duration = 2.001; // of the stimulus, in ms

        HeartConfig::Instance()->SetOutputDirectory("ElectrodesSwitchOffAtCorrectTime");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetSimulationDuration(5.0);  //ms

        //////////////////////////////////////////////////////
        // solve with printing_dt = 0.01
        //////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.001, 0.01, 0.01);  //ms

        BidomainWithBathProblem<2> bidomain_problem1( &cell_factory );
        TetrahedralMesh<2,2>* p_mesh1 = Load2dMeshAndSetCircularTissue<TetrahedralMesh<2,2> >(
           "mesh/test/data/2D_0_to_1mm_400_elements", 0.05, 0.05, 0.02);

        HeartConfig::Instance()->SetElectrodeParameters(false, 0, boundary_flux, start_time, duration);

        bidomain_problem1.SetMesh(p_mesh1);
        bidomain_problem1.PrintOutput(false);
        bidomain_problem1.Initialise();

        TS_ASSERT_THROWS_THIS(bidomain_problem1.Solve(), "Additional times are now deprecated.  Use only to check whether the given times are met: "
                                                          "e.g. Electrode events should only happen on printing steps.");
        delete p_mesh1;
    }

    void TestBidomainWithBathCanOutputVariables()
    {
        HeartConfig::Instance()->SetSimulationDuration(0.01);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_10_elements_with_two_attributes");
        HeartConfig::Instance()->SetOutputDirectory("BidomainBathOutputVariables");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer();

        std::vector<std::string> output_variables;
        output_variables.push_back("cytosolic_calcium_concentration");
        HeartConfig::Instance()->SetOutputVariables(output_variables);

        std::set<unsigned> tissue_ids;
        tissue_ids.insert(0); // Same as default value defined in HeartConfig
        std::set<unsigned> bath_ids;
        bath_ids.insert(1);
        HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        BidomainWithBathProblem<1> bidomain_problem(&cell_factory);

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        FileFinder calcium_results("BidomainBathOutputVariables/output/BidomainLR91_1d_cytosolic_calcium_concentration.dat",
                                   RelativeTo::ChasteTestOutput);

        TS_ASSERT_EQUALS(calcium_results.IsFile(), true);

        FileFinder reference_results("heart/test/data/BidomainBathOutputVariables/BidomainLR91_1d_cytosolic_calcium_concentration.dat",
                                     RelativeTo::ChasteSourceRoot);

        NumericFileComparison comparer(calcium_results, reference_results);
        TS_ASSERT_EQUALS(comparer.CompareFiles(), true);
    }
};

#endif /*TESTBIDOMAINWITHBATH_HPP_*/
