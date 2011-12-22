/*

Copyright (C) Fujitsu Laboratories of Europe, 2009-2010

*/

/*

Copyright (C) University of Oxford, 2005-2011

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



#ifndef TESTADAPTIVEBIDOMAINPROBLEMNIGHTLY_HPP_
#define TESTADAPTIVEBIDOMAINPROBLEMNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>
#include <fstream>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#ifdef CHASTE_ADAPTIVITY
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkTriangle.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetAttributes.h>
#include <vtkCellDerivatives.h>
#include <vtkDataSetToUnstructuredGridFilter.h>
#include <vtkCellDataToPointData.h>
#include <vtkExtractVectorComponents.h>
#include <vtkPointData.h>
#endif //CHASTE_ADAPTIVITY

#include "AdaptiveBidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "Hdf5DataReader.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "HeartEventHandler.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "HeartConfig.hpp"

class TestAdaptiveBidomainProblemNightly : public CxxTest::TestSuite
{
public:

    void tearDown()
    {
        HeartConfig::Reset();
    }

    void TestBidomain3dLong() throw (Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.1);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));
        HeartConfig::Instance()->SetSimulationDuration(4.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_0_to_1mm_6000_elements");
        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveBidomainProblem");
        HeartConfig::Instance()->SetOutputFilenamePrefix("coarse_slab");

        double target_error = 1.0;
        double sigma = 0.01;
        double max_edge_length = 0.04;
        double min_edge_length = 0.005;
        double gradation = 1.3;
        unsigned max_num_nodes = 1000;
        unsigned num_adaptive_sweeps = 5;

        HeartConfig::Instance()->SetAdaptivityParameters( target_error, sigma, max_edge_length, min_edge_length,
                                                          gradation, max_num_nodes, num_adaptive_sweeps );

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> bidomain_cell_factory(-600.0*1000);

        AdaptiveBidomainProblem bidomain_problem( &bidomain_cell_factory );
        bidomain_problem.PrintOutput(false);
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();

        Vec voltage=bidomain_problem.GetSolution();
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(voltage);

        /*
         * Test the top right node against the right one in the non-adaptive case,
         * comparing voltage, and then test all the nodes on the right hand
         * face of the cube against the top right one, comparing voltage.
         */
        bool need_initialisation = true;
        double probe_voltage=0;

        need_initialisation = true;

        // Test the RHF of the mesh
        for (AbstractTetrahedralMesh<3,3>::NodeIterator it = bidomain_problem.rGetMesh().GetNodeIteratorBegin();
             it != bidomain_problem.rGetMesh().GetNodeIteratorEnd();
             ++it)
        {
            if (it->GetPoint()[0] == 0.1)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.1cm is the other end of the mesh and where we want to
                //       to test the value of the nodes

                if (need_initialisation)
                {
                    probe_voltage = voltage_replicated[2*it->GetIndex()];
                    probe_voltage = 20.35;        // <- This is the voltage from the non-adaptive corner node
                    need_initialisation = false;
                }
                else
                {
                    // the voltage at the end face varies a little because
                    // of drift due to the orientation of the tets in the mesh,
                    // hence the tolerance of 0.2
                    TS_ASSERT_DELTA(voltage_replicated[2*it->GetIndex()], probe_voltage, 0.25);
                }

                // if a 1D simulation is run for 4ms on the 0_1mm_10elements mesh
                // the result at the end node is 20.0755
                TS_ASSERT_DELTA(voltage_replicated[2*it->GetIndex()], 20.0755, 1.3);
            }
        }

        TS_ASSERT_DELTA( bidomain_problem.mpAdaptiveMesh->GetNumNodes(), 75U, 35 );
        TS_ASSERT_DELTA( bidomain_problem.mpAdaptiveMesh->GetNumElements(), 225U, 75 );

#endif
    }

    void TestBidomain3dFromVtu() throw (Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.1);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));
        HeartConfig::Instance()->SetSimulationDuration(2.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/adaptivity/coarse_slab0020.vtu");
        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveBidomainProblem");
        HeartConfig::Instance()->SetOutputFilenamePrefix("coarse_slab_restore");

        double target_error = 1.0;
        double sigma = 0.01;
        double max_edge_length = 0.04;
        double min_edge_length = 0.005;
        double gradation = 1.3;
        unsigned max_num_nodes = 1000;
        unsigned num_adaptive_sweeps = 5;

        HeartConfig::Instance()->SetAdaptivityParameters( target_error, sigma, max_edge_length, min_edge_length,
                                                          gradation, max_num_nodes, num_adaptive_sweeps );

        ZeroStimulusCellFactory<CellLuoRudy1991FromCellML, 3> bidomain_cell_factory;

        AdaptiveBidomainProblem bidomain_problem( &bidomain_cell_factory );
        bidomain_problem.PrintOutput(false);
        bidomain_problem.LoadSimulationFromVtuFile();
        bidomain_problem.Solve();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();

        Vec voltage=bidomain_problem.GetSolution();
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(voltage);

        /*
         * Test the top right node against the right one in the non-adaptive case,
         * comparing voltage, and then test all the nodes on the right hand
         * face of the cube against the top right one, comparing voltage.
         */
        bool need_initialisation = true;
        double probe_voltage=0;

        need_initialisation = true;

        // Test the RHF of the mesh
        for (AbstractTetrahedralMesh<3,3>::NodeIterator it = bidomain_problem.rGetMesh().GetNodeIteratorBegin();
             it != bidomain_problem.rGetMesh().GetNodeIteratorEnd();
             ++it)
        {
            if (it->GetPoint()[0] == 0.1)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.1cm is the other end of the mesh and where we want to
                //       to test the value of the nodes

                if (need_initialisation)
                {
                    probe_voltage = voltage_replicated[2*it->GetIndex()];
                    probe_voltage = 20.35;        // <- This is the voltage from the non-adaptive corner node
                    need_initialisation = false;
                }
                else
                {
                    // the voltage at the end face varies a little because
                    // of drift due to the orientation of the tets in the mesh,
                    // hence the tolerance of 0.2
                    TS_ASSERT_DELTA(voltage_replicated[2*it->GetIndex()], probe_voltage, 0.25);
                }

                // if a 1D simulation is run for 4ms on the 0_1mm_10elements mesh
                // the result at the end node is 20.0755
                TS_ASSERT_DELTA(voltage_replicated[2*it->GetIndex()], 20.0755, 1.3);
            }
        }

        TS_ASSERT_DELTA( bidomain_problem.mpAdaptiveMesh->GetNumNodes(), 75U, 35 );
        TS_ASSERT_DELTA( bidomain_problem.mpAdaptiveMesh->GetNumElements(), 225U, 75 );

#endif
    }

    void TestBidomain3dLongWithNeumannBoundaryConditions() throw (Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.1);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));
        HeartConfig::Instance()->SetSimulationDuration(4.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_0_to_1mm_6000_elements");
        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveBidomainProblem");
        HeartConfig::Instance()->SetOutputFilenamePrefix("coarse_slab_neumann");

        double target_error = 1.0;
        double sigma = 0.01;
        double max_edge_length = 0.04;
        double min_edge_length = 0.005;
        double gradation = 1.3;
        unsigned max_num_nodes = 1000;
        unsigned num_adaptive_sweeps = 5;

        HeartConfig::Instance()->SetAdaptivityParameters( target_error, sigma, max_edge_length, min_edge_length,
                                                          gradation, max_num_nodes, num_adaptive_sweeps );

        ZeroStimulusCellFactory<CellLuoRudy1991FromCellML, 3> bidomain_cell_factory;

        AdaptiveBidomainProblem bidomain_problem( &bidomain_cell_factory );
        bidomain_problem.UseNeumannBoundaryCondition();
        bidomain_problem.SetNeumannStimulusMagnitudeAndDuration(4000.0, 0.5);
//        bidomain_problem.SetAdaptCriterion(0.8, 1e-6);

        bidomain_problem.PrintOutput(false);

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
        Vec voltage=bidomain_problem.GetSolution();
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(voltage);

        /*
         * Test the top right node against the right one in the non-adaptive case,
         * comparing voltage, and then test all the nodes on the right hand
         * face of the cube against the top right one, comparing voltage.
         */
        bool need_initialisation = true;
        double probe_voltage=0;

        need_initialisation = true;

        // Test the RHF of the mesh
        for (unsigned i = 0; i < bidomain_problem.rGetMesh().GetNumNodes(); i++)
        {
            if (bidomain_problem.rGetMesh().GetNode(i)->GetPoint()[0] == 0.1)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.1cm is the other end of the mesh and where we want to
                //       to test the value of the nodes

                if (need_initialisation)
                {
                    probe_voltage = 19.88;        // <- This is the voltage from the non-adaptive corner node
                    need_initialisation = false;
                }
                else
                {
                    // the voltage at the end face varies a little because
                    // of drift due to the orientation of the tets in the mesh,
                    // hence the tolerance of 0.2
                    TS_ASSERT_DELTA(voltage_replicated[2*i], probe_voltage, 0.25);
                }

                // if a non-adaptive simulation is run for 4ms on the 0_1mm_10elements mesh
                // the result at the end node is 18.9436
                TS_ASSERT_DELTA(voltage_replicated[2*i], 18.9436, 1.3);
            }
        }

        TS_ASSERT_DELTA( bidomain_problem.mpAdaptiveMesh->GetNumNodes(), 71U, 15 );
        TS_ASSERT_DELTA( bidomain_problem.mpAdaptiveMesh->GetNumElements(), 220U, 100 );
#endif
    }

    void TestBidomain3dLongWithNoAdapt() throw (Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.1);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));
        HeartConfig::Instance()->SetSimulationDuration(4.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_0_to_1mm_6000_elements");
        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveBidomainProblem");
        HeartConfig::Instance()->SetOutputFilenamePrefix("coarse_slab_neumann_no_adapt");

        ZeroStimulusCellFactory<CellLuoRudy1991FromCellML, 3> bidomain_cell_factory;

        AdaptiveBidomainProblem bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.DoNotAdaptMesh();
        bidomain_problem.UseNeumannBoundaryCondition();
        bidomain_problem.SetNeumannStimulusMagnitudeAndDuration(4000.0, 0.5);
        bidomain_problem.PrintOutput(false);

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
        Vec voltage=bidomain_problem.GetSolution();
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(voltage);

        /*
         * Test the top right node against the right one in the non-adaptive case,
         * comparing voltage, and then test all the nodes on the right hand
         * face of the cube against the top right one, comparing voltage.
         */
        bool need_initialisation = true;
        double probe_voltage=0;

        need_initialisation = true;

        // Test the RHF of the mesh
        for (unsigned i = 0; i < bidomain_problem.rGetMesh().GetNumNodes(); i++)
        {
            if (bidomain_problem.rGetMesh().GetNode(i)->GetPoint()[0] == 0.1)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.1cm is the other end of the mesh and where we want to
                //       to test the value of the nodes

                if (need_initialisation)
                {
                    probe_voltage = 19.88;
                    need_initialisation = false;
                }
                else
                {
                    // the voltage at the end face varies a little because
                    // of drift due to the orientation of the tets in the mesh,
                    // hence the tolerance of 0.2
                    TS_ASSERT_DELTA(voltage_replicated[2*i], probe_voltage, 0.25);
                }

                // if a non-adaptive simulation is run for 4ms on the 0_1mm_10elements mesh
                // the result at the end node is 18.9436
                TS_ASSERT_DELTA(voltage_replicated[2*i], 18.9436, 1.3);
            }
        }

        TS_ASSERT_EQUALS( bidomain_problem.mpAdaptiveMesh->GetNumNodes(), 1331U );
        TS_ASSERT_EQUALS( bidomain_problem.mpAdaptiveMesh->GetNumElements(), 6000U );
#endif
    }

};

#endif /*TESTADAPTIVEBIDOMAINPROBLEMNIGHTLY_HPP_*/
