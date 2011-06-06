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
#ifndef TESTFORCESNOTFORRELEASE_HPP_
#define TESTFORCESNOTFORRELEASE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "BuskeInteractionForce.hpp"
#include "BuskeCompressionForce.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"

class TestForcesNotForRelease : public AbstractCellBasedTestSuite
{
public:

    void noTestBuskeInteractionForceMethods() throw (Exception)
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        unsigned num_cells_depth = 1;
        unsigned num_cells_width = 2;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);
        cell_population.Update();

        // Create force
        BuskeInteractionForce<2> buske_interaction_force;

        // Test set/get methods
        TS_ASSERT_EQUALS(buske_interaction_force.GetUseCutOffLength(), false);
        TS_ASSERT_DELTA(buske_interaction_force.GetCutOffLength(), DBL_MAX, 1e-6);
        TS_ASSERT_DELTA(buske_interaction_force.GetAdhesionEnergyParameter(), 200, 1e-6);
        TS_ASSERT_DELTA(buske_interaction_force.GetDeformationEnergyParameter(), 4.0/3.0, 1e-6);

        buske_interaction_force.SetCutOffLength(1.5);
        buske_interaction_force.SetAdhesionEnergyParameter(1.0);
        buske_interaction_force.SetDeformationEnergyParameter(1.0);

        TS_ASSERT_EQUALS(buske_interaction_force.GetUseCutOffLength(), true);
        TS_ASSERT_DELTA(buske_interaction_force.GetCutOffLength(), 1.5, 1e-6);
        TS_ASSERT_DELTA(buske_interaction_force.GetAdhesionEnergyParameter(), 1.0, 1e-6);
        TS_ASSERT_DELTA(buske_interaction_force.GetDeformationEnergyParameter(), 1.0, 1e-6);

        buske_interaction_force.SetCutOffLength(DBL_MAX);
        buske_interaction_force.SetAdhesionEnergyParameter(200);
        buske_interaction_force.SetDeformationEnergyParameter(4.0/3.0);

        // Test node force calculation

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

		for (unsigned i=0; i<40; i++)
		{
			// Move nodes close together
			double seperation = 4.0 - (double)i/10.0;
			cell_population.GetNode(1)->rGetModifiableLocation()[0]=seperation;

			// Reset the vector of node forces
			node_forces.clear();
			for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
			{
				 node_forces.push_back(zero_vector<double>(2));
			}

			buske_interaction_force.AddForceContribution(node_forces, cell_population);

			// Test forces on nodes
			double analytical_force_magnitude = 100.0*M_PI*seperation;
			if (seperation < 2.0)
			{
				analytical_force_magnitude += 3.0/4.0*pow((2.0-seperation),1.5)*sqrt(0.5);
			}

			TS_ASSERT_DELTA(node_forces[0][0], -analytical_force_magnitude, 1e-4);
			TS_ASSERT_DELTA(node_forces[0][1], 0.0, 1e-4);

			TS_ASSERT_DELTA(node_forces[1][0], analytical_force_magnitude, 1e-4);
			TS_ASSERT_DELTA(node_forces[1][1], 0.0, 1e-4);
		}
    }

    void TestBuskeCompressionForceMethods() throw (Exception)
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        unsigned num_cells_depth = 1;
        unsigned num_cells_width = 2;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);
        cell_population.Update();

        // Create force
        BuskeCompressionForce<2> buske_compression_force;

        // Test set/get methods
        TS_ASSERT_DELTA(buske_compression_force.GetCompressionEnergyParameter(), 1.0, 1e-6);
        buske_compression_force.SetCompressionEnergyParameter(15.0);

        TS_ASSERT_DELTA(buske_compression_force.GetCompressionEnergyParameter(), 15.0, 1e-6);
        buske_compression_force.SetCompressionEnergyParameter(1.0);

        // Test node force calculation

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

		for (unsigned i=0; i<40; i++)
		{
			// Move nodes close together
			double seperation = 4.0 - (double)i/10.0;
			cell_population.GetNode(1)->rGetModifiableLocation()[0]=seperation;

			// Reset the vector of node forces
			node_forces.clear();
			for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
			{
				 node_forces.push_back(zero_vector<double>(2));
			}

			buske_compression_force.AddForceContribution(node_forces, cell_population);

			// Test forces on nodes
			double analytical_force_magnitude = 0.0;
			if (seperation < 2.0)
			{
				analytical_force_magnitude = M_PI/6.0/5.0*(5.0-M_PI/3.0*(4.0-pow((1.0-seperation/2.0),2.0)*(2.0-seperation/2.0)))
													      *(3.0/4.0*pow(seperation,2.0)-4.0*seperation+5.0);
			}

			TS_ASSERT_DELTA(node_forces[0][0], -analytical_force_magnitude, 1e-4);
			TS_ASSERT_DELTA(node_forces[0][1], 0.0, 1e-4);

			TS_ASSERT_DELTA(node_forces[1][0], analytical_force_magnitude, 1e-4);
			TS_ASSERT_DELTA(node_forces[1][1], 0.0, 1e-4);
		}
    }


    void TestBuskeCompressionForceWithMultipleCells() throw (Exception)
       {
           SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
           // Create triangle mesh
           std::vector<Node<2>*> nodes;
		   nodes.push_back(new Node<2>(0, true,  0.0, 0.0));
		   nodes.push_back(new Node<2>(1, true, 1.0, -1.0));
		   nodes.push_back(new Node<2>(2, true, 1.0, 1.0));

		   NodesOnlyMesh<2> mesh;
		   mesh.ConstructNodesWithoutMesh(nodes);

           // Create cells
           std::vector<CellPtr> cells;
           CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
           cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

           // Create a node-based cell population
           NodeBasedCellPopulation<2> cell_population(mesh, cells);
           cell_population.SetMechanicsCutOffLength(1.5);
           cell_population.Update();

           // Create force
           BuskeCompressionForce<2> buske_compression_force;

           // Initialise a vector of node forces
           std::vector<c_vector<double, 2> > node_forces;
           node_forces.reserve(cell_population.GetNumNodes());

   			for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
   			{
   				 node_forces.push_back(zero_vector<double>(2));
   			}

   			buske_compression_force.AddForceContribution(node_forces, cell_population);

   			// Test forces on nodes

   			// This node should only move in the x-direction
   			TS_ASSERT_DELTA(node_forces[0][0], -0.1302, 1e-4);
   			TS_ASSERT_DELTA(node_forces[0][1], 0.0, 1e-4);
   			TS_ASSERT_DELTA(node_forces[1][0], 0.0578, 1e-4);
   			TS_ASSERT_DELTA(node_forces[1][1], -0.0578, 1e-4);
   			TS_ASSERT_DELTA(node_forces[2][0], 0.0578, 1e-4);
   			TS_ASSERT_DELTA(node_forces[2][1], 0.0578, 1e-4);
       }


    void TestForceOutputParameters()
    {
        std::string output_directory = "TestNotForReleaseForcesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with BuskeInteractionForce
        BuskeInteractionForce<2> buske_interaction_force;
        buske_interaction_force.SetCutOffLength(1.5);
        TS_ASSERT_EQUALS(buske_interaction_force.GetIdentifier(), "BuskeInteractionForce-2");

        out_stream buske_interaction_force_parameter_file = output_file_handler.OpenOutputFile("buske_results.parameters");
        buske_interaction_force.OutputForceParameters(buske_interaction_force_parameter_file);
        buske_interaction_force_parameter_file->close();

        std::string buske_interaction_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + buske_interaction_force_results_dir + "buske_results.parameters notforrelease_cell_based/test/data/TestForcesNotForRelease/buske_results.parameters").c_str()), 0);

        // Test with BuskeCompressionForce
        BuskeCompressionForce<2> buske_compression_force;
        TS_ASSERT_EQUALS(buske_compression_force.GetIdentifier(), "BuskeCompressionForce-2");

        out_stream buske_compression_force_parameter_file = output_file_handler.OpenOutputFile("buske_compression_results.parameters");
        buske_compression_force.OutputForceParameters(buske_compression_force_parameter_file);
        buske_compression_force_parameter_file->close();

        std::string buske_compression_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + buske_compression_force_results_dir + "buske_compression_results.parameters notforrelease_cell_based/test/data/TestForcesNotForRelease/buske_compression_results.parameters").c_str()), 0);
    }

    void TestBuskeInteractionForceArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "buske_interaction_force_system.arch";

        {
            // Create a NodesOnlyMesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
            MutableMesh<2,2> generating_mesh;
            generating_mesh.ConstructFromMeshReader(mesh_reader);

            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(generating_mesh);

            // Set up SimulationTime
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

            // Create a cell population
            std::vector<CellPtr> cells;
            boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
                p_model->SetCellProliferativeType(STEM);

                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetBirthTime(-50.0);
                cells.push_back(p_cell);
            }

            NodeBasedCellPopulation<2> cell_population(mesh, cells);

            // Create a force object and set member variables
            BuskeInteractionForce<2> buske_interaction_force;
            buske_interaction_force.SetCutOffLength(1.7);
            buske_interaction_force.SetAdhesionEnergyParameter(12.0);
            buske_interaction_force.SetDeformationEnergyParameter(13.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize force object via pointer
            AbstractForce<2>* const p_buske_interaction_force = &buske_interaction_force;

            output_arch << p_buske_interaction_force;
        }

        {
            ArchiveLocationInfo::SetMeshPathname("mesh/test/data", "square_2_elements");

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<2>* p_buske_interaction_force;

            // Restore force object from the archive
            input_arch >> p_buske_interaction_force;

            // Test member variables
            TS_ASSERT_EQUALS(dynamic_cast<BuskeInteractionForce<2>*>(p_buske_interaction_force)->GetUseCutOffLength(), true);
            TS_ASSERT_DELTA(dynamic_cast<BuskeInteractionForce<2>*>(p_buske_interaction_force)->GetCutOffLength(), 1.7, 1e-6);
            TS_ASSERT_DELTA(dynamic_cast<BuskeInteractionForce<2>*>(p_buske_interaction_force)->GetAdhesionEnergyParameter(), 12.0, 1e-6);
            TS_ASSERT_DELTA(dynamic_cast<BuskeInteractionForce<2>*>(p_buske_interaction_force)->GetDeformationEnergyParameter(), 13.0, 1e-6);

            // Tidy up
            delete p_buske_interaction_force;
        }
    }

    void TestBuskeCompressionForceArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "buske_compression_force_system.arch";

        {
            // Create a NodesOnlyMesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
            MutableMesh<2,2> generating_mesh;
            generating_mesh.ConstructFromMeshReader(mesh_reader);

            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(generating_mesh);

            // Set up SimulationTime
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

            // Create a cell population
            std::vector<CellPtr> cells;
            boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
                p_model->SetCellProliferativeType(STEM);

                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetBirthTime(-50.0);
                cells.push_back(p_cell);
            }

            NodeBasedCellPopulation<2> cell_population(mesh, cells);

            // Create a force object and set member variables
            BuskeCompressionForce<2> buske_compression_force;
            buske_compression_force.SetCompressionEnergyParameter(14.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize force object via pointer
            AbstractForce<2>* const p_buske_compression_force = &buske_compression_force;

            output_arch << p_buske_compression_force;
        }

        {
            ArchiveLocationInfo::SetMeshPathname("mesh/test/data", "square_2_elements");

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<2>* p_buske_compression_force;

            // Restore force object from the archive
            input_arch >> p_buske_compression_force;

            // Test member variable
            TS_ASSERT_DELTA(dynamic_cast<BuskeCompressionForce<2>*>(p_buske_compression_force)->GetCompressionEnergyParameter(), 14.0, 1e-6);

            // Tidy up
            delete p_buske_compression_force;
        }
    }
};
#endif /*TESTFORCESNOTFORRELEASE_HPP_*/
