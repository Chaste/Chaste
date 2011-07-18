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
#ifndef TESTPOTTSUPDATERULES_HPP_
#define TESTPOTTSUPDATERULES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractPottsUpdateRule.hpp"
#include "VolumeConstraintUpdateRule.hpp"
#include "AdhesionUpdateRule.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CellwiseData.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"

#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"


class TestPottsUpdateRules : public AbstractCellBasedTestSuite
{
public:

    void TestVolumeConstraintUpdateRuleMethods() throw (Exception)
    {
//		unsigned cells_across = 7;
//		unsigned cells_up = 5;
//		unsigned thickness_of_ghost_layer = 0;
//
//		SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
//
//		HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
//		MutableMesh<2,2>* p_mesh = generator.GetMesh();
//		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
//
//		// Create cells
//		std::vector<CellPtr> cells;
//		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//		cells_generator.GenerateBasic(cells, location_indices.size(), location_indices);
//
//		boost::shared_ptr<AbstractCellProperty> p_label(new CellLabel);
//		for (unsigned i=0; i<cells.size(); i++)
//		{
//			cells[i]->SetBirthTime(-10);
//			cells[i]->AddCellProperty(p_label);
//		}
//
//		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
//
//		// Set up cellwise data and associate it with the cell population
//		CellwiseData<2>* p_data = CellwiseData<2>::Instance();
//		p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
//		p_data->SetCellPopulation(&cell_population);
//
//		for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
//		{
//			double x = p_mesh->GetNode(i)->rGetLocation()[0];
//			p_data->SetValue(x/50.0, p_mesh->GetNode(i)->GetIndex());
//		}
//
//		ChemotacticForce<2> chemotactic_force;
//
//		// Initialise a vector of new node forces
//		std::vector<c_vector<double, 2> > node_forces;
//		node_forces.reserve(cell_population.GetNumNodes());
//
//		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
//		{
//			 node_forces.push_back(zero_vector<double>(2));
//		}
//		chemotactic_force.AddForceContribution(node_forces, cell_population);
//
//		for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
//			 cell_iter != cell_population.End();
//			 ++cell_iter)
//		{
//			unsigned index = cell_population.GetLocationIndexUsingCell(*cell_iter);
//			double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
//			double c = x/50;
//			double norm_grad_c = 1.0/50.0;
//			double force_magnitude = chemotactic_force.GetChemotacticForceMagnitude(c, norm_grad_c);
//
//			// Fc = force_magnitude*(1,0), Fspring = 0
//			TS_ASSERT_DELTA(node_forces[index][0], force_magnitude, 1e-4);
//			TS_ASSERT_DELTA(node_forces[index][1], 0.0, 1e-4);
//		}
	}

	void TestVolumeConstraintUpdateRuleArchiving() throw (Exception)
	{
//		// Set up
//		OutputFileHandler handler("archive", false);    // don't erase contents of folder
//		std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemotaxis_spring_system.arch";
//
//		{
//			// Create mesh
//			TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
//			MutableMesh<2,2> mesh;
//			mesh.ConstructFromMeshReader(mesh_reader);
//
//			// SimulationTime is usually set up by a CellBasedSimulation
//			SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
//
//			// Create cells
//			std::vector<CellPtr> cells;
//			CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//			cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
//			for (unsigned i=0; i<cells.size(); i++)
//			{
//				cells[i]->SetBirthTime(-50);
//			}
//
//			// Create cell population
//			MeshBasedCellPopulation<2> cell_population(mesh, cells);
//
//			// Create force
//			ChemotacticForce<2> chemotactic_force;
//
//			// Serialize force via pointer
//			std::ofstream ofs(archive_filename.c_str());
//			boost::archive::text_oarchive output_arch(ofs);
//
//			ChemotacticForce<2>* const p_chemotactic_force = &chemotactic_force;
//			output_arch << p_chemotactic_force;
//		}
//
//		{
//			ArchiveLocationInfo::SetMeshPathname("mesh/test/data/", "square_2_elements");
//
//			// Create an input archive
//			std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
//			boost::archive::text_iarchive input_arch(ifs);
//
//			// Restore force from the archive
//			ChemotacticForce<2>* p_chemotactic_force;
//			input_arch >> p_chemotactic_force;
//
//			///\todo test something here, for example that member variables have been correctly archived
//
//			// Tidy up
//			delete p_chemotactic_force;
//		}
	}

	void TestUpdateRuleOutputUpdateRuleInfo()
	{
		std::string output_directory = "TestPottsUpdateRulesOutputParameters";
		OutputFileHandler output_file_handler(output_directory, false);

		// Test with VolumeConstraintUpdateRule
		VolumeConstraintUpdateRule<2> volume_constraint;
		//volume_constraint.SetSomething(1.5);
		TS_ASSERT_EQUALS(volume_constraint.GetIdentifier(), "VolumeConstraintUpdateRule-2");

		out_stream volume_constraint_parameter_file = output_file_handler.OpenOutputFile("volume_constraint_results.parameters");
		volume_constraint.OutputUpdateRuleInfo(volume_constraint_parameter_file);
		volume_constraint_parameter_file->close();

		std::string volume_constraint_results_dir = output_file_handler.GetOutputDirectoryFullPath();
		TS_ASSERT_EQUALS(system(("diff " + volume_constraint_results_dir + "volume_constraint_results.parameters notforrelease_cell_based/test/data/TestUpdateRules/volume_constraint_results.parameters").c_str()), 0);

		// Test with VolumeConstraintUpdateRule
		AdhesionUpdateRule<2> adhesion_update;
		//adhesion_update.SetSomething(1.5);
		TS_ASSERT_EQUALS(adhesion_update.GetIdentifier(), "AdhesionUpdateRule-2");

		out_stream adhesion_update_parameter_file = output_file_handler.OpenOutputFile("adhesion_update_results.parameters");
		adhesion_update.OutputUpdateRuleInfo(adhesion_update_parameter_file);
		adhesion_update_parameter_file->close();

		std::string adhesion_update_results_dir = output_file_handler.GetOutputDirectoryFullPath();
		TS_ASSERT_EQUALS(system(("diff " + adhesion_update_results_dir + "adhesion_update_results.parameters notforrelease_cell_based/test/data/TestUpdateRules/adhesion_update_results.parameters").c_str()), 0);
	}


    void TestIncompatibleUpdateRules() throw (Exception)
    {
        // Create a NodeBasedCellPopulation
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 10;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double x = (double)(i);
            double y = (double)(i);
            nodes.push_back(new Node<2>(i, true, x, y));
        }

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, num_nodes);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Test that VolumeConstraintUpdateRule throws the correct exception
        VolumeConstraintUpdateRule<2> volume_constraint;
        TS_ASSERT_THROWS_THIS(volume_constraint.EvaluateHamiltonianContribution(0u, 1u, cell_population),
                "VolumeConstraintUpdateRule is to be used with a PottsBasedCellPopulation only");

        // Test that AdhesionUpdateRule throws the correct exception
        AdhesionUpdateRule<2> adhesion_update;
        TS_ASSERT_THROWS_THIS(adhesion_update.EvaluateHamiltonianContribution(0u, 1u, cell_population),
                "AdhesionUpdateRule is to be used with a PottsBasedCellPopulation only");

        // When the node-only mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif /*TESTPOTTSUPDATERULES_HPP_*/
