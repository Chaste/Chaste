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
#include "PottsMeshGenerator.hpp"

#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"


class TestPottsUpdateRules : public AbstractCellBasedTestSuite
{
public:

    void TestVolumeConstraintUpdateRuleMethods() throw (Exception)
    {
		// Create a simple 2D PottsMesh with 2 elements
		PottsMeshGenerator generator(4, 4, 1, 2, 2, 2);
		PottsMesh<2>* p_mesh = generator.GetMesh();

		// Create cells
		std::vector<CellPtr> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

		// Create cell population
		PottsBasedCellPopulation cell_population(*p_mesh, cells);

		// Create an update law system
		VolumeConstraintUpdateRule<2> volume_constraint;

		// Test get/set methods
		TS_ASSERT_DELTA(volume_constraint.GetDeformationEnergyParameter(), 0.5, 1e-12);
		TS_ASSERT_DELTA(volume_constraint.GetMatureCellTargetVolume(), 16.0, 1e-12);

		volume_constraint.SetDeformationEnergyParameter(0.5);
		volume_constraint.SetMatureCellTargetVolume(0.6);

		TS_ASSERT_DELTA(volume_constraint.GetDeformationEnergyParameter(), 0.5, 1e-12);
		TS_ASSERT_DELTA(volume_constraint.GetMatureCellTargetVolume(), 0.6, 1e-12);

		volume_constraint.SetDeformationEnergyParameter(1.0);
		volume_constraint.SetMatureCellTargetVolume(16.0);

		//TODO Add tests of other methods
	}

	void TestVolumeConstraintUpdateRuleArchiving() throw (Exception)
	{
		// TODO
	}

    void TestAdhesionUpdateRuleMethods() throw (Exception)
    {
    	// Create a simple 2D PottsMesh with 2 elements
		PottsMeshGenerator generator(4, 4, 1, 2, 2, 2);
		PottsMesh<2>* p_mesh = generator.GetMesh();

		// Create cells
		std::vector<CellPtr> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

		// Create cell population
		PottsBasedCellPopulation cell_population(*p_mesh, cells);

		// Create an update law system
		AdhesionUpdateRule<2> adhesion_update;

		// Test get/set methods
	 	TS_ASSERT_DELTA(adhesion_update.GetCellCellAdhesionEnergyParameter(), 0.1, 1e-12);
		TS_ASSERT_DELTA(adhesion_update.GetCellBoundaryAdhesionEnergyParameter(), 0.2, 1e-12);

		adhesion_update.SetCellCellAdhesionEnergyParameter(0.5);
		adhesion_update.SetCellBoundaryAdhesionEnergyParameter(0.6);

     	TS_ASSERT_DELTA(adhesion_update.GetCellCellAdhesionEnergyParameter(), 0.5, 1e-12);
		TS_ASSERT_DELTA(adhesion_update.GetCellBoundaryAdhesionEnergyParameter(), 0.6, 1e-12);

		adhesion_update.SetCellCellAdhesionEnergyParameter(0.1);
		adhesion_update.SetCellBoundaryAdhesionEnergyParameter(0.2);

		//TODO Add tests of other methods
    }

	void TestAdhesionUpdateRuleArchiving() throw (Exception)
	{
		// TODO
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
