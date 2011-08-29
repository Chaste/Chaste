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

#ifndef TESTNODEBASEDCELLPOPULATIONWITHBUSKEUPDATE_HPP_
#define TESTNODEBASEDCELLPOPULATIONWITHBUSKEUPDATE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"

#include "TrianglesMeshReader.hpp"
#include "NodesOnlyMesh.hpp"

// The header file below must be included in any TEST(!) that uses Petsc
#include "PetscSetupAndFinalize.hpp"

class TestNodeBasedCellPopulationWithBuskeUpdate : public AbstractCellBasedTestSuite
{
public:
	void TestMethods()
	{
		// Create a simple mesh
		TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
		MutableMesh<2,2> generating_mesh;
		generating_mesh.ConstructFromMeshReader(mesh_reader);

		// Convert this to a NodesOnlyMesh
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(generating_mesh);

		std::vector<CellPtr> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

		// Create a cell population, with no ghost nodes at the moment
		NodeBasedCellPopulationWithBuskeUpdate<2> cell_population(mesh, cells);

		TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "NodeBasedCellPopulationWithBuskeUpdate-2");

		// Test NodeBasedCellPopulationWithBuskeUpdate::UpdateNodeLocations()

		// Make up some forces
		std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
		std::vector<c_vector<double, 2> > forces_on_nodes(cell_population.GetNumNodes());

		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
		{
			old_posns[i][0] = cell_population.GetNode(i)->rGetLocation()[0];
			old_posns[i][1] = cell_population.GetNode(i)->rGetLocation()[1];

			forces_on_nodes[i][0] = i*0.01;
			forces_on_nodes[i][1] = 2*i*0.01;
		}

		// Call method
		double time_step = 0.01;
		cell_population.UpdateNodeLocations(forces_on_nodes, time_step);

		// Check that node locations were correctly updated
		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
		{
			TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetLocation()[0], old_posns[i][0] +   i*0.01*0.01, 1e-9);
			TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetLocation()[1], old_posns[i][1] + 2*i*0.01*0.01, 1e-9);
		}

		// Test NodeBasedCellPopulationWithBuskeUpdate::OutputCellPopulationParameters()
		std::string output_directory = "TestNodeBasedCellPopulationWithBuskeUpdate";
		OutputFileHandler output_file_handler(output_directory, false);

		// Test that the cell population parameters are output correctly
		out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

		// Write cell population parameters to file
		cell_population.OutputCellPopulationParameters(parameter_file);
		parameter_file->close();

		// Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
		TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.parameters notforrelease_cell_based/test/data/TestNodeBasedCellPopulationWithBuskeUpdate/results.parameters").c_str()), 0);
	}

    void TestArchivingCellPopulation() throw (Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "node_based_cell_population.arch";
        ArchiveLocationInfo::SetMeshFilename("node_based_cell_population_mesh");

        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create a simple mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
            TetrahedralMesh<2,2> generating_mesh;
            generating_mesh.ConstructFromMeshReader(mesh_reader);

            // Convert this to a NodesOnlyMesh
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(generating_mesh);

            // Create cells
            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

            // Create a cell population
            NodeBasedCellPopulationWithBuskeUpdate<2>* const p_cell_population = new NodeBasedCellPopulationWithBuskeUpdate<2>(mesh, cells);

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // loop over them to run to time 0.0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                cell_iter != p_cell_population->End();
                ++cell_iter)
            {
                cell_iter->ReadyToDivide();
            }

            p_cell_population->SetMechanicsCutOffLength(1.5);

            // Create an output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write the cell population to the archive
            (*p_arch) << static_cast<const SimulationTime&>(*p_simulation_time);
            (*p_arch) << p_cell_population;

            // Tidy up
            SimulationTime::Destroy();
            delete p_cell_population;
        }

        {
            // Need to set up time
            unsigned num_steps = 10;

            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            NodeBasedCellPopulationWithBuskeUpdate<2>* p_cell_population;

            // Restore the cell population
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) >> *p_simulation_time;
            (*p_arch) >> p_cell_population;

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0;

            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(), (double)(counter), 1e-7);
                counter++;
            }

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);

            // Check the cell population has been restored
            TS_ASSERT_EQUALS(p_cell_population->rGetCells().size(), 5u);

            TS_ASSERT_DELTA(p_cell_population->GetMechanicsCutOffLength(), 1.5, 1e-6);

            // Check number of nodes
            TS_ASSERT_EQUALS(p_cell_population->GetNumNodes(), 5u);

            // Check some node positions
            TS_ASSERT_EQUALS(p_cell_population->GetNode(3)->GetIndex(), 3u);
            TS_ASSERT_EQUALS(p_cell_population->GetNode(4)->GetIndex(), 4u);

            TS_ASSERT_DELTA(p_cell_population->GetNode(3)->rGetLocation()[0], 0.0, 1e-9);
            TS_ASSERT_DELTA(p_cell_population->GetNode(3)->rGetLocation()[1], 1.0, 1e-9);
            TS_ASSERT_DELTA(p_cell_population->GetNode(4)->rGetLocation()[0], 0.5, 1e-9);
            TS_ASSERT_DELTA(p_cell_population->GetNode(4)->rGetLocation()[1], 0.5, 1e-9);

            // Check the member variables have been restored
            TS_ASSERT_DELTA(p_cell_population->GetMechanicsCutOffLength(), 1.5, 1e-9);

            // Tidy up
            delete p_cell_population;
        }
    }
};

#endif /*TESTNODEBASEDCELLPOPULATIONWITHBUSKEUPDATE_HPP_*/
