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

#ifndef TESTNODEBASEDCELLPOPULATIONWITHBUSKEUPDATE_HPP_
#define TESTNODEBASEDCELLPOPULATIONWITHBUSKEUPDATE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

#include "TrianglesMeshReader.hpp"
#include "NodesOnlyMesh.hpp"
#include "FileComparison.hpp"

// The header file below must be included in any TEST(!) that uses PETSc
#include "PetscSetupAndFinalize.hpp"

class TestNodeBasedCellPopulationWithBuskeUpdate : public AbstractCellBasedTestSuite
{
public:
    void TestMethods()
    {
        EXIT_IF_PARALLEL;    // NodeBasedCellPopulationWithBuskeUpdate doesn't work in parallel.

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.2);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population, with no ghost nodes at the moment
        NodeBasedCellPopulationWithBuskeUpdate<2> cell_population(mesh, cells);
        cell_population.Update();

        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "NodeBasedCellPopulationWithBuskeUpdate-2");

        // Test NodeBasedCellPopulationWithBuskeUpdate::UpdateNodeLocations()

        // Make up some forces
        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());

        for (NodesOnlyMesh<2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            unsigned i = mesh.SolveNodeMapping(node_index);

            c_vector<double, 2> force;
            old_posns[i][0] = node_iter->rGetLocation()[0];
            old_posns[i][1] = node_iter->rGetLocation()[1];

            force[0] = i*0.01;
            force[1] = 2*i*0.01;

            node_iter->ClearAppliedForce();
            node_iter->AddAppliedForceContribution(force);
        }

        // Call method
        double time_step = 0.01;

        cell_population.UpdateNodeLocations(time_step);

        // Check that node locations were correctly updated
        for (NodesOnlyMesh<2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            unsigned i = mesh.SolveNodeMapping(node_iter->GetIndex());
            TS_ASSERT_DELTA(node_iter->rGetLocation()[0], old_posns[i][0] +   i*0.01*0.01, 1e-9);
            TS_ASSERT_DELTA(node_iter->rGetLocation()[1], old_posns[i][1] + 2*i*0.01*0.01, 1e-9);
        }

        // Test NodeBasedCellPopulationWithBuskeUpdate::OutputCellPopulationParameters()
        std::string output_directory = "TestNodeBasedCellPopulationWithBuskeUpdate";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

        // Write cell population parameters to file
        cell_population.OutputCellPopulationParameters(parameter_file);
        parameter_file->close();

        // Compare the generated file in test output with a reference copy in the source code.
        FileFinder generated = output_file_handler.FindFile("results.parameters");
        FileFinder reference("cell_based/test/data/TestNodeBasedCellPopulationWithBuskeUpdate/results.parameters",
                RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated, reference);
        TS_ASSERT(comparer.CompareFiles());
    }

    void TestArchivingCellPopulation()
    {
        EXIT_IF_PARALLEL;    // Cell-based archiving doesn't yet work in parallel.

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
            mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

            // Create cells
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
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

            // Create an output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write the cell population to the archive
            (*p_arch) << static_cast<const SimulationTime&>(*p_simulation_time);
            (*p_arch) << p_cell_population;

            // Avoid memory leak
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
