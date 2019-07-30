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

#ifndef TESTCELLKILLERS_HPP_
#define TESTCELLKILLERS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellLabel.hpp"
#include "CellsGenerator.hpp"
#include "TargetedCellKiller.hpp"
#include "RandomCellKiller.hpp"
#include "ApoptoticCellKiller.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "IsolatedLabelledCellKiller.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "TrianglesMeshReader.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"
#include "SmartPointers.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * This class contains tests for methods on classes
 * inheriting from AbstractCellKiller.
 */
class TestCellKillers : public AbstractCellBasedTestSuite
{
public:

    void TestTargetedCellKiller()
    {
        // Set up singleton classes
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(32.0, 32);

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Get a reference to the cells held in cell population
        std::list<CellPtr>& r_cells = cell_population.rGetCells();

        // Create cell killer
        TargetedCellKiller<2> single_cell_killer(&cell_population, 1u);

        TS_ASSERT_EQUALS(single_cell_killer.GetIdentifier(), "TargetedCellKiller-2");

        // Check that some of the vector of cells reach apotosis
        single_cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        std::set<double> old_locations;

        std::list<CellPtr>::iterator cell_it = r_cells.begin();
        TS_ASSERT(!(*cell_it)->IsDead());
        ++cell_it;
        TS_ASSERT((*cell_it)->IsDead());
        ++cell_it;

        while (cell_it != r_cells.end())
        {
            TS_ASSERT(!(*cell_it)->IsDead());
            ++cell_it;
        }

        // Store 'locations' of cells which are not dead
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            if (!(*cell_iter)->IsDead())
            {
                Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
                c_vector<double, 2> location;
                location = p_node->rGetLocation();
                old_locations.insert(location[0] + location[1]*1000);
            }
        }

        // Remove dead cells
        cell_population.RemoveDeadCells();

        // Check that dead cells are removed from the mesh
        std::set< double > new_locations;
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS((*cell_iter)->IsDead(), false);
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location;
            location = p_node->rGetLocation();
            new_locations.insert(location[0] + location[1]*1000);
        }

        TS_ASSERT(new_locations == old_locations);
    }

    void TestRandomCellKiller()
    {
        // Set up singleton classes
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(32.0, 32);

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        double death_time = p_simulation_time->GetTime() + cells[0]->GetApoptosisTime();

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Get a reference to the cells held in cell population
        std::list<CellPtr>& r_cells = cell_population.rGetCells();

        // Check for bad probabilities being passed in
        TS_ASSERT_THROWS_THIS(RandomCellKiller<2> random_cell_killer(&cell_population, -0.1),
                             "Probability of death must be between zero and one");

        TS_ASSERT_THROWS_THIS(RandomCellKiller<2> random_cell_killer(&cell_population,  1.1),
                             "Probability of death must be between zero and one");

        // Create cell killer
        RandomCellKiller<2> random_cell_killer(&cell_population, 0.05);

        TS_ASSERT_EQUALS(random_cell_killer.GetIdentifier(), "RandomCellKiller-2");

        // Check that a single cell reaches apoptosis
        unsigned max_tries = 0;
        while (!(*r_cells.begin())->HasApoptosisBegun() && max_tries<99)
        {
            random_cell_killer.CheckAndLabelSingleCellForApoptosis(*r_cells.begin());
            max_tries++;
        }
        TS_ASSERT_DIFFERS(max_tries, 99u);
        TS_ASSERT_DIFFERS(max_tries, 0u);

        // Check that some of the vector of cells reach apotosis
        random_cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        std::set<double> old_locations;

        bool apoptosis_cell_found = false;
        std::list<CellPtr>::iterator cell_it = r_cells.begin();
        ++cell_it;
        while (cell_it != r_cells.end() && !apoptosis_cell_found)
        {
            if ((*cell_it)->HasApoptosisBegun())
            {
                apoptosis_cell_found = true;
            }
            ++cell_it;
        }

        TS_ASSERT_EQUALS(apoptosis_cell_found, true);

        // Increment time to a time after cell death
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_DELTA(p_simulation_time->GetTime(), 1.0, 1e-3);
        TS_ASSERT_DELTA(death_time, 0.25, 1e-3);
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(death_time+1.0, 1);
        p_simulation_time->IncrementTimeOneStep();

        // Store 'locations' of cells which are not dead
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
            cell_iter != r_cells.end();
            ++cell_iter)
        {
            if (!(*cell_iter)->IsDead())
            {
                Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
                c_vector<double, 2> location;
                location = p_node->rGetLocation();
                old_locations.insert(location[0] + location[1]*1000);
            }
        }

        // Remove dead cells
        cell_population.RemoveDeadCells();

        // Check that dead cells are removed from the mesh
        std::set< double > new_locations;
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
            cell_iter != r_cells.end();
            ++cell_iter)
        {
            TS_ASSERT_EQUALS((*cell_iter)->IsDead(), false);
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location;
            location = p_node->rGetLocation();
            new_locations.insert(location[0] + location[1]*1000);
        }

        TS_ASSERT(new_locations == old_locations);
    }

    void TestApoptoticCellKiller()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double end_time = 1.0;
        unsigned num_timesteps = 100*(unsigned)end_time; // ensure the time step is not too small
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        double lo_oxygen_concentration = 0.0;

        // Set the oxygen level for the cells
        for (unsigned index=0; index < cells.size(); index++)
        {
            cells[index]->GetCellData()->SetItem("oxygen", lo_oxygen_concentration);
        }

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        ApoptoticCellKiller<2> bad_cell_killer(&cell_population);

        // Get a reference to the cells held in cell population
        std::list<CellPtr>& r_cells = cell_population.rGetCells();

        // Reset each cell to have a StemCellProliferativeType
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        }

        ApoptoticCellKiller<2> oxygen_based_cell_killer(&cell_population);
        TS_ASSERT_EQUALS(oxygen_based_cell_killer.GetIdentifier(), "ApoptoticCellKiller-2");
        oxygen_based_cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Check that a single cell reaches apoptosis
        TS_ASSERT_EQUALS((*r_cells.begin())->HasApoptosisBegun(), false);

        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());
        (*r_cells.begin())->AddCellProperty(p_apoptotic_state);
        oxygen_based_cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        TS_ASSERT_EQUALS((*r_cells.begin())->HasApoptosisBegun(), true);

        // Increment time to a time after death
        p_simulation_time->IncrementTimeOneStep();

        // Store 'locations' of cells which are not dead
        std::set< double > old_locations;
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            if (!(*cell_iter)->IsDead())
            {
                Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
                c_vector<double, 2> location;
                location = p_node->rGetLocation();
                old_locations.insert(location[0] + location[1]*1000);
            }
        }

        // Remove the dead cell
        cell_population.RemoveDeadCells();

        // Check that dead cells are removed from the mesh
        std::set< double > new_locations;
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS((*cell_iter)->IsDead(), false);
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location;
            location = p_node->rGetLocation();
            new_locations.insert(location[0] + location[1]*1000);
        }

        TS_ASSERT(new_locations == old_locations);
    }

    void TestPlaneBasedCellKillerIn1d()
    {
        // Create 1D mesh
        unsigned num_cells = 14;
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(num_cells-1);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<1> cell_population(mesh, cells);

        // Create cell killer and kill cells
        c_vector<double, 1> point;
        point(0) = 10.0;
        PlaneBasedCellKiller<1> cell_killer(&cell_population, point, unit_vector<double>(1,0)); // x<10
        cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractCellPopulation<1>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            if (x > point(0))
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
            }
        }

        // Check that dead cells were correctly removed
        cell_population.RemoveDeadCells();

        for (AbstractCellPopulation<1>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            TS_ASSERT_LESS_THAN_EQUALS(x, point(0));
        }
    }

    void TestPlaneBasedCellKillerIn2d()
    {
        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(-0.25,-0.25);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Create cell killer and kill cells
        PlaneBasedCellKiller<2> cell_killer(&cell_population, zero_vector<double>(2), unit_vector<double>(2,1)); // y<0
        cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
            if (y > 0.0)
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
            }
        }

        cell_population.RemoveDeadCells();

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
            TS_ASSERT_LESS_THAN_EQUALS(y, 0.0);
        }
    }

    void TestPlaneBasedCellKillerIn3d()
    {
        // Create 3D mesh
        MutableMesh<3,3> mesh;
        mesh.ConstructCuboid(4, 5, 6);
        mesh.Translate(-2.0,-2.0, -2.0);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        // Create cell killer
        PlaneBasedCellKiller<3> cell_killer(&cell_population,  zero_vector<double>(3), unit_vector<double>(3,2)); // z<0
        cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double z = cell_population.GetLocationOfCellCentre(*cell_iter)[2];
            if (z > 0.0)
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
            }
        }

        cell_population.RemoveDeadCells();

        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double z = cell_population.GetLocationOfCellCentre(*cell_iter)[2];
            TS_ASSERT_LESS_THAN_EQUALS(z, 0.0);
        }
    }

    void TestIsolatedLabelledCellKiller()
    {
        // Create a non-vertex based cell population
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> non_vertex_mesh;
        non_vertex_mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> non_vertex_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> non_vertex_cells_generator;
        non_vertex_cells_generator.GenerateBasic(non_vertex_cells, non_vertex_mesh.GetNumNodes());

        MeshBasedCellPopulation<2> non_vertex_cell_population(non_vertex_mesh, non_vertex_cells);

        // Test that the correct exception is thrown when trying to create cell killer without a vertex-based cell population
        TS_ASSERT_THROWS_THIS(IsolatedLabelledCellKiller<2> cell_killer(&non_vertex_cell_population),
                              "IsolatedLabelledCellKiller only works with a VertexBasedCellPopulation.");

        // Create a vertex-based cell population
        HoneycombVertexMeshGenerator generator(4, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Label cells 0 and 1 (which are neighbours) and 3 (which is isolated
        boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
        cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
        cell_population.GetCellUsingLocationIndex(1)->AddCellProperty(p_label);
        cell_population.GetCellUsingLocationIndex(3)->AddCellProperty(p_label);
        TS_ASSERT_EQUALS(cell_population.GetCellPropertyRegistry()->Get<CellLabel>()->GetCellCount(), 3u);

        // Create cell killer
        IsolatedLabelledCellKiller<2> cell_killer(&cell_population);

        // No cells should yet have been killed
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
        }

        // Kill any isolated labelled cells
        cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Cell 3 should have been killed
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            bool should_be_dead = (cell_population.GetLocationIndexUsingCell(*cell_iter) == 3);
            TS_ASSERT_EQUALS(cell_iter->IsDead(), should_be_dead);
        }

        // Remove the dead cell
        cell_population.RemoveDeadCells();

        // Remove the label from cell 1
        cell_population.GetCellUsingLocationIndex(1)->RemoveCellProperty<CellLabel>();
        TS_ASSERT_EQUALS(cell_population.GetCellPropertyRegistry()->Get<CellLabel>()->GetCellCount(), 1u);

        // Again kill any isolated labelled cells
        cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Although isolated, cell 0 should not have been killed, since it is the only labelled cell
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
        }
    }

    void TestArchivingOfTargetedCellKiller()
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "single_cell_killer.arch";

        {
            // Create an output archive
             TargetedCellKiller<2> cell_killer(NULL, 1u);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            TargetedCellKiller<2>* const p_cell_killer = &cell_killer;
            output_arch << p_cell_killer;

            TS_ASSERT_EQUALS(p_cell_killer->GetTargetIndex(), 1u);
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            TargetedCellKiller<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the Target Cell correctly
            TS_ASSERT_EQUALS(p_cell_killer->GetTargetIndex(), 1u);

            delete p_cell_killer;
       }
    }

    void TestArchivingOfRandomCellKiller()
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "random_killer.arch";

        {
            // Create an output archive
            RandomCellKiller<2> cell_killer(NULL, 0.134);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            RandomCellKiller<2>* const p_cell_killer = &cell_killer;
            output_arch << p_cell_killer;

            TS_ASSERT_DELTA(p_cell_killer->GetDeathProbabilityInAnHour(), 0.134, 1e-9);
       }

       {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            RandomCellKiller<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the probability correctly
            TS_ASSERT_DELTA(p_cell_killer->GetDeathProbabilityInAnHour(), 0.134, 1e-9);

            delete p_cell_killer;
        }
    }

    void TestArchivingOfApoptoticCellKiller()
    {
        // Set up
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "oxygen_based_killer.arch";
        {
            // Create an output archive
            ApoptoticCellKiller<2> cell_killer(NULL);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            ApoptoticCellKiller<2>* const p_cell_killer = &cell_killer;

            output_arch << p_cell_killer;
       }

       {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            ApoptoticCellKiller<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            TS_ASSERT(p_cell_killer != NULL);

            // Tidy up
            delete p_cell_killer;
        }
    }

    void TestArchivingOfPlaneBasedCellKiller()
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "region_based_killer.arch";

        {
            // Create an output archive

            PlaneBasedCellKiller<2> cell_killer(NULL, zero_vector<double>(2), unit_vector<double>(2,1));

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            PlaneBasedCellKiller<2>* const p_cell_killer = &cell_killer;
            output_arch << p_cell_killer;

            TS_ASSERT_EQUALS(p_cell_killer->rGetPointOnPlane()[0], 0.0);
            TS_ASSERT_EQUALS(p_cell_killer->rGetPointOnPlane()[1], 0.0);
            TS_ASSERT_EQUALS(p_cell_killer->rGetNormalToPlane()[0], 0.0);
            TS_ASSERT_EQUALS(p_cell_killer->rGetNormalToPlane()[1], 1.0);
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            PlaneBasedCellKiller<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the region properties correctly
            TS_ASSERT_EQUALS(p_cell_killer->rGetPointOnPlane()[0], 0.0);
            TS_ASSERT_EQUALS(p_cell_killer->rGetPointOnPlane()[1], 0.0);
            TS_ASSERT_EQUALS(p_cell_killer->rGetNormalToPlane()[0], 0.0);
            TS_ASSERT_EQUALS(p_cell_killer->rGetNormalToPlane()[1], 1.0);

            delete p_cell_killer;
        }
    }

    void TestArchivingOfIsolatedLabelledCellKiller()
    {
        // Set up singleton classes
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "isolated_killer.arch";
        ArchiveLocationInfo::SetMeshFilename("vertex_mesh_isolated");

        {
            HoneycombVertexMeshGenerator generator(4,4);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            IsolatedLabelledCellKiller<2> cell_killer(&cell_population);

            // Create an output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Serialize via pointer
            IsolatedLabelledCellKiller<2>* const p_cell_killer = &cell_killer;
            (*p_arch) << p_cell_killer;
        }

        {
            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            IsolatedLabelledCellKiller<2>* p_cell_killer;

            // Restore from the archive
            (*p_arch) >> p_cell_killer;

            TS_ASSERT(p_cell_killer != NULL);
            TS_ASSERT(p_cell_killer->GetCellPopulation() != NULL);

            // Tidy up
            delete p_cell_killer->mpCellPopulation;
            delete p_cell_killer;
        }
    }

    void TestCellKillersOutputParameters()
    {
        std::string output_directory = "TestCellKillersOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with TargetedCellKiller
        TargetedCellKiller<2> targeted_cell_killer(NULL, 1u);
        TS_ASSERT_EQUALS(targeted_cell_killer.GetIdentifier(), "TargetedCellKiller-2");

        out_stream targeted_cell_killer_parameter_file = output_file_handler.OpenOutputFile("targeted_results.parameters");
        targeted_cell_killer.OutputCellKillerParameters(targeted_cell_killer_parameter_file);
        targeted_cell_killer_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code
            FileFinder generated = output_file_handler.FindFile("targeted_results.parameters");
            FileFinder reference("cell_based/test/data/TestCellKillers/targeted_results.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with RandomCellKiller
        RandomCellKiller<2> random_cell_killer(NULL, 0.01);
        TS_ASSERT_EQUALS(random_cell_killer.GetIdentifier(), "RandomCellKiller-2");

        out_stream random_cell_killer_parameter_file = output_file_handler.OpenOutputFile("random_results.parameters");
        random_cell_killer.OutputCellKillerParameters(random_cell_killer_parameter_file);
        random_cell_killer_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code.
            FileFinder generated = output_file_handler.FindFile("random_results.parameters");
            FileFinder reference("cell_based/test/data/TestCellKillers/random_results.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with ApoptoticCellKiller
        ApoptoticCellKiller<2> apop_cell_killer(NULL);
        TS_ASSERT_EQUALS(apop_cell_killer.GetIdentifier(), "ApoptoticCellKiller-2");

        out_stream apop_cell_killer_parameter_file = output_file_handler.OpenOutputFile("apop_results.parameters");
        apop_cell_killer.OutputCellKillerParameters(apop_cell_killer_parameter_file);
        apop_cell_killer_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code
            FileFinder generated = output_file_handler.FindFile("apop_results.parameters");
            FileFinder reference("cell_based/test/data/TestCellKillers/apop_results.parameters", RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with PlaneBasedCellKiller
        PlaneBasedCellKiller<2> region_cell_killer(NULL, zero_vector<double>(2), unit_vector<double>(2,1)); // y<0;
        TS_ASSERT_EQUALS(region_cell_killer.GetIdentifier(), "PlaneBasedCellKiller-2");

        out_stream region_cell_killer_parameter_file = output_file_handler.OpenOutputFile("region_results.parameters");
        region_cell_killer.OutputCellKillerParameters(region_cell_killer_parameter_file);
        region_cell_killer_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code
            FileFinder generated = output_file_handler.FindFile("region_results.parameters");
            FileFinder reference("cell_based/test/data/TestCellKillers/region_results.parameters", RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with IsolatedLabelledCellKiller
        HoneycombVertexMeshGenerator generator(4,4);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        IsolatedLabelledCellKiller<2> isolated_labelled_cell_killer(&cell_population);
        TS_ASSERT_EQUALS(isolated_labelled_cell_killer.GetIdentifier(), "IsolatedLabelledCellKiller-2");

        out_stream isolated_labelled_cell_killer_parameter_file = output_file_handler.OpenOutputFile("isolated_labelled_results.parameters");
        isolated_labelled_cell_killer.OutputCellKillerParameters(isolated_labelled_cell_killer_parameter_file);
        isolated_labelled_cell_killer_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code
            FileFinder generated = output_file_handler.FindFile("isolated_labelled_results.parameters");
            FileFinder reference("cell_based/test/data/TestCellKillers/isolated_labelled_results.parameters", RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif /*TESTCELLKILLERS_HPP_*/
