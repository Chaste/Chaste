/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef TESTCELLKILLERSNOTFORRELEASE_HPP_
#define TESTCELLKILLERSNOTFORRELEASE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "FixedDurationGenerationBasedCellCycleModelCellsGenerator.hpp"
#include "RadialSloughingCellKiller.hpp"
#include "OxygenBasedCellKiller.hpp"
#include "CellwiseData.hpp"
#include "TrianglesMeshReader.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/**
 * This class contains tests for methods on classes
 * inheriting from AbstractCellKiller that are not
 * yet ready for release.
 */
class TestCellKillersNotForRelease : public AbstractCellBasedTestSuite
{
public:

    void TestRadialSloughingCellKiller() throw(Exception)
    {
        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(-0.5,-0.5);

        // Get centre of mesh (we know it's at the origin, really)
        c_vector<double,2> centre(2);
        centre[0] = 0.0;
        centre[1] = 0.0;
        for (unsigned i=0; i< mesh.GetNumNodes(); i++)
        {
            centre += mesh.GetNode(i)->rGetLocation();
        }
        centre = centre/mesh.GetNumNodes();

        // Choose radius of cell killer
        double radius = 0.4;

        // Create cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Create cell killer and kill cells
        RadialSloughingCellKiller radial_cell_killer(&tissue, centre, radius);
        radial_cell_killer.TestAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            double r = norm_2(tissue.GetLocationOfCellCentre(*cell_iter) - centre);

            if (r > radius)
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
            }
        }

        // Now get rid of dead cells
        tissue.RemoveDeadCells();

        // Check that we are correctly left with cells inside the circle of death
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            double r = norm_2(tissue.GetLocationOfCellCentre(*cell_iter) - centre);
            TS_ASSERT_LESS_THAN_EQUALS(r, radius);
        }
    }

    void TestOxygenBasedCellKiller() throw(Exception)
    {
        // Set up
        TissueConfig::Instance()->SetHepaOneParameters();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double end_time = 1.0;
        unsigned num_timesteps = 100*(unsigned)end_time; // ensure the time step is not too small
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Before we can do anything with the cell killer, we need to set up CellwiseData
        std::vector<double> oxygen_concentration;

        // Set the oxygen concentration to be zero
        oxygen_concentration.push_back(0.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        OxygenBasedCellKiller<2> bad_cell_killer(&tissue);

        // Get a reference to the cells held in tissue
        std::list<TissueCell>& r_cells = tissue.rGetCells();

        // Reset cell types to STEM
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            cell_iter->SetCellProliferativeType(STEM);
        }

        TS_ASSERT_THROWS_NOTHING(OxygenBasedCellKiller<2> oxygen_based_cell_killer(&tissue));

        OxygenBasedCellKiller<2> oxygen_based_cell_killer(&tissue);

        TS_ASSERT_THROWS_NOTHING(oxygen_based_cell_killer.TestAndLabelSingleCellForApoptosis(*r_cells.begin()));

        // Check that a single cell reaches apoptosis
        TS_ASSERT(!r_cells.begin()->HasApoptosisBegun());
        r_cells.begin()->SetCellProliferativeType(APOPTOTIC);
        oxygen_based_cell_killer.TestAndLabelSingleCellForApoptosis(*r_cells.begin());

        TS_ASSERT(r_cells.begin()->HasApoptosisBegun());

        // Increment time to a time after death
        p_simulation_time->IncrementTimeOneStep();

        // Store 'locations' of cells which are not dead
        std::set< double > old_locations;
        for (std::list<TissueCell>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            if (!cell_iter->IsDead())
            {
                Node<2>* p_node = tissue.GetNodeCorrespondingToCell(*cell_iter);
                c_vector<double, 2> location = p_node->rGetLocation();
                old_locations.insert(location[0] + location[1]*1000);
            }
        }

        // Remove the dead cell
        tissue.RemoveDeadCells();

        // Check that dead cells are removed from the mesh
        std::set< double > new_locations;
        for (std::list<TissueCell>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            TS_ASSERT(!cell_iter->IsDead());
            Node<2>* p_node = tissue.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location = p_node->rGetLocation();
            new_locations.insert(location[0] + location[1]*1000);
        }

        TS_ASSERT(new_locations == old_locations);
        CellwiseData<2>::Destroy();
    }

    void TestArchivingOfRadialSloughingCellKiller() throw (Exception)
    {
        // Set up
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "radial_killer.arch";

        c_vector<double,2> centre(2);
        centre[0] = 0.1;
        centre[1] = 0.2;

        double radius = 0.4;

        {
            // Create an output archive
            RadialSloughingCellKiller cell_killer(NULL, centre, radius);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            RadialSloughingCellKiller* const p_cell_killer = &cell_killer;
            output_arch << p_cell_killer;

            TS_ASSERT_DELTA(p_cell_killer->GetCentre()[0], 0.1, 1e-9);
            TS_ASSERT_DELTA(p_cell_killer->GetCentre()[1], 0.2, 1e-9);
            TS_ASSERT_DELTA(p_cell_killer->GetRadius(), 0.4, 1e-9);
        }

        // Change centre and radius prior to restoring the cell killer
        centre[0] = 0.0;
        centre[1] = 0.0;
        radius = 0.0;

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            RadialSloughingCellKiller* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the sloughing properties correctly
            TS_ASSERT_DELTA(p_cell_killer->GetCentre()[0], 0.1, 1e-9);
            TS_ASSERT_DELTA(p_cell_killer->GetCentre()[1], 0.2, 1e-9);
            TS_ASSERT_DELTA(p_cell_killer->GetRadius(), 0.4, 1e-9);

            delete p_cell_killer;
        }
    }


    void TestArchivingOfOxygenBasedCellKiller() throw (Exception)
    {
        // Set up
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "oxygen_based_killer.arch";

        {
            // Create an output archive
            OxygenBasedCellKiller<2> cell_killer(NULL);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            OxygenBasedCellKiller<2>* const p_cell_killer = &cell_killer;

            p_cell_killer->SetHypoxicConcentration(0.3);

            output_arch << p_cell_killer;

            TS_ASSERT_DELTA(p_cell_killer->GetHypoxicConcentration(), 0.3, 1e-5);
       }

       {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            OxygenBasedCellKiller<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the sloughing properties correctly
            TS_ASSERT_DELTA(p_cell_killer->GetHypoxicConcentration(), 0.3, 1e-5);

            delete p_cell_killer;
        }
    }
};

#endif /*TESTCELLKILLERSNOTFORRELEASE_HPP_*/
