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

#ifndef TESTSLOUGHINGCELLKILLERS_HPP_
#define TESTSLOUGHINGCELLKILLERS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "SloughingCellKiller.hpp"
#include "RadialSloughingCellKiller.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "TrianglesMeshReader.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSloughingCellKillers : public AbstractCellBasedTestSuite
{
public:

    void TestSloughingCellKillerTopAndSides()
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
        SloughingCellKiller<2> sloughing_cell_killer(&cell_population, 0.5, true, 0.5);

        sloughing_cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        TS_ASSERT_EQUALS(sloughing_cell_killer.GetIdentifier(), "SloughingCellKiller-2");

        // Check that cells were labelled for death correctly
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

            if ((x<0) || (x>0.5) || (y>0.5))
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
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

            TS_ASSERT_LESS_THAN_EQUALS(x, 0.5);
            TS_ASSERT_LESS_THAN_EQUALS(y, 0.5);
        }
    }

    void TestSloughingCellKillerTopOnly()
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
        SloughingCellKiller<2> sloughing_cell_killer(&cell_population, 0.5);

        sloughing_cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
            if (y > 0.5)
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
            TS_ASSERT_LESS_THAN_EQUALS(y, 0.5);
        }
    }

    void TestSloughingCellKillerIn1d()
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

        // Set the crypt length so that 2 cells should be sloughed off
        double crypt_length = 12.5;

        // Create cell killer and kill cells
        SloughingCellKiller<1> sloughing_cell_killer(&cell_population, crypt_length);
        sloughing_cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractCellPopulation<1>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            if (x > crypt_length)
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
            TS_ASSERT_LESS_THAN_EQUALS(x, crypt_length);
        }
    }

    void TestSloughingCellKillerIn3d()
    {
        // Create 3D mesh
        MutableMesh<3,3> mesh;
        mesh.ConstructCuboid(4, 5, 6);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        // Create cell killer
        SloughingCellKiller<3> sloughing_cell_killer(&cell_population, 1.0); // number is irrelevent as long as its positive.

        // Check that an exception is thrown, as this method is not yet implemented in 3D
        TS_ASSERT_THROWS_THIS(sloughing_cell_killer.CheckAndLabelCellsForApoptosisOrDeath(), "SloughingCellKiller is not yet implemented in 3D");
    }

    void TestArchivingOfSloughingCellKiller()
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sloughing_killer.arch";

        {
            // Create an output archive
            SloughingCellKiller<2> cell_killer(NULL, 10.0, true, 5.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            SloughingCellKiller<2>* const p_cell_killer = &cell_killer;
            output_arch << p_cell_killer;

            TS_ASSERT_EQUALS(p_cell_killer->GetSloughSides(), true);
            TS_ASSERT_DELTA(p_cell_killer->GetSloughHeight(), 10.0, 1e-9);
            TS_ASSERT_DELTA(p_cell_killer->GetSloughWidth(), 5.0, 1e-9);
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SloughingCellKiller<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the sloughing properties correctly
            TS_ASSERT_EQUALS(p_cell_killer->GetSloughSides(), true);
            TS_ASSERT_DELTA(p_cell_killer->GetSloughHeight(), 10.0, 1e-9);
            TS_ASSERT_DELTA(p_cell_killer->GetSloughWidth(), 5.0, 1e-9);

            delete p_cell_killer;
        }
    }

    void TestRadialSloughingCellKillerMethods()
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
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            centre += mesh.GetNode(i)->rGetLocation();
        }
        centre = centre/mesh.GetNumNodes();

        // Choose radius of cell killer
        double radius = 0.4;

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Create cell killer and kill cells
        RadialSloughingCellKiller radial_cell_killer(&cell_population, centre, radius);
        radial_cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double r = norm_2(cell_population.GetLocationOfCellCentre(*cell_iter) - centre);

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
        cell_population.RemoveDeadCells();

        // Check that we are correctly left with cells inside the circle of death
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double r = norm_2(cell_population.GetLocationOfCellCentre(*cell_iter) - centre);
            TS_ASSERT_LESS_THAN_EQUALS(r, radius);
        }
    }

    void TestArchivingOfRadialSloughingCellKiller()
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

    void TestCellKillersOutputParameters()
   {
       std::string output_directory = "TestSloughingCellKillersOutputParameters";
       OutputFileHandler output_file_handler(output_directory, false);

       // Test with SloughingCellKiller
       SloughingCellKiller<2> sloughing_cell_killer(NULL, 20.0, true, 10.0);
       TS_ASSERT_EQUALS(sloughing_cell_killer.GetIdentifier(), "SloughingCellKiller-2");

       out_stream sloughing_cell_killer_parameter_file = output_file_handler.OpenOutputFile("sloughing_results.parameters");
       sloughing_cell_killer.OutputCellKillerParameters(sloughing_cell_killer_parameter_file);
       sloughing_cell_killer_parameter_file->close();

       std::string sloughing_cell_killer_results_dir = output_file_handler.GetOutputDirectoryFullPath();
       FileComparison( sloughing_cell_killer_results_dir + "sloughing_results.parameters", "crypt/test/data/TestSloughingCellKillers/sloughing_results.parameters").CompareFiles();

       // Test with RadialSloughingCellKiller
       c_vector<double,2> centre(2);
       centre[0] = 0.1;
       centre[1] = 0.2;
       double radius = 0.4;
       RadialSloughingCellKiller radial_cell_killer(NULL, centre, radius);
       TS_ASSERT_EQUALS(radial_cell_killer.GetIdentifier(), "RadialSloughingCellKiller");

       out_stream radial_cell_killer_parameter_file = output_file_handler.OpenOutputFile("radial_results.parameters");
       radial_cell_killer.OutputCellKillerParameters(radial_cell_killer_parameter_file);
       radial_cell_killer_parameter_file->close();

       std::string radial_cell_killer_results_dir = output_file_handler.GetOutputDirectoryFullPath();
       FileComparison( radial_cell_killer_results_dir + "radial_results.parameters", "crypt/test/data/TestSloughingCellKillers/radial_results.parameters").CompareFiles();
   }
};

#endif /*TESTSLOUGHINGCELLKILLER_HPP_*/
