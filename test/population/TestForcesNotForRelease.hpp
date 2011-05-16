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
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"

class TestForcesNotForRelease : public AbstractCellBasedTestSuite
{
public:

    void TestBuskeInteractionForceMethods() throw (Exception)
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        // Create a simple mesh
		unsigned num_cells_depth = 5;
		unsigned num_cells_width = 5;
		HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
		TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

		// Create cells
		std::vector<CellPtr> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

		// Create a node-based cell population
		NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
		cell_population.SetMechanicsCutOffLength(1.5);
		cell_population.Update();


        // Create force
        BuskeInteractionForce<2> buske_force;

        // Test set/get method
        TS_ASSERT_EQUALS(buske_force.GetUseCutOffLength(), false);
        TS_ASSERT_DELTA(buske_force.GetCutOffLength(), DBL_MAX, 1e-6);

        buske_force.SetCutOffLength(1.5);

        TS_ASSERT_EQUALS(buske_force.GetUseCutOffLength(), true);
        TS_ASSERT_DELTA(buske_force.GetCutOffLength(), 1.5, 1e-6);

        // Reset cut off length
        buske_force.SetCutOffLength(DBL_MAX);
        /*
         ************************************************************************
         ************************************************************************
         *  Test node force calculation
         ************************************************************************
         ************************************************************************
         */

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }

        buske_force.AddForceContribution(node_forces, cell_population);

        // Test forces on non-ghost nodes
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);

            TS_ASSERT_DELTA(node_forces[node_index][0], 0.0, 1e-4);
            TS_ASSERT_DELTA(node_forces[node_index][1], 0.0, 1e-4);
        }

    }

    void TestForceOutputParameters()
    {
		std::string output_directory = "TestForcesOutputParameters";
		OutputFileHandler output_file_handler(output_directory, false);

		// Test with BuskeInteractionForce
		BuskeInteractionForce<2> buske_force;
		buske_force.SetCutOffLength(1.5);
		TS_ASSERT_EQUALS(buske_force.GetIdentifier(), "BuskeInteractionForce-2");

		out_stream buske_force_parameter_file = output_file_handler.OpenOutputFile("buske_results.parameters");
		buske_force.OutputForceParameters(buske_force_parameter_file);
		buske_force_parameter_file->close();

		std::string buske_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
		TS_ASSERT_EQUALS(system(("diff " + buske_force_results_dir + "buske_results.parameters notforrelease_cell_based/test/data/TestForcesNotForRelease/buske_results.parameters").c_str()), 0);
    }

    void TestBuskeInteractionForceArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "buske_force_system.arch";

        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");

            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

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
            BuskeInteractionForce<2> buske_force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            BuskeInteractionForce<2>* const p_buske_force = &buske_force;

            ///\todo set member variables (#1568)

            output_arch << p_buske_force;
        }

        {
            ArchiveLocationInfo::SetMeshPathname("mesh/test/data", "square_2_elements");

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            BuskeInteractionForce<2>* p_buske_force;

            // Restore from the archive
            input_arch >> p_buske_force;

            ///\todo Test the member data (#1568)

            // Tidy up
            delete p_buske_force;
        }
    }

};
#endif /*TESTFORCESNOTFORRELEASE_HPP_*/
