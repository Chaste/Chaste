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

#ifndef TESTCRYPTSTATISTICS_HPP_
#define TESTCRYPTSTATISTICS_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CryptProjectionForce.hpp"
#include "CryptProjectionStatistics.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "OffLatticeSimulation.hpp"
#include "RadialSloughingCellKiller.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "WntConcentration.hpp"

#include "FakePetscSetup.hpp"

class TestCryptProjectionStatistics : public AbstractCellBasedTestSuite
{

public:
    void TestGetSection()
    {
        double a = 0.2;
        double b = 2.0;
        WntConcentration<2>::Instance()->SetCryptProjectionParameterA(a);
        WntConcentration<2>::Instance()->SetCryptProjectionParameterB(b);

        int num_cells_depth = 20;
        int num_cells_width = 20;
        unsigned thickness_of_ghost_layer = 3;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer);
        MutableMesh<2, 2>* p_mesh = generator.GetMesh();

        double crypt_length = (double)num_cells_depth * sqrt(3.0) / 2.0;

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        ChasteCuboid<2> bounding_box = p_mesh->CalculateBoundingBox();
        double width_of_mesh = (num_cells_width / (num_cells_width + 2.0 * thickness_of_ghost_layer)) * (bounding_box.GetWidth(0));
        double height_of_mesh = (num_cells_depth / (num_cells_depth + 2.0 * thickness_of_ghost_layer)) * (bounding_box.GetWidth(1));

        p_mesh->Translate(-width_of_mesh / 2, -height_of_mesh / 2);

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        for (unsigned i = 0; i < location_indices.size(); i++)
        {
            SimpleWntCellCycleModel* p_model = new SimpleWntCellCycleModel;
            p_model->SetDimension(2);
            p_model->SetWntStemThreshold(0.95);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_transit_type);

            p_cell->InitialiseCellCycleModel();
            double birth_time = -RandomNumberGenerator::Instance()->ranf() * (p_model->GetTransitCellG1Duration()
                                                                              + p_model->GetSG2MDuration());
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Make a cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        double crypt_radius = pow(crypt_length / a, 1.0 / b);

        // Set up the Wnt gradient
        WntConcentration<2>::Instance()->SetType(RADIAL);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        CryptProjectionStatistics statistics(crypt);

        std::vector<CellPtr> test_section = statistics.GetCryptSection(M_PI / 2.0);

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section.size(), 10u);

        unsigned expected_indices[10] = { 350, 377, 402, 429, 454, 481, 506, 533, 558, 585 };

        for (unsigned i = 0; i < test_section.size(); i++)
        {
            TS_ASSERT_EQUALS(crypt.GetLocationIndexUsingCell(test_section[i]), expected_indices[i]);
        }

        // Make a cell-based simulation
        OffLatticeSimulation<2> crypt_projection_simulator(crypt, false, false);

        // Create a force law and pass it to the OffLatticeSimulation
        MAKE_PTR(CryptProjectionForce, p_crypt_projection_force);
        crypt_projection_simulator.AddForce(p_crypt_projection_force);

        // Create a radial cell killer and pass it in to the cell-based simulation
        c_vector<double, 2> centre = zero_vector<double>(2);

        MAKE_PTR_ARGS(RadialSloughingCellKiller, p_killer, (&crypt, centre, crypt_radius));
        crypt_projection_simulator.AddCellKiller(p_killer);

        // Set up the simulation
        crypt_projection_simulator.SetOutputDirectory("CryptProjectionStatistics");
        crypt_projection_simulator.SetEndTime(0.25);
        crypt_projection_simulator.Solve();

        statistics.LabelSPhaseCells();

        std::vector<CellPtr> test_section2 = statistics.GetCryptSection();
        std::vector<bool> labelled_cells = statistics.AreCryptSectionCellsLabelled(test_section2);

        TS_ASSERT_EQUALS(test_section2.size(), labelled_cells.size());
        TS_ASSERT_EQUALS(test_section2.size(), 13u);

        crypt_projection_simulator.SetEndTime(0.3);
        crypt_projection_simulator.Solve();

        // Tidy up
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTCRYPTSTATISTICS_HPP_*/
