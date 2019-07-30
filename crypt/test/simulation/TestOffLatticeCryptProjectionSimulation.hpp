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

#ifndef TESTOFFLATTICECRYPTPROJECTIONSIMULATION_HPP_
#define TESTOFFLATTICECRYPTPROJECTIONSIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>

#include "CheckpointArchiveTypes.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "CryptCellsGenerator.hpp"
#include "CryptProjectionForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NumericFileComparison.hpp"
#include "OffLatticeSimulation.hpp"
#include "RadialSloughingCellKiller.hpp"
#include "RandomCellKiller.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellCycleModelProteinConcentrationsWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"

// Cell population writers
#include "CellPopulationAreaWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestOffLatticeCryptProjectionSimulation : public AbstractCellBasedWithTimingsTestSuite
{
public:
    void TestOutputStatistics()
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2, 2>* p_mesh = generator.GetMesh();

        double crypt_length = (double)num_cells_depth * sqrt(3.0) / 2.0;

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticWntCellCycleModel> cell_generator;
        cell_generator.Generate(cells, p_mesh, location_indices, true);

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddPopulationWriter<CellPopulationAreaWriter>();
        cell_population.AddCellWriter<CellCycleModelProteinConcentrationsWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up Wnt Gradient
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        TS_ASSERT_EQUALS(simulator.GetIdentifier(), "OffLatticeSimulation-2-2");
        simulator.SetOutputDirectory("OffLatticeSimulationWritingProteins");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        TS_ASSERT_DELTA(simulator.GetDt(), 1.0 / 120.0, 1e-12);

        // Run cell-based simulation
        TS_ASSERT_EQUALS(simulator.GetOutputDirectory(), "OffLatticeSimulationWritingProteins");
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        OutputFileHandler handler("OffLatticeSimulationWritingProteins", false);

        std::string protein_conc_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/proteinconcentrations.dat";
        NumericFileComparison comp_protein_conc(protein_conc_file, "crypt/test/data/OffLatticeSimulationWritingProteins/proteinconcentrations.dat");
        TS_ASSERT(comp_protein_conc.CompareFiles(1e-2));

        std::string cell_cycle_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellcyclephases.dat";
        NumericFileComparison comp_cell_cycle(cell_cycle_file, "crypt/test/data/OffLatticeSimulationWritingProteins/cellcyclephases.dat");
        TS_ASSERT(comp_cell_cycle.CompareFiles(1e-2));

        std::string cell_ages_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellages.dat";
        NumericFileComparison comp_cell_ages(cell_ages_file, "crypt/test/data/OffLatticeSimulationWritingProteins/cellages.dat");
        TS_ASSERT(comp_cell_ages.CompareFiles(1e-2));

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    /**
     * Test a cell-based simulation with a non-Meineke spring system.
     *
     * This test consists of a standard crypt projection model simulation with a
     * radial sloughing cell killer, a crypt projection cell-cycle model that
     * depends on a radial Wnt gradient, and the crypt projection model spring
     * system, and store the results for use in later archiving tests.
     */
    void TestOffLatticeSimulationWithCryptProjectionSpringSystem()
    {
        double a = 0.2;
        double b = 2.0;
        WntConcentration<2>::Instance()->SetCryptProjectionParameterA(a);
        WntConcentration<2>::Instance()->SetCryptProjectionParameterB(b);

        // Set up mesh
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

        // To start off with, set up all cells to have TransitCellProliferativeType
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        for (unsigned i = 0; i < location_indices.size(); i++)
        {
            SimpleWntCellCycleModel* p_model = new SimpleWntCellCycleModel();
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

        // Set up the Wnt gradient
        WntConcentration<2>::Instance()->SetType(RADIAL);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Make a cell-based simulation
        OffLatticeSimulation<2> crypt_projection_simulator(crypt, false, false);

        // Create a force law and pass it to the simulation
        MAKE_PTR(CryptProjectionForce, p_force);
        crypt_projection_simulator.AddForce(p_force);

        // Create a radial cell killer and pass it in to the cell-based simulation
        c_vector<double, 2> centre = zero_vector<double>(2);
        double crypt_radius = pow(crypt_length / a, 1.0 / b);

        MAKE_PTR_ARGS(RadialSloughingCellKiller, p_killer, (&crypt, centre, crypt_radius));
        crypt_projection_simulator.AddCellKiller(p_killer);

        // Set up the simulation
        crypt_projection_simulator.SetOutputDirectory("CryptProjectionSimulation");
        crypt_projection_simulator.SetEndTime(0.25);

        // Run the simulation
        TS_ASSERT_THROWS_NOTHING(crypt_projection_simulator.Solve());

        // These cells just divided and have been gradually moving apart.
        // This is happening around (4, -1).  The exact spring length is slightly
        // compiler/architecture dependent.
        // These results are from time 0.25.
        std::vector<double> node_a_location = crypt_projection_simulator.GetNodeLocation(301);
        std::vector<double> node_b_location = crypt_projection_simulator.GetNodeLocation(504);

        c_vector<double, 2> distance_between;
        distance_between(0) = node_a_location[0] - node_b_location[0];
        distance_between(1) = node_a_location[1] - node_b_location[1];
        // Note that this distance varies based on the quality of the original honeycomb mesh,
        // the precision of the machine and the optimisation level
        TS_ASSERT_DELTA(norm_2(distance_between), 6.1701, 4.0e-3);

        // Test the Wnt concentration result
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(crypt.GetCellUsingLocationIndex(257)), 0.7777, 2e-3);
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(crypt.GetCellUsingLocationIndex(503)), 0.8598, 1e-3);
        // Tidy up
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTOFFLATTICECRYPTPROJECTIONSIMULATION_HPP_*/
