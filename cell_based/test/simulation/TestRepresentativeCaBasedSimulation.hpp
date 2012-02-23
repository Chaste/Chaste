/*

Copyright (c) 2005-2012, University of Oxford.
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
#ifndef TESTREPRESENTATIVECABASEDONLATTICESIMULATION_HPP_
#define TESTREPRESENTATIVECABASEDONLATTICESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "OnLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "CaBasedCellPopulation.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "RandomCellKiller.hpp"
#include "NumericFileComparison.hpp"
#include "CellBasedEventHandler.hpp"

/**
 * This class consists of a single test - a 2D Cellular Automaton-based cell population
 * simulation. It has 9 cells that will divide in a square in the middle,
 * and the rest are empty sites.
 *
 * This test is used for profiling, to establish the run time
 * variation as the code is developed.
 */
class TestRepresentativeLatticeBasedCellBasedSimulation : public AbstractCellBasedTestSuite
{
public:

    void TestRepresentativeCaBasedSimulationForProfiling() throw (Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(20, 20, true); // 21*21 nodes

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 9);

        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->GetCellCycleModel()->SetCellProliferativeType(TRANSIT);
            dynamic_cast<FixedDurationGenerationBasedCellCycleModel*>(cells[i]->GetCellCycleModel())->SetMaxTransitGenerations(UINT_MAX);
            real_node_indices.push_back(21 * (9 + i/3) + 9 + i%3);
         }

        // Create a cell population
        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestRepresentativeCaBasedSimulationForProfiling");
        simulator.SetEndTime(50);

        // Run simulation
        simulator.Solve();

        RandomNumberGenerator::Destroy();
    }
};

#endif /* TESTREPRESENTATIVECABASEDONLATTICESIMULATION_HPP_ */
