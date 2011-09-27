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
#ifndef TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_
#define TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "CryptCellsGenerator.hpp"
#include "OnLatticeSimulation.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "SloughingCellKiller.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PottsMeshGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "Warnings.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"

/**
 * This class consists of a single test - a 2D Potts-based cell population
 * simulation of 100 cells with differential adhesion and no birth or death.
 *
 * This test is used for profiling, to establish the run time
 * variation as the code is developed.
 */
class TestRepresentativePottsBasedOnLatticeSimulationForProfiling : public AbstractCellBasedTestSuite
{
public:

    void TestPottsMonolayerCellSorting() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(60, 10, 4, 60, 10, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);

        // Make this pointer first as if we move it after creating the cell population the label numbers aren't tracked
        MAKE_PTR(CellLabel, p_label);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.SetOutputCellMutationStates(true); // So outputs the labeled cells

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (RandomNumberGenerator::Instance()->ranf() < 0.5)
            {
                (*cell_iter)->AddCellProperty(p_label);
            }
        }

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestRepresentativePottsBasedSimulationForProfiling");
        simulator.SetDt(0.1);
        simulator.SetEndTime(10);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(DifferentialAdhesionPottsUpdateRule<2>, p_differential_adhesion_update_rule);
        p_differential_adhesion_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.11);
        p_differential_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02);
        p_differential_adhesion_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16);
        simulator.AddPottsUpdateRule(p_differential_adhesion_update_rule);

        // Run simulation
        simulator.Solve();
    }
};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/
