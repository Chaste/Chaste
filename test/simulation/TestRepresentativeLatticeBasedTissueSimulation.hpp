/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef TESTREPRESENTATIVELATTICEBASEDTISSUESIMULATION_HPP_
#define TESTREPRESENTATIVELATTICEBASEDTISSUESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "TissueSimulationArchiver.hpp"

#include "LatticeBasedTissueSimulation.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "LatticeBasedTissue.hpp"
#include "DiffusionUpdateRule.hpp"
#include "AdvectionUpdateRule.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "RandomCellKiller.hpp"
#include "NumericFileComparison.hpp"
#include "CellBasedEventHandler.hpp"

/**
 * This class consists of a single test - a 2D lattice-based tissue
 * simulation.  It has 9 cells that will divide in a square in the middle,
 * and the rest are empty sites.
 *
 * This test is used for profiling, to establish the run time
 * variation as the code is developed.
 */
class TestRepresentativeLatticeBasedTissueSimulation : public AbstractCellBasedTestSuite
{
public:

    void TestRepresentativeLatticeBasedTissueSimulationForProfiling() throw (Exception)
    {
        // Create mesh
         TetrahedralMesh<2,2> mesh;
         mesh.ConstructRectangularMesh(20, 20, true); // 21*21 nodes

         // Create cell which keeps dividing
         TissueConfig::Instance()->SetMaxTransitGenerations(UINT_MAX);
         boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

         std::vector<TissueCell> cells;
         std::vector<unsigned> real_node_indices;
          for(unsigned i = 0; i < 9; i++)
         {
             FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
             p_model->SetCellProliferativeType(TRANSIT);
             double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                  (TissueConfig::Instance()->GetStemCellG1Duration()
                                      + TissueConfig::Instance()->GetSG2MDuration() );

             TissueCell cell(p_state, p_model);
        	 cell.SetBirthTime(birth_time);

        	 cells.push_back(cell);
        	 real_node_indices.push_back(21 * (9 + i/3) + 9 + i%3);
         }

         // Create a tissue
         LatticeBasedTissue<2> tissue(mesh, cells, real_node_indices);

         // Create an empty UpdateRule system so only movement from cell birth
         std::vector<AbstractUpdateRule<2>* > update_rule_collection;

         // Set up tissue simulation
         LatticeBasedTissueSimulation<2> simulator(tissue, update_rule_collection);
         simulator.SetOutputDirectory("TestRepresentativeLatticeBasedSimulationForProfiling");
         simulator.SetEndTime(50);

         // Run simulation
         simulator.Solve();

         RandomNumberGenerator::Destroy();
    }
};

#endif /* TESTREPRESENTATIVELATTICEBASEDTISSUESIMULATION_HPP_ */
