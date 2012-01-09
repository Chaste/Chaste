/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef TESTVERTEXBASEDCELLPOPULATIONFORCRYPT_HPP_
#define TESTVERTEXBASEDCELLPOPULATIONFORCRYPT_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "WntCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"

class TestVertexBasedCellPopulationForCrypt : public AbstractCellBasedTestSuite
{
public:

    /**
     * Test that post-#878, WntConcentration copes with a VertexBasedCellPopulation.
     * \todo When vertex-based cell population code is added to cell_based folder, move this
     *       test to TestWntConcentration.hpp
     */
    void TestWntConcentrationWithVertexBasedCellPopulation() throw(Exception)
    {
        // Make some nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 2.0, -1.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 1.0));
        nodes.push_back(new Node<2>(2, true, -2.0, 1.0));
        nodes.push_back(new Node<2>(3, true, -2.0, -1.0));
        nodes.push_back(new Node<2>(4, true, 0.0, 2.0));

        // Make a rectangular element out of nodes 0,1,2,3
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);

        // Make a triangular element out of nodes 1,4,2
        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[1]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[2]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_2));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<vertex_mesh.GetNumElements(); i++)
        {
            WntCellCycleModel* p_cell_cycle_model = new WntCellCycleModel;
            p_cell_cycle_model->SetDimension(2);
            p_cell_cycle_model->SetCellProliferativeType(DIFFERENTIATED);

            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            double birth_time = 0.0 - i;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set the top of this cell_population, for the purposes of computing the WntConcentration
        double crypt_length = 4.0;

        // Set up an instance of the WntConcentration singleton object
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);

        // Check that the singleton can be set up
        p_wnt->SetType(LINEAR);
        p_wnt->SetCellPopulation(cell_population);
        p_wnt->SetCryptLength(crypt_length);

        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);

        // Check that the singleton can be destroyed then recreated
        WntConcentration<2>::Destroy();
        WntConcentration<2>::Instance()->SetType(NONE);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        TS_ASSERT_EQUALS(WntConcentration<2>::Instance()->IsWntSetUp(), false); // not fully set up now it is a NONE type

        WntConcentration<2>::Destroy();
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true); // set up again

        double wnt_at_cell0 = p_wnt->GetWntLevel(cell_population.GetCellUsingLocationIndex(0));
        double wnt_at_cell1 = p_wnt->GetWntLevel(cell_population.GetCellUsingLocationIndex(1));

        // We have set the top of the cell population to be 4, so the WntConcentration should decrease linearly
        // up the cell_population, from one at height 0 to zero at height 4.

        // Cell 0 has centre of mass (0,0)
        TS_ASSERT_DELTA(wnt_at_cell0, 1.0, 1e-4);

        // Cell 1 has centre of mass (0, 4/3)
        TS_ASSERT_DELTA(wnt_at_cell1, 2.0/3.0, 1e-4);
    }
};

#endif /*TESTVERTEXBASEDCELLPOPULATIONFORCRYPT_HPP_*/
