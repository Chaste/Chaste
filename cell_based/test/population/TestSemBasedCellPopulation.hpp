/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef TESTSEMBASEDCELLPOPULATION_HPP_
#define TESTSEMBASEDCELLPOPULATION_HPP_


#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "SemBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "ApcOneHitCellMutationState.hpp"

#include "AbstractCellBasedTestSuite.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSemBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    static std::vector<Node<2>*> CreateNodes()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1u, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(2u, false, 0.05, 0.0));
        return nodes;
    }

    static std::vector<SemElement<2>*> CreateElements(std::vector<Node<2>*>& rNodes)
    {
        std::vector<SemElement<2>*> elements;

        std::vector<Node<2>*> element_0_nodes;
        element_0_nodes.push_back(rNodes[0]);
        element_0_nodes.push_back(rNodes[1]);
        elements.push_back(new SemElement<2>(0u, element_0_nodes));

        std::vector<Node<2>*> element_1_nodes;
        element_1_nodes.push_back(rNodes[1]);
        element_1_nodes.push_back(rNodes[2]);
        elements.push_back(new SemElement<2>(1u, element_1_nodes));

        return elements;
    }

    struct TwoElementSemMesh
    {
        std::vector<Node<2>*> nodes;
        std::vector<SemElement<2>*> elements;
        SemMesh<2> mesh;

        TwoElementSemMesh()
            : nodes(CreateNodes()),
              elements(CreateElements(nodes)),
              mesh(nodes, elements)
        {
        }
    };

    static std::vector<CellPtr> CreateCells(unsigned numCells)
    {
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, numCells);
        return cells;
    }

public:

    void TestConstructorRespectsLocationIndices()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        CellPtr p_cell_0 = cells[0];
        CellPtr p_cell_1 = cells[1];
        std::vector<unsigned> location_indices;
        location_indices.push_back(1u);
        location_indices.push_back(0u);

        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells, false, true, location_indices);

        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_cell_0), 1u);
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_cell_1), 0u);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(1u), p_cell_0);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(0u), p_cell_1);
    }

    void TestConstructorRejectsInvalidLocationIndices()
    {
        {
            TwoElementSemMesh fixture;
            std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
            std::vector<unsigned> location_indices;
            location_indices.push_back(0u);
            location_indices.push_back(0u);

            TS_ASSERT_THROWS_THIS(SemBasedCellPopulation<2> cell_population(fixture.mesh, cells, false, true, location_indices),
                                  "A SemElement location index is assigned to more than one cell");
        }

        {
            TwoElementSemMesh fixture;
            std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
            std::vector<unsigned> location_indices;
            location_indices.push_back(0u);
            location_indices.push_back(2u);

            TS_ASSERT_THROWS_THIS(SemBasedCellPopulation<2> cell_population(fixture.mesh, cells, false, true, location_indices),
                                  "A supplied location index does not correspond to a SemElement");
        }
    }

    void TestConstructorRejectsWrongNumberOfCells()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(1u);

        TS_ASSERT_THROWS_THIS(SemBasedCellPopulation<2> cell_population(fixture.mesh, cells),
                              "There must be precisely one CellPtr for each SemElement");
    }

    void TestBasicAccessorsWidthAndDefaultTimeStep()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        TS_ASSERT_EQUALS(&(cell_population.rGetMesh()), &(fixture.mesh));
        TS_ASSERT_EQUALS(cell_population.GetElement(0u), fixture.mesh.GetElement(0u));
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(cell_population.GetNode(2u), fixture.mesh.GetNode(2u));
        TS_ASSERT_DELTA(cell_population.GetWidth(0u), 0.05, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetWidth(1u), 1.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetDefaultTimeStep(), 0.002, 1e-12);

        ChastePoint<2> new_location(0.25, 0.25);
        cell_population.SetNode(2u, new_location);
        TS_ASSERT_DELTA(cell_population.GetNode(2u)->rGetLocation()[0], 0.25, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(2u)->rGetLocation()[1], 0.25, 1e-12);
    }

    void TestUpdatePopulatesNeighbouringNodeAndLocationIndices()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        c_vector<double, 4> domain;
        domain[0] = -0.1;
        domain[1] =  0.2;
        domain[2] = -0.1;
        domain[3] =  1.1;
        fixture.mesh.SetUpBoxCollection(0.2, domain);

        cell_population.Update(false);

        std::set<unsigned> neighbouring_node_indices = cell_population.GetNeighbouringNodeIndices(0u);
        TS_ASSERT_EQUALS(neighbouring_node_indices.size(), 1u);
        TS_ASSERT_EQUALS(neighbouring_node_indices.count(2u), 1u);

        std::set<unsigned> neighbouring_location_indices = cell_population.GetNeighbouringLocationIndices(cell_population.GetCellUsingLocationIndex(0u));
        TS_ASSERT_EQUALS(neighbouring_location_indices.size(), 1u);
        TS_ASSERT_EQUALS(neighbouring_location_indices.count(1u), 1u);
    }

    void TestDampingConstantUsesContainingElementCells()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);
        cell_population.SetDampingConstantNormal(2.0);
        cell_population.SetDampingConstantMutant(8.0);

        boost::shared_ptr<AbstractCellProperty> p_apc1(cell_population.GetCellPropertyRegistry()->Get<ApcOneHitCellMutationState>());
        cell_population.GetCellUsingLocationIndex(1u)->SetMutationState(p_apc1);

        TS_ASSERT_DELTA(cell_population.GetDampingConstant(0u), 2.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetDampingConstant(2u), 8.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetDampingConstant(1u), 5.0, 1e-12);
    }

    void TestRemoveDeadCellsDeletesElementsAndMappings()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        CellPtr p_dead_cell = cell_population.GetCellUsingLocationIndex(1u);
        p_dead_cell->Kill();

        TS_ASSERT(!cell_population.GetElement(1u)->IsDeleted());
        TS_ASSERT(!cell_population.IsCellAssociatedWithADeletedLocation(p_dead_cell));

        TS_ASSERT_EQUALS(cell_population.RemoveDeadCells(), 1u);
        TS_ASSERT(cell_population.GetElement(1u)->IsDeleted());
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1u);
        TS_ASSERT_EQUALS(fixture.mesh.GetNode(2u)->rGetContainingElementIndices().count(1u), 0u);
    }

    void TestAddCellIsExplicitlyUnsupported()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);
        std::vector<CellPtr> new_cells = CreateCells(1u);

        TS_ASSERT_THROWS_THIS(cell_population.AddCell(new_cells[0], cell_population.GetCellUsingLocationIndex(0u)),
                              "SemBasedCellPopulation does not support AddCell() because SEM element division is not implemented");
    }
};

#endif /*TESTSEMBASEDCELLPOPULATION_HPP_*/
