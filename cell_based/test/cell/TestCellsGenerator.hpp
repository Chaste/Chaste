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

#ifndef TESTCELLSGENERATOR_HPP_
#define TESTCELLSGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "TrianglesMeshReader.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"

class TestCellsGenerator : public AbstractCellBasedTestSuite
{
public:

    void TestGenerateBasicWithFixedDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_DELTA(cells[i]->GetBirthTime(), -(double)(i), 1e-9);
            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetDimension(), 2u);
        }

        // Test with extra input argument
        std::vector<unsigned> location_indices;
        location_indices.push_back(2);
        location_indices.push_back(7);
        location_indices.push_back(9);

        std::vector<CellPtr> cells2;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator2;
        cells_generator2.GenerateBasic(cells2, 3, location_indices);

        TS_ASSERT_EQUALS(cells2.size(), 3u);
        TS_ASSERT_DELTA(cells2[0]->GetBirthTime(), -2.0, 1e-4);
        TS_ASSERT_DELTA(cells2[1]->GetBirthTime(), -7.0, 1e-4);
        TS_ASSERT_DELTA(cells2[2]->GetBirthTime(), -9.0, 1e-4);
    }

    void TestGenerateGivenLocationIndicesWithFixedDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        // Use a mesh generator to generate some location indices corresponding to real cells
        HoneycombMeshGenerator mesh_generator(6, 7, 2);
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells again with basic
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        TS_ASSERT_THROWS_THIS(cells_generator.GenerateBasic(cells, 83511u, location_indices),
                              "The size of the locationIndices vector must match the required number of output cells");
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        // Test that cells were generated correctly
        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_DELTA(cells[i]->GetBirthTime(), -(double)(location_indices[i]), 1e-9);
            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetDimension(), 2u);
        }
    }

    void TestGenerateBasicRandomWithFixedDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), TRANSIT);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            // Should lie between -24 and 0
            TS_ASSERT_LESS_THAN_EQUALS(cells[i]->GetBirthTime(), 0.0);
            TS_ASSERT_LESS_THAN_EQUALS(-24.0, cells[i]->GetBirthTime());
            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetDimension(), 2u);
            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetCellProliferativeType(),TRANSIT);
        }

        // Test exact random numbers as test re-seeds random number generator.
        TS_ASSERT_DELTA(cells[0]->GetBirthTime(), -4.7325, 1e-4);
        TS_ASSERT_DELTA(cells[1]->GetBirthTime(), -9.5812, 1e-4);
        TS_ASSERT_DELTA(cells[2]->GetBirthTime(), -2.3706, 1e-4);
    }

    void TestGenerateBasicRandomWithFixedDurationGenerationBasedCellCycleModelandVertexCells() throw(Exception)
    {
        // Create mesh
        HoneycombVertexMeshGenerator mesh_generator(2, 2);
        VertexMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), TRANSIT);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumElements());

        for (unsigned i=0; i<cells.size(); i++)
        {
            // Shold lie between -24 and 0
            TS_ASSERT_LESS_THAN_EQUALS(cells[i]->GetBirthTime(), 0.0);
            TS_ASSERT_LESS_THAN_EQUALS(-24.0, cells[i]->GetBirthTime());
            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetDimension(), 2u);
            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetCellProliferativeType(),TRANSIT);
        }

        // Test exact random numbers as test re-seeds random number generator.
        TS_ASSERT_DELTA(cells[0]->GetBirthTime(), -4.7325, 1e-4);
        TS_ASSERT_DELTA(cells[1]->GetBirthTime(), -9.5812, 1e-4);
        TS_ASSERT_DELTA(cells[2]->GetBirthTime(), -2.3706, 1e-4);
    }
};

#endif /*TESTCELLSGENERATOR_HPP_*/
