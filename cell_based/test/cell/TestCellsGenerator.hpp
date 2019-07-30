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

#ifndef TESTCELLSGENERATOR_HPP_
#define TESTCELLSGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "TrianglesMeshReader.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestCellsGenerator : public AbstractCellBasedTestSuite
{
public:

    void TestGenerateBasicWithFixedG1GenerationalCellCycleModel()
    {
        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
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
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator2;
        cells_generator2.GenerateBasic(cells2, 3, location_indices);

        TS_ASSERT_EQUALS(cells2.size(), 3u);
        TS_ASSERT_DELTA(cells2[0]->GetBirthTime(), -2.0, 1e-4);
        TS_ASSERT_DELTA(cells2[1]->GetBirthTime(), -7.0, 1e-4);
        TS_ASSERT_DELTA(cells2[2]->GetBirthTime(), -9.0, 1e-4);
    }

    void TestGenerateGivenLocationIndicesWithFixedG1GenerationalCellCycleModel()
    {
        EXIT_IF_PARALLEL;
        // Use a mesh generator to generate some location indices corresponding to real cells
        HoneycombMeshGenerator mesh_generator(6, 7, 2);
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
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

    void TestGenerateGivenLocationIndicesWithSpecifiedCellProliferativeType()
    {
        EXIT_IF_PARALLEL;
        // Use a mesh generator to generate some location indices corresponding to real cells
        HoneycombMeshGenerator mesh_generator(6, 7, 2);
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;

        cells_generator.GenerateGivenLocationIndices(cells,
                                                     location_indices,
                                                     CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

        // Test that cells were generated correctly
        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_DELTA(cells[i]->GetBirthTime(), -(double)(location_indices[i]), 1e-9);
            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetDimension(), 2u);
            TS_ASSERT_EQUALS(cells[i]->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>(), true);
        }
    }

    void TestGenerateBasicRandomWithNoSpecifiedProliferativeCellType()
    {
        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetDimension(), 2u);
            TS_ASSERT_EQUALS(cells[i]->GetCellProliferativeType()->IsType<StemCellProliferativeType>(), true);
            // Should lie between -24 and 0
            double birth_time=cells[i]->GetBirthTime();
            ///\todo Breaks Intel 10? TS_ASSERT_LESS_THAN_EQUALS(birth_time, 0.0);
            TS_ASSERT_LESS_THAN_EQUALS(-24.0, birth_time);
        }
    }

    void TestGenerateBasicRandomWithFixedG1GenerationalCellCycleModel()
    {
        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            // Should lie between -24 and 0
            TS_ASSERT_LESS_THAN_EQUALS(cells[i]->GetBirthTime(), 0.0);
            TS_ASSERT_LESS_THAN_EQUALS(-24.0, cells[i]->GetBirthTime());
            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetDimension(), 2u);
            TS_ASSERT_EQUALS(cells[i]->GetCellProliferativeType(), p_transit_type);
        }

        // Test exact random numbers as test re-seeds random number generator.
        TS_ASSERT_DELTA(cells[0]->GetBirthTime(), -7.1141, 1e-4);
        TS_ASSERT_DELTA(cells[1]->GetBirthTime(), -10.1311, 1e-4);
        TS_ASSERT_DELTA(cells[2]->GetBirthTime(), -10.2953, 1e-4);
    }

    void TestGenerateBasicRandomWithFixedG1GenerationalCellCycleModelandVertexCells()
    {
        EXIT_IF_PARALLEL;
        // Create mesh
        HoneycombVertexMeshGenerator mesh_generator(2, 2);
        VertexMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumElements());

        for (unsigned i=0; i<cells.size(); i++)
        {
            // Shold lie between -24 and 0
            TS_ASSERT_LESS_THAN_EQUALS(cells[i]->GetBirthTime(), 0.0);
            TS_ASSERT_LESS_THAN_EQUALS(-24.0, cells[i]->GetBirthTime());
            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetDimension(), 2u);
            TS_ASSERT_EQUALS(cells[i]->GetCellProliferativeType(), p_transit_type);
        }

        // Test exact random numbers as test re-seeds random number generator.
        TS_ASSERT_DELTA(cells[0]->GetBirthTime(), -7.1141, 1e-4);
        TS_ASSERT_DELTA(cells[1]->GetBirthTime(), -10.1311, 1e-4);
        TS_ASSERT_DELTA(cells[2]->GetBirthTime(), -10.2953, 1e-4);
    }
};

#endif /*TESTCELLSGENERATOR_HPP_*/
