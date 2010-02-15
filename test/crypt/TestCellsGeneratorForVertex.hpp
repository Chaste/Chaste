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
#ifndef TESTCELLSGENERATORFORVERTEX_HPP_
#define TESTCELLSGENERATORFORVERTEX_HPP_

#include <cxxtest/TestSuite.h>

#include "FixedDurationGenerationBasedCellCycleModelCellsGeneratorForVertex.hpp"
#include "SimpleWntCellCycleModelCellsGeneratorForVertex.hpp"
#include "StochasticDurationGenerationBasedCellCycleModelCellsGeneratorForVertex.hpp"
//#include "StochasticWntCellCycleModelCellsGenerator.hpp"
//#include "TysonNovakCellCycleModelCellsGenerator.hpp"
//#include "WntCellCycleModelCellsGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "TrianglesMeshReader.hpp"

#include "AbstractCellBasedTestSuite.hpp"

/**
 * This class contains tests for methods on classes
 * inheriting from AbstractCellsGenerator.
 */
class TestCellsGeneratorForVertex : public AbstractCellBasedTestSuite
{
public:

//    void TestFixedDurationGenerationBasedCellCycleModelCellsGeneratorGenerateBasic() throw(Exception)
//    {
//        // Create mesh
//        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
//        TetrahedralMesh<2,2> mesh;
//        mesh.ConstructFromMeshReader(mesh_reader);
//
//        // Create mesh
//        unsigned crypt_width = 2;
//        unsigned crypt_height = 2;
//        HoneycombVertexMeshGenerator mesh_generator(crypt_width, crypt_height, true);
//        Cylindrical2dVertexMesh* p_mesh = mesh_generator.GetCylindricalMesh();
//
//        // Create cells
//        std::vector<TissueCell> cells;
//        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
//        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
//
//        // Test that cells were generated correctly
//        TS_ASSERT_EQUALS(cells.size(), mesh.GetNumNodes());
//
//        for (unsigned i=0; i<cells.size(); i++)
//        {
//            TS_ASSERT_DELTA(cells[i].GetBirthTime(), -(double)(i), 1e-9);
//        }
//    }
    
    void TestFixedAndStochasticDurationGenerationBasedCellCycleModelCellsGeneratorForVertexGenerateForVertexCrypt() throw(Exception)
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        HoneycombVertexMeshGenerator mesh_generator(crypt_width, crypt_height, true);
        Cylindrical2dVertexMesh* p_mesh = mesh_generator.GetCylindricalMesh();
        
        
        double y0 = 1.0;
        double y1 = 2.0;
        double y2 = 3.0;
        double y3 = 4.0;

		// Create cells
        std::vector<TissueCell> fixed_cells, stochastic_cells;
        FixedDurationGenerationBasedCellCycleModelCellsGeneratorForVertex<2> fixed_cells_generator;
        fixed_cells_generator.GenerateForVertexCrypt(fixed_cells, *p_mesh, std::vector<unsigned>(), true, y0, y1, y2,y3 );

		StochasticDurationGenerationBasedCellCycleModelCellsGeneratorForVertex<2> stochastic_cells_generator;
        stochastic_cells_generator.GenerateForVertexCrypt(stochastic_cells, *p_mesh, std::vector<unsigned>(), true, y0, y1, y2,y3 );


        TS_ASSERT_EQUALS(fixed_cells.size(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(stochastic_cells.size(), p_mesh->GetNumElements());

        // Test that cells were generated correctly
        for (unsigned i=0; i<fixed_cells.size(); i++)
        {
            double height = p_mesh->GetCentroidOfElement(i)[1];
            unsigned fixed_generation = static_cast<FixedDurationGenerationBasedCellCycleModel*>(fixed_cells[i].GetCellCycleModel())->GetGeneration();
            unsigned stochastic_generation = static_cast<StochasticDurationGenerationBasedCellCycleModel*>(stochastic_cells[i].GetCellCycleModel())->GetGeneration();

            if (height <= y0)
            {
                TS_ASSERT_EQUALS(fixed_generation, 0u);
                TS_ASSERT_EQUALS(stochastic_generation, 0u);
            }
            else if (height < y1)
            {
                TS_ASSERT_EQUALS(fixed_generation, 1u);
                TS_ASSERT_EQUALS(stochastic_generation, 1u);
            }
            else if (height < y2)
            {
                TS_ASSERT_EQUALS(fixed_generation, 2u);
                TS_ASSERT_EQUALS(stochastic_generation, 2u);
            }
            else if (height < y3)
            {
                TS_ASSERT_EQUALS(fixed_generation, 3u);
                TS_ASSERT_EQUALS(stochastic_generation, 3u);
            }
            else
            {
                TS_ASSERT_EQUALS(fixed_generation, 4u);
                TS_ASSERT_EQUALS(stochastic_generation, 4u);
            }
        }
    }
    

    void TestSimpleWntCellCycleModelCellsGeneratorForVertexGenerateForVertexCrypt() throw(Exception)
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        HoneycombVertexMeshGenerator mesh_generator(crypt_width, crypt_height, true);
        Cylindrical2dVertexMesh* p_mesh = mesh_generator.GetCylindricalMesh();
        
        // Create cells
        std::vector<TissueCell> cells;
        SimpleWntCellCycleModelCellsGeneratorForVertex<2> cells_generator;
        cells_generator.GenerateForVertexCrypt(cells, *p_mesh, std::vector<unsigned>(), false);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumElements());
        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_DELTA(cells[i].GetBirthTime(), 0.0, 1e-9);
        }
    }
};

#endif /*TESTCELLSGENERATORFORVERTEX_HPP_*/
