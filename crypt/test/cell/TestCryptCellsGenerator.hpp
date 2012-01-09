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

#ifndef TESTCRYPTCELLSGENERATOR_HPP_
#define TESTCRYPTCELLSGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "CryptCellsGenerator.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TrianglesMeshReader.hpp"
#include "AbstractCellBasedTestSuite.hpp"

// For archiving?
#include "WildTypeCellMutationState.hpp"

class TestCryptCellsGenerator : public AbstractCellBasedTestSuite
{
public:

    void TestCryptCellsGeneratorWithFixedDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;

        double y0 = 0.2;
        double y1 = 1.0;
        double y2 = 2.0;
        double y3 = 3.0;

        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, true, y0, y1, y2, y3);

        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        // Test that cells were generated correctly
        for (unsigned i=0; i<cells.size(); i++)
        {
            double height = p_mesh->GetNode(i)->rGetLocation()[1];
            unsigned generation = static_cast<FixedDurationGenerationBasedCellCycleModel*>(cells[i]->GetCellCycleModel())->GetGeneration();

            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetDimension(), 2u);

            if (height <= y0)
            {
                TS_ASSERT_EQUALS(generation, 0u);
                TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetCellProliferativeType(), STEM);
            }
            else if (height < y1)
            {
                TS_ASSERT_EQUALS(generation, 1u);
                TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);
            }
            else if (height < y2)
            {
                TS_ASSERT_EQUALS(generation, 2u);
                TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);
            }
            else if (height < y3)
            {
                TS_ASSERT_EQUALS(generation, 3u);
                TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);
            }
            else
            {
                TS_ASSERT_EQUALS(generation, 4u);
                TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetCellProliferativeType(), DIFFERENTIATED);
            }
        }
    }

    void TestCryptCellsGeneratorWithStochasticDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, false);

        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        double y0 = 0.3;
        double y1 = 2.0;
        double y2 = 3.0;
        double y3 = 4.0;

        // Test that cells were generated correctly
        for (unsigned i=0; i<cells.size(); i++)
        {
            double height = p_mesh->GetNode(i)->rGetLocation()[1];
            unsigned generation = static_cast<StochasticDurationGenerationBasedCellCycleModel*>(cells[i]->GetCellCycleModel())->GetGeneration();

            if (height <= y0)
            {
                TS_ASSERT_EQUALS(generation, 0u);
            }
            else if (height < y1)
            {
                TS_ASSERT_EQUALS(generation, 1u);
            }
            else if (height < y2)
            {
                TS_ASSERT_EQUALS(generation, 2u);
            }
            else if (height < y3)
            {
                TS_ASSERT_EQUALS(generation, 3u);
            }
            else
            {
                TS_ASSERT_EQUALS(generation, 4u);
            }

            TS_ASSERT_DELTA(cells[i]->GetBirthTime(), 0.0, 1e-9);
        }

        // Create cells again with basic
        std::vector<CellPtr> new_cells;
        generator.GenerateBasic(new_cells, p_mesh->GetNumNodes());
        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(new_cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<new_cells.size(); i++)
        {
            TS_ASSERT_DELTA(new_cells[i]->GetBirthTime(), -(double)(i), 1e-9);
        }
    }

    void TestCryptCellsGeneratorWithTysonNovakCellCycleModel() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<TysonNovakCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, true);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        // Create cells again with basic
        std::vector<CellPtr> new_cells;
        generator.GenerateBasic(new_cells, p_mesh->GetNumNodes());

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(new_cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<new_cells.size(); i++)
        {
            TS_ASSERT_DELTA(new_cells[i]->GetBirthTime(), -(double)(i), 1e-9);
        }
    }

    void TestCryptCellsGeneratorWithWntCellCycleModel() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<WntCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, false);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_DELTA(cells[i]->GetBirthTime(), 0.0, 1e-9);
        }

        // Create cells again with basic
        std::vector<CellPtr> new_cells;
        generator.GenerateBasic(new_cells, p_mesh->GetNumNodes());
        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(new_cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<new_cells.size(); i++)
        {
            TS_ASSERT_DELTA(new_cells[i]->GetBirthTime(), -(double)(i), 1e-9);
        }
    }

    void TestCryptCellsGeneratorWithSimpleWntCellCycleModel() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, false);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_DELTA(cells[i]->GetBirthTime(), 0.0, 1e-9);
        }
    }

    void TestCryptCellsGeneratorWithStochasticWntCellCycleModel() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticWntCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, false);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_DELTA(cells[i]->GetBirthTime(), 0.0, 1e-9);
        }

        // Create cells again with basic
        std::vector<CellPtr> new_cells;
        generator.GenerateBasic(new_cells, p_mesh->GetNumNodes());
        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(new_cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<new_cells.size(); i++)
        {
            TS_ASSERT_DELTA(new_cells[i]->GetBirthTime(), -(double)(i), 1e-9);
        }
    }

    void TestCryptCellsGeneratorWithStochasticDurationGenerationBasedCellCycleModelAndVertexMesh() throw(Exception)
      {
          // Create mesh
          unsigned crypt_width = 4;
          unsigned crypt_height = 6;
          CylindricalHoneycombVertexMeshGenerator mesh_generator(crypt_width, crypt_height);
          Cylindrical2dVertexMesh* p_mesh = mesh_generator.GetCylindricalMesh();

          double y0 = 1.0;
          double y1 = 2.0;
          double y2 = 3.0;
          double y3 = 4.0;

          // Create cells
          std::vector<CellPtr> fixed_cells, stochastic_cells;
          CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> fixed_cells_generator;
          fixed_cells_generator.Generate(fixed_cells, p_mesh, std::vector<unsigned>(), true, y0, y1, y2, y3, true);

          CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> stochastic_cells_generator;
          stochastic_cells_generator.Generate(stochastic_cells, p_mesh, std::vector<unsigned>(), true, y0, y1, y2, y3, true);

          TS_ASSERT_EQUALS(fixed_cells.size(), p_mesh->GetNumElements());
          TS_ASSERT_EQUALS(stochastic_cells.size(), p_mesh->GetNumElements());

          // Test that cells were generated correctly
          for (unsigned i=0; i<fixed_cells.size(); i++)
          {
              double height = p_mesh->GetCentroidOfElement(i)[1];
              unsigned fixed_generation = static_cast<FixedDurationGenerationBasedCellCycleModel*>(fixed_cells[i]->GetCellCycleModel())->GetGeneration();
              unsigned stochastic_generation = static_cast<StochasticDurationGenerationBasedCellCycleModel*>(stochastic_cells[i]->GetCellCycleModel())->GetGeneration();

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

      void TestCryptCellsGeneratorWithSimpleWntCellCycleModelAndVertexMesh() throw(Exception)
      {
          // Create mesh
          unsigned crypt_width = 4;
          unsigned crypt_height = 6;
          CylindricalHoneycombVertexMeshGenerator mesh_generator(crypt_width, crypt_height);
          Cylindrical2dVertexMesh* p_mesh = mesh_generator.GetCylindricalMesh();

          // Create cells
          std::vector<CellPtr> cells;
          CryptCellsGenerator<SimpleWntCellCycleModel> cells_generator;
          cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, true);

          // Test that the correct number cells was generated
          TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumElements());
      }
};

#endif /*TESTCRYPTCELLSGENERATOR_HPP_*/
