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

// For archiving?
#include "WildTypeCellMutationState.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

class TestCryptCellsGenerator : public AbstractCellBasedTestSuite
{
public:

    void TestCryptCellsGeneratorWithFixedG1GenerationalCellCycleModel()
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

        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, true, y0, y1, y2, y3);

        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        // Test that cells were generated correctly
        for (unsigned i=0; i<cells.size(); i++)
        {
            double height = p_mesh->GetNode(i)->rGetLocation()[1];
            unsigned generation = static_cast<FixedG1GenerationalCellCycleModel*>(cells[i]->GetCellCycleModel())->GetGeneration();

            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetDimension(), 2u);

            if (height <= y0)
            {
                TS_ASSERT_EQUALS(generation, 0u);
                TS_ASSERT_EQUALS(cells[i]->GetCellProliferativeType()->IsType<StemCellProliferativeType>(), true);
            }
            else if (height < y1)
            {
                TS_ASSERT_EQUALS(generation, 1u);
                TS_ASSERT_EQUALS(cells[i]->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);
            }
            else if (height < y2)
            {
                TS_ASSERT_EQUALS(generation, 2u);
                TS_ASSERT_EQUALS(cells[i]->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);
            }
            else if (height < y3)
            {
                TS_ASSERT_EQUALS(generation, 3u);
                TS_ASSERT_EQUALS(cells[i]->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(generation, 4u);
                TS_ASSERT_EQUALS(cells[i]->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>(), true);
            }
        }
    }

    void TestCryptCellsGeneratorWithUniformG1GenerationalCellCycleModel()
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<UniformG1GenerationalCellCycleModel> generator;
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
            unsigned generation = static_cast<UniformG1GenerationalCellCycleModel*>(cells[i]->GetCellCycleModel())->GetGeneration();

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

    void TestCryptCellsGeneratorWithTysonNovakCellCycleModel()
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

    void TestCryptCellsGeneratorWithWntCellCycleModel()
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

    void TestCryptCellsGeneratorWithSimpleWntCellCycleModel()
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

    void TestCryptCellsGeneratorWithStochasticWntCellCycleModel()
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

    void TestCryptCellsGeneratorWithUniformG1GenerationalCellCycleModelAndVertexMesh()
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
          CryptCellsGenerator<FixedG1GenerationalCellCycleModel> fixed_cells_generator;
          fixed_cells_generator.Generate(fixed_cells, p_mesh, std::vector<unsigned>(), true, y0, y1, y2, y3, true);

          CryptCellsGenerator<UniformG1GenerationalCellCycleModel> stochastic_cells_generator;
          stochastic_cells_generator.Generate(stochastic_cells, p_mesh, std::vector<unsigned>(), true, y0, y1, y2, y3, true);

          TS_ASSERT_EQUALS(fixed_cells.size(), p_mesh->GetNumElements());
          TS_ASSERT_EQUALS(stochastic_cells.size(), p_mesh->GetNumElements());

          // Test that cells were generated correctly
          for (unsigned i=0; i<fixed_cells.size(); i++)
          {
              double height = p_mesh->GetCentroidOfElement(i)[1];
              unsigned fixed_generation = static_cast<FixedG1GenerationalCellCycleModel*>(fixed_cells[i]->GetCellCycleModel())->GetGeneration();
              unsigned stochastic_generation = static_cast<UniformG1GenerationalCellCycleModel*>(stochastic_cells[i]->GetCellCycleModel())->GetGeneration();

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

      void TestCryptCellsGeneratorWithSimpleWntCellCycleModelAndVertexMesh()
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
