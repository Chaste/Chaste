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

#ifndef TESTCELLWISEDATAGRADIENT_HPP_
#define TESTCELLWISEDATAGRADIENT_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include <fstream>

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CellwiseDataGradient.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "TrianglesMeshReader.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * This class contains tests for methods on the class CellwiseDataGradient.
 */
class TestCellwiseDataGradient : public AbstractCellBasedTestSuite
{

public:

    void TestCellwiseDataGradientVerySmallMesh()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create a cell population
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up data: C(x,y) = x^2
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(mesh.GetNode(i)->GetIndex());
            p_cell->GetCellData()->SetItem("x^2", x*x);
        }

        CellwiseDataGradient<2> gradient;
        gradient.SetupGradients(cell_population, "x^2");

        // With the algorithm being used, the numerical gradient is (1,0)
        // for each of the nodes
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(gradient.rGetGradient(i)(0), 1.0, 1e-9);
            TS_ASSERT_DELTA(gradient.rGetGradient(i)(1), 0.0, 1e-9);
        }
    }

    void TestCellwiseDataGradientFineMesh()
    {
        // Create a mesh: [0,2]x[0,2]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create a cell population
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        //////////////////////////////////
        // C(x,y) = const
        //////////////////////////////////
        cell_population.SetDataOnAllCells("const", 1.0);

        CellwiseDataGradient<2> gradient;
        gradient.SetupGradients(cell_population, "const");

        // Check gradient
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(gradient.rGetGradient(i)(0), 0.0, 1e-9);
            TS_ASSERT_DELTA(gradient.rGetGradient(i)(1), 0.0, 1e-9);
        }

        //////////////////////////////////
        // Combined setup for
        // C(x,y) = x-y
        // and
        // C(x,y) = x^2 - y^2
        //////////////////////////////////
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(mesh.GetNode(i)->GetIndex());
            p_cell->GetCellData()->SetItem("x-y",       x-y);
            p_cell->GetCellData()->SetItem("x^2 - y^2", x*x - y*y);
        }

        // Check gradient
        gradient.SetupGradients(cell_population, "x-y");
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(gradient.rGetGradient(i)(0),  1.0, 1e-9);
            TS_ASSERT_DELTA(gradient.rGetGradient(i)(1), -1.0, 1e-9);
        }

        // Check gradient - here there is some numerical error
        gradient.SetupGradients(cell_population, "x^2 - y^2");
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];

            double tol = 0.3;
            if (x==0 || x==2 || y==0 || y==2) //ie on boundary
            {
                tol = 0.6;
            }

            TS_ASSERT_DELTA(gradient.rGetGradient(i)(0),  2*x, tol);
            TS_ASSERT_DELTA(gradient.rGetGradient(i)(1), -2*y, tol);
        }
    }

///\todo uncomment this test or remove it
//    void TestCellwiseDataGradientWithGhostNodes()
//    {
//        // Create a mesh: [0,2]x[0,2]
//        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
//        MutableMesh<2,2> mesh;
//        mesh.ConstructFromMeshReader(mesh_reader);
//
//        // Set boundary nodes to be ghost nodes, interior nodes to be cells
//        std::vector<unsigned> cell_location_indices;
//        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
//        {
//            if (!(mesh.GetNode(i)->IsBoundaryNode()))
//            {
//                cell_location_indices.push_back(i);
//            }
//        }
//
//        // Set up cells
//        std::vector<CellPtr> cells;
//        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
//        cells_generator.GenerateBasic(cells, cell_location_indices.size());
//
//        // Create a cell population
//        MeshBasedCellPopulationWithGhostNodes<2> cell_population(mesh, cells, cell_location_indices);
//
//        // Create an instance of CellwiseData and associate it with the cell population
//        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
//        p_data->SetPopulationAndNumVars(&cell_population, 1);
//
//        //////////////////////////////////
//        // C(x,y) = x^2 - y^2
//        //////////////////////////////////
//        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
//        {
//            double x = mesh.GetNode(i)->rGetLocation()[0];
//            double y = mesh.GetNode(i)->rGetLocation()[1];
//            if (mesh.GetNode(i)->IsBoundaryNode())
//            {
//                p_data->SetValue(DBL_MAX, mesh.GetNode(i)->GetIndex());
//            }
//            else
//            {
//                p_data->SetValue(x*x - y*y, mesh.GetNode(i)->GetIndex());
//            }
//        }
//
//        // Check gradient - here there is some numerical error
//
//        // The corner nodes are special because they have no adjacent real elements
//        CellwiseDataGradient<2> gradient;
//        gradient.SetupGradients(cell_population);
//        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
//        {
//            double x = mesh.GetNode(i)->rGetLocation()[0];
//            double y = mesh.GetNode(i)->rGetLocation()[1];
//
//            if (!(mesh.GetNode(i)->IsBoundaryNode())) // i.e. not ghost
//            {
//                int x_corner = 0;
//
//                // Work out if on left or right
//                if (x == 0.03125)
//                {
//                    x_corner = -1;
//                }
//                if (x == 1.96875)
//                {
//                    x_corner = 1;
//                }
//                int y_corner=0;
//
//                // Work out if on top or bottom
//                if (y == 0.03125)
//                {
//                    y_corner = -1;
//                }
//                if (y == 1.96875)
//                {
//                    y_corner = 1;
//                }
//
//                switch (x_corner*y_corner)
//                {
//                    case 1: // bottom left or top right
//                        TS_ASSERT_DELTA(gradient.rGetGradient(i)(0),  0.0, 1e-9);
//                        TS_ASSERT_DELTA(gradient.rGetGradient(i)(1),  0.0, 1e-9);
//                        break;
//                    case -1: // bottom right or top left
//                        TS_ASSERT_DELTA(gradient.rGetGradient(i)(0),  2.0, 1e-9);
//                        TS_ASSERT_DELTA(gradient.rGetGradient(i)(1), -2.0, 1e-9);
//                        break;
//                    case 0: // otherwise
//                        TS_ASSERT_DELTA(gradient.rGetGradient(i)(0),  2*x, 0.3);
//                        TS_ASSERT_DELTA(gradient.rGetGradient(i)(1), -2*y, 0.3);
//                        break;
//                    default:
//                        break;
//                }
//            }
//        }
//
//        CellwiseData<2>::Destroy();
//    }
};

#endif /*TESTCELLWISEDATAGRADIENT_HPP_*/
