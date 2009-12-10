/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef TESTCELLWISEDATAGRADIENT_HPP_
#define TESTCELLWISEDATAGRADIENT_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "MeshBasedTissueWithGhostNodes.hpp"
#include "CellwiseDataGradient.hpp"
#include "FixedDurationGenerationBasedCellCycleModelCellsGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "TrianglesMeshReader.hpp"

/**
 * This class contains tests for methods on the class CellwiseData.
 */
class TestCellwiseDataGradient : public AbstractCellBasedTestSuite
{

public:

    void TestCellwiseDataGradientVerySmallMesh() throw(Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create a tissue
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        MeshBasedTissue<2> tissue(mesh,cells);

        // Set up data: C(x,y) = x^2
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1);
        p_data->SetTissue(tissue);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            p_data->SetValue(x*x, mesh.GetNode(i));
        }

        CellwiseDataGradient<2> gradient;
        gradient.SetupGradients();

        // With the algorithm being used, the numerical gradient is (1,0)
        // for each of the nodes
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(0), 1.0, 1e-9);
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(1), 0.0, 1e-9);
        }

        CellwiseData<2>::Destroy();
    }


    void TestCellwiseDataGradientFineMesh() throw(Exception)
    {
        // Create a mesh: [0,2]x[0,2]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create a tissue
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        MeshBasedTissue<2> tissue(mesh,cells);

        //////////////////////////////////
        // C(x,y) = const
        //////////////////////////////////
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1);
        p_data->SetTissue(tissue);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, mesh.GetNode(i));
        }

        CellwiseDataGradient<2> gradient;
        gradient.SetupGradients();

        // Check gradient
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(0), 0.0, 1e-9);
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(1), 0.0, 1e-9);
        }

        //////////////////////////////////
        // C(x,y) = x-y
        //////////////////////////////////
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            p_data->SetValue(x-y, mesh.GetNode(i));
        }

        // Check gradient
        gradient.SetupGradients();
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(0),  1.0, 1e-9);
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(1), -1.0, 1e-9);
        }

        //////////////////////////////////
        // C(x,y) = x^2 - y^2
        //////////////////////////////////
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            p_data->SetValue(x*x - y*y, mesh.GetNode(i));
        }

        // Check gradient - here there is some numerical error
        gradient.SetupGradients();
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];

            double tol = 0.3;
            if (x==0 || x==2 || y==0 || y==2) //ie on boundary
            {
                tol = 0.6;
            }

            TS_ASSERT_DELTA( gradient.rGetGradient(i)(0),  2*x, tol);
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(1), -2*y, tol);
        }

        CellwiseData<2>::Destroy();
    }

    void TestCellwiseDataGradientWithGhostNodes() throw(Exception)
    {
        // Create a mesh: [0,2]x[0,2]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set boundary nodes to be ghost nodes, interior nodes to be cells
        std::vector<unsigned> cell_location_indices;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if ( !(mesh.GetNode(i)->IsBoundaryNode()) )
            {
                cell_location_indices.push_back(i);
            }
        }

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, cell_location_indices.size());

        // Create a tissue
        MeshBasedTissueWithGhostNodes<2> tissue(mesh, cells, cell_location_indices);

        // Create an instance of CellwiseData and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1);
        p_data->SetTissue(tissue);

        //////////////////////////////////
        // C(x,y) = x^2 - y^2
        //////////////////////////////////
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            if (mesh.GetNode(i)->IsBoundaryNode())
            {
                p_data->SetValue(DBL_MAX, mesh.GetNode(i));
            }
            else
            {
                p_data->SetValue(x*x - y*y, mesh.GetNode(i));
            }
        }

        // Check gradient - here there is some numerical error

        // The corner nodes are special because they have no adjacent real elements
        CellwiseDataGradient<2> gradient;
        gradient.SetupGradients();
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];

            if ( !mesh.GetNode(i)->IsBoundaryNode() ) // ie not ghost
            {
                int x_corner=0;

                // Work out if on left or right
                if (x==0.03125)
                {
                    x_corner = -1;
                }
                if (x==1.96875)
                {
                    x_corner = 1;
                }
                int y_corner=0;

                // Work out if on top or bottom
                if (y==0.03125)
                {
                    y_corner = -1;
                }
                if (y==1.96875)
                {
                    y_corner = 1;
                }

                switch (x_corner*y_corner)
                {
                    case 1: // bottom left or top right
                        TS_ASSERT_DELTA( gradient.rGetGradient(i)(0),  0.0, 1e-9);
                        TS_ASSERT_DELTA( gradient.rGetGradient(i)(1),  0.0, 1e-9);
                        break;
                    case -1: // bottom right or top left
                        TS_ASSERT_DELTA( gradient.rGetGradient(i)(0),  2.0, 1e-9);
                        TS_ASSERT_DELTA( gradient.rGetGradient(i)(1), -2.0, 1e-9);
                        break;
                    case 0: // otherwise
                        TS_ASSERT_DELTA( gradient.rGetGradient(i)(0),  2*x, 0.3);
                        TS_ASSERT_DELTA( gradient.rGetGradient(i)(1), -2*y, 0.3);
                        break;
                    default:
                        break;
                }
            }
        }

        CellwiseData<2>::Destroy();
    }

};

#endif /*TESTCELLWISEDATAGRADIENT_HPP_*/
