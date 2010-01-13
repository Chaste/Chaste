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
#ifndef TESTNUTRIENTPDES_HPP_
#define TESTNUTRIENTPDES_HPP_

#include <cxxtest/TestSuite.h>

#include "MeshBasedTissue.hpp"
#include "SimpleNutrientPde.hpp"
#include "CellwiseNutrientSinkPde.hpp"
#include "AveragedSinksPde.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModelCellsGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"

class TestNutrientPdes : public AbstractCellBasedTestSuite
{
public:

    void TestSimpleNutrientPde()
    {
        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0, false);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up PDE
        SimpleNutrientPde<2> pde(1.0);

        // Test Compute source term
        ChastePoint<2> unused_point;
        double value_at_elem_0 = pde.ComputeConstantInUSourceTerm(unused_point);
        double value_at_elem_1 = pde.ComputeLinearInUCoeffInSourceTerm(unused_point,p_mesh->GetElement(0));

        TS_ASSERT_DELTA(value_at_elem_0, 0.0, 1e-6);
        TS_ASSERT_DELTA(value_at_elem_1, -1.0, 1e-6);
    }

    void TestCellwiseNutrientSinkPde()
    {
        // Set up tissue
        HoneycombMeshGenerator generator(5, 5, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Make one cell apoptotic
        cells[0].SetCellProliferativeType(APOPTOTIC);

        MeshBasedTissue<2> tissue(*p_mesh, cells);

        // Set up PDE
        CellwiseNutrientSinkPde<2> pde(tissue, 1.0);

        // Test compute source terms
        ChastePoint<2> unused_point;
        double constant_in_u_source_term = pde.ComputeConstantInUSourceTerm(unused_point);

        TS_ASSERT_DELTA(constant_in_u_source_term, 0.0, 1e-6);

        Node<2>* p_node_0 = tissue.GetNodeCorrespondingToCell(tissue.rGetCellUsingLocationIndex(0));
        Node<2>* p_node_1 = tissue.GetNodeCorrespondingToCell(tissue.rGetCellUsingLocationIndex(1));

        double source_term_at_node_0 = pde.ComputeLinearInUCoeffInSourceTermAtNode(*p_node_0);
        double source_term_at_node_1 = pde.ComputeLinearInUCoeffInSourceTermAtNode(*p_node_1);

        TS_ASSERT_DELTA(source_term_at_node_0, 0.0, 1e-6);
        TS_ASSERT_DELTA(source_term_at_node_1, -1.0, 1e-6);

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    void TestAveragedSinksPde()
    {
        // Set up tissue
        HoneycombMeshGenerator generator(5, 5, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        // Create a coarse mesh - element 1 contains all the cells,
        // element 0 contains none
        TetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructRectangularMesh(1,1);
        coarse_mesh.Scale(100,100);
        //Top right
        TS_ASSERT_DELTA(coarse_mesh.GetElement(0)->CalculateCentroid()[0], 200.0/3.0, 0.1);
        TS_ASSERT_DELTA(coarse_mesh.GetElement(0)->CalculateCentroid()[1], 200.0/3.0, 0.1);
        //Bottom left
        TS_ASSERT_DELTA(coarse_mesh.GetElement(1)->CalculateCentroid()[0], 100.0/3.0, 0.1);
        TS_ASSERT_DELTA(coarse_mesh.GetElement(1)->CalculateCentroid()[1], 100.0/3.0, 0.1);
        // Set up PDE
        AveragedSinksPde<2> pde(tissue, -1.0);
        pde.SetupSourceTerms(coarse_mesh);

        // Test Compute source term
        ChastePoint<2> unused_point;
        double value_at_elem_0 = pde.ComputeLinearInUCoeffInSourceTerm(unused_point,coarse_mesh.GetElement(0));
        double value_at_elem_1 = pde.ComputeLinearInUCoeffInSourceTerm(unused_point,coarse_mesh.GetElement(1));

        TS_ASSERT_DELTA(value_at_elem_0, 0.0, 1e-6);
        c_matrix <double, 2, 2> jacobian;
        double det;
        coarse_mesh.GetElement(1)->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(value_at_elem_1, -(tissue.GetNumRealCells()/coarse_mesh.GetElement(1)->GetVolume(det)), 1e-6);

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTNUTRIENTPDES_HPP_*/
