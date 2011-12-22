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

#ifndef TESTPDESFOROFFLATTICESIMULATIONS_HPP_
#define TESTPDESFOROFFLATTICESIMULATIONS_HPP_

#include <cxxtest/TestSuite.h>

#include "MeshBasedCellPopulation.hpp"
#include "SimpleUniformSourcePde.hpp"
#include "CellwiseSourcePde.hpp"
#include "AveragedSourcePde.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

class TestPdesForOffLatticeSimulations : public AbstractCellBasedTestSuite
{
public:

    void TestSimpleUniformSourcePde()
    {
        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up PDE
        SimpleUniformSourcePde<2> pde(-1.0);

        // Test Compute source term
        ChastePoint<2> unused_point;
        double value_at_elem_0 = pde.ComputeConstantInUSourceTerm(unused_point, NULL);
        double value_at_elem_1 = pde.ComputeLinearInUCoeffInSourceTerm(unused_point,p_mesh->GetElement(0));

        TS_ASSERT_DELTA(value_at_elem_0, 0.0, 1e-6);
        TS_ASSERT_DELTA(value_at_elem_1, -1.0, 1e-6);
    }

    void TestCellwiseSourcePde()
    {
        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Make one cell apoptotic
        MAKE_PTR(ApoptoticCellProperty, p_apoptotic_state);
        cells[0]->AddCellProperty(p_apoptotic_state);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up PDE
        CellwiseSourcePde<2> pde(cell_population, -1.0);

        // Test compute source terms
        ChastePoint<2> unused_point;
        double constant_in_u_source_term = pde.ComputeConstantInUSourceTerm(unused_point, NULL);

        TS_ASSERT_DELTA(constant_in_u_source_term, 0.0, 1e-6);

        Node<2>* p_node_0 = cell_population.GetNodeCorrespondingToCell(cell_population.GetCellUsingLocationIndex(0));
        Node<2>* p_node_1 = cell_population.GetNodeCorrespondingToCell(cell_population.GetCellUsingLocationIndex(1));

        double source_term_at_node_0 = pde.ComputeLinearInUCoeffInSourceTermAtNode(*p_node_0);
        double source_term_at_node_1 = pde.ComputeLinearInUCoeffInSourceTermAtNode(*p_node_1);

        TS_ASSERT_DELTA(source_term_at_node_0, 0.0, 1e-6);
        TS_ASSERT_DELTA(source_term_at_node_1, -1.0, 1e-6);

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    void TestAveragedSourcePde()
    {
        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a coarse mesh - element 1 contains all the cells,
        // element 0 contains none
        TetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructRegularSlabMesh(100.0, 100.0,100.0);

        //Top right
        TS_ASSERT_DELTA(coarse_mesh.GetElement(0)->CalculateCentroid()[0], 200.0/3.0, 0.1);
        TS_ASSERT_DELTA(coarse_mesh.GetElement(0)->CalculateCentroid()[1], 200.0/3.0, 0.1);
        //Bottom left
        TS_ASSERT_DELTA(coarse_mesh.GetElement(1)->CalculateCentroid()[0], 100.0/3.0, 0.1);
        TS_ASSERT_DELTA(coarse_mesh.GetElement(1)->CalculateCentroid()[1], 100.0/3.0, 0.1);
        // Set up PDE
        AveragedSourcePde<2> pde(cell_population, -1.0);
        pde.SetupSourceTerms(coarse_mesh);

        // Test Compute source term
        ChastePoint<2> unused_point;
        double value_at_elem_0 = pde.ComputeLinearInUCoeffInSourceTerm(unused_point,coarse_mesh.GetElement(0));
        double value_at_elem_1 = pde.ComputeLinearInUCoeffInSourceTerm(unused_point,coarse_mesh.GetElement(1));

        TS_ASSERT_DELTA(value_at_elem_0, 0.0, 1e-6);
        c_matrix <double, 2, 2> jacobian;
        double det;
        coarse_mesh.GetElement(1)->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(value_at_elem_1, -(cell_population.GetNumRealCells()/coarse_mesh.GetElement(1)->GetVolume(det)), 1e-6);

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTPDESFOROFFLATTICESIMULATIONS_HPP_*/
