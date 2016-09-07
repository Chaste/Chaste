/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTPDESFOROFFLATTICESIMULATIONS_HPP_
#define TESTPDESFOROFFLATTICESIMULATIONS_HPP_

#include <cxxtest/TestSuite.h>

#include "MeshBasedCellPopulation.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "ApoptoticCellProperty.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

///\todo move into cell_based/test/cell_based_pde and rename test suite (#2687)
class TestPdesForOffLatticeSimulations : public AbstractCellBasedTestSuite
{
public:

    void TestUniformSourceEllipticPde()
    {
        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up PDE
        UniformSourceEllipticPde<2> pde(-1.0);

        // Test Compute source term
        ChastePoint<2> unused_point;
        double value_at_elem_0 = pde.ComputeConstantInUSourceTerm(unused_point, NULL);
        double value_at_elem_1 = pde.ComputeLinearInUCoeffInSourceTerm(unused_point,p_mesh->GetElement(0));

        TS_ASSERT_DELTA(value_at_elem_0, 0.0, 1e-6);
        TS_ASSERT_DELTA(value_at_elem_1, -1.0, 1e-6);
    }

    void TestCellwiseSourceEllipticPde()
    {
        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Make one cell apoptotic
        MAKE_PTR(ApoptoticCellProperty, p_apoptotic_state);
        cells[0]->AddCellProperty(p_apoptotic_state);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up PDE
        CellwiseSourceEllipticPde<2> pde(cell_population, -1.0);

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

    void TestAveragedSourceEllipticPde()
    {
        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a coarse mesh - element 1 contains all the cells, element 0 contains none
        TetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructRegularSlabMesh(100.0, 100.0,100.0);

        // Top right
        TS_ASSERT_DELTA(coarse_mesh.GetElement(0)->CalculateCentroid()[0], 200.0/3.0, 0.1);
        TS_ASSERT_DELTA(coarse_mesh.GetElement(0)->CalculateCentroid()[1], 200.0/3.0, 0.1);

        // Bottom left
        TS_ASSERT_DELTA(coarse_mesh.GetElement(1)->CalculateCentroid()[0], 100.0/3.0, 0.1);
        TS_ASSERT_DELTA(coarse_mesh.GetElement(1)->CalculateCentroid()[1], 100.0/3.0, 0.1);

        // Set up PDE
        AveragedSourceEllipticPde<2> pde(cell_population, -1.0);
        pde.SetupSourceTerms(coarse_mesh);

        // Test compute source term
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
