/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef TESTIMMERSEDBOUNDARYCELLPOPULATION_HPP_
#define TESTIMMERSEDBOUNDARYCELLPOPULATION_HPP_

// Needed for the test environment
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

// Includes from trunk
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "ArchiveOpener.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellDivisionLocationsWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellLabel.hpp"
#include "CellLabelWriter.hpp"
#include "CellLocationIndexWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "CellVolumesWriter.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "FixedVertexBasedDivisionRule.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "ImmersedBoundaryEnumerations.hpp"
#include "OffLatticeSimulation.hpp"
#include "ShortAxisImmersedBoundaryDivisionRule.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"

// Includes from Immersed Boundary
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryLinearMembraneForce.hpp"
#include "ImmersedBoundaryLinearInteractionForce.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundaryCellPopulation : public AbstractCellBasedTestSuite
{
public:

    void TestGetAndSetMethods()
    {
        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        // Test that GetDampingConstant() returns the correct value
        for (unsigned node_index = 0; node_index < cell_population.GetNumNodes(); ++node_index)
        {
            TS_ASSERT_DELTA(cell_population.GetDampingConstant(node_index), 0.0, 1e-6);
        }

        // Test that GetIntrinsicSpacing() returns the correct value
        TS_ASSERT_DELTA(cell_population.GetIntrinsicSpacing(), 0.01, 1e-6);

        // Test GetInteractionDistance() and SetInteractionDistance() work correctly
        TS_ASSERT_DELTA(cell_population.GetInteractionDistance(), 0.0139, 1e-4);

        cell_population.SetInteractionDistance(0.1234);
        TS_ASSERT_DELTA(cell_population.GetInteractionDistance(), 0.1234, 1e-6);

        // Test DoesPopulationHaveActiveSources() and SetIfPopulationHasActiveSources() work correctly
        TS_ASSERT_EQUALS(cell_population.DoesPopulationHaveActiveSources(), false);

        cell_population.SetIfPopulationHasActiveSources(true);
        TS_ASSERT_EQUALS(cell_population.DoesPopulationHaveActiveSources(), true);

        TS_ASSERT_EQUALS(cell_population.GetReMeshFrequency(), UINT_MAX);
        cell_population.SetReMeshFrequency(5u);
        TS_ASSERT_EQUALS(cell_population.GetReMeshFrequency(), 5u);

        TS_ASSERT_EQUALS(cell_population.mOutputNodeRegionToVtk, false);
        cell_population.SetOutputNodeRegionToVtk(true);
        TS_ASSERT_EQUALS(cell_population.mOutputNodeRegionToVtk, true);

        auto division_rule = boost::shared_ptr<ShortAxisImmersedBoundaryDivisionRule<2>>(new ShortAxisImmersedBoundaryDivisionRule<2>());
        cell_population.SetImmersedBoundaryDivisionRule(division_rule);
        TS_ASSERT_EQUALS(cell_population.GetImmersedBoundaryDivisionRule().get(), division_rule.get());
    }

    void TestMeshMethods()
    {
        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        // Test that GetNumElements() and GetNumNodes() return the correct values
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 5u);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 583u);

        // Test that rGetMesh() returns the mesh correctly
        ImmersedBoundaryMesh<2, 2>& r_mesh = cell_population.rGetMesh();
        TS_ASSERT_EQUALS(r_mesh.GetNumNodes(), 583u);
        TS_ASSERT_EQUALS(r_mesh.GetNumGridPtsX(), 128u);
        TS_ASSERT_EQUALS(r_mesh.GetNumGridPtsY(), 128u);
        TS_ASSERT_DELTA(r_mesh.GetCharacteristicNodeSpacing(), 0.0115, 1e-4);
        TS_ASSERT_DELTA(r_mesh.GetSpacingRatio(), 1.4805, 1e-4);

        // Test that GetElement() returns an element correctly
        ImmersedBoundaryElement<2, 2>* p_element_0 = cell_population.GetElement(0);
        TS_ASSERT_EQUALS(p_element_0->GetNumNodes(), 100u);

        // Test that GetElement() returns the lamina correctly
        ImmersedBoundaryElement<1, 2>* p_lamina_0 = cell_population.GetLamina(0);
        TS_ASSERT_EQUALS(p_lamina_0->GetNumNodes(), 83u);

        // Test that GetElementCorrespondingToCell() returns an element correctly
        CellPtr p_cell_0 = *(cell_population.Begin());
        ImmersedBoundaryElement<2, 2>* p_element_0_again = cell_population.GetElementCorrespondingToCell(p_cell_0);
        TS_ASSERT_EQUALS(p_element_0_again, p_element_0);

        // Test that GetNode() returns a node correctly
        Node<2>* p_node_3 = cell_population.GetNode(3);
        TS_ASSERT_EQUALS(p_node_3->IsBoundaryNode(), true);
        TS_ASSERT_DELTA(p_node_3->rGetLocation()[0], 0.0421, 1e-4);
        TS_ASSERT_DELTA(p_node_3->rGetLocation()[1], 0.2800, 1e-4);

        // Test GetLocationOfCellCentre() returns the correct location
        c_vector<double, 2> cell_0_centre = cell_population.GetLocationOfCellCentre(p_cell_0);
        TS_ASSERT_DELTA(cell_0_centre[0], 0.1950, 1e-6);
        TS_ASSERT_DELTA(cell_0_centre[1], 0.505641, 1e-6);

        // Test GetWidth() returns the correct values
        TS_ASSERT_DELTA(cell_population.GetWidth(0), 0.9891, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetWidth(1), 0.4481, 1e-4);

        // Test GetVolumeOfCell() returns the correct value
        CellPtr p_cell_1 = *(++(cell_population.Begin()));
        double cell_1_volume = cell_population.GetVolumeOfCell(p_cell_1);
        TS_ASSERT_DELTA(cell_1_volume, 0.0774, 1e-4);

        // Test SetNode() works correctly
        c_vector<double, 2> new_location;
        new_location[0] = 0.0;
        new_location[1] = 0.0;
        new_location = cell_population.GetNode(0)->rGetLocation();
        new_location[0] += 1e-2;
        new_location[1] += 1e-2;
        ChastePoint<2> new_location_point(new_location);
        cell_population.SetNode(0, new_location_point);

        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[0], new_location[0], 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[1], new_location[1], 1e-12);

        // Test GetNumLaminas
        TS_ASSERT_EQUALS(cell_population.GetNumLaminas(), 1);

        // Test adding ndoe from population - doesn't do anything
        TS_ASSERT_EQUALS(cell_population.AddNode(nullptr), 0);

        TetrahedralMesh<2, 2>* p_tet_mesh = cell_population.GetTetrahedralMeshForPdeModifier();
        TS_ASSERT_DIFFERS(p_tet_mesh, nullptr);

        ImmersedBoundaryMesh<3, 3> ib_mesh_3d;
        std::vector<CellPtr> cells_3d;
        ImmersedBoundaryCellPopulation<3> cell_population_3d(ib_mesh_3d, cells_3d);
        TS_ASSERT_THROWS_CONTAINS(cell_population_3d.GetTetrahedralMeshForPdeModifier(), "only implemented in 2D");

        // Test adding new cell
        boost::shared_ptr<AbstractCellProperty> p_wildtype(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        UniformCellCycleModel* p_model = new UniformCellCycleModel();
        CellPtr p_cell(new Cell(p_wildtype, p_model));
        TS_ASSERT_THROWS_NOTHING(cell_population.AddCell(p_cell, *(cell_population.rGetCells().begin())));

        TS_ASSERT(cell_population.IsPdeNodeAssociatedWithNonApoptoticCell(0));
        TS_ASSERT(!cell_population.IsCellOnBoundary(p_cell));

        // Tidy up
        delete p_tet_mesh;
    }

    void TestStepSizeException()
    {
        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetThrowsStepSizeException(true);

        c_vector<double, 2> displacement;
        displacement[0] = 0.8;
        displacement[1] = 0.8;
        TS_ASSERT_THROWS_ANYTHING(cell_population.CheckForStepSizeException(0, displacement, 0.1));
    }

    void TestValidateException()
    {
        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 1, p_diff_type);

        TS_ASSERT_THROWS_CONTAINS(ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells), "does not appear to have a cell associated");
    }

    void TestOverlyLargeDisplacements()
    {
        { // UpdateNodeLocations() coverage > 10x

            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);

            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            auto& r_velocity_field = mesh.rGetModifiable2dVelocityGrids();
            for (unsigned dim = 0; dim < 2; ++dim)
            {
                for (unsigned x = 0; x < 10; ++x)
                {
                    for (unsigned y = 0; y < 10; ++y)
                    {
                        r_velocity_field[dim][x][y] = 10000.0;
                    }
                }
            }

            mesh.SetCharacteristicNodeSpacing(0.000001);
            cell_population.SetReMeshFrequency(1);
            TS_ASSERT_THROWS_CONTAINS(cell_population.UpdateNodeLocations(0.1), "10x Characteristic");
        }

        { // UpdateNodeLocations() coverage

            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);

            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            auto& r_velocity_field = mesh.rGetModifiable2dVelocityGrids();
            for (unsigned dim = 0; dim < 2; ++dim)
            {
                for (unsigned x = 0; x < 10; ++x)
                {
                    for (unsigned y = 0; y < 10; ++y)
                    {
                        r_velocity_field[dim][x][y] = 0.1;
                    }
                }
            }

            mesh.SetCharacteristicNodeSpacing(0.01);
            cell_population.SetReMeshFrequency(1);
            TS_ASSERT_THROWS_NOTHING(cell_population.UpdateNodeLocations(0.1));
        }

        { // UpdateNodeLocations() coverage > 10x with fluid sources

            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            // Position all nodes near the middle
            nodes.push_back(new Node<2>(0, true, 0.58, 0.58));
            nodes.push_back(new Node<2>(1, true, 0.59, 0.59));
            nodes.push_back(new Node<2>(2, true, 0.58, 0.59));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));
            auto source = std::make_shared<FluidSource<2>>(0, 0.01, 0.01);
            elems.back()->SetFluidSource(source);

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);

            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);
            cell_population.SetIfPopulationHasActiveSources(true);

            auto& r_velocity_field = mesh.rGetModifiable2dVelocityGrids();
            for (unsigned dim = 0; dim < 2; ++dim)
            {
                r_velocity_field[dim][0][0] = 0.1;
                r_velocity_field[dim][0][1] = 0.1;
                r_velocity_field[dim][1][0] = 0.1;
                r_velocity_field[dim][1][1] = 0.1;
            }

            mesh.SetCharacteristicNodeSpacing(0.000001);
            cell_population.SetReMeshFrequency(1);
            TS_ASSERT_THROWS_CONTAINS(cell_population.UpdateNodeLocations(0.1), "Sources are moving more than 10x Characteristic");
        }

        { // UpdateNodeLocations() coverage with fluid sources

            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            // Position all nodes in top right corner
            nodes.push_back(new Node<2>(0, true, 0.58, 0.58));
            nodes.push_back(new Node<2>(1, true, 0.59, 0.59));
            nodes.push_back(new Node<2>(2, true, 0.58, 0.59));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));
            // Position fluid source in bottom left corner
            auto source = std::make_shared<FluidSource<2>>(0, 0.01, 0.01);
            elems.back()->SetFluidSource(source);

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);

            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);
            cell_population.SetIfPopulationHasActiveSources(true);

            auto& r_velocity_field = mesh.rGetModifiable2dVelocityGrids();
            for (unsigned dim = 0; dim < 2; ++dim)
            {
                r_velocity_field[dim][0][0] = 0.1;
                r_velocity_field[dim][0][1] = 0.1;
                r_velocity_field[dim][1][0] = 0.1;
                r_velocity_field[dim][1][1] = 0.1;
            }

            mesh.SetCharacteristicNodeSpacing(0.01);
            cell_population.SetReMeshFrequency(1);
            TS_ASSERT_THROWS_NOTHING(cell_population.UpdateNodeLocations(0.1));

            // Remesh coverage
            SimulationTime::Instance()->IncrementTimeOneStep();
            cell_population.UpdateNodeLocations(0.1);
        }
    }

    void TestWritersWithImmersedBoundaryCellPopulation()
    {
        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        boost::shared_ptr<AbstractCellProperty> p_stem(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_transit(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_diff(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_wildtype(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        std::vector<CellPtr> cells;
        for (unsigned elem_index = 0; elem_index < p_mesh->GetNumElements(); ++elem_index)
        {
            UniformCellCycleModel* p_model = new UniformCellCycleModel();

            CellPtr p_cell(new Cell(p_wildtype, p_model));
            if (elem_index % 3 == 0)
            {
                p_cell->SetCellProliferativeType(p_stem);
            }
            else if (elem_index % 3 == 1)
            {
                p_cell->SetCellProliferativeType(p_transit);
            }
            else
            {
                p_cell->SetCellProliferativeType(p_diff);
            }

            double birth_time = 0.0 - elem_index;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.InitialiseCells();
        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "ImmersedBoundaryCellPopulation-2");

        // Allocate some cells to have a different cell mutation state, cell label or apoptotic cell property
        cell_population.GetCellPropertyRegistry()->Get<WildTypeCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_apc1(cell_population.GetCellPropertyRegistry()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc2(cell_population.GetCellPropertyRegistry()->Get<ApcTwoHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat1(cell_population.GetCellPropertyRegistry()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(cell_population.GetCellPropertyRegistry()->Get<ApoptoticCellProperty>());
        boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());

        cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
        cell_population.GetCellUsingLocationIndex(1)->SetMutationState(p_apc1);
        cell_population.GetCellUsingLocationIndex(2)->SetMutationState(p_apc2);;
        cell_population.GetCellUsingLocationIndex(3)->SetMutationState(p_bcat1);;
        cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(p_apoptotic_state);;
        cell_population.SetCellAncestorsToLocationIndices();

        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        using WriterType = CellDivisionLocationsWriter<2, 2>;
        MAKE_PTR(WriterType, cdlw);
        cell_population.AddCellPopulationEventWriter(cdlw);

        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellLabelWriter>();
        cell_population.AddCellWriter<CellLocationIndexWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Coverage of writing CellData to VTK
        for (auto cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("var0", 0.0);
            cell_iter->GetCellData()->SetItem("var1", 3.0);
            cell_iter->GetCellData()->SetItem("target area", 0.1);
        }

        // Set up node regions
        for (auto iter = p_mesh->GetNodeIteratorBegin();
             iter != p_mesh->GetNodeIteratorEnd();
             ++iter)
        {
            iter->SetRegion(LAMINA_REGION);
        }

        std::string output_directory = "TestImmersedBoundaryPopulationWriters";
        OutputFileHandler output_file_handler(output_directory, false);

        cell_population.OpenWritersFiles(output_file_handler);
        cell_population.WriteResultsToFiles(output_directory);
        cell_population.SetOutputNodeRegionToVtk(true);

        SimulationTime::Instance()->IncrementTimeOneStep();
        cell_population.Update();
        cell_population.WriteResultsToFiles(output_directory);

        cell_population.CloseWritersFiles();

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.viznodes").CompareFiles();
        FileComparison(results_dir + "results.vizelements", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizelements").CompareFiles();
        FileComparison(results_dir + "cellages.dat", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/cellages.dat").CompareFiles();
        FileComparison(results_dir + "results.vizancestors", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizancestors").CompareFiles();
        FileComparison(results_dir + "loggedcell.dat", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/loggedcell.dat").CompareFiles();
        FileComparison(results_dir + "results.vizlabels", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizlabels").CompareFiles();
        FileComparison(results_dir + "results.vizlocationindices", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizlocationindices").CompareFiles();
        FileComparison(results_dir + "results.vizmutationstates", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizmutationstates").CompareFiles();
        FileComparison(results_dir + "results.vizcellphases", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizcellphases").CompareFiles();
        FileComparison(results_dir + "results.vizcelltypes", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizcelltypes").CompareFiles();
        FileComparison(results_dir + "cellareas.dat", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/cellareas.dat").CompareFiles();
        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/cellmutationstates.dat").CompareFiles();
        FileComparison(results_dir + "celltypes.dat", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/celltypes.dat").CompareFiles();

        // Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");
        cell_population.OutputCellPopulationParameters(parameter_file);
        parameter_file->close();

        FileComparison( results_dir + "results.parameters", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.parameters").CompareFiles();

#ifdef CHASTE_VTK
        // Test that VTK writer has produced some files

        // Initial condition file
        FileFinder vtk_file(results_dir + "results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file.Exists());

        // Final file
        FileFinder vtk_file2(results_dir + "results_1.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file2.Exists());

        // PVD file
        FileFinder vtk_file3(results_dir + "results.pvd", RelativeTo::Absolute);
        TS_ASSERT(vtk_file3.Exists());
 #endif //CHASTE_VTK
    }

    void TestArchiving()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "immersed_boundary_cell_population_2d.arch";

        // The following line is required because the loading of a cell population
        // is usually called by the method CellBasedSimulation::Load()
        ArchiveLocationInfo::SetMeshFilename("immersed_boundary_mesh_2d");

        // Create mesh
        ImmersedBoundaryMeshReader<2, 2> mesh_reader("mesh/test/data/ib_mesh_2d");
        ImmersedBoundaryMesh<2, 2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Archive cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create an Immersed Boundary cell population object
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumElements());

            // Create cell population
            AbstractCellPopulation<2>* const p_cell_population = new ImmersedBoundaryCellPopulation<2>(mesh, cells);

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // loop over them to run to time 0.0;
            for (auto cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                cell_iter->ReadyToDivide();
            }

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            auto ib_pop = dynamic_cast<ImmersedBoundaryCellPopulation<2>*>(p_cell_population);
            ib_pop->SetReMeshFrequency(4);
            ib_pop->SetInteractionDistance(0.3);
            ib_pop->SetIfPopulationHasActiveSources(true);

            // Archive the cell population
            (*p_arch) << static_cast<const SimulationTime&>(*p_simulation_time);
            (*p_arch) << p_cell_population;

            // Tidy up
            SimulationTime::Destroy();
            delete p_cell_population;
        }

        // Restore cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            AbstractCellPopulation<2>* p_cell_population;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore the cell population
            (*p_arch) >> *p_simulation_time;
            (*p_arch) >> p_cell_population;

            ImmersedBoundaryCellPopulation<2>* p_static_population = static_cast<ImmersedBoundaryCellPopulation<2>*>(p_cell_population);
            TS_ASSERT_EQUALS(p_static_population->GetReMeshFrequency(), 4);
            TS_ASSERT_EQUALS(p_static_population->GetInteractionDistance(), 0.3);
            TS_ASSERT_EQUALS(p_static_population->DoesPopulationHaveActiveSources(), true);

            // Tidy up
            delete p_cell_population;
        }
    }

    void TestRemoveDeadCells()
    {
        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        TS_ASSERT_EQUALS(cell_population.GetNumAllCells(), 5);
        auto cell = cell_population.GetCellUsingLocationIndex(2);
        cell->Kill();
        cell_population.RemoveDeadCells();
        cell_population.rGetMesh().ReMesh();
        TS_ASSERT_EQUALS(cell_population.GetNumAllCells(), 4);
    }

    void TestGetCellDataItemAtPdeNode()
    {
        {
            // Create an immersed boundary cell population object
            ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
            ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

            for (auto& p_cell : cell_population.rGetCells())
            {
                p_cell->GetCellData()->SetItem("cell data", 0.2);
            }

            std::string str = "cell data";
            TS_ASSERT_DELTA(cell_population.GetCellDataItemAtPdeNode(0, str, true, 0.1), 0.1, 1e-9);
            TS_ASSERT_DELTA(cell_population.GetCellDataItemAtPdeNode(p_mesh->GetNumNodes() + 1, str, true, 0.1), 0.2, 1e-9);
        }

        {
            // Create a small mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));
            elems.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes));
            elems.push_back(new ImmersedBoundaryElement<2, 2>(2, nodes));

            std::vector<ImmersedBoundaryElement<1, 2>*> lams;
            lams.push_back(new ImmersedBoundaryElement<1, 2>(0, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(1, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(2, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, lams);
            auto p_mesh = &mesh;
            p_mesh->SetNumGridPtsXAndY(32);

            // Create a minimal cell population
            std::vector<CellPtr> cells;
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
            for (auto& p_cell : cell_population.rGetCells())
            {
                p_cell->GetCellData()->SetItem("cell data", 0.2);
            }

            std::string str = "cell data";
            TS_ASSERT_DELTA(cell_population.GetCellDataItemAtPdeNode(0, str, false, 0.1), 0.2, 0.0001);
        }
    }

    void TestPdeNonApoptoticCell()
    {
        /*
         * Test overridden IsPdeNodeAssociatedWithNonApoptoticCell() method,
         * which returns whether a node, specified by its index in a tetrahedral
         * mesh for use with a PDE modifier, is associated with a non-apoptotic
         * cell.
         *
         * \todo this seems an odd test, since there is not actually a PDE
         * modifier present.
         */

        // Create a small mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
        nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
        nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

        std::vector<ImmersedBoundaryElement<2, 2>*> elems;
        elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));
        elems.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes));
        elems.push_back(new ImmersedBoundaryElement<2, 2>(2, nodes));

        std::vector<ImmersedBoundaryElement<1, 2>*> lams;
        lams.push_back(new ImmersedBoundaryElement<1, 2>(0, nodes));
        lams.push_back(new ImmersedBoundaryElement<1, 2>(1, nodes));
        lams.push_back(new ImmersedBoundaryElement<1, 2>(2, nodes));

        ImmersedBoundaryMesh<2,2> mesh(nodes, elems, lams);
        auto p_mesh = &mesh;
        p_mesh->SetNumGridPtsXAndY(32);

        // Create a minimal cell population
        std::vector<CellPtr> cells;
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        for (auto& p_cell : cell_population.rGetCells())
        {
            p_cell->GetCellData()->SetItem("cell data", 0.2);
        }
        MAKE_PTR(ApoptoticCellProperty, cellProperty);
        cell_population.rGetCells().front()->AddCellProperty(cellProperty);

        std::string str = "cell data";

        TS_ASSERT_EQUALS(cell_population.IsPdeNodeAssociatedWithNonApoptoticCell(0), false);
        TS_ASSERT_EQUALS(cell_population.IsPdeNodeAssociatedWithNonApoptoticCell(4), false);
    }

    void TestGetNeighbouringNodeIndices()
    {
        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        std::set<unsigned> expected_neighbours {};
        TS_ASSERT_EQUALS(cell_population.GetNeighbouringNodeIndices(1), expected_neighbours);
    }
};

#endif /*TESTIMMERSEDBOUNDARYCELLPOPULATION_HPP_*/
