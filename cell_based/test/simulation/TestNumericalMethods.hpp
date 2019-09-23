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

#ifndef TESTNUMERICALMETHODS_HPP_
#define TESTNUMERICALMETHODS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NodeBasedCellPopulationWithParticles.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "CellAncestor.hpp"
#include "CellId.hpp"
#include "CellPropertyRegistry.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "PopulationTestingForce.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "Warnings.hpp"


#include "PetscSetupAndFinalize.hpp"

class TestNumericalMethods : public AbstractCellBasedTestSuite
{
public:

    void TestMethodsAndExceptions()
    {

        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

        {
            unsigned cells_across = 7;
            unsigned cells_up = 5;

            HoneycombMeshGenerator generator(cells_across, cells_up, 0);
            MutableMesh<2,2>* p_mesh = generator.GetMesh();

            // Create cells
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

            // Create cell population
            MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Create a force collection
            std::vector<boost::shared_ptr<AbstractForce<2,2> > > force_collection;
            MAKE_PTR(PopulationTestingForce<2>, p_test_force);
            force_collection.push_back(p_test_force);

            // Create Numerical method
            ForwardEulerNumericalMethod<2> numerical_method;

            numerical_method.SetCellPopulation(&cell_population);
            numerical_method.SetForceCollection(&force_collection);
            numerical_method.SetUseAdaptiveTimestep(true);

            std::vector<c_vector<double, 2> > saved_locations;

            saved_locations = numerical_method.SaveCurrentLocations();

            TS_ASSERT_EQUALS(saved_locations.size(), cell_population.GetNumRealCells());

        }

        // This tests the exceptions for Node based with Buske Update
        {
            // Create a simple mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
            MutableMesh<2,2> generating_mesh;
            generating_mesh.ConstructFromMeshReader(mesh_reader);

            // Convert this to a NodesOnlyMesh
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(generating_mesh, 1.2);

            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

            // Create a cell population, with no ghost nodes at the moment
            NodeBasedCellPopulationWithBuskeUpdate<2> cell_population(mesh, cells);


            // Create Numerical method
            ForwardEulerNumericalMethod<2> numerical_method;

            numerical_method.SetCellPopulation(&cell_population);

            TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Non-Euler steppers are not yet implemented for NodeBasedCellPopulationWithBuskeUpdate");
            Warnings::QuietDestroy();

        }
    }

    void TestUpdateAllNodePositionsWithMeshBased()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population, with no ghost nodes at the moment
        MeshBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetDampingConstantNormal(1.1);

        // Create a force collection
        std::vector<boost::shared_ptr<AbstractForce<2,2> > > force_collection;
        MAKE_PTR(PopulationTestingForce<2>, p_test_force);
        force_collection.push_back(p_test_force);

        // Create numerical method for testing
        MAKE_PTR(ForwardEulerNumericalMethod<2>, p_fe_method);

        double dt = 0.01;

        p_fe_method->SetCellPopulation(&cell_population);
        p_fe_method->SetForceCollection(&force_collection);

        // Save starting positions
        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
        for (unsigned j=0; j<cell_population.GetNumNodes(); j++)
        {
            old_posns[j][0] = cell_population.GetNode(j)->rGetLocation()[0];
            old_posns[j][1] = cell_population.GetNode(j)->rGetLocation()[1];
        }

        // Update positions and check the answer
        p_fe_method->UpdateAllNodePositions(dt);

        for (unsigned j=0; j<cell_population.GetNumNodes(); j++)
        {
            c_vector<double, 2> actualLocation = cell_population.GetNode(j)->rGetLocation();

            double damping =  cell_population.GetDampingConstant(j);
            c_vector<double, 2> expectedLocation;
            expectedLocation = p_test_force->GetExpectedOneStepLocationFE(j, damping, old_posns[j], dt);

            TS_ASSERT_DELTA(norm_2(actualLocation - expectedLocation), 0, 1e-6);
        }
    }

    void TestUpdateAllNodePositionsWithMeshBasedWithGhosts()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

        HoneycombMeshGenerator generator(3, 3, 1);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.SetDampingConstantNormal(1.1);

        // Create a force collection
        std::vector<boost::shared_ptr<AbstractForce<2,2> > > force_collection;
        MAKE_PTR(PopulationTestingForce<2>, p_test_force);
        force_collection.push_back(p_test_force);

        // Create numerical methods for testing
        MAKE_PTR(ForwardEulerNumericalMethod<2>, p_fe_method);

        double dt = 0.01;

        p_fe_method->SetCellPopulation(&cell_population);
        p_fe_method->SetForceCollection(&force_collection);

        // Save starting positions
        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
        for (unsigned j=0; j<cell_population.GetNumNodes(); j++)
        {
            old_posns[j][0] = cell_population.GetNode(j)->rGetLocation()[0];
            old_posns[j][1] = cell_population.GetNode(j)->rGetLocation()[1];
        }

        // Update positions
        p_fe_method->UpdateAllNodePositions(dt);

        //Check the answer (for cell associated nodes only)
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            int j = cell_population.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, 2> actualLocation = cell_population.GetNode(j)->rGetLocation();

            double damping =  cell_population.GetDampingConstant(j);
            c_vector<double, 2> expectedLocation;
            expectedLocation = p_test_force->GetExpectedOneStepLocationFE(j, damping, old_posns[j], dt);
            TS_ASSERT_DELTA(norm_2(actualLocation - expectedLocation), 0, 1e-9);
        }
    }

    void TestUpdateAllNodePositionsWithNodeBased()
    {
        EXIT_IF_PARALLEL;    // This test doesn't work in parallel.

        HoneycombMeshGenerator generator(3, 3, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        MAKE_PTR(NodesOnlyMesh<2>, p_mesh);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.1);

        // Create a force collection
        std::vector<boost::shared_ptr<AbstractForce<2,2> > > force_collection;
        MAKE_PTR(PopulationTestingForce<2>, p_test_force);
        force_collection.push_back(p_test_force);

        // Create numerical method for testing
        MAKE_PTR(ForwardEulerNumericalMethod<2>, p_fe_method);

        double dt = 0.01;

        p_fe_method->SetCellPopulation(&cell_population);
        p_fe_method->SetForceCollection(&force_collection);

        // Save starting positions
        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
        for (unsigned j=0; j<cell_population.GetNumNodes(); j++)
        {
            old_posns[j][0] = cell_population.GetNode(j)->rGetLocation()[0];
            old_posns[j][1] = cell_population.GetNode(j)->rGetLocation()[1];
        }

        // Update positions and check the answer
        p_fe_method->UpdateAllNodePositions(dt);

        for (unsigned j=0; j<cell_population.GetNumNodes(); j++)
        {
            c_vector<double, 2> actualLocation = cell_population.GetNode(j)->rGetLocation();

            double damping =  cell_population.GetDampingConstant(j);
            c_vector<double, 2> expectedLocation;
            expectedLocation = p_test_force->GetExpectedOneStepLocationFE(j, damping, old_posns[j], dt);

            TS_ASSERT_DELTA(norm_2(actualLocation - expectedLocation), 0, 1e-12);
        }
    }

    void TestUpdateAllNodePositionsWithNodeBasedWithBuskeUpdate()
    {
        EXIT_IF_PARALLEL;    // This test doesn't work in parallel.

        HoneycombMeshGenerator generator(3, 3, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        MAKE_PTR(NodesOnlyMesh<2>, p_mesh);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        NodeBasedCellPopulationWithBuskeUpdate<2> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.1);

        // Create a force collection
        std::vector<boost::shared_ptr<AbstractForce<2,2> > > force_collection;
        MAKE_PTR(PopulationTestingForce<2>, p_test_force);
        force_collection.push_back(p_test_force);

        // Create numerical method for testing
        MAKE_PTR(ForwardEulerNumericalMethod<2>, p_fe_method);

        double dt = 0.01;

        p_fe_method->SetCellPopulation(&cell_population);
        p_fe_method->SetForceCollection(&force_collection);

        // Save starting positions
        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
        for (unsigned j=0; j<cell_population.GetNumNodes(); j++)
        {
            old_posns[j][0] = cell_population.GetNode(j)->rGetLocation()[0];
            old_posns[j][1] = cell_population.GetNode(j)->rGetLocation()[1];
        }

        // Update positions and check the answer
        // Currently this throws an error as not set up correctly as it is in a simulation #2087
        TS_ASSERT_THROWS_THIS(p_fe_method->UpdateAllNodePositions(dt),"You must provide a rowPreallocation argument for a large sparse system");

        // for (unsigned j=0; j<cell_population.GetNumNodes(); j++)
        // {
        //     c_vector<double, 2> actualLocation = cell_population.GetNode(j)->rGetLocation();

        //     double damping =  cell_population.GetDampingConstant(j);
        //     c_vector<double, 2> expectedLocation;
        //     expectedLocation = p_test_force->GetExpectedOneStepLocationFE(j, damping, old_posns[j], dt);

        //     TS_ASSERT_DELTA(norm_2(actualLocation - expectedLocation), 0, 1e-12);
        // }
    }

    void TestUpdateAllNodePositionsWithNodeBasedWithParticles()
    {
        EXIT_IF_PARALLEL;    // This test doesn't work in parallel.

        HoneycombMeshGenerator generator(3, 3, 1);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        MAKE_PTR(NodesOnlyMesh<2>, p_mesh);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());


        NodeBasedCellPopulationWithParticles<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.SetDampingConstantNormal(1.1);

        // Create a force collection
        std::vector<boost::shared_ptr<AbstractForce<2,2> > > force_collection;
        MAKE_PTR(PopulationTestingForce<2>, p_test_force);
        force_collection.push_back(p_test_force);

        // Create numerical method for testing
        MAKE_PTR(ForwardEulerNumericalMethod<2>, p_fe_method);

        double dt = 0.01;

        p_fe_method->SetCellPopulation(&cell_population);
        p_fe_method->SetForceCollection(&force_collection);

        // Save starting positions
        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
        for (unsigned j=0; j<cell_population.GetNumNodes(); j++)
        {
            old_posns[j][0] = cell_population.GetNode(j)->rGetLocation()[0];
            old_posns[j][1] = cell_population.GetNode(j)->rGetLocation()[1];
        }

        // Update positions and check the answer
        p_fe_method->UpdateAllNodePositions(dt);

        for (unsigned j=0; j<cell_population.GetNumNodes(); j++)
        {
            c_vector<double, 2> actualLocation = cell_population.GetNode(j)->rGetLocation();

            double damping =  cell_population.GetDampingConstant(j);
            c_vector<double, 2> expectedLocation;
            expectedLocation = p_test_force->GetExpectedOneStepLocationFE(j, damping, old_posns[j], dt);

            TS_ASSERT_DELTA(norm_2(actualLocation - expectedLocation), 0, 1e-12);
        }
    }

    void TestUpdateAllNodePositionsWithVertexBased()
    {
        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(5, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Impose a larger cell rearrangement threshold so that motion is uninhibited (see #1376)
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create a cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.1);

        // Create a force collection
        std::vector<boost::shared_ptr<AbstractForce<2,2> > > force_collection;
        MAKE_PTR(PopulationTestingForce<2>, p_test_force);
        force_collection.push_back(p_test_force);

        // Create numerical methods for testing
        MAKE_PTR(ForwardEulerNumericalMethod<2>, p_fe_method);

        double dt = 0.01;

        p_fe_method->SetCellPopulation(&cell_population);
        p_fe_method->SetForceCollection(&force_collection);

        // Save starting positions
        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
        for (unsigned j=0; j<cell_population.GetNumNodes(); j++)
        {
            old_posns[j][0] = cell_population.GetNode(j)->rGetLocation()[0];
            old_posns[j][1] = cell_population.GetNode(j)->rGetLocation()[1];
        }

        // Update positions and check the answer
        p_fe_method->UpdateAllNodePositions(dt);

        for (unsigned j=0; j<cell_population.GetNumNodes(); j++)
        {
            c_vector<double, 2> actualLocation = cell_population.GetNode(j)->rGetLocation();

            double damping =  cell_population.GetDampingConstant(j);
            c_vector<double, 2> expectedLocation;
            expectedLocation = p_test_force->GetExpectedOneStepLocationFE(j, damping, old_posns[j], dt);

            TS_ASSERT_DELTA(norm_2(actualLocation - expectedLocation), 0, 1e-12);
        }
    }

    void TestSettingAndGettingFlags()
    {
        // Create numerical methods for testing
        MAKE_PTR(ForwardEulerNumericalMethod<2>, p_fe_method);

        // mUseUpdateNodeLocation should default to false
        TS_ASSERT(!(p_fe_method->GetUseUpdateNodeLocation()));

        // Set mUseUpdateNodeLocation to true and check
        p_fe_method->SetUseUpdateNodeLocation(true);
        TS_ASSERT(p_fe_method->GetUseUpdateNodeLocation());
    }
};

#endif /*TESTMESHBASEDCELLPOPULATION_HPP_*/
