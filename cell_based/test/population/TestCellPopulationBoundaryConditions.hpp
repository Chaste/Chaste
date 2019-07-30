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

#ifndef TESTCELLPOPULATIONBOUNDARYCONDITIONS_HPP_
#define TESTCELLPOPULATIONBOUNDARYCONDITIONS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellsGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "SphereGeometryBoundaryCondition.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "TrianglesMeshReader.hpp"
#include "WildTypeCellMutationState.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "ArchiveOpener.hpp"
#include "ArchiveLocationInfo.hpp"
#include "FileComparison.hpp"
#include "SmartPointers.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * This class contains tests for methods on classes inheriting from AbstractCellPopulationBoundaryCondition.
 */
class TestCellPopulationBoundaryConditions : public AbstractCellBasedTestSuite
{
public:

    void TestPlaneBoundaryConditionWithNodeBasedCellPopulation()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

        // Create mesh
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up cell population boundary condition
        c_vector<double,2> point = zero_vector<double>(2);
        point(0) = 2.0;
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        PlaneBoundaryCondition<2> boundary_condition(&cell_population, point, normal);

        TS_ASSERT_EQUALS(boundary_condition.GetIdentifier(), "PlaneBoundaryCondition-2-2");

        // Impose boundary condition
        std::map<Node<2>*, c_vector<double,2> > old_locations;
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            old_locations[p_node] = p_node->rGetLocation();
        }

        boundary_condition.ImposeBoundaryCondition(old_locations);

        // Test that all nodes satisfy the boundary condition
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location;
            location = p_node->rGetLocation();
            if (old_locations[p_node][0] < 2.0)
            {
                TS_ASSERT_DELTA(2.0, location[0], 1e-6);
                TS_ASSERT_DELTA(location[1], old_locations[p_node][1], 1e-6);
            }
            else
            {
                TS_ASSERT_DELTA(location[0], old_locations[p_node][0], 1e-6);
                TS_ASSERT_DELTA(location[1], old_locations[p_node][1], 1e-6);
            }
        }

        // Test VerifyBoundaryCondition() method
        TS_ASSERT_EQUALS(boundary_condition.VerifyBoundaryCondition(), true);
    }

    void TestPlaneBoundaryConditionWithNodeBasedCellPopulationAndJiggledNodes()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

        // Create mesh
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up cell population boundary condition
        c_vector<double,2> point = zero_vector<double>(2);
        point(0) = 2.0;
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        PlaneBoundaryCondition<2> boundary_condition(&cell_population, point, normal);

        // Turn on the jiggle nodes
        TS_ASSERT(!(boundary_condition.GetUseJiggledNodesOnPlane()));
        boundary_condition.SetUseJiggledNodesOnPlane(true);
        TS_ASSERT(boundary_condition.GetUseJiggledNodesOnPlane());

        // Impose boundary condition
        std::map<Node<2>*, c_vector<double,2> > old_locations;
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            old_locations[p_node] = p_node->rGetLocation();
        }

        boundary_condition.ImposeBoundaryCondition(old_locations);

        // Test that all nodes satisfy the boundary condition. Note these are random so depend on the seed
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetNode(0)->rGetLocation()[0],2.00005,1e-5);
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetNode(0)->rGetLocation()[1],0.0,1e-5);

        TS_ASSERT_DELTA(cell_population.rGetMesh().GetNode(1)->rGetLocation()[0],2.00006,1e-5);
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetNode(1)->rGetLocation()[1],0.0,1e-5);

        TS_ASSERT_DELTA(cell_population.rGetMesh().GetNode(2)->rGetLocation()[0],2.00007,1e-5);
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetNode(2)->rGetLocation()[1],sqrt(3.0)/2.0,1e-5);

        TS_ASSERT_DELTA(cell_population.rGetMesh().GetNode(3)->rGetLocation()[0],2.00008,1e-5);
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetNode(3)->rGetLocation()[1],sqrt(3.0)/2.0,1e-5);

        // Test VerifyBoundaryCondition() method
        TS_ASSERT_EQUALS(boundary_condition.VerifyBoundaryCondition(), true);
    }

    void TestPlaneBoundaryConditionWithVertexBasedCellPopulation()
    {
        // Create mesh
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // There are four cells and their x-extents are (0, 1), (0.5, 1.5), (1, 2), (1.5, 2.5)
        // Imposing x>0.75 will squash two of these cells, but not so much that they degenerate to zero area
        // Set up cell population boundary condition x>0.75
        double x_boundary = 0.75;
        c_vector<double,2> point = zero_vector<double>(2);
        point(0) = x_boundary;
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        PlaneBoundaryCondition<2> boundary_condition(&cell_population, point, normal);

        TS_ASSERT_EQUALS(boundary_condition.GetIdentifier(), "PlaneBoundaryCondition-2-2");

        // Impose boundary condition
        std::map<Node<2>*, c_vector<double,2> > old_locations;

        for (MutableVertexMesh<2,2>::NodeIterator node_iter = cell_population.rGetMesh().GetNodeIteratorBegin();
             node_iter != cell_population.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            old_locations[&(*node_iter)] = node_iter->rGetLocation();
        }

        boundary_condition.ImposeBoundaryCondition(old_locations);

        // Test that all nodes satisfy the boundary condition
        for (MutableVertexMesh<2,2>::NodeIterator node_iter = cell_population.rGetMesh().GetNodeIteratorBegin();
             node_iter != cell_population.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            c_vector<double, 2> location;
            location = node_iter->rGetLocation();
            if (old_locations[&(*node_iter)][0] < x_boundary)
            {
                TS_ASSERT_DELTA(x_boundary, location[0], 1e-6);
                TS_ASSERT_DELTA(location[1], old_locations[&(*node_iter)][1], 1e-6);
            }
            else
            {
                TS_ASSERT_DELTA(location[0], old_locations[&(*node_iter)][0], 1e-6);
                TS_ASSERT_DELTA(location[1], old_locations[&(*node_iter)][1], 1e-6);
            }
        }

        // Test VerifyBoundaryCondition() method
        TS_ASSERT_EQUALS(boundary_condition.VerifyBoundaryCondition(), true);

        // Now test with the jiggled nodes for coverage
        x_boundary = 0.8;
        point(0) = x_boundary;
        PlaneBoundaryCondition<2> boundary_condition_2(&cell_population, point, normal);

        TS_ASSERT(!(boundary_condition_2.GetUseJiggledNodesOnPlane()));
        boundary_condition_2.SetUseJiggledNodesOnPlane(true);
        TS_ASSERT(boundary_condition_2.GetUseJiggledNodesOnPlane());

        boundary_condition_2.ImposeBoundaryCondition(old_locations);

        // Test that all nodes satisfy the jiggled boundary condition
        for (MutableVertexMesh<2,2>::NodeIterator node_iter = cell_population.rGetMesh().GetNodeIteratorBegin();
             node_iter != cell_population.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            c_vector<double, 2> location;
            location = node_iter->rGetLocation();
            if (old_locations[&(*node_iter)][0] < x_boundary)
            {
                TS_ASSERT_LESS_THAN(x_boundary, location[0]);// note strict inequality
                TS_ASSERT_DELTA(location[1], old_locations[&(*node_iter)][1], 1e-6);
            }
            else
            {
                TS_ASSERT_DELTA(location[0], old_locations[&(*node_iter)][0], 1e-6);
                TS_ASSERT_DELTA(location[1], old_locations[&(*node_iter)][1], 1e-6);
            }
        }
    }

    void TestPlaneBoundaryConditionExceptions()
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<2> potts_cell_population(*p_mesh, cells);

        // Attempt to set up cell population boundary condition
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = 1.0;

        PlaneBoundaryCondition<2> plane_boundary_condition(&potts_cell_population, point, normal);
        std::map<Node<2>*, c_vector<double, 2> > old_locations;
        TS_ASSERT_THROWS_THIS(plane_boundary_condition.ImposeBoundaryCondition(old_locations),
            "PlaneBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");

        // Test the correct exception is thrown in 1D
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0, true,  0.0));
        nodes.push_back(new Node<1>(1, false, 1.0));
        NodesOnlyMesh<1> nodes_only_mesh;
        nodes_only_mesh.ConstructNodesWithoutMesh(nodes, 1.5);
        std::vector<CellPtr> node_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator_1d;
        cells_generator_1d.GenerateBasicRandom(node_based_cells, nodes_only_mesh.GetNumNodes(), p_diff_type);
        NodeBasedCellPopulation<1> node_based_cell_population(nodes_only_mesh, node_based_cells);

        c_vector<double,1> point_1d = zero_vector<double>(1);
        c_vector<double,1> normal_1d = zero_vector<double>(1);
        normal_1d(0) = 1.0;
        PlaneBoundaryCondition<1> plane_bc_1d(&node_based_cell_population, point_1d, normal_1d);
        TS_ASSERT_THROWS_THIS(plane_bc_1d.VerifyBoundaryCondition(),
            "PlaneBoundaryCondition is not implemented in 1D");

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestSphereGeometryBoundaryCondition()
    {
        // We first test that the correct exception is thrown in 1D
        TrianglesMeshReader<1,1> mesh_reader_1d("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> generating_mesh_1d;
        generating_mesh_1d.ConstructFromMeshReader(mesh_reader_1d);

        NodesOnlyMesh<1> mesh_1d;
        mesh_1d.ConstructNodesWithoutMesh(generating_mesh_1d, 1.5);

        std::vector<CellPtr> cells_1d;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 1> cells_generator_1d;
        cells_generator_1d.GenerateBasic(cells_1d, mesh_1d.GetNumNodes());

        NodeBasedCellPopulation<1> population_1d(mesh_1d, cells_1d);

        c_vector<double,1> centre_1d = zero_vector<double>(1);

        TS_ASSERT_THROWS_THIS(SphereGeometryBoundaryCondition<1> bc_1d(&population_1d, centre_1d, 1.0),
            "This boundary condition is not implemented in 1D.");

        // Next we test that the correct exception is thrown if not using a NodeBasedCellPopulation
        TrianglesMeshReader<2,2> mesh_reader_2d("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh_2d;
        mesh_2d.ConstructFromMeshReader(mesh_reader_2d);

        std::vector<CellPtr> cells_2d;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator_2d;
        cells_generator_2d.GenerateBasic(cells_2d, mesh_2d.GetNumNodes());

        MeshBasedCellPopulation<2> population_2d(mesh_2d, cells_2d);

        c_vector<double,2> centre_2d = zero_vector<double>(2);

        TS_ASSERT_THROWS_THIS(SphereGeometryBoundaryCondition<2> bc_2d(&population_2d, centre_2d, 1.0),
            "A NodeBasedCellPopulation must be used with this boundary condition object.");

        // We now test the methods of this class
        TrianglesMeshReader<3,3> mesh_reader_3d("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> generating_mesh_3d;
        generating_mesh_3d.ConstructFromMeshReader(mesh_reader_3d);

        NodesOnlyMesh<3> mesh_3d;
        mesh_3d.ConstructNodesWithoutMesh(generating_mesh_3d, 1.5);

        std::vector<CellPtr> cells_3d;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator_3d;
        cells_generator_3d.GenerateBasic(cells_3d, mesh_3d.GetNumNodes());

        NodeBasedCellPopulation<3> population_3d(mesh_3d, cells_3d);

        c_vector<double,3> centre_3d = zero_vector<double>(3);
        centre_3d(0) = 0.1; //The centre is at (0.1, 0.0, 0.0).  This is because there is no node at the centre

        SphereGeometryBoundaryCondition<3> bc_3d(&population_3d, centre_3d, 0.4, 1e-4);

        // Test that member variables were initialised correctly
        TS_ASSERT_DELTA(bc_3d.rGetCentreOfSphere()[0], centre_3d(0), 1e-4);
        TS_ASSERT_DELTA(bc_3d.rGetCentreOfSphere()[1], centre_3d(1), 1e-4);
        TS_ASSERT_DELTA(bc_3d.GetRadiusOfSphere(), 0.4, 1e-4);

        // If there are any cells on this process, this should return false
        if (mesh_3d.GetNumNodes() > 0)
        {
            TS_ASSERT_EQUALS(bc_3d.VerifyBoundaryCondition(), false);
        }

        // Store the location of each node prior to imposing the boundary condition
        std::map<Node<3>*, c_vector<double,3> > old_locations;

        for (std::list<CellPtr>::iterator cell_iter = population_3d.rGetCells().begin();
             cell_iter != population_3d.rGetCells().end();
             ++cell_iter)
        {
            c_vector<double,3> location = population_3d.GetLocationOfCellCentre(*cell_iter);
            unsigned index = population_3d.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = mesh_3d.GetNode(index);
            old_locations[p_node] = location;
        }

        bc_3d.ImposeBoundaryCondition(old_locations);

        for (std::list<CellPtr>::iterator cell_iter = population_3d.rGetCells().begin();
             cell_iter != population_3d.rGetCells().end();
             ++cell_iter)
        {
            c_vector<double,3> new_direction = population_3d.GetLocationOfCellCentre(*cell_iter)- centre_3d;

            //Check it's at the right distance from the centre
            double new_distance = norm_2(new_direction);
            TS_ASSERT_DELTA(new_distance, 0.4, 1e-10);

            //Check that the direction is correct i.e that the new direction is parallel to the old direction
            new_direction /= new_distance;
            unsigned index = population_3d.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = mesh_3d.GetNode(index);
            c_vector<double,3> old_direction =  old_locations[p_node] - centre_3d;
            old_direction /= norm_2(old_direction);
            TS_ASSERT_DELTA(norm_inf(new_direction-old_direction), 0.0, 1e-10);
        }

        // Test that the boundary condition was imposed correctly
        TS_ASSERT_EQUALS(bc_3d.VerifyBoundaryCondition(), true);
    }

    void TestArchivingOfPlaneBoundaryCondition()
    {
        EXIT_IF_PARALLEL;    // We cannot archive parallel cell based simulations yet.
        // Set up singleton classes

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> population(mesh, cells);

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "PlaneBoundaryCondition.arch";
        ArchiveLocationInfo::SetMeshFilename("PlaneBoundaryCondition");

        {
            // Create an output archive
            PlaneBoundaryCondition<2> boundary_condition(&population, zero_vector<double>(2), unit_vector<double>(2,1));

            TS_ASSERT_DELTA(boundary_condition.rGetPointOnPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(boundary_condition.rGetPointOnPlane()[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(boundary_condition.rGetNormalToPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(boundary_condition.rGetNormalToPlane()[1], 1.0, 1e-6);

            // Create an output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Serialize via pointer
            AbstractCellPopulationBoundaryCondition<2>* const p_boundary_condition = &boundary_condition;
            (*p_arch) << p_boundary_condition;
        }

        {
            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCellPopulationBoundaryCondition<2,2>* p_boundary_condition;

            // Restore from the archive
            (*p_arch) >> p_boundary_condition;

            // Test we have restored the plane geometry correctly
            TS_ASSERT_DELTA(static_cast<PlaneBoundaryCondition<2>*>(p_boundary_condition)->rGetPointOnPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(static_cast<PlaneBoundaryCondition<2>*>(p_boundary_condition)->rGetPointOnPlane()[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(static_cast<PlaneBoundaryCondition<2>*>(p_boundary_condition)->rGetNormalToPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(static_cast<PlaneBoundaryCondition<2>*>(p_boundary_condition)->rGetNormalToPlane()[1], 1.0, 1e-6);

            // Tidy up
            delete p_boundary_condition->mpCellPopulation;
            delete p_boundary_condition;
       }
    }

    void TestArchivingOfSphereGeometryBoundaryCondition()
    {
        EXIT_IF_PARALLEL;    // We cannot archive parallel cell based simulations yet.

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> population(mesh, cells);

        c_vector<double,2> centre = zero_vector<double>(2);
        centre(0) = 0.5;
        centre(1) = 0.7;

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "SphereGeometryBoundaryCondition.arch";
        ArchiveLocationInfo::SetMeshFilename("SphereGeometryBoundaryCondition");

        {
            SphereGeometryBoundaryCondition<2> bc(&population, centre, 0.56, 1e-3);

            // Create an output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Serialize via pointer
            AbstractCellPopulationBoundaryCondition<2>* const p_bc = &bc;
            (*p_arch) << p_bc;
        }

        {
            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCellPopulationBoundaryCondition<2>* p_bc;

            // Restore from the archive
            (*p_arch) >> p_bc;

            // Test we have restored the object correctly
            TS_ASSERT_DELTA(static_cast<SphereGeometryBoundaryCondition<2>*>(p_bc)->rGetCentreOfSphere()[0], 0.5, 1e-6);
            TS_ASSERT_DELTA(static_cast<SphereGeometryBoundaryCondition<2>*>(p_bc)->rGetCentreOfSphere()[1], 0.7, 1e-6);
            TS_ASSERT_DELTA(static_cast<SphereGeometryBoundaryCondition<2>*>(p_bc)->GetRadiusOfSphere(), 0.56, 1e-6);

            // Tidy up
            delete p_bc->mpCellPopulation;
            delete p_bc;
        }
    }

    void TestCellBoundaryConditionsOutputParameters()
    {
        EXIT_IF_PARALLEL;
        std::string output_directory = "TestCellBoundaryConditionsOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> population(mesh, cells);

        // Test with PlaneBoundaryCondition
        PlaneBoundaryCondition<2> plane_boundary_condition(&population, zero_vector<double>(2), unit_vector<double>(2,1));
        TS_ASSERT_EQUALS(plane_boundary_condition.GetIdentifier(), "PlaneBoundaryCondition-2-2");

        out_stream plane_boundary_condition_parameter_file = output_file_handler.OpenOutputFile("plane_results.parameters");
        plane_boundary_condition.OutputCellPopulationBoundaryConditionParameters(plane_boundary_condition_parameter_file);
        plane_boundary_condition_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code.
            FileFinder generated = output_file_handler.FindFile("plane_results.parameters");
            FileFinder reference("cell_based/test/data/TestCellBoundaryConditionsOutputParameters/plane_results.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with SphereGeometryBoundaryCondition
        c_vector<double,2> centre = zero_vector<double>(2);

        SphereGeometryBoundaryCondition<2> sphere_boundary_condition(&population, centre, 0.56, 1e-3);
        TS_ASSERT_EQUALS(sphere_boundary_condition.GetIdentifier(), "SphereGeometryBoundaryCondition-2");

        out_stream sphere_boundary_condition_parameter_file = output_file_handler.OpenOutputFile("sphere_results.parameters");
        sphere_boundary_condition.OutputCellPopulationBoundaryConditionParameters(sphere_boundary_condition_parameter_file);
        sphere_boundary_condition_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code.
            FileFinder generated = output_file_handler.FindFile("sphere_results.parameters");
            FileFinder reference("cell_based/test/data/TestCellBoundaryConditionsOutputParameters/sphere_results.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test OutputCellPopulationBoundaryConditionInfo() method
        out_stream plane_boundary_condition_info_file = output_file_handler.OpenOutputFile("plane_results.info");
        plane_boundary_condition.OutputCellPopulationBoundaryConditionInfo(plane_boundary_condition_info_file);
        plane_boundary_condition_info_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code.
            FileFinder generated = output_file_handler.FindFile("plane_results.info");
            FileFinder reference("cell_based/test/data/TestCellBoundaryConditionsOutputParameters/plane_results.info",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif /*TESTCELLPOPULATIONBOUNDARYCONDITIONS_HPP_*/
