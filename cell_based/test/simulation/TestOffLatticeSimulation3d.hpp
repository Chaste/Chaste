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

#ifndef TESTOFFLATTICESIMULATION3D_HPP_
#define TESTOFFLATTICESIMULATION3D_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "TrianglesMeshReader.hpp"
#include "OffLatticeSimulation.hpp"
#include "TrianglesMeshWriter.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SmartPointers.hpp"
#include "CellVolumesWriter.hpp"

// Cell population writers
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "VoronoiDataWriter.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestOffLatticeSimulation3d : public AbstractCellBasedWithTimingsTestSuite
{
private:
    double mLocationGhosts;
    double mLocationWithoutGhosts;

    MutableMesh<3,3>* Make3dMesh(unsigned width=3, unsigned height=3, unsigned depth=3)
    {
        MutableMesh<3,3>* p_mesh = new MutableMesh<3,3>;
        p_mesh->ConstructCuboid(width, height, depth);
        TrianglesMeshWriter<3,3> mesh_writer("", "3dSpringMesh");
        mesh_writer.WriteFilesUsingMesh(*p_mesh);

        return p_mesh;
    }

public:

    void TestDoCellBirth()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            dynamic_cast<FixedG1GenerationalCellCycleModel*>(cells[i]->GetCellCycleModel())->SetGeneration(0);
            cells[i]->SetBirthTime(0.0);
        }
        cells[50]->SetBirthTime(-50.0);

        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);

        unsigned num_births = simulator.DoCellBirth();

        TS_ASSERT_EQUALS(num_births, 1u);
    }

    void TestBirthOccursDuringSolve()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");

        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            dynamic_cast<FixedG1GenerationalCellCycleModel*>(cells[i]->GetCellCycleModel())->SetGeneration(0);
            cells[i]->SetBirthTime(0.0);
        }
        cells[mesh.GetNumNodes()-1]->SetBirthTime(-50.0);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);

        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);

        TrianglesMeshWriter<3,3> mesh_writer1("Test3DCellBirth", "StartMesh");
        mesh_writer1.WriteFilesUsingMesh(mesh);

        simulator.SetOutputDirectory("Test3DCellBirth");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Start with a single 3D tetrahedron, add one node get three tetrahedral elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);

        TrianglesMeshWriter<3,3> mesh_writer2("Test3DCellBirth", "EndMesh", false);
        mesh_writer2.WriteFilesUsingMesh(mesh);
    }

    void TestSolveMethodSpheroidSimulation3D()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<3,3> mesh_writer("TestSolveMethodSpheroidSimulation3DMesh", "StartMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            dynamic_cast<FixedG1GenerationalCellCycleModel*>(cells[i]->GetCellCycleModel())->SetGeneration(0);
        }

        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestSolveMethodSpheroidSimulation3D");

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Test SetSamplingTimestepMultiple method
        TS_ASSERT_EQUALS(simulator.mSamplingTimestepMultiple, 1u);
        simulator.SetSamplingTimestepMultiple(2);
        TS_ASSERT_EQUALS(simulator.mSamplingTimestepMultiple, 2u);

        // Uncommenting this line calls an error in accessing nodes in the vertex elements #
        //cell_population.AddPopulationWriter<VoronoiDataWriter>();

        simulator.SetEndTime(0.1);
        simulator.Solve();

        TrianglesMeshWriter<3,3> mesh_writer2("TestSolveMethodSpheroidSimulation3DMesh", "EndMesh", false);
        mesh_writer2.WriteFilesUsingMesh(mesh);
    }

    void TestGhostNodesSpheroidSimulation3DandSave()
    {
        unsigned width = 3;
        unsigned height = 3;
        unsigned depth = 3;

        MutableMesh<3,3>* p_mesh = Make3dMesh(width, height, depth);
        TrianglesMeshWriter<3,3> mesh_writer("TestGhostNodesSpheroidSimulation3D", "StartMesh");
        mesh_writer.WriteFilesUsingMesh(*p_mesh);

        c_vector<double, 3> spheroid_centre;
        spheroid_centre[0] = 0.5*((double) width);
        spheroid_centre[1] = 0.5*((double) height);
        spheroid_centre[2] = 0.5*((double) depth);

        // Set up cells by iterating through the mesh nodes
        unsigned num_nodes = p_mesh->GetNumAllNodes();

        std::vector<CellPtr> cells;
        std::vector<CellPtr> cells2;
        std::vector<unsigned> location_indices;

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        for (unsigned i=0; i<num_nodes; i++)
        {
            c_vector<double, 3> node_location = p_mesh->GetNode(i)->rGetLocation();

            unsigned min_spatial_dimension;
            if (width <= height && width <= depth)
            {
                min_spatial_dimension = width;
            }
            else
            {
                if (height <= depth)
                {
                    min_spatial_dimension = height;
                }
                else
                {
                    min_spatial_dimension = depth;
                }
            }

            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            p_model->SetGeneration(0);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->SetBirthTime(-RandomNumberGenerator::Instance()->ranf()*
                                 (p_model->GetStemCellG1Duration() + p_model->GetSG2MDuration()));

            cells2.push_back(p_cell);

            if (norm_2(node_location - spheroid_centre) <= 0.5*sqrt(3.0)*1.01*((double) min_spatial_dimension)/3.0)
            {
                location_indices.push_back(i);
                cells.push_back(p_cell);
            }
        }

        TS_ASSERT_EQUALS(location_indices.size(), cells.size());
        TS_ASSERT_LESS_THAN(location_indices.size(), num_nodes);
        TS_ASSERT_EQUALS(location_indices.size(), 8u);

        // Test Save() with a MeshBasedCellPopulationWithGhostNodes
        MeshBasedCellPopulationWithGhostNodes<3> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestGhostNodesSpheroidSimulation3D");
        simulator.SetEndTime(0.1);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.0);
        simulator.AddForce(p_linear_force);

        simulator.Solve();
        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);

        // To check consistency with for test below
        mLocationGhosts = p_mesh->GetNode(23)->rGetLocation()[2];

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        // Test Save() with a MeshBasedCellPopulation - one cell born during this
        MeshBasedCellPopulationWithGhostNodes<3> cell_population2(*p_mesh, cells2);
        cell_population2.AddPopulationWriter<VoronoiDataWriter>();
        cell_population2.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population2.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population2.AddCellWriter<CellVolumesWriter>();

        OffLatticeSimulation<3> simulator2(cell_population2);
        simulator2.SetOutputDirectory("TestGhostNodesSpheroidSimulation3DNoGhosts");
        simulator2.SetEndTime(0.1);

        // Pass force_law to the simulation
        simulator2.AddForce(p_linear_force);

        simulator2.Solve();

        // For coverage check that mUpdateCellPopulation is archived
        simulator2.SetUpdateCellPopulationRule(false);

        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator2);

        // To check consistency with for test below
        mLocationWithoutGhosts = p_mesh->GetNode(23)->rGetLocation()[2];
        delete p_mesh;
    }

    void TestLoadOf3DSimulation()
    {
        {
            // With ghost nodes - 56 ghosts 8 real cells
            OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("TestGhostNodesSpheroidSimulation3D", 0.1);
            unsigned num_cells = p_simulator->rGetCellPopulation().GetNumRealCells();

            TS_ASSERT_EQUALS(num_cells, 8u);
            TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.1, 1e-9);
            TS_ASSERT_DELTA(p_simulator->rGetCellPopulation().GetNode(23)->rGetLocation()[2],
                            mLocationGhosts, 1e-6);

            delete p_simulator;
        }

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        {
            // Without ghost nodes - all 64 are real cells
            OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("TestGhostNodesSpheroidSimulation3DNoGhosts", 0.1);
            unsigned num_cells = p_simulator->rGetCellPopulation().GetNumRealCells();

            TS_ASSERT_EQUALS(num_cells, 64u);
            TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.1, 1e-9);
            TS_ASSERT_DELTA(p_simulator->rGetCellPopulation().GetNode(23)->rGetLocation()[2],
                            mLocationWithoutGhosts, 1e-6);
            TS_ASSERT(!p_simulator->GetUpdateCellPopulationRule());

            delete p_simulator;
        }
    }
};

#endif /*TESTOFFLATTICESIMULATION3D_HPP_*/
