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

#ifndef TESTCRYPTSIMULATION2DWITHVERTEXBASEDCELLPOPULATION_HPP_
#define TESTCRYPTSIMULATION2DWITHVERTEXBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based or crypt headers
#include "CellBasedSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "VertexCryptBoundaryForce.hpp"
#include "PopulationTestingForce.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "CryptCellsGenerator.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "Warnings.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "CellAncestorWriter.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestCryptSimulation2dWithVertexBasedCellPopulation : public AbstractCellBasedWithTimingsTestSuite
{
private:

    /**
     * Compare two meshes to see if they are 'the same'.  Doesn't check everything,
     * but is fairly thorough.  Used for testing serialization.
     */
    template<unsigned DIM>
    void CompareMeshes(VertexMesh<DIM,DIM>* pMesh1, VertexMesh<DIM,DIM>* pMesh2)
    {
        TS_ASSERT_EQUALS(pMesh1->GetNumNodes(), pMesh2->GetNumNodes());

        for (unsigned i=0; i<pMesh1->GetNumNodes(); i++)
        {
            Node<DIM>* p_node1 = pMesh1->GetNode(i);
            Node<DIM>* p_node2 = pMesh2->GetNode(i);

            TS_ASSERT_EQUALS(p_node1->IsDeleted(), p_node2->IsDeleted());
            TS_ASSERT_EQUALS(p_node1->GetIndex(), p_node2->GetIndex());

            TS_ASSERT_EQUALS(p_node1->IsBoundaryNode(), p_node2->IsBoundaryNode());

            for (unsigned j=0; j<DIM; j++)
            {
                TS_ASSERT_DELTA(p_node1->rGetLocation()[j], p_node2->rGetLocation()[j], 1e-4);
            }
        }

        TS_ASSERT_EQUALS(pMesh1->GetNumElements(), pMesh2->GetNumElements());

        for (typename VertexMesh<DIM,DIM>::VertexElementIterator iter = pMesh1->GetElementIteratorBegin();
             iter != pMesh1->GetElementIteratorEnd();
             ++iter)
        {
            unsigned elem_index = iter->GetIndex();
            VertexElement<DIM,DIM>* p_elt2 = pMesh2->GetElement(elem_index);
            TS_ASSERT_EQUALS(iter->GetNumNodes(), p_elt2->GetNumNodes());

            for (unsigned j=0; j<iter->GetNumNodes(); j++)
            {
                TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(j), p_elt2->GetNodeGlobalIndex(j));
            }
        }
    }

public:

    void TestBoundaryConditionsAtCryptBase()
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Impose a larger cell rearrangement threshold so that motion is uninhibited (see #1376)
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(cell_population);

        // Add a simple testing force
        bool positionDependentForce = false;
        MAKE_PTR_ARGS(PopulationTestingForce<2>, p_force,(positionDependentForce));
        simulator.AddForce(p_force);

        // Save old node locations
        std::vector<c_vector<double, 2> > old_node_locations(p_mesh->GetNumNodes());
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            old_node_locations[i][0] = p_mesh->GetNode(i)->rGetLocation()[0];
            old_node_locations[i][1] = p_mesh->GetNode(i)->rGetLocation()[1];
        }

        double dt = 0.01;
        simulator.SetDt(dt);
        simulator.SetupSolve();

        simulator.UpdateCellLocationsAndTopology();

        for (unsigned node_index=0; node_index<simulator.rGetCellPopulation().GetNumNodes(); node_index++)
        {
            c_vector<double, 2> node_location;
            node_location = simulator.rGetCellPopulation().GetNode(node_index)->rGetLocation();

            AbstractOffLatticeCellPopulation<2,2>* p_offLattice_pop = dynamic_cast<AbstractOffLatticeCellPopulation<2,2>* >(&(simulator.rGetCellPopulation()));
            double damping = p_offLattice_pop->GetDampingConstant(node_index);
            c_vector<double, 2> expected_location = p_force->GetExpectedOneStepLocationFE(node_index, damping, old_node_locations[node_index], dt);

            TS_ASSERT_DELTA(node_location[0], expected_location[0], 1e-9);

            if (old_node_locations[node_index][1] > 0.0)
            {
                TS_ASSERT_DELTA(node_location[1], expected_location[1], 1e-9);
            }
            else
            {
                TS_ASSERT_DELTA(node_location[1], old_node_locations[node_index][1], 1e-9);
            }
        }

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        Warnings::QuietDestroy();
    }

    void TestUsingJiggledBottomSurface()
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create cell population
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);

        simulator.SetOutputDirectory("VertexCrypt2DJiggledBottomCells");
        simulator.SetEndTime(0.01);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.UseJiggledBottomCells();

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Move the first node (which should be on y=0) down a bit
        TS_ASSERT_DELTA(crypt.GetNode(0)->rGetLocation()[1], 0.0, 1e-6);

        // Move the node (can't use the iterator for this as it is const)
        crypt.rGetMesh().GetNode(0)->rGetModifiableLocation()[1] = -1.0;
        TS_ASSERT_LESS_THAN(crypt.GetNode(0)->rGetLocation()[1], 0.0);

        // The time step should have been modified in the constructor
        TS_ASSERT_DELTA(simulator.GetDt(), 0.002, 1e-12);

        // Run simulation
        simulator.Solve();

        // The node should have been pulled up, but not above y=0. However it should
        // then been moved to above y=0 by the jiggling
        TS_ASSERT_LESS_THAN(0.0, crypt.GetNode(0)->rGetLocation()[1]);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    /**
     * Test that a short crypt simulation without cell birth runs without throwing any errors.
     */
    void TestCryptWithNoBirth()
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells, all differentiated
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 0.0, 0.0, 0.0, 0.0);

        // Create cell population
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetEndTime(0.1);
        simulator.SetSamplingTimestepMultiple(50);

        simulator.SetOutputDirectory("TestVertexCryptWithNoBirth");

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        // No cell killer

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    /**
     * Test that a short crypt simulation, in which cell birth occurs,
     * runs without throwing any errors.
     */
    void TestCryptWithBirth()
    {
        double crypt_length = 5.0;

        // Create mesh
        CylindricalHoneycombVertexMeshGenerator generator(4, 6);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells: the bottom row have StemCellProliferativeType and the rest have DifferentiatedCellProliferativeType
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 0.8, 0.8, 0.8, 0.8);

        // Cell 1 should divide at time t=0.05
        cells[0]->SetBirthTime(-23.95);

        // Cells 2-4 should divide later
        cells[1]->SetBirthTime(-23.0);
        cells[2]->SetBirthTime(-22.0);
        cells[3]->SetBirthTime(-21.0);

        // Create cell population
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(0.1);
        simulator.SetOutputDirectory("TestVertexCryptWithBirth");

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Make crypt shorter for sloughing
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 2u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "A Cell is removed without performing a T2 swap. This could leave a void in the mesh.");
        Warnings::QuietDestroy();
    }

    /**
     * Commented test of a long crypt simulation. Used to generate attachment
     * VertexSimulation.mpeg on #1095.
     */
    void noTestCryptSimulationLong()
    {
        double crypt_length = 20.0;

        // Create mesh
        unsigned crypt_width = 10;
        unsigned crypt_height = 20;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height, true);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<UniformG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create cell population
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(10);
        simulator.SetOutputDirectory("TestVertexCryptLong");

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Make crypt shorter for sloughing
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    /**
     * Set up and briefly solve a vertex crypt simulation in which
     * cell proliferation is Wnt-based, to check that WntConcentration
     * doesn't throw a wobbly.
     */
    void TestShortWntBasedCryptSimulation()
    {
        double crypt_length = 10.0;

        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height, true);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create cell population
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Set up Wnt gradient
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(0.1);
        simulator.SetOutputDirectory("TestShortWntBasedCryptSimulation");

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Make crypt shorter for sloughing
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    /** Longer Wnt based simulation
     */
    void noTestWntBasedCryptSimulationLong()
    {
        double crypt_length = 20.0;

        // Create mesh
        unsigned crypt_width = 10;
        unsigned crypt_height = 20;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height, true);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create cell population
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Set up Wnt gradient
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(50.0);
        simulator.SetOutputDirectory("TestLongWntBasedVertexCryptSimulation");

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Make crypt shorter for sloughing
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    // Test a crypt simulation with a boundary force on the crypt base.
    void TestCryptSimulationWithBoundaryForce()
    {
        double crypt_length = 6.0;

        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells: the bottom row have StemCellProliferativeType and the rest have DifferentiatedCellProliferativeType
        std::vector<CellPtr> cells;
        CryptCellsGenerator<UniformG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 0.8, 0.8, 0.8, 0.8);

        // Create cell population
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetSamplingTimestepMultiple(50);
        double end_time = 0.1;
        simulator.SetEndTime(end_time);
        simulator.SetOutputDirectory("TestVertexCryptWithBoundaryForce");

        // Create a force laws and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);
        MAKE_PTR_ARGS(VertexCryptBoundaryForce<2>, p_boundary_force_law, (150));
        simulator.AddForce(p_boundary_force_law);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Make crypt shorter for sloughing
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Coverage
        CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);
        CryptSimulation2d* p_simulator;
        p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load("TestVertexCryptWithBoundaryForce", end_time);
        delete p_simulator;

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    /**
     * Test that archiving a crypt simulation correctly archives its mesh.
     */
    void TestMeshSurvivesSaveLoad()
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<UniformG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create cell population
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("VertexCrypt2DArchive");
        simulator.SetEndTime(0.1);

        // Create a force laws and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        /*
         * Memory leak (unconditional jump) without the following line. The
         * archiver assumes that a Solve has been called and simulation time
         * has been set up properly. In this test it hasn't, so we need this
         * to avoid a memory leak.
         */
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.1, 100);

        // Save
        CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);

        // Load
        CryptSimulation2d* p_simulator;
        p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load("VertexCrypt2DArchive", 0.0);

        // Create an identical mesh for comparison purposes
        Cylindrical2dVertexMesh* p_mesh2 = generator.GetCylindricalMesh();

        // Compare meshes
        VertexMesh<2,2>& r_mesh = (static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation())))->rGetMesh();
        CompareMeshes(p_mesh2, &r_mesh);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        Warnings::QuietDestroy();

        // Tidy up
        delete p_simulator;
    }

    void TestStandardResultForArchivingTestsBelow()
    {
        double crypt_length = 22.0;

        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<UniformG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create cell population
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // We have a Wnt Gradient - but not Wnt dependent cells
        // so that the test runs quickly, but we test archiving of it!
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("VertexCrypt2DPeriodicStandardResult");
        simulator.SetEndTime(0.25);

        // Create a force laws and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Test the locations of a few nodes
        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 0.0641, 1e-4);

        std::vector<double> node_5_location = simulator.GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 1.0, 1e-3);
        TS_ASSERT_DELTA(node_5_location[1], 0.0641, 1e-4);

        // Test the Wnt concentration result
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(crypt.GetCellUsingLocationIndex(2)), 0.9769, 1e-4);
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(crypt.GetCellUsingLocationIndex(3)), 0.9769, 1e-4);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    void TestVertexCryptSimulation2DParameterOutput()
    {
        double crypt_length = 22.0;

        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<UniformG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create cell population
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetEndTime(1.5);

        // Create a force laws and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        std::string output_directory = "TestVertexCryptSimulation2dOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);
        out_stream parameter_file = output_file_handler.OpenOutputFile("vertex_crypt_sim_2d_results.parameters");
        simulator.OutputSimulationParameters(parameter_file);
        parameter_file->close();

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
        FileComparison( results_dir + "vertex_crypt_sim_2d_results.parameters", "crypt/test/data/TestVertexCryptSimulationOutputParameters/vertex_crypt_sim_2d_results.parameters").CompareFiles();

        ///\todo check output of simulator.OutputSimulationSetup()
    }

    // Testing Save
    void TestSave()
    {
        double crypt_length = 22.0;

        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<UniformG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create cell population
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("VertexCrypt2DPeriodicSaveAndLoad");

        // Our full end time is 0.25, here we run until 0.1 then load and run more below.
        simulator.SetEndTime(0.1);

        // Create a force laws and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Save the results
        CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    // Testing Load (based on previous two tests)
    void TestLoad()
    {
        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 0.2
        CryptSimulation2d* p_simulator1;

        WntConcentration<2>::Instance();   // Make sure there is no existing Wnt Gradient before load
        WntConcentration<2>::Destroy();

        p_simulator1 = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load("VertexCrypt2DPeriodicSaveAndLoad", 0.1);
        p_simulator1->SetEndTime(0.2);
        p_simulator1->Solve();

        // Get mesh
        MutableVertexMesh<2,2>& r_mesh1 = (static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator1->rGetCellPopulation())))->rGetMesh();

        // Save then reload, compare meshes either side
        CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(p_simulator1);

        CryptSimulation2d* p_simulator2 = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load("VertexCrypt2DPeriodicSaveAndLoad", 0.2);
        MutableVertexMesh<2,2>& r_mesh2 = (static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator2->rGetCellPopulation())))->rGetMesh();

        CompareMeshes(&r_mesh1, &r_mesh2);

        // Run a bit further...
        p_simulator2->SetEndTime(0.25);

        // Run simulation
        p_simulator2->Solve();

        // Test the locations of a few nodes
        std::vector<double> node_4_location = p_simulator2->GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 0.0641, 1e-4);

        std::vector<double> node_5_location = p_simulator2->GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 1.0, 1e-3);
        TS_ASSERT_DELTA(node_5_location[1], 0.0641, 1e-4);

        // Test Wnt concentration was set up correctly
        TS_ASSERT_EQUALS(WntConcentration<2>::Instance()->IsWntSetUp(), true);

        // Test the Wnt concentration result
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(p_simulator2->rGetCellPopulation().GetCellUsingLocationIndex(2)), 0.9769, 1e-4);
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(p_simulator2->rGetCellPopulation().GetCellUsingLocationIndex(3)), 0.9769, 1e-4);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();

        // Tidy up
        delete p_simulator1;
        delete p_simulator2;
        WntConcentration<2>::Destroy();
    }

    void TestWriteBetaCateninAndAncestors()
    {
        // Create mesh
        unsigned crypt_width = 6;
        unsigned crypt_height = 4;
        double crypt_length = crypt_height*(sqrt(3.0)/2);
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<VanLeeuwen2009WntSwatCellCycleModelHypothesisOne> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create crypt
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Set crypt to output cell types and cell ancestors
        crypt.AddCellWriter<CellAncestorWriter>();

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("VertexCryptBetaCatenin");
        simulator.SetEndTime(0.1);
        simulator.SetBottomCellAncestors();

        // Create a force laws and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Check writing of beta-catenin data
        OutputFileHandler handler("VertexCryptBetaCatenin", false);
        std::string ancestor_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizancestors";
        std::string results_setup_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizsetup";

        FileComparison( ancestor_results_file, "crypt/test/data/VertexCryptBetaCatenin/results.vizancestors").CompareFiles();
        FileComparison( results_setup_file, "crypt/test/data/VertexCryptBetaCatenin/results.vizsetup").CompareFiles();

        // Tidy up
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTCRYPTSIMULATION2DWITHVERTEXBASEDCELLPOPULATION_HPP_*/
