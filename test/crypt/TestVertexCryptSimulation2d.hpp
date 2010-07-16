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
#ifndef TESTVERTEXCRYPTSIMULATION2D_HPP_
#define TESTVERTEXCRYPTSIMULATION2D_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "TissueSimulationArchiver.hpp"

#include "VertexCryptSimulation2d.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "NagaiHondaForce.hpp"
#include "VertexCryptBoundaryForce.hpp"
#include "CryptCellsGenerator.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "Warnings.hpp"

class TestVertexCryptSimulation2d : public AbstractCellBasedTestSuite
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

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void TestBoundaryConditionsAtCryptBase() throw (Exception)
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Impose a larger cell rearrangement threshold so that motion is uninhibited (see #1376)
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(tissue, force_collection);

        std::vector<c_vector<double, 2> > old_node_locations(p_mesh->GetNumNodes());
        std::vector<c_vector<double, 2> > forces(p_mesh->GetNumNodes());

        // Make up some forces
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            old_node_locations[i][0] = p_mesh->GetNode(i)->rGetLocation()[0];
            old_node_locations[i][1] = p_mesh->GetNode(i)->rGetLocation()[1];

            forces[i][0] = i*0.01;
            forces[i][1] = 2*i*0.01;
       }

        simulator.SetDt(0.01);
        simulator.UpdateNodePositions(forces);

        for (unsigned node_index=0; node_index<simulator.rGetTissue().GetNumNodes(); node_index++)
        {
            c_vector<double, 2> node_location = simulator.rGetTissue().GetNode(node_index)->rGetLocation();

            TS_ASSERT_DELTA(node_location[0], old_node_locations[node_index][0] + node_index*0.01*0.01, 1e-9);

            if (old_node_locations[node_index][1] > 0.0)
            {
                TS_ASSERT_DELTA(node_location[1], old_node_locations[node_index][1] + 2*node_index*0.01*0.01, 1e-9);
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
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create tissue
        VertexBasedTissue<2> crypt(*p_mesh, cells);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);

        simulator.SetOutputDirectory("VertexCrypt2DJiggledBottomCells");
        simulator.SetEndTime(0.01);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.UseJiggledBottomCells();

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
        ///\todo this should be only one warning (see #1394)
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 136u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearangementThreshold this could cause elements to become inverted the motion has been restricted: - To avoid these warnings use a smaller timestep");
        Warnings::QuietDestroy();
    }

    /**
     * Test that a short crypt simulation without cell birth runs without throwing any errors.
     */
    void TestCryptWithNoBirth() throw (Exception)
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells, all differentiated
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 0.0, 0.0, 0.0, 0.0);

        // Create tissue
        VertexBasedTissue<2> crypt(*p_mesh, cells);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetEndTime(0.1);
        simulator.SetSamplingTimestepMultiple(50);

        simulator.SetOutputDirectory("TestVertexCryptWithNoBirth");

        // No cell killer

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test Warnings
        ///\todo this should be only one warning (see #1394)
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 972u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearangementThreshold this could cause elements to become inverted the motion has been restricted: - To avoid these warnings use a smaller timestep");
        Warnings::QuietDestroy();
    }

    /**
     * Test that a short crypt simulation, in which cell birth occurs,
     * runs without throwing any errors.
     */
    void TestCryptWithBirth() throw (Exception)
    {
        // Create mesh
        CylindricalHoneycombVertexMeshGenerator generator(4, 6);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells, bottom row are STEM rest are DIFFERENTIATED.
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 0.8,0.8,0.8,0.8);

        // Cell 1 should divide at time t=0.05
        cells[0]->SetBirthTime(-23.95);
        // Cell 2-4 should divide later
        cells[1]->SetBirthTime(-23.0);
        cells[2]->SetBirthTime(-22.0);
        cells[3]->SetBirthTime(-21.0);

        // Create tissue
        VertexBasedTissue<2> crypt(*p_mesh, cells);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(0.1);
        simulator.SetOutputDirectory("TestVertexCryptWithBirth");

        // Make crypt shorter for sloughing
        TissueConfig::Instance()->SetCryptLength(5.0);
        SloughingCellKiller<2> sloughing_cell_killer(&crypt);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test Warnings
        ///\todo this should be only one warning (see #1394)
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 971u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Vertices are moving more than half the CellRearangementThreshold this could cause elements to become inverted the motion has been restricted: - To avoid these warnings use a smaller timestep");
        Warnings::QuietDestroy();
    }

    /**
     * Commented test of a long crypt simulation. Used to generate attachment
     * VertexSimulation.mpeg on #1095.
     */
    void noTestCryptSimulationLong() throw (Exception)
    {
        // Create mesh
        unsigned crypt_width = 10;
        unsigned crypt_height = 20;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height, true);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create tissue
        VertexBasedTissue<2> crypt(*p_mesh, cells);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(10);
        simulator.SetOutputDirectory("TestVertexCryptLong");

        // Make crypt shorter for sloughing
        TissueConfig::Instance()->SetCryptLength(20.0);
        SloughingCellKiller<2> sloughing_cell_killer(&crypt);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test Warnings
        ///\todo this should be only one warning (see #1394)
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u); // Number not known
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearangementThreshold this could cause elements to become inverted the motion has been restricted: - To avoid these warnings use a smaller timestep");
        Warnings::QuietDestroy();
    }

    /**
     * Set up and briefly solve a vertex crypt simulation in which
     * cell proliferation is Wnt-based, to check that WntConcentration
     * doesn't throw a wobbly.
     */
    void TestShortWntBasedCryptSimulation() throw (Exception)
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height, true);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create tissue
        VertexBasedTissue<2> crypt(*p_mesh, cells);

        // Set up Wnt gradient
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(0.1);
        simulator.SetOutputDirectory("TestShortWntBasedCryptSimulation");

        // Make crypt shorter for sloughing
        TissueConfig::Instance()->SetCryptLength(10.0);
        SloughingCellKiller<2> sloughing_cell_killer(&crypt);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test Warnings
        ///\todo this should be only one warning (see #1394)
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 761u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearangementThreshold this could cause elements to become inverted the motion has been restricted: - To avoid these warnings use a smaller timestep");
        Warnings::QuietDestroy();

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    /** Longer Wnt based simulation
     */
    void noTestWntBasedCryptSimulationLong() throw (Exception)
    {
        // Create mesh
        unsigned crypt_width = 10;
        unsigned crypt_height = 20;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height, true);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create tissue
        VertexBasedTissue<2> crypt(*p_mesh, cells);

        // Set up Wnt gradient
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(50.0);
        simulator.SetOutputDirectory("TestLongWntBasedVertexCryptSimulation");

        // Make crypt shorter for sloughing
        TissueConfig::Instance()->SetCryptLength(20.0);
        SloughingCellKiller<2> sloughing_cell_killer(&crypt);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test Warnings
        ///\todo this should be only one warning (see #1394)
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u); // Number not known.
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearangementThreshold this could cause elements to become inverted the motion has been restricted: - To avoid these warnings use a smaller timestep");
        Warnings::QuietDestroy();

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    // Test a crypt simulation with a boundary force on the crypt base.
    void TestCryptSimulationWithBoundaryForce() throw (Exception)
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells, bottom row are STEM rest are DIFFERENTIATED.
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 0.8, 0.8, 0.8, 0.8);

        // Create tissue
        VertexBasedTissue<2> crypt(*p_mesh, cells);

        // Create boundary force law
        VertexCryptBoundaryForce<2> boundary_force_law(150);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;

        force_collection.push_back(&boundary_force_law);
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetSamplingTimestepMultiple(50);
        double end_time = 0.1;
        simulator.SetEndTime(end_time);
        simulator.SetOutputDirectory("TestVertexCryptWithBoundaryForce");

        // Make crypt shorter for sloughing
        TissueConfig::Instance()->SetCryptLength(6.0);
        SloughingCellKiller<2> sloughing_cell_killer(&crypt);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Coverage
        TissueSimulationArchiver<2, VertexCryptSimulation2d>::Save(&simulator);
        VertexCryptSimulation2d* p_simulator;
        p_simulator = TissueSimulationArchiver<2, VertexCryptSimulation2d>::Load("TestVertexCryptWithBoundaryForce", end_time);
        delete p_simulator;

        // Test Warnings
        ///\todo this should be only one warning (see #1394)
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 968u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearangementThreshold this could cause elements to become inverted the motion has been restricted: - To avoid these warnings use a smaller timestep");
        Warnings::QuietDestroy();
    }

    /**
     * Test that archiving a crypt simulation correctly archives its mesh.
     */
    void TestMeshSurvivesSaveLoad() throw (Exception)
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create tissue
        VertexBasedTissue<2> crypt(*p_mesh, cells);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("VertexCrypt2DArchive");
        simulator.SetEndTime(0.1);

        /*
         * Memory leak (unconditional jump) without the following line. The
         * archiver assumes that a Solve has been called and simulation time
         * has been set up properly. In this test it hasn't, so we need this
         * to avoid a memory leak.
         */
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.1, 100);

        // Save
        TissueSimulationArchiver<2, VertexCryptSimulation2d>::Save(&simulator);

        // Load
        VertexCryptSimulation2d* p_simulator;
        p_simulator = TissueSimulationArchiver<2, VertexCryptSimulation2d>::Load("VertexCrypt2DArchive", 0.0);

        // Create an identical mesh for comparison purposes
        Cylindrical2dVertexMesh* p_mesh2 = generator.GetCylindricalMesh();

        // Compare meshes
        VertexMesh<2,2>& r_mesh = (static_cast<VertexBasedTissue<2>*>(&(p_simulator->rGetTissue())))->rGetMesh();
        CompareMeshes(p_mesh2, &r_mesh);

        //Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        Warnings::QuietDestroy();

        // Tidy up
        delete p_simulator;
    }

    void TestStandardResultForArchivingTestsBelow() throw (Exception)
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create tissue
        VertexBasedTissue<2> crypt(*p_mesh, cells);

        // We have a Wnt Gradient - but not Wnt dependent cells
        // so that the test runs quickly, but we test archiving of it!
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("VertexCrypt2DPeriodicStandardResult");
        simulator.SetEndTime(0.25);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller<2> sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();

        // Test the locations of a few nodes
        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 3.9992, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 0.1757, 1e-4);

        std::vector<double> node_5_location = simulator.GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.9755, 1e-3);
        TS_ASSERT_DELTA(node_5_location[1], 0.1715, 1e-4);

        // Test the Wnt concentration result
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(crypt.GetCellUsingLocationIndex(2)), 0.9748, 1e-4);
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(crypt.GetCellUsingLocationIndex(3)), 0.9736, 1e-4);

        // Test Warnings
        ///\todo this should be only one warning (see #1394)
        //TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1178u); // This gives 1177 on bob!! see #1376
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Vertices are moving more than half the CellRearangementThreshold this could cause elements to become inverted the motion has been restricted: - To avoid these warnings use a smaller timestep");
        Warnings::QuietDestroy();

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    // Testing Save
    void TestSave() throw (Exception)
    {
        // Create mesh
        unsigned crypt_width = 4;
        unsigned crypt_height = 6;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create tissue
        VertexBasedTissue<2> crypt(*p_mesh, cells);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("VertexCrypt2DPeriodicSaveAndLoad");

        // Our full end time is 0.25, here we run until 0.1 then load and run more below.
        simulator.SetEndTime(0.1);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller<2> sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();

        // Save the results
        TissueSimulationArchiver<2, VertexCryptSimulation2d>::Save(&simulator);

        // Test Warnings
        ///\todo this should be only one warning (see #1394)
        //TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 905u); // This gives  902 on bob!! see #1394
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Vertices are moving more than half the CellRearangementThreshold this could cause elements to become inverted the motion has been restricted: - To avoid these warnings use a smaller timestep");
        Warnings::QuietDestroy();

        // Tidy up
        WntConcentration<2>::Destroy();
    }


    // Testing Load (based on previous two tests)
    void TestLoad() throw (Exception)
    {
        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 0.2
        VertexCryptSimulation2d* p_simulator1;

        WntConcentration<2>::Instance();   // Make sure there is no existing Wnt Gradient before load
        WntConcentration<2>::Destroy();

        p_simulator1 = TissueSimulationArchiver<2, VertexCryptSimulation2d>::Load("VertexCrypt2DPeriodicSaveAndLoad", 0.1);
        p_simulator1->SetEndTime(0.2);
        p_simulator1->Solve();

        // Get mesh
        MutableVertexMesh<2,2>& r_mesh1 = (static_cast<VertexBasedTissue<2>*>(&(p_simulator1->rGetTissue())))->rGetMesh();

        // Save then reload, compare meshes either side
        TissueSimulationArchiver<2, VertexCryptSimulation2d>::Save(p_simulator1);

        VertexCryptSimulation2d* p_simulator2 = TissueSimulationArchiver<2, VertexCryptSimulation2d>::Load("VertexCrypt2DPeriodicSaveAndLoad", 0.2);
        MutableVertexMesh<2,2>& r_mesh2 = (static_cast<VertexBasedTissue<2>*>(&(p_simulator2->rGetTissue())))->rGetMesh();

        CompareMeshes(&r_mesh1, &r_mesh2);

        // Run a bit further...
        p_simulator2->SetEndTime(0.25);

        // Run simulation
        p_simulator2->Solve();

        // Test the locations of a few nodes
        std::vector<double> node_4_location = p_simulator2->GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 3.9992, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 0.1757, 1e-4);

        std::vector<double> node_5_location = p_simulator2->GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.9755, 1e-3);
        TS_ASSERT_DELTA(node_5_location[1], 0.1715, 1e-4);

        // Test Wnt concentration was set up correctly
        TS_ASSERT_EQUALS(WntConcentration<2>::Instance()->IsWntSetUp(), true);

        // Test the Wnt concentration result
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(p_simulator2->rGetTissue().GetCellUsingLocationIndex(2)), 0.9748, 1e-4);
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(p_simulator2->rGetTissue().GetCellUsingLocationIndex(3)), 0.9736, 1e-4);

        // Test Warnings
        ///\todo this should be only one warning (see #1394)
        //TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 273u); // This gives  275 on bob!! see #1376
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Vertices are moving more than half the CellRearangementThreshold this could cause elements to become inverted the motion has been restricted: - To avoid these warnings use a smaller timestep");
        Warnings::QuietDestroy();

        // Tidy up
        delete p_simulator1;
        delete p_simulator2;
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTVERTEXCRYPTSIMULATION2D_HPP_*/
