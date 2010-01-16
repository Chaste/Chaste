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
#ifndef TESTTISSUESIMULATIONWITHVERTEXBASEDTISSUE_HPP_
#define TESTTISSUESIMULATIONWITHVERTEXBASEDTISSUE_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "TissueSimulationArchiver.hpp"

#include "TissueSimulation.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "VertexBasedTissue.hpp"
#include "NagaiHondaForce.hpp"
#include "AbstractCellKiller.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexMeshWriter.hpp"

/**
 * Simple cell killer which at the fisrt timestep kills any cell
 * whose corresponding location index is a given number.
 *
 * For testing purposes.
 */
class TargetedCellKiller : public AbstractCellKiller<2>
{
private :

    unsigned mTargetIndex;
    bool mBloodLust;

public :
    TargetedCellKiller(AbstractTissue<2>* pTissue, unsigned targetIndex)
        : AbstractCellKiller<2>(pTissue),
          mTargetIndex(targetIndex),
          mBloodLust(true)
    {
    }

    virtual void TestAndLabelCellsForApoptosisOrDeath()
    {
        if ( !mBloodLust || mpTissue->GetNumRealCells()==0 || mpTissue->GetNumRealCells()<mTargetIndex+1)
        {
            return;
        }
        mpTissue->rGetCellUsingLocationIndex(mTargetIndex).Kill();
        mBloodLust = false;
    }
};


class TestTissueSimulationWithVertexBasedTissue : public AbstractCellBasedTestSuite
{
private:

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

    void TestSolveThrowsNothing() throw (Exception)
    {
        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(6, 6);
        VertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement. Give each cell
        // a birth time of -elem_index, so its age is elem_index
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            TissueCell cell(DIFFERENTIATED, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - elem_index;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        // Create a force system
        NagaiHondaForce<2> force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("TestSolveThrowsNothing");
        simulator.SetEndTime(0.1);

        TS_ASSERT_DELTA(simulator.GetDt(), 0.002, 1e-12);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }


    void TestSingleCellRelaxation() throw (Exception)
    {
        // Construct a 2D vertex mesh consisting of a single element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 20;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        double cell_swap_threshold = 0.01;
        double edge_division_threshold = 2.0;
        VertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

        // Set up cells, one for each VertexElement. Give each cell
        // a birth time of 0
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<mesh.GetNumElements(); elem_index++)
        {
            CellProliferativeType cell_type = DIFFERENTIATED;
            double birth_time = -1.0;

            TissueCell cell(cell_type, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);

        // Create a force system
        NagaiHondaForce<2> force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("TestSingleCellRelaxation");
        simulator.SetEndTime(1.0);

        // Run simulation
        simulator.Solve();

        // Test relaxes to circle (can be more stringent with more nodes and more time)
        TS_ASSERT_DELTA(tissue.rGetMesh().GetAreaOfElement(0), 1.0, 0.05);
        TS_ASSERT_DELTA(tissue.rGetMesh().GetPerimeterOfElement(0), 3.5449077, 0.1);
    }

    /*
     * In this test edges can divide if they become too long.
     * \todo This isn't dealt with by ReMesh() so will fall
     * over at later times when a T1 swap occurs (see #1110)
     */
    void TestSimpleVertexMonolayerWithCellBirth() throw (Exception)
    {
        // Create a simple 2D VertexMesh with only one cell
        HoneycombVertexMeshGenerator generator(1, 1);
        VertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetEdgeDivisionThreshold(0.5);

        // Set up cell.
        std::vector<TissueCell> cells;
        CellProliferativeType cell_type = TRANSIT;
        double birth_time = -20.0; // Divides Straight Away

        TissueCell cell(cell_type, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
        cell.SetBirthTime(birth_time);
        cells.push_back(cell);

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        unsigned old_num_nodes = tissue.GetNumNodes();
        unsigned old_num_elements = tissue.GetNumElements();
        unsigned old_num_cells = tissue.GetNumRealCells();

        // Create a force system
        NagaiHondaForce<2> force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("TestSimpleVertexMonolayerWithCellBirth");
        simulator.SetEndTime(1.0);

        // Run simulation
        simulator.Solve();

        // Check that cell divided successfully
        unsigned new_num_nodes = simulator.rGetTissue().GetNumNodes();
        unsigned new_num_elements = (static_cast<VertexBasedTissue<2>*>(&(simulator.rGetTissue())))->GetNumElements();
        unsigned new_num_cells = simulator.rGetTissue().GetNumRealCells();

        TS_ASSERT_EQUALS(new_num_nodes, old_num_nodes+9); // as division of element is longer than threshold so is divided
        TS_ASSERT_EQUALS(new_num_elements, old_num_elements+1);
        TS_ASSERT_EQUALS(new_num_cells, old_num_cells+1);
    }

    void TestVertexMonolayerWithCellBirth() throw (Exception)
    {
        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(5, 5);
        VertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement. Give each cell
        // a birth time of -elem_index, so its age is elem_index
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            CellProliferativeType cell_type = DIFFERENTIATED;
            double birth_time = 0.0 - elem_index;

            // Cell should divide at time t=0.5
            if (elem_index==12) //(elem_index==6 ||elem_index==7 || elem_index==8 || elem_index==11 ||elem_index==12 || elem_index==13 ||elem_index==16 ||elem_index==17 || elem_index==18)
            {
                 cell_type = STEM;
                birth_time = -23.5;
            }

            TissueCell cell(cell_type, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        unsigned old_num_nodes = tissue.GetNumNodes();
        unsigned old_num_elements = tissue.GetNumElements();
        unsigned old_num_cells = tissue.GetNumRealCells();

        // Create a force system
        NagaiHondaForce<2> force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("TestVertexMonolayerWithCellBirth");
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(1.0);

        // Run simulation
        simulator.Solve();

        // Check that cell 12 divided successfully
        unsigned new_num_nodes = simulator.rGetTissue().GetNumNodes();
        unsigned new_num_elements = (static_cast<VertexBasedTissue<2>*>(&(simulator.rGetTissue())))->GetNumElements();
        unsigned new_num_cells = simulator.rGetTissue().GetNumRealCells();

        TS_ASSERT_EQUALS(new_num_nodes, old_num_nodes+2); // as division of element is longer than threshold so is divided
        TS_ASSERT_EQUALS(new_num_elements, old_num_elements+1);
        TS_ASSERT_EQUALS(new_num_cells, old_num_cells+1);
    }

    void noTestVertexMonolayerLong() throw (Exception)
    {
        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(3, 3);
        VertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement.
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            CellProliferativeType cell_type = TRANSIT;

            double birth_time = - RandomNumberGenerator::Instance()->ranf()*
                                 ( TissueConfig::Instance()->GetTransitCellG1Duration()
                                    + TissueConfig::Instance()->GetSG2MDuration() );


            TissueCell cell(cell_type, HEALTHY, new StochasticDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        // Create a force system
        NagaiHondaForce<2> force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("TestVertexMonolayerLong");
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(30.0); // at 35.6 a void forms and need to deal with this!!

        // Adjust Max Generations so cells keep proliferating
        TissueConfig::Instance()->SetMaxTransitGenerations(8u);

        // Run simulation
        simulator.Solve();
    }


    void TestVertexMonolayerWithCellDeath() throw (Exception)
    {
        /*
         * We don't want apoptosing cells to be labelled as dead after a certain time in
         * vertex simulations, so set the apoptosis time to something large.
         */

        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(5, 5);
        VertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.2);
        p_mesh->SetT2Threshold(sqrt(3.0)/1000.0); // so T2Swaps once it becomes a triangle

        // Set up cells, one for each VertexElement. Give each cell
        // a birth time of -elem_index, so its age is elem_index
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            double birth_time = 0.0 - elem_index;
            TissueCell cell(DIFFERENTIATED, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);

            if (elem_index==18)
            {
                cell.StartApoptosis(false);
            }

            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        unsigned old_num_nodes = tissue.GetNumNodes();
        unsigned old_num_elements = tissue.GetNumElements();
        unsigned old_num_cells = tissue.GetNumRealCells();

        // Create a force system
        NagaiHondaForce<2> force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("TestVertexMonolayerWithCellDeath");
        simulator.SetEndTime(4.5); // Any longer and cell needs to T2 Swap \todo implement T2 Swaps

        // Longer appoptosis time so cells shrink over a longer time
        TissueConfig::Instance()->SetApoptosisTime(1.5);

        // Create a cell killer and pass in to simulation (note we must account for element index changes following each kill)
        TargetedCellKiller cell0_killer(&tissue, 0);    // element on the bottom boundary
        TargetedCellKiller cell2_killer(&tissue, 2);    // element in the interior
        TargetedCellKiller cell12_killer(&tissue, 12);  // element on the corner boundary

        simulator.AddCellKiller(&cell0_killer);
        simulator.AddCellKiller(&cell2_killer);
        simulator.AddCellKiller(&cell12_killer);

        // Run simulation
        simulator.Solve();

        // Check that cell 8 has died
        unsigned new_num_nodes = simulator.rGetTissue().GetNumNodes();
        unsigned new_num_elements = (static_cast<VertexBasedTissue<2>*>(&(simulator.rGetTissue())))->GetNumElements();
        unsigned new_num_cells = simulator.rGetTissue().GetNumRealCells();

        TS_ASSERT_EQUALS(new_num_nodes, old_num_nodes-6);
        TS_ASSERT_EQUALS(new_num_elements, old_num_elements-4);
        TS_ASSERT_EQUALS(new_num_cells, old_num_cells-3); //\todo this should match with the elements
    }

    void TestSingleCellRelaxationAndApoptosis() throw (Exception)
    {
        // Construct a 2D vertex mesh consisting of a single element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 20;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        VertexMesh<2,2> mesh(nodes, elements);
        mesh.SetCellRearrangementThreshold(0.1);

        // Set up cells, one for each VertexElement
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<mesh.GetNumElements(); elem_index++)
        {
            CellProliferativeType cell_type = DIFFERENTIATED;
            double birth_time = -1.0;

            TissueCell cell(cell_type, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);

        // Create a force system
        NagaiHondaForce<2> force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("TestVertexSingleCellApoptosis");
        simulator.SetEndTime(2.0);

        // Run simulation
        simulator.Solve();

        // Test relaxes to circle (can be more stringent with more nodes and more time)
        TS_ASSERT_DELTA(tissue.rGetMesh().GetAreaOfElement(0), 1.0, 1e-1);
        TS_ASSERT_DELTA(tissue.rGetMesh().GetPerimeterOfElement(0), 3.5449077, 1e-1);

        TissueCell& r_cell = simulator.rGetTissue().rGetCellUsingLocationIndex(0);
        r_cell.StartApoptosis(false);

        simulator.SetEndTime(2.25); // any longer and cell target area is zero but element can't be removed as its the only one.

        // Run simulation
        simulator.Solve();

        TS_ASSERT_DELTA(tissue.rGetMesh().GetAreaOfElement(0), 0.0068, 1e-4);
        TS_ASSERT_DELTA(tissue.rGetMesh().GetPerimeterOfElement(0), 0.3783, 1e-3);
    }

    /**
     * Test archiving of a TissueSimulation that uses a VertexBasedTissue.
     */
    void TestArchiving() throw (Exception)
    {
        // Set end time
        double end_time = 0.1;

        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(6, 6);
        VertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement. Give each cell
        // a birth time of -elem_index, so its age is elem_index
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            TissueCell cell(DIFFERENTIATED, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - elem_index;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        // Create a force system
        NagaiHondaForce<2> force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("TestTissueSimulationWithVertexBasedTissueSaveAndLoad");
        simulator.SetEndTime(end_time);

        // Run simulation
        simulator.Solve();

        TissueSimulationArchiver<2, TissueSimulation<2> >::Save(&simulator);

        TissueSimulation<2>* p_simulator
            = TissueSimulationArchiver<2, TissueSimulation<2> >::Load("TestTissueSimulationWithVertexBasedTissueSaveAndLoad", end_time);

        p_simulator->SetEndTime(0.2);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(p_simulator->Solve());

        /// \todo add further tests (see #821 and #862)

        // Tidy up
        delete p_simulator;
    }
};

#endif /*TESTTISSUESIMULATIONWITHVERTEXBASEDTISSUE_HPP_*/
