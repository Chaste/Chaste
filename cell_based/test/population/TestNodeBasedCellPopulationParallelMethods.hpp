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

#ifndef TESTNODEBASEDCELLPOPULATIONPARALLELMETHODS_HPP_
#define TESTNODEBASEDCELLPOPULATIONPARALLELMETHODS_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "TrianglesMeshReader.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "CellPropertyRegistry.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellVolumesWriter.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestNodeBasedCellPopulationParallelMethods : public AbstractCellBasedTestSuite
{
private:

    NodesOnlyMesh<3>* mpNodesOnlyMesh;

    NodeBasedCellPopulation<3>* mpNodeBasedCellPopulation;

    void setUp()
    {
        AbstractCellBasedTestSuite::setUp();

        std::vector<Node<3>* > nodes;
        for (unsigned i=0; i<PetscTools::GetNumProcs(); i++)
        {
            nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.5+(double)i));
        }

        mpNodesOnlyMesh = new NodesOnlyMesh<3>;
        mpNodesOnlyMesh->ConstructNodesWithoutMesh(nodes, 1.0);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mpNodesOnlyMesh->GetNumNodes());

        mpNodeBasedCellPopulation = new NodeBasedCellPopulation<3>(*mpNodesOnlyMesh, cells);

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void tearDown()
    {
        delete mpNodeBasedCellPopulation;
        delete mpNodesOnlyMesh;

        AbstractCellBasedTestSuite::tearDown();
    }
public:

    void TestGetCellAndNodePair()
    {
        unsigned node_index = mpNodesOnlyMesh->GetNodeIteratorBegin()->GetIndex();

        std::pair<CellPtr, Node<3>* > pair = mpNodeBasedCellPopulation->GetCellNodePair(node_index);

        TS_ASSERT_EQUALS(pair.second->GetIndex(), node_index);

        CellPtr p_returned_cell = pair.first;
        TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->GetLocationIndexUsingCell(p_returned_cell), node_index);
    }

    void TestAddNodeAndCellsToSend()
    {
        unsigned index_of_node_to_send = mpNodesOnlyMesh->GetNodeIteratorBegin()->GetIndex();
        mpNodeBasedCellPopulation->AddNodeAndCellToSendRight(index_of_node_to_send);
        mpNodeBasedCellPopulation->AddNodeAndCellToSendLeft(index_of_node_to_send);

        TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mCellsToSendRight.size(), 1u);
        TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mCellsToSendLeft.size(), 1u);

        unsigned node_right_index = (*mpNodeBasedCellPopulation->mCellsToSendRight.begin()).second->GetIndex();
        TS_ASSERT_EQUALS(node_right_index, index_of_node_to_send);

        unsigned node_left_index = (*mpNodeBasedCellPopulation->mCellsToSendLeft.begin()).second->GetIndex();
        TS_ASSERT_EQUALS(node_left_index, index_of_node_to_send);
    }

    void TestSendAndReceiveCells()
    {
        unsigned index_of_node_to_send = mpNodesOnlyMesh->GetNodeIteratorBegin()->GetIndex();;
        mpNodeBasedCellPopulation->AddNodeAndCellToSendRight(index_of_node_to_send);
        mpNodeBasedCellPopulation->AddNodeAndCellToSendLeft(index_of_node_to_send);

        TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mCellCommunicationTag, 123u);

        TS_ASSERT(!(mpNodeBasedCellPopulation->mpCellsRecvRight));
        TS_ASSERT(!(mpNodeBasedCellPopulation->mpCellsRecvLeft));

        mpNodeBasedCellPopulation->SendCellsToNeighbourProcesses();

        if (!PetscTools::AmTopMost())
        {
            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mpCellsRecvRight->size(), 1u);

            unsigned index = (*mpNodeBasedCellPopulation->mpCellsRecvRight->begin()).second->GetIndex();
            TS_ASSERT_EQUALS(index, PetscTools::GetMyRank() + 1);
        }
        if (!PetscTools::AmMaster())
        {
            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mpCellsRecvLeft->size(), 1u);

            unsigned index = (*mpNodeBasedCellPopulation->mpCellsRecvLeft->begin()).second->GetIndex();
            TS_ASSERT_EQUALS(index, PetscTools::GetMyRank() - 1);
        }
    }

    void TestSendAndReceiveCellsNonBlocking()
    {
        unsigned index_of_node_to_send = mpNodesOnlyMesh->GetNodeIteratorBegin()->GetIndex();;
        mpNodeBasedCellPopulation->AddNodeAndCellToSendRight(index_of_node_to_send);
        mpNodeBasedCellPopulation->AddNodeAndCellToSendLeft(index_of_node_to_send);

        TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mCellCommunicationTag, 123u);

        TS_ASSERT(!(mpNodeBasedCellPopulation->mpCellsRecvRight));
        TS_ASSERT(!(mpNodeBasedCellPopulation->mpCellsRecvLeft));

        mpNodeBasedCellPopulation->NonBlockingSendCellsToNeighbourProcesses();

        mpNodeBasedCellPopulation->GetReceivedCells();

        if (!PetscTools::AmTopMost())
        {
            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mpCellsRecvRight->size(), 1u);

            unsigned index = (*mpNodeBasedCellPopulation->mpCellsRecvRight->begin()).second->GetIndex();
            TS_ASSERT_EQUALS(index, PetscTools::GetMyRank() + 1);
        }
        if (!PetscTools::AmMaster())
        {
            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mpCellsRecvLeft->size(), 1u);

            unsigned index = (*mpNodeBasedCellPopulation->mpCellsRecvLeft->begin()).second->GetIndex();
            TS_ASSERT_EQUALS(index, PetscTools::GetMyRank() - 1);
        }
    }

    void TestUpdateCellProcessLocation()
    {
        if (PetscTools::GetNumProcs() > 1)
        {
            if (PetscTools::AmMaster())
            {
                // Move node to the next location.
                c_vector<double, 3> new_location = zero_vector<double>(3);
                new_location[2] = 1.6;
                ChastePoint<3> point(new_location);
                mpNodesOnlyMesh->GetNode(0)->SetPoint(point);
            }
            if (PetscTools::GetMyRank() == 1)
            {
                // Move node to the next location.
                c_vector<double, 3> new_location = zero_vector<double>(3);
                new_location[2] = 0.5;
                ChastePoint<3> point(new_location);
                mpNodesOnlyMesh->GetNode(1)->SetPoint(point);
            }
            mpNodeBasedCellPopulation->UpdateCellProcessLocation();

            if (PetscTools::AmMaster())
            {
                TS_ASSERT_EQUALS(mpNodesOnlyMesh->GetNumNodes(), 1u);
                TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->GetNumRealCells(), 1u);
            }
            if (PetscTools::GetMyRank() == 1)
            {
                TS_ASSERT_EQUALS(mpNodesOnlyMesh->GetNumNodes(), 1u);
                TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->GetNumRealCells(), 1u);

                AbstractMesh<3,3>::NodeIterator node_iter = mpNodesOnlyMesh->GetNodeIteratorBegin();
                TS_ASSERT_DELTA(node_iter->rGetLocation()[2], 1.6, 1e-4);
            }
        }
    }

    void TestRefreshHaloCells()
    {
        // Set up the halo boxes and nodes.
        mpNodeBasedCellPopulation->Update();

        // Send and receive halo nodes.
        mpNodeBasedCellPopulation->RefreshHaloCells();

        mpNodeBasedCellPopulation->AddReceivedHaloCells();

        if (!PetscTools::AmMaster() && !PetscTools::AmTopMost())
        {
           TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mHaloCells.size(), 2u);
           TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mHaloCellLocationMap[mpNodeBasedCellPopulation->mHaloCells[0]], PetscTools::GetMyRank() - 1);
           TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mHaloCellLocationMap[mpNodeBasedCellPopulation->mHaloCells[1]], PetscTools::GetMyRank() + 1);
        }
        else if (!PetscTools::AmMaster() || !PetscTools::AmTopMost())
        {
           TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mHaloCells.size(), 1u);
        }
    }

    void TestUpdateWithLoadBalanceDoesntThrow()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        mpNodeBasedCellPopulation->SetLoadBalanceMesh(true);

        mpNodeBasedCellPopulation->SetLoadBalanceFrequency(50);

        TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mLoadBalanceFrequency, 50u);

        TS_ASSERT_THROWS_NOTHING(mpNodeBasedCellPopulation->Update());
    }

    void TestGetCellUsingLocationIndexWithHaloCell()
    {
        boost::shared_ptr<Node<3> > p_node(new Node<3>(10, false, 0.0, 0.0, 0.0));

        // Create a cell.
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_type);
        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

        CellPtr p_cell(new Cell(p_state, p_model));

        mpNodeBasedCellPopulation->AddHaloCell(p_cell, p_node);

        TS_ASSERT_THROWS_NOTHING(mpNodeBasedCellPopulation->GetCellUsingLocationIndex(10));
        TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->GetCellUsingLocationIndex(10), p_cell);
    }

    void TestNodeBasedCellPopulationOutputInParallel()
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        std::vector<Node<2>* > nodes;
        for (unsigned i=0; i<5; i++)
        {
            nodes.push_back(new Node<2>(i, false, 0.0, 0.75*i));
        }

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(0.0);
        }

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        node_based_cell_population.Update(); // so cell neighbours are calculated

        TS_ASSERT(!node_based_cell_population.HasWriter<CellMutationStatesCountWriter>());
        node_based_cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        TS_ASSERT(node_based_cell_population.HasWriter<CellMutationStatesCountWriter>());

        TS_ASSERT(!node_based_cell_population.HasWriter<CellIdWriter>());
        node_based_cell_population.AddCellWriter<CellIdWriter>();
        TS_ASSERT(node_based_cell_population.HasWriter<CellIdWriter>());

        TS_ASSERT(!node_based_cell_population.HasWriter<CellProliferativeTypesCountWriter>());
        node_based_cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        TS_ASSERT(node_based_cell_population.HasWriter<CellProliferativeTypesCountWriter>());

        TS_ASSERT(!node_based_cell_population.HasWriter<CellProliferativePhasesCountWriter>());
        node_based_cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        TS_ASSERT(node_based_cell_population.HasWriter<CellProliferativePhasesCountWriter>());

        TS_ASSERT(!node_based_cell_population.HasWriter<CellProliferativePhasesWriter>());
        node_based_cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        TS_ASSERT(node_based_cell_population.HasWriter<CellProliferativePhasesWriter>());

        TS_ASSERT(!node_based_cell_population.HasWriter<CellAgesWriter>());
        node_based_cell_population.AddCellWriter<CellAgesWriter>();
        TS_ASSERT(node_based_cell_population.HasWriter<CellAgesWriter>());

        TS_ASSERT(!node_based_cell_population.HasWriter<CellVolumesWriter>());
        node_based_cell_population.AddCellWriter<CellVolumesWriter>();
        TS_ASSERT(node_based_cell_population.HasWriter<CellVolumesWriter>());

        TS_ASSERT(!node_based_cell_population.HasWriter<CellAncestorWriter>());
        node_based_cell_population.AddCellWriter<CellAncestorWriter>();
        TS_ASSERT(node_based_cell_population.HasWriter<CellAncestorWriter>());

        // Test set methods
        std::string output_directory = "TestNodeBasedCellPopulationWritersParallel";
        OutputFileHandler output_file_handler(output_directory, false);

        node_based_cell_population.OpenWritersFiles(output_file_handler);

        // Write out the files here
        node_based_cell_population.WriteResultsToFiles(output_directory);
        node_based_cell_population.CloseWritersFiles();

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestNodeBasedCellPopulationWritersParallel/results.viznodes").CompareFiles();
        FileComparison(results_dir + "results.vizcelltypes", "cell_based/test/data/TestNodeBasedCellPopulationWritersParallel/results.vizcelltypes").CompareFiles();
        FileComparison(results_dir + "results.vizancestors", "cell_based/test/data/TestNodeBasedCellPopulationWritersParallel/results.vizancestors").CompareFiles();
        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestNodeBasedCellPopulationWritersParallel/cellmutationstates.dat").CompareFiles();

        if (PetscTools::IsSequential())
        {
            // Cell ages and volumes file differs because it writes the global index of the node to file, which is different depending on how many processes there are.
            FileComparison(results_dir + "cellages.dat", "cell_based/test/data/TestNodeBasedCellPopulationWritersParallel/cellages.dat").CompareFiles();
            FileComparison(results_dir + "cellareas.dat", "cell_based/test/data/TestNodeBasedCellPopulationWritersParallel/cellareas.dat").CompareFiles();
            FileComparison(results_dir + "cellareas.dat", "cell_based/test/data/TestNodeBasedCellPopulationWritersParallel/cellareas.dat").CompareFiles();
        }

#ifdef CHASTE_VTK
        if (PetscTools::IsParallel())
        {
            // Meta-file links to parallel files (which link to the fragments)
            FileComparison(results_dir + "results.pvd", "cell_based/test/data/TestNodeBasedCellPopulationWritersParallel/results_from_parallel.pvd").CompareFiles();
        }
        else
        {
            FileComparison(results_dir + "results.pvd", "cell_based/test/data/TestNodeBasedCellPopulationWritersParallel/results.pvd").CompareFiles();
        }
#endif //CHASTE_VTK
        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif /*TESTNODEBASEDCELLPOPULATIONPARALLELMETHODS_HPP_*/
