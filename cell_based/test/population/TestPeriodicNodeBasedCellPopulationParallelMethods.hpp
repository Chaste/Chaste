/*

Copyright (c) 2005-2022, University of Oxford.
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

#ifndef TESTPERIODICNODEBASEDCELLPOPULATIONPARALLELMETHODS_HPP_
#define TESTPERIODICNODEBASEDCELLPOPULATIONPARALLELMETHODS_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellLabel.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PeriodicNodesOnlyMesh.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "TrianglesMeshReader.hpp"
#include "WildTypeCellMutationState.hpp"

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

class TestPeriodicNodeBasedCellPopulationParallelMethods : public AbstractCellBasedTestSuite
{
private:
    PeriodicNodesOnlyMesh<3>* mpPeriodicNodesOnlyMesh;

    NodeBasedCellPopulation<3>* mpNodeBasedCellPopulation;

    void setUp()
    {
        AbstractCellBasedTestSuite::setUp();

        if (PetscTools::GetNumProcs() > 2)
        {
            std::vector<Node<3>*> nodes;
            unsigned num_processors = PetscTools::GetNumProcs();
            for (unsigned i = 0; i < num_processors; i++)
            {
                nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.5 + (double)i));
            }

            c_vector<double, 3> periodic_width = zero_vector<double>(3);
            periodic_width[2] = (double)num_processors;

            mpPeriodicNodesOnlyMesh = new PeriodicNodesOnlyMesh<3>(periodic_width);
            mpPeriodicNodesOnlyMesh->ConstructNodesWithoutMesh(nodes, 1.0); // Cutoff length of means we need at least 3 processors

            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
            cells_generator.GenerateBasic(cells, mpPeriodicNodesOnlyMesh->GetNumNodes());

            mpNodeBasedCellPopulation = new NodeBasedCellPopulation<3>(*mpPeriodicNodesOnlyMesh, cells);

            for (unsigned i = 0; i < nodes.size(); i++)
            {
                delete nodes[i];
            }
        }
    }

    void tearDown()
    {
        if (PetscTools::GetNumProcs() > 2)
        {
            delete mpNodeBasedCellPopulation;
            delete mpPeriodicNodesOnlyMesh;
        }
        AbstractCellBasedTestSuite::tearDown();
    }

public:
    void TestSendAndReceiveCells()
    {
        if (PetscTools::GetNumProcs() > 2)
        {
            unsigned index_of_node_to_send = mpPeriodicNodesOnlyMesh->GetNodeIteratorBegin()->GetIndex();
            ;
            mpNodeBasedCellPopulation->AddNodeAndCellToSendRight(index_of_node_to_send);
            mpNodeBasedCellPopulation->AddNodeAndCellToSendLeft(index_of_node_to_send);

            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mCellCommunicationTag, 123u);

            TS_ASSERT(!(mpNodeBasedCellPopulation->mpCellsRecvRight));
            TS_ASSERT(!(mpNodeBasedCellPopulation->mpCellsRecvLeft));

            mpNodeBasedCellPopulation->SendCellsToNeighbourProcesses();

            // Check the nodes transferred correctly
            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mpCellsRecvRight->size(), 1u);
            unsigned index = (*mpNodeBasedCellPopulation->mpCellsRecvRight->begin()).second->GetIndex();
            if (PetscTools::AmTopMost())
            {
                TS_ASSERT_EQUALS(index, 0u);
            }
            else
            {
                TS_ASSERT_EQUALS(index, PetscTools::GetMyRank() + 1u);
            }

            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mpCellsRecvLeft->size(), 1u);
            index = (*mpNodeBasedCellPopulation->mpCellsRecvLeft->begin()).second->GetIndex();
            if (PetscTools::AmMaster())
            {
                TS_ASSERT_EQUALS(index, PetscTools::GetNumProcs() - 1u);
            }
            else
            {
                TS_ASSERT_EQUALS(index, PetscTools::GetMyRank() - 1u);
            }
        }
    }

    void TestSendAndReceiveCellsNonBlocking()
    {
        if (PetscTools::GetNumProcs() > 2)
        {
            unsigned index_of_node_to_send = mpPeriodicNodesOnlyMesh->GetNodeIteratorBegin()->GetIndex();
            ;
            mpNodeBasedCellPopulation->AddNodeAndCellToSendRight(index_of_node_to_send);
            mpNodeBasedCellPopulation->AddNodeAndCellToSendLeft(index_of_node_to_send);

            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mCellCommunicationTag, 123u);

            TS_ASSERT(!(mpNodeBasedCellPopulation->mpCellsRecvRight));
            TS_ASSERT(!(mpNodeBasedCellPopulation->mpCellsRecvLeft));

            mpNodeBasedCellPopulation->NonBlockingSendCellsToNeighbourProcesses();

            mpNodeBasedCellPopulation->GetReceivedCells();

            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mpCellsRecvRight->size(), 1u);
            unsigned index = (*mpNodeBasedCellPopulation->mpCellsRecvRight->begin()).second->GetIndex();
            if (PetscTools::AmTopMost())
            {
                TS_ASSERT_EQUALS(index, 0u);
            }
            else
            {
                TS_ASSERT_EQUALS(index, PetscTools::GetMyRank() + 1u);
            }

            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mpCellsRecvLeft->size(), 1u);
            index = (*mpNodeBasedCellPopulation->mpCellsRecvLeft->begin()).second->GetIndex();
            if (PetscTools::AmMaster())
            {
                TS_ASSERT_EQUALS(index, PetscTools::GetNumProcs() - 1u);
            }
            else
            {
                TS_ASSERT_EQUALS(index, PetscTools::GetMyRank() - 1u);
            }
        }
    }

    void TestUpdateCellProcessLocation()
    {
        if (PetscTools::GetNumProcs() > 2)
        {
            if (PetscTools::AmMaster())
            {
                // Move node to the next location.
                c_vector<double, 3> new_location = zero_vector<double>(3);
                new_location[2] = (double)PetscTools::GetNumProcs() - 0.4;
                ChastePoint<3> point(new_location);
                mpPeriodicNodesOnlyMesh->GetNode(0)->SetPoint(point);
            }
            if (PetscTools::GetMyRank() == (PetscTools::GetNumProcs() - 1))
            {
                // Move node to the next location.
                c_vector<double, 3> new_location = zero_vector<double>(3);
                new_location[2] = 0.5;
                ChastePoint<3> point(new_location);
                mpPeriodicNodesOnlyMesh->GetNode(PetscTools::GetNumProcs() - 1)->SetPoint(point);
            }
            mpNodeBasedCellPopulation->UpdateCellProcessLocation();

            if (PetscTools::AmMaster())
            {
                TS_ASSERT_EQUALS(mpPeriodicNodesOnlyMesh->GetNumNodes(), 1u);
                TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->GetNumRealCells(), 1u);
            }
            if (PetscTools::GetMyRank() == (PetscTools::GetNumProcs() - 1))
            {
                TS_ASSERT_EQUALS(mpPeriodicNodesOnlyMesh->GetNumNodes(), 1u);
                TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->GetNumRealCells(), 1u);

                AbstractMesh<3, 3>::NodeIterator node_iter = mpPeriodicNodesOnlyMesh->GetNodeIteratorBegin();
                TS_ASSERT_DELTA(node_iter->rGetLocation()[2], (double)PetscTools::GetNumProcs() - 0.4, 1e-4);
            }
        }
    }

    void TestRefreshHaloCells()
    {
        if (PetscTools::GetNumProcs() > 2)
        {
            // Set up the halo boxes and nodes.
            mpNodeBasedCellPopulation->Update();

            // Send and receive halo nodes.
            mpNodeBasedCellPopulation->RefreshHaloCells();

            mpNodeBasedCellPopulation->AddReceivedHaloCells();

            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mHaloCells.size(), 2u);
            if (PetscTools::AmMaster())
            {
                TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mHaloCellLocationMap[mpNodeBasedCellPopulation->mHaloCells[0]], PetscTools::GetNumProcs() - 1u);
            }
            else
            {
                TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mHaloCellLocationMap[mpNodeBasedCellPopulation->mHaloCells[0]], PetscTools::GetMyRank() - 1u);
            }
            if (PetscTools::AmTopMost())
            {
                TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mHaloCellLocationMap[mpNodeBasedCellPopulation->mHaloCells[1]], 0u);
            }
            else
            {
                TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mHaloCellLocationMap[mpNodeBasedCellPopulation->mHaloCells[1]], PetscTools::GetMyRank() + 1u);
            }
        }
    }
};

#endif /*TESTPERIODICNODEBASEDCELLPOPULATIONPARALLELMETHODS_HPP_*/
