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

#ifndef TESTMUTABLEVERTEXEDGES_HPP_
#define TESTMUTABLEVERTEXEDGES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "VertexMeshReader.hpp"
#include "VertexMeshWriter.hpp"
#include "MutableVertexMesh.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "ArchiveOpener.hpp"
#include "Edge.hpp"
#include <HoneycombVertexMeshGenerator.hpp>

class TestMutableVertexEdges : public CxxTest::TestSuite
{
public:

    void TestEdgesCreation()
    {
        const unsigned ELEMENT_DIM = 2;
        const unsigned SPACE_DIM = 2;

        unsigned exampleEdgeIndex = 1;

        // Empty edge should be invalid
        Edge<2>* testEdge = new Edge<2>(exampleEdgeIndex);
        TS_ASSERT(!testEdge->IsEdgeValid());

        TS_ASSERT_EQUALS(exampleEdgeIndex, testEdge->GetIndex());

        // Add nodes to the edge, check validity and centre position
        Node<SPACE_DIM> node0(0, ChastePoint<SPACE_DIM>(0,0,0));
        Node<SPACE_DIM> node1(1, ChastePoint<SPACE_DIM>(1,1,1));

        testEdge->SetNodes(&node0, &node1);
        auto edgeCentreLocation = testEdge->rGetCentreLocation();
        TS_ASSERT(testEdge->IsEdgeValid());
        TS_ASSERT(edgeCentreLocation[0] == 0.5 && edgeCentreLocation[1] == 0.5);

        // Generate two hex vertex elements
        std::vector<Node<SPACE_DIM>*> nodes0, nodes1, allnodes;
        nodes0.push_back(new Node<SPACE_DIM>(0));
        nodes0.push_back(new Node<SPACE_DIM>(1));
        nodes0.push_back(new Node<SPACE_DIM>(2));
        nodes0.push_back(new Node<SPACE_DIM>(3));
        nodes0.push_back(new Node<SPACE_DIM>(4));
        nodes0.push_back(new Node<SPACE_DIM>(5));

        nodes1.push_back(nodes0[0]);
        nodes1.push_back(nodes0[1]);
        nodes1.push_back(new Node<SPACE_DIM>(6));
        nodes1.push_back(new Node<SPACE_DIM>(7));
        nodes1.push_back(new Node<SPACE_DIM>(8));
        nodes1.push_back(new Node<SPACE_DIM>(9));

        allnodes.insert(allnodes.end(), nodes0.begin(), nodes0.end());
        allnodes.insert(allnodes.end(), nodes1.begin(), nodes1.end());

        std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> elements;
        elements.push_back(new VertexElement<ELEMENT_DIM,SPACE_DIM>(0, nodes0));
        elements.push_back(new VertexElement<ELEMENT_DIM,SPACE_DIM>(1, nodes1));

        // Generate a mesh which will automatically build the edges in the constructor
        VertexMesh<ELEMENT_DIM, SPACE_DIM>* mesh = new VertexMesh<ELEMENT_DIM, SPACE_DIM>(allnodes, elements);
        for (unsigned i = 0; i < mesh->GetNumElements(); i++)
        {
            auto element = mesh->GetElement(i);
            TS_ASSERT(element->CheckEdgesAreValid());
        }

        // Also test with honeycomb mesh (MutableVertexMesh)
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* honeycombMesh = generator.GetMesh();
        TS_ASSERT_EQUALS(honeycombMesh->GetNumEdges(), 19);
        for (unsigned i = 0; i < honeycombMesh->GetNumElements(); i++)
        {
            auto element = honeycombMesh->GetElement(i);
            TS_ASSERT(element->CheckEdgesAreValid());
        }

    }

    void TestT1Swap()
    {
        // Also test with honeycomb mesh (MutableVertexMesh)
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* honeycombMesh = generator.GetMesh();
        TS_ASSERT_EQUALS(honeycombMesh->GetNumEdges(), 19);
        for (unsigned i = 0; i < honeycombMesh->GetNumElements(); i++)
        {
            auto element = honeycombMesh->GetElement(i);
            TS_ASSERT(element->CheckEdgesAreValid());
        }

        // Perform a test T1 swap on a shared edge
        int numSwapsPerformed = 0;
        for (unsigned i = 0 ;i < honeycombMesh->GetNumEdges(); i++)
        {
            auto edge = honeycombMesh->GetEdge(i);
            if (edge->GetNumElements() > 1)
            {
                honeycombMesh->IdentifySwapType(edge->GetNode(0), edge->GetNode(1));
                honeycombMesh->RemoveDeletedNodes();
                numSwapsPerformed++;
            }
        }

        // Check edges again
        for (unsigned i = 0; i < honeycombMesh->GetNumElements(); i++)
        {
            auto element = honeycombMesh->GetElement(i);
            TS_ASSERT(element->CheckEdgesAreValid());
        }

        // Count edge operations
        printf("Num swaps performd: %i with num edge operations: %lu \n",numSwapsPerformed, honeycombMesh->GetEdgeOperations().size());

        TS_ASSERT(honeycombMesh->GetEdgeOperations().size() == 8);
        for (auto edgeOperation: honeycombMesh->GetEdgeOperations())
        {
            TS_ASSERT(edgeOperation->GetOperation() != EDGE_OPERATION_DIVIDE);
        }
    }

    void TestMeshElementDivide()
    {
        // Test cell division with honeycomb mesh
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* honeycombMesh = generator.GetMesh();
        std::vector<unsigned> edgeIdBeforeDivide;

        auto element0 = honeycombMesh->GetElement(0);
        printf("Node ids: ");
        for (unsigned i = 0; i < element0->GetNumNodes(); i++)
        {
            printf(" %i ", element0->GetNode(i)->GetIndex());
        }
        printf("\n");
        printf("Edge ids: ");
        for (unsigned i = 0; i < element0->GetNumEdges(); i++)
        {
            printf(" %i ", element0->GetEdge(i)->GetIndex());
            edgeIdBeforeDivide.push_back(element0->GetEdge(i)->GetIndex());
        }
        printf("\n");

        auto newElemIndex = honeycombMesh->DivideElementAlongShortAxis(honeycombMesh->GetElement(0));
        auto newElem = honeycombMesh->GetElement(newElemIndex);

        printf("Redivide node ids: ");
        for (unsigned i = 0; i < element0->GetNumNodes(); i++)
        {
            printf(" %i ", element0->GetNode(i)->GetIndex());
        }
        printf("\n");
        printf("Redivide edge ids: ");
        for (unsigned i = 0; i < element0->GetNumEdges(); i++)
        {
            printf(" %i ", element0->GetEdge(i)->GetIndex());
        }
        printf("\n");

        printf("New Element Node ids: ");
        for (unsigned i = 0; i < newElem->GetNumNodes(); i++)
        {
            printf(" %i ", newElem->GetNode(i)->GetIndex());
        }
        printf("\n");
        printf("New element edge ids: ");
        for (unsigned i = 0; i < newElem->GetNumEdges(); i++)
        {
            printf(" %i ", newElem->GetEdge(i)->GetIndex());
        }
        printf("\n");

        // Check edges again
        for ( unsigned i = 0; i < honeycombMesh->GetNumElements(); i++)
        {
            auto element = honeycombMesh->GetElement(i);
            TS_ASSERT(element->CheckEdgesAreValid());
        }

        printf("Num edge operations: %lu \n", honeycombMesh->GetEdgeOperations().size());
        auto edgeOperations = honeycombMesh->GetEdgeOperations();

        // Only one divide operation
        TS_ASSERT(edgeOperations.size() == 1);
        TS_ASSERT(edgeOperations[0]->GetOperation() == EDGE_OPERATION_DIVIDE);

        // Test edge remapping data for the old cell
        {
            auto edgeRemap = edgeOperations[0]->GetNewEdges();
            auto elem = element0;

            TS_ASSERT(edgeRemap->GetEdgesMapping().size() == elem->GetNumEdges());
            TS_ASSERT(edgeRemap->GetEdgesStatus().size() == elem->GetNumEdges());

            int numRemapEdges = 0;
            int numSplitEdges = 0;
            int numNewEdges = 0;
            for (unsigned i = 0 ; i < edgeRemap->GetEdgesStatus().size(); i++)
            {
                switch (edgeRemap->GetEdgesStatus()[i])
                {
                    case 0:
                    {
                        auto remapIndex = edgeRemap->GetEdgesMapping()[i];
                        auto oldEdgeId = edgeIdBeforeDivide[remapIndex];
                        auto newEdgeId = elem->GetEdge(i)->GetIndex();
                        TS_ASSERT(oldEdgeId == newEdgeId);
                        numRemapEdges++;
                        break;
                    }
                    case 1:
                    {
                        TS_ASSERT(edgeRemap->GetEdgesMapping()[i] > -1);
                        numSplitEdges++;
                        break;
                    }
                    case 2:
                    {
                        TS_ASSERT(edgeRemap->GetEdgesMapping()[i] < 0);
                        numNewEdges++;
                        break;
                    }
                    default:
                    {
                        TS_FAIL("Edge status invalid");
                        break;
                    }
                }
            }

            TS_ASSERT_EQUALS(numRemapEdges, 2u);
            TS_ASSERT_EQUALS(numSplitEdges, 2u);
            TS_ASSERT_EQUALS(numNewEdges, 1u);
        }

        // Test edge remapping data for the new cell
        {
            auto edgeRemap = edgeOperations[0]->GetNewEdges2();
            auto elem = newElem;

            TS_ASSERT(edgeRemap->GetEdgesMapping().size() == elem->GetNumEdges());
            TS_ASSERT(edgeRemap->GetEdgesStatus().size() == elem->GetNumEdges());

            int numRemapEdges = 0;
            int numSplitEdges = 0;
            int numNewEdges = 0;
            for (unsigned i = 0; i < edgeRemap->GetEdgesStatus().size(); i++)
            {
                switch (edgeRemap->GetEdgesStatus()[i])
                {
                    case 0:
                    {
                        auto remapIndex = edgeRemap->GetEdgesMapping()[i];
                        auto oldEdgeId = edgeIdBeforeDivide[remapIndex];
                        auto newEdgeId = elem->GetEdge(i)->GetIndex();
                        TS_ASSERT(oldEdgeId == newEdgeId);
                        numRemapEdges++;
                        break;
                    }
                    case 1:
                    {
                        TS_ASSERT(edgeRemap->GetEdgesMapping()[i] > -1);
                        numSplitEdges++;
                        break;
                    }
                    case 2:
                    {
                        TS_ASSERT(edgeRemap->GetEdgesMapping()[i] < 0);
                        numNewEdges++;
                        break;
                    }
                    default:
                    {
                        TS_FAIL("Edge status invalid");
                        break;
                    }
                }
            }

            TS_ASSERT_EQUALS(numRemapEdges, 2u);
            TS_ASSERT_EQUALS(numSplitEdges, 2u);
            TS_ASSERT_EQUALS(numNewEdges, 1u);
        }
    }
};

#endif //TESTMUTABLEVERTEXEDGES_HPP_
