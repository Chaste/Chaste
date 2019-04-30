//
// Created by twin on 25/01/19.
//

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
        for( unsigned i = 0; i < mesh->GetNumElements(); i++)
        {
            auto element = mesh->GetElement(i);
            TS_ASSERT(element->CheckEdgesAreValid());
        }

        // Also test with honeycomb mesh (MutableVertexMesh)
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* honeycombMesh = generator.GetMesh();
        TS_ASSERT_EQUALS(honeycombMesh->GetNumEdges(), 19);
        for( unsigned i = 0; i < honeycombMesh->GetNumElements(); i++)
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
        for( unsigned i = 0; i < honeycombMesh->GetNumElements(); i++)
        {
            auto element = honeycombMesh->GetElement(i);
            TS_ASSERT(element->CheckEdgesAreValid());
        }


        // Perform a test T1 swap on a shared edge
        for(unsigned i = 0 ;i < honeycombMesh->GetNumEdges(); i++)
        {
            auto edge = honeycombMesh->GetEdge(i);
            if(edge->GetNumElements() > 1)
            {
                honeycombMesh->IdentifySwapType(edge->GetNode(0), edge->GetNode(1));
                honeycombMesh->RemoveDeletedNodes();
            }

        }

        // Check edges again
        for( unsigned i = 0; i < honeycombMesh->GetNumElements(); i++)
        {
            auto element = honeycombMesh->GetElement(i);
            TS_ASSERT(element->CheckEdgesAreValid());

        }


    }

//    void TestMeshElementDivide()
//    {
//        // Test cell division with honeycomb mesh
//        HoneycombVertexMeshGenerator generator(2, 2);
//        MutableVertexMesh<2,2>* honeycombMesh = generator.GetMesh();
//
//        auto element0 = honeycombMesh->GetElement(0);
//        printf("Node ids: ");
//        for(unsigned i = 0; i < element0->GetNumNodes(); i++)
//        {
//            printf(" %i ", element0->GetNode(i)->GetIndex());
//        }
//        printf("\n");
//        printf("Edge ids: ");
//        for(unsigned i = 0; i < element0->GetNumEdges(); i++)
//        {
//            printf(" %i ", element0->GetEdge(i)->GetIndex());
//        }
//        printf("\n");
//
//
//
//        auto newElemIndex = honeycombMesh->DivideElementAlongShortAxis(honeycombMesh->GetElement(0));
//        auto newElem = honeycombMesh->GetElement(newElemIndex);
//
//        printf("Redivide node ids: ");
//        for(unsigned i = 0; i < element0->GetNumNodes(); i++)
//        {
//            printf(" %i ", element0->GetNode(i)->GetIndex());
//        }
//        printf("\n");
//        printf("Redivide edge ids: ");
//        for(unsigned i = 0; i < element0->GetNumEdges(); i++)
//        {
//            printf(" %i ", element0->GetEdge(i)->GetIndex());
//        }
//        printf("\n");
//
//        printf("New Element Node ids: ");
//        for(unsigned i = 0; i < newElem->GetNumNodes(); i++)
//        {
//            printf(" %i ", newElem->GetNode(i)->GetIndex());
//        }
//        printf("\n");
//        printf("New element edge ids: ");
//        for(unsigned i = 0; i < newElem->GetNumEdges(); i++)
//        {
//            printf(" %i ", newElem->GetEdge(i)->GetIndex());
//        }
//        printf("\n");
//
//        // Check edges again
//        for( unsigned i = 0; i < honeycombMesh->GetNumElements(); i++)
//        {
//            auto element = honeycombMesh->GetElement(i);
//            TS_ASSERT(element->CheckEdgesAreValid());
//
//        }
//
//    }

    void TestMeshElementDivideShort()
    {
        const unsigned ELEMENT_DIM = 2;
        const unsigned SPACE_DIM = 2;


        // Generate two hex vertex elements
        std::vector<Node<SPACE_DIM>*> nodes0;
        nodes0.push_back(new Node<SPACE_DIM>(0, ChastePoint<2>(-0.5, 1)));
        nodes0.push_back(new Node<SPACE_DIM>(1, ChastePoint<2>(-0.5, 0)));
        nodes0.push_back(new Node<SPACE_DIM>(2, ChastePoint<2>(-0.5, -1)));
        nodes0.push_back(new Node<SPACE_DIM>(3, ChastePoint<2>(0.5, -1)));
        nodes0.push_back(new Node<SPACE_DIM>(4, ChastePoint<2>(0.5, 0)));
        nodes0.push_back(new Node<SPACE_DIM>(5, ChastePoint<2>(0.5, 1)));

        std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> elements;
        elements.push_back(new VertexElement<ELEMENT_DIM,SPACE_DIM>(0, nodes0));

        MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>* mesh = new MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>(nodes0, elements);

        mesh->DivideElementAlongShortAxis(mesh->GetElement(0));
    }

};

#endif //TESTMUTABLEVERTEXEDGES_HPP_
