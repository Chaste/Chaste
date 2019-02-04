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

        //Empty edge should be invalid
        Edge<2>* emptyEdge = new Edge<2>(exampleEdgeIndex);
        TS_ASSERT(!emptyEdge->IsEdgeValid());

        TS_ASSERT_EQUALS(exampleEdgeIndex, emptyEdge->GetIndex());

        //Generate two hex vertex elements
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

        //Generate a mesh which will automatically build the edges in the constructor
        VertexMesh<ELEMENT_DIM, SPACE_DIM>* mesh = new VertexMesh<ELEMENT_DIM, SPACE_DIM>(allnodes, elements);

        TS_ASSERT_EQUALS(mesh->GetNumEdges(), 11);

        for( unsigned i = 0 ; i < mesh->GetNumEdges(); i++)
        {
            TS_ASSERT(mesh->GetEdge(i)->IsEdgeValid());
        }





        //Create a mesh, creating edges
//        HoneycombVertexMeshGenerator generator(5, 5);
//        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();


    }

};

#endif //TESTMUTABLEVERTEXEDGES_HPP_
