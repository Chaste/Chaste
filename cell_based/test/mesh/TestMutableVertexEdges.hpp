//
// Created by twin on 25/01/19.
//

#ifndef TESTMUTABLEVERTEXEDGES_HPP_
#define TESTMUTABLEVERTEXEDGES_HPP_

#include "Edge.hpp"
#include <HoneycombVertexMeshGenerator.hpp>
#include "AbstractCellBasedTestSuite.hpp"

class TestMutableVertexEdges : public AbstractCellBasedTestSuite
{
public:

    void TestEdges()
    {
        const unsigned SPACE_DIM = 2;

        //Generate a hexagonal mesh
        std::vector<Node<SPACE_DIM>*> nodes;
        for(int i = 0 ;i < 6; i++){
            nodes.push_back(new Node<SPACE_DIM>(i));
        }

        //Generate the edges for the nodes
        for(int i = 0 ;i < 6; i++){

            unsigned i_next = (i + 1) % 6;

            //Create an edge
            Edge<2>* edge = new Edge<2>(1);
            edge->SetNodes(nodes[i], nodes[i_next]);

        }




        //Create a mesh, creating edges
//        HoneycombVertexMeshGenerator generator(5, 5);
//        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();


    }

};

#endif //TESTMUTABLEVERTEXEDGES_HPP_
