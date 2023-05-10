/*

Copyright (c) 2005-2023, University of Oxford.
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


#ifndef _TESTEDGE_HPP_
#define _TESTEDGE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "MutableVertexMesh.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "ArchiveOpener.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestEdge : public CxxTest::TestSuite
{
public:

    /**
     * Check whether edge building works as intended
     */
    void TestEdgeInitialisation()
    {
        const unsigned ELEMENT_DIM = 2;
        const unsigned SPACE_DIM = 2;

        unsigned exampleEdgeIndex = 1;

        Edge<SPACE_DIM>* testEdge = new Edge<SPACE_DIM>(exampleEdgeIndex);
        TS_ASSERT(testEdge->GetNumNodes()==0);

        TS_ASSERT_EQUALS(exampleEdgeIndex, testEdge->GetIndex());

        // Add nodes to the edge, check validity and centre position
        Node<SPACE_DIM> node0(0, ChastePoint<SPACE_DIM>(0,0,0));
        Node<SPACE_DIM> node1(1, ChastePoint<SPACE_DIM>(1,1,1));

        testEdge->SetNodes(&node0, &node1);
        auto edgeCentreLocation = testEdge->rGetCentreLocation();
        TS_ASSERT(testEdge->GetNumNodes()==2);
        TS_ASSERT(testEdge->GetNode(0)==&node0 ||testEdge->GetNode(0) == &node1);
        TS_ASSERT(testEdge->GetNode(1)==&node0 ||testEdge->GetNode(1) == &node1);
        TS_ASSERT(testEdge->GetNode(0)!=testEdge->GetNode(1));
        TS_ASSERT(edgeCentreLocation[0] == 0.5 && edgeCentreLocation[1] == 0.5);

        // Generate two square vertex elements
        std::vector<Node<SPACE_DIM>*> nodes0, nodes1, allnodes;
        nodes0.push_back(new Node<SPACE_DIM>(0, false, 0.0, 0.0));
        nodes0.push_back(new Node<SPACE_DIM>(1, false, 1.0, 0.0));
        nodes0.push_back(new Node<SPACE_DIM>(2, false, 1.0, 1.0));
        nodes0.push_back(new Node<SPACE_DIM>(3, false, 0.0, 1.0));

        nodes1.push_back(nodes0[1]);
        nodes1.push_back(nodes0[0]);
        nodes1.push_back(new Node<SPACE_DIM>(4, false, 0.0, -1.0));
        nodes1.push_back(new Node<SPACE_DIM>(5, false, 1.0, -1.0));

        allnodes.insert(allnodes.end(), nodes0.begin(), nodes0.end());
        allnodes.insert(allnodes.end(), nodes1.begin(), nodes1.end());

        std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> elements;
        elements.push_back(new VertexElement<ELEMENT_DIM,SPACE_DIM>(0, nodes0));
        elements.push_back(new VertexElement<ELEMENT_DIM,SPACE_DIM>(1, nodes1));

        // Generate a mesh which will automatically build the edges in the constructor
        VertexMesh<ELEMENT_DIM, SPACE_DIM>* mesh = new VertexMesh<ELEMENT_DIM, SPACE_DIM>(allnodes, elements);
        EdgeHelper<SPACE_DIM> edge_helper = mesh->GetEdgeHelper();
        //There are two elements in our mesh
        //We test Edge class methods here
        for (unsigned int i=0; i<2; i++)
        {
            VertexElement<ELEMENT_DIM, SPACE_DIM>* element = elements[i];
            const unsigned int n_edges = element->GetNumEdges();
            for (unsigned int index = 0; index<n_edges; index++)
            {
                Edge<SPACE_DIM>* p_edge = element->GetEdge(index);
                TS_ASSERT(p_edge->GetNumNodes()==2);
                TS_ASSERT(p_edge->GetNode(0)!=p_edge->GetNode(1));
                TS_ASSERT(p_edge->GetNumElements()>0&&p_edge->GetNumElements()<=2);
                TS_ASSERT(element->ContainsEdge(p_edge));
            }
        }
        for (unsigned int i=0; i<mesh->GetNumEdges(); i++)
        {
            Edge<SPACE_DIM>* p_edge = mesh->GetEdge(i);
            TS_ASSERT(edge_helper.GetEdge(i)==p_edge);
            TS_ASSERT(elements[0]->ContainsEdge(p_edge)||elements[1]->ContainsEdge(p_edge));
        }

        //For coverage
        {
            const VertexMesh<ELEMENT_DIM, SPACE_DIM>* mesh_const = new VertexMesh<ELEMENT_DIM, SPACE_DIM>(allnodes, elements);
            const EdgeHelper<SPACE_DIM> edge_helper_for_const = mesh_const->GetEdgeHelper();
            for (unsigned int i=0; i<mesh_const->GetNumEdges(); i++)
            {
                Edge<SPACE_DIM>* p_edge = mesh_const->GetEdge(i);
                TS_ASSERT(edge_helper_for_const.GetEdge(i)==p_edge);
                TS_ASSERT(edge_helper_for_const[i]==p_edge);
                TS_ASSERT(elements[0]->ContainsEdge(p_edge)||elements[1]->ContainsEdge(p_edge));
            }
        }

        // Also test constructors in honeycomb mesh (MutableVertexMesh)
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* honeycombMesh = generator.GetMesh();
        TS_ASSERT_EQUALS(honeycombMesh->GetNumEdges(), 19);
    }

    /**
     * Check whether edge building works as intended
     */
    void TestEdgeInitialisation1d()
    {
        const unsigned SPACE_DIM = 2;

        unsigned exampleEdgeIndex = 1;

        Edge<SPACE_DIM>* testEdge = new Edge<SPACE_DIM>(exampleEdgeIndex);
        TS_ASSERT(testEdge->GetNumNodes()==0);

        TS_ASSERT_EQUALS(exampleEdgeIndex, testEdge->GetIndex());

        // Add nodes to the edge, check validity and centre position
        Node<SPACE_DIM> node0(0, ChastePoint<SPACE_DIM>(0));
        Node<SPACE_DIM> node1(1, ChastePoint<SPACE_DIM>(1));

        testEdge->SetNodes(&node0, &node1);
        auto edgeCentreLocation = testEdge->rGetCentreLocation();
        TS_ASSERT(testEdge->GetNumNodes()==2);
        TS_ASSERT(testEdge->GetNode(0)==&node0 ||testEdge->GetNode(0) == &node1);
        TS_ASSERT(testEdge->GetNode(1)==&node0 ||testEdge->GetNode(1) == &node1);
        TS_ASSERT(testEdge->GetNode(0)!=testEdge->GetNode(1));

        // Generate two line vertex elements
        std::vector<Node<SPACE_DIM>*> nodes0, nodes1, allnodes;
        nodes0.push_back(new Node<SPACE_DIM>(0, false, 0.0));
        nodes0.push_back(new Node<SPACE_DIM>(1, false, 1.0));

        nodes1.push_back(nodes0[0]);
        nodes1.push_back(new Node<SPACE_DIM>(2, false, -1.0));

        allnodes.insert(allnodes.end(), nodes0.begin(), nodes0.end());
        allnodes.insert(allnodes.end(), nodes1.begin(), nodes1.end());

        std::vector<VertexElement<1,SPACE_DIM>*> elements;
        elements.push_back(new VertexElement<1,SPACE_DIM>(0, nodes0));
        elements.push_back(new VertexElement<1,SPACE_DIM>(1, nodes1));

        // Generate a mesh which will automatically build the edges in the constructor
        VertexMesh<1, SPACE_DIM>* mesh = new VertexMesh<1, SPACE_DIM>(allnodes, elements);
        EdgeHelper<SPACE_DIM> edge_helper = mesh->GetEdgeHelper();
        //There are two elements in our mesh
        //We test Edge class methods here
        for (unsigned int i=0; i<2; i++)
        {
            VertexElement<1, SPACE_DIM>* element = elements[i];
            //MutableElement<1, SPACE_DIM>* element2 = elements[i];
            //VertexElement<1, SPACE_DIM>* element = static_cast<VertexElement<1, SPACE_DIM>*>(element2);
            const unsigned int n_edges = element->GetNumEdges();
            for (unsigned int index = 0; index<n_edges; index++)
            {
                Edge<SPACE_DIM>* p_edge = element->GetEdge(index);
                TS_ASSERT(p_edge->GetNumNodes()==2);
                TS_ASSERT(p_edge->GetNode(0)!=p_edge->GetNode(1));
                TS_ASSERT(p_edge->GetNumElements()>0&&p_edge->GetNumElements()<=2);
                TS_ASSERT(element->GetLocalEdgeIndex(p_edge) > 0 && element->GetLocalEdgeIndex(p_edge) < 2);
                TS_ASSERT(element->GetEdgeGlobalIndex(index) <2);
                TS_ASSERT(element->ContainsEdge(p_edge));
            }
        }
        for (unsigned int i=0; i<mesh->GetNumEdges(); i++)
        {
            Edge<SPACE_DIM>* p_edge = mesh->GetEdge(i);
            TS_ASSERT(edge_helper.GetEdge(i)==p_edge);
            TS_ASSERT(elements[0]->ContainsEdge(p_edge)||elements[1]->ContainsEdge(p_edge));
        }

        //For coverage
        {
            const VertexMesh<1, SPACE_DIM>* mesh_const = new VertexMesh<1, SPACE_DIM>(allnodes, elements);
            const EdgeHelper<SPACE_DIM> edge_helper_for_const = mesh_const->GetEdgeHelper();
            for (unsigned int i=0; i<mesh_const->GetNumEdges(); i++)
            {
                Edge<SPACE_DIM>* p_edge = mesh_const->GetEdge(i);
                TS_ASSERT(edge_helper_for_const.GetEdge(i)==p_edge);
                TS_ASSERT(edge_helper_for_const[i]==p_edge);
                TS_ASSERT(elements[0]->ContainsEdge(p_edge)||elements[1]->ContainsEdge(p_edge));
            }
        }

    }

    void TestEdgeProperties()
    {
        const unsigned ELEMENT_DIM = 2;
        const unsigned SPACE_DIM = 2;

        unsigned exampleEdgeIndex = 1;

        Edge<SPACE_DIM>* testEdge = new Edge<SPACE_DIM>(exampleEdgeIndex);
        TS_ASSERT(testEdge->GetNumNodes()==0);

        TS_ASSERT_EQUALS(exampleEdgeIndex, testEdge->GetIndex());

        // Add nodes to the edge, check validity and centre position
        Node<SPACE_DIM> node0(0, ChastePoint<SPACE_DIM>(0,0,0));
        Node<SPACE_DIM> node1(1, ChastePoint<SPACE_DIM>(1,1,1));

       testEdge->SetNodes(&node1, &node0);
       auto edge_node_indices = testEdge->GetMapIndex();
       TS_ASSERT_EQUALS(edge_node_indices.first, 0u);
       TS_ASSERT_EQUALS(edge_node_indices.second, 1u);

       node0.AddEdge(exampleEdgeIndex);
       node1.AddEdge(exampleEdgeIndex);
       testEdge->MarkAsDeleted();
       TS_ASSERT(testEdge->IsDeleted());

       TS_ASSERT_THROWS_THIS(node0.RemoveEdge(0u), "Tried to remove an edge index which was not in the set");

       std::set<unsigned int> edge_indices;
       edge_indices.insert(0);
       edge_indices.insert(1);
       node0.SetEdgeIndices(edge_indices);
       TS_ASSERT_EQUALS(node0.GetEdgeIndices(), edge_indices);

       // Generate two square vertex elements
       std::vector<Node<SPACE_DIM>*> nodes0, nodes1, allnodes;
       nodes0.push_back(new Node<SPACE_DIM>(0, false, 0.0, 0.0));
       nodes0.push_back(new Node<SPACE_DIM>(1, false, 1.0, 0.0));
       nodes0.push_back(new Node<SPACE_DIM>(2, false, 1.0, 1.0));
       nodes0.push_back(new Node<SPACE_DIM>(3, false, 0.0, 1.0));

       nodes1.push_back(nodes0[1]);
       nodes1.push_back(nodes0[0]);
       nodes1.push_back(new Node<SPACE_DIM>(4, false, 0.0, -1.0));
       nodes1.push_back(new Node<SPACE_DIM>(5, false, 1.0, -1.0));

        allnodes.insert(allnodes.end(), nodes0.begin(), nodes0.end());
        allnodes.insert(allnodes.end(), nodes1.begin(), nodes1.end());

        std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> elements;
        elements.push_back(new VertexElement<ELEMENT_DIM,SPACE_DIM>(0, nodes0));
        elements.push_back(new VertexElement<ELEMENT_DIM,SPACE_DIM>(1, nodes1));

        // Generate a mesh which will automatically build the edges in the constructor
        VertexMesh<ELEMENT_DIM, SPACE_DIM>* mesh = new VertexMesh<ELEMENT_DIM, SPACE_DIM>(allnodes, elements);
        Edge<SPACE_DIM>* edge0 = mesh->GetElement(0)->GetEdge(0);
        Edge<SPACE_DIM>* edge1 = mesh->GetElement(0)->GetEdge(1);
        TS_ASSERT(!edge0->IsBoundaryEdge());
        TS_ASSERT(edge1->IsBoundaryEdge());

        std::set<unsigned> elements_contain_edge = edge0->GetNeighbouringElementIndices();
        std::set<unsigned> expected_elements;
        expected_elements.insert(0);
        expected_elements.insert(1);
        TS_ASSERT_EQUALS(elements_contain_edge, expected_elements);
    }

    void TestArchiveEdge()
    {
        OutputFileHandler handler("TestEdge", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "edge.arch";

        {
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            Edge<2>* testEdge = new Edge<2>(0);
            testEdge->AddElement(0);
            testEdge->AddElement(1);

            Node<2> node0(0, ChastePoint<2>(0,0.0,0.0));
            Node<2> node1(1, ChastePoint<2>(1,1.0,1.0));

            testEdge->SetNodes(&node0, &node1);
            // Write the edge to file
            output_arch << testEdge;
            delete testEdge;
        }

        {
            // Restore the nodes
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            Edge<2>* testEdge;
            input_arch >> testEdge;

            Node<2>* node0 = testEdge->GetNode(0);
            Node<2>* node1 = testEdge->GetNode(1);

            TS_ASSERT_EQUALS(node0->GetIndex(), 0u);
            TS_ASSERT_EQUALS(node1->GetIndex(), 1u);

            TS_ASSERT_EQUALS(node0->rGetLocation()[0], 0.0);
            TS_ASSERT_EQUALS(node0->rGetLocation()[1], 0.0);
            TS_ASSERT_EQUALS(node1->rGetLocation()[0], 1.0);
            TS_ASSERT_EQUALS(node1->rGetLocation()[1], 1.0);

            std::set<unsigned int> arch_elements = testEdge->GetOtherElements(5);
            std::set<unsigned int> check_set;
            check_set.insert(0);
            check_set.insert(1);

            TS_ASSERT_EQUALS(arch_elements, check_set);
            TS_ASSERT_EQUALS(testEdge->GetIndex(), 0u);
            delete testEdge;
        }
    }

    void TestEdgeOperations()
    {
        //Testing EdgeOperation and EdgeMapInfo class methods
        EdgeRemapInfo* info = new EdgeRemapInfo;
        info->SetSplitProportions(std::vector<double>(2, 0.5));
        EdgeOperation edge_oper (EDGE_OPERATION_ADD, 2, info, false);

        TS_ASSERT_EQUALS(edge_oper.GetOperation(), EDGE_OPERATION_ADD);
        TS_ASSERT_EQUALS(edge_oper.GetElementIndex(), 2u);
        TS_ASSERT_EQUALS(edge_oper.IsElementIndexRemapped(), false);

        edge_oper.SetElementIndex(1);
        edge_oper.SetElementIndex2(5);
        TS_ASSERT_EQUALS(edge_oper.GetElementIndex(), 1u);
        TS_ASSERT_EQUALS(edge_oper.GetElementIndex2(), 5u);

        TS_ASSERT_EQUALS(info, edge_oper.GetRemapInfo());
        TS_ASSERT_EQUALS(edge_oper.GetRemapInfo()->GetSplitProportions(), std::vector<double>(2,0.5));

        EdgeRemapInfo* info_1 = new EdgeRemapInfo(std::vector<long>(4, 7), std::vector<unsigned> (4, 8));
        EdgeRemapInfo* info_2 = new EdgeRemapInfo(std::vector<long>(3, 5), std::vector<unsigned> (3, 4));
        EdgeOperation edge_oper2 (0, 1, info_1, info_2);
        TS_ASSERT_EQUALS(edge_oper2.GetRemapInfo()->GetEdgesMapping(), std::vector<long>(4, 7));
        TS_ASSERT_EQUALS(edge_oper2.GetRemapInfo()->GetEdgesStatus(), std::vector<unsigned>(4, 8));
        TS_ASSERT_EQUALS(edge_oper2.GetRemapInfo2()->GetEdgesMapping(), std::vector<long>(3, 5));
        TS_ASSERT_EQUALS(edge_oper2.GetRemapInfo2()->GetEdgesStatus(), std::vector<unsigned>(3, 4));
    }

    void TestEdgeOperationArchive()
    {
        OutputFileHandler handler("TestEdgeOperation", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "edgeOperation.arch";

        {
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            EdgeRemapInfo* info = new EdgeRemapInfo;
            info->SetSplitProportions(std::vector<double>(2, 0.5));

            EdgeRemapInfo* info_2 = new EdgeRemapInfo(std::vector<long>(3, 5), std::vector<unsigned> (3, 4));
            EdgeOperation* edge_oper = new EdgeOperation(2, 10, info, info_2);

            // Write the edge to file
            output_arch << edge_oper;
            delete edge_oper;
        }

        {
            // Restore the nodes
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            EdgeOperation* edge_oper;
            input_arch >> edge_oper;

            TS_ASSERT_EQUALS(edge_oper->GetOperation(), EDGE_OPERATION_DIVIDE);
            TS_ASSERT_EQUALS(edge_oper->GetElementIndex(), 2u);
            TS_ASSERT_EQUALS(edge_oper->GetElementIndex2(), 10u);
            TS_ASSERT_EQUALS(edge_oper->IsElementIndexRemapped(), false);

            TS_ASSERT_EQUALS(edge_oper->GetRemapInfo()->GetSplitProportions(), std::vector<double>(2,0.5));
            TS_ASSERT_EQUALS(edge_oper->GetRemapInfo2()->GetEdgesMapping(), std::vector<long>(3, 5));
            TS_ASSERT_EQUALS(edge_oper->GetRemapInfo2()->GetEdgesStatus(), std::vector<unsigned>(3, 4));

            delete edge_oper;
        }
    }
};

#endif //_TESTEDGE_HPP_
