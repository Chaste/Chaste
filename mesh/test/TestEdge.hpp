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

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/shared_ptr.hpp>

#include "ArchiveOpener.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestEdge : public CxxTest::TestSuite
{
public:
    /**
     * Check whether edge building works as intended
     */
    void TestEdgeInitialisation()
    {
        const unsigned example_edge_index = 1;
        auto p_test_edge = std::make_unique<Edge<2> >(example_edge_index);

        TS_ASSERT_EQUALS(p_test_edge->GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(example_edge_index, p_test_edge->GetIndex());

        // Add nodes to the edge, check validity and centre position
        Node<2> node0(0, ChastePoint<2>(0, 0, 0));
        Node<2> node1(1, ChastePoint<2>(1, 1, 1));

        p_test_edge->SetNodes(&node0, &node1);
        auto edgeCentreLocation = p_test_edge->rGetCentreLocation();
        TS_ASSERT_EQUALS(p_test_edge->GetNumNodes(), 2u);
        TS_ASSERT(p_test_edge->GetNode(0) == &node0 || p_test_edge->GetNode(0) == &node1);
        TS_ASSERT(p_test_edge->GetNode(1) == &node0 || p_test_edge->GetNode(1) == &node1);
        TS_ASSERT_DIFFERS(p_test_edge->GetNode(0), p_test_edge->GetNode(1));
        TS_ASSERT(edgeCentreLocation[0] == 0.5 && edgeCentreLocation[1] == 0.5);

        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.0, -1.0));
        nodes.push_back(new Node<2>(5, false, 1.0, -1.0));

        // Generate two line vertex elements
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[5]);

        std::vector<VertexElement<2, 2>*> elements;
        elements.push_back(new VertexElement<2, 2>(0, nodes_elem_0));
        elements.push_back(new VertexElement<2, 2>(1, nodes_elem_1));

        // Generate a mesh which will automatically build the edges in the constructor
        auto p_mesh = std::make_unique<VertexMesh<2, 2> >(nodes, elements);
        const EdgeHelper<2>& edge_helper = p_mesh->rGetEdgeHelper();
        // There are two elements in our mesh
        // We test Edge class methods here
        for (unsigned i = 0; i < 2; i++)
        {
            VertexElement<2, 2>* p_element = elements[i];
            const unsigned num_edges = p_element->GetNumEdges();
            for (unsigned index = 0; index < num_edges; index++)
            {
                Edge<2>* p_edge = p_element->GetEdge(index);
                TS_ASSERT_EQUALS(p_edge->GetNumNodes(), 2u);
                TS_ASSERT_DIFFERS(p_edge->GetNode(0), p_edge->GetNode(1));
                TS_ASSERT(p_edge->GetNumElements() > 0 && p_edge->GetNumElements() <= 2);
                TS_ASSERT_EQUALS(p_element->ContainsEdge(p_edge), true);
            }
        }
        for (unsigned i = 0; i < p_mesh->GetNumEdges(); i++)
        {
            Edge<2>* p_edge = p_mesh->GetEdge(i);
            TS_ASSERT_EQUALS(edge_helper.GetEdge(i), p_edge);
            TS_ASSERT(elements[0]->ContainsEdge(p_edge) || elements[1]->ContainsEdge(p_edge));
        }

        for (unsigned i = 0; i < 2; i++)
        {
            VertexElement<2, 2>* p_element = elements[i];
            p_element->ClearEdges();
            TS_ASSERT_EQUALS(p_element->GetNumEdges(), 0u);
        }

        // Also test constructors in honeycomb mesh (MutableVertexMesh)
        HoneycombVertexMeshGenerator generator(2, 2);
        boost::shared_ptr<MutableVertexMesh<2, 2> > honeycombMesh = generator.GetMesh();
        TS_ASSERT_EQUALS(honeycombMesh->GetNumEdges(), 19u);
    }

    /**
     * Check whether edge building works as intended
     */
    void TestEdgeInitialisation1d()
    {
        unsigned example_edge_index = 1;
        auto p_test_edge = std::make_unique<Edge<2> >(example_edge_index);

        TS_ASSERT_EQUALS(p_test_edge->GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(example_edge_index, p_test_edge->GetIndex());

        // Add nodes to the edge, check validity and centre position
        Node<2> node0(0, ChastePoint<2>(0));
        Node<2> node1(1, ChastePoint<2>(1));

        p_test_edge->SetNodes(&node0, &node1);
        auto edgeCentreLocation = p_test_edge->rGetCentreLocation();
        TS_ASSERT_EQUALS(p_test_edge->GetNumNodes(), 2u);
        TS_ASSERT(p_test_edge->GetNode(0) == &node0 || p_test_edge->GetNode(0) == &node1);
        TS_ASSERT(p_test_edge->GetNode(1) == &node0 || p_test_edge->GetNode(1) == &node1);
        TS_ASSERT_DIFFERS(p_test_edge->GetNode(0), p_test_edge->GetNode(1));

        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0));
        nodes.push_back(new Node<2>(2, false, -1.0));

        // Generate two line vertex elements
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);

        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[1]);
        nodes_elem_2.push_back(nodes[2]);

        std::vector<VertexElement<1, 2>*> elements;
        elements.push_back(new VertexElement<1, 2>(0, nodes_elem_1));
        elements.push_back(new VertexElement<1, 2>(1, nodes_elem_2));

        // Generate a mesh which will automatically build the edges in the constructor
        auto p_mesh = std::make_unique<VertexMesh<1, 2> >(nodes, elements);
        const EdgeHelper<2>& edge_helper = p_mesh->rGetEdgeHelper();
        // There are two elements in our mesh
        // We test Edge class methods here
        for (unsigned i = 0; i < 2; i++)
        {
            VertexElement<1, 2>* p_element = elements[i];
            const unsigned num_edges = p_element->GetNumEdges();
            for (unsigned index = 0; index < num_edges; index++)
            {
                Edge<2>* p_edge = p_element->GetEdge(index);
                TS_ASSERT_EQUALS(p_edge->GetNumNodes(), 2u);
                TS_ASSERT_DIFFERS(p_edge->GetNode(0), p_edge->GetNode(1));
                TS_ASSERT(p_edge->GetNumElements() > 0 && p_edge->GetNumElements() <= 2);
                TS_ASSERT(p_element->GetLocalEdgeIndex(p_edge) > 0 && p_element->GetLocalEdgeIndex(p_edge) < 2);
                TS_ASSERT_LESS_THAN(p_element->GetEdgeGlobalIndex(index), 2);

                unsigned edge_local_index = p_element->GetLocalEdgeIndex(p_edge);
                TS_ASSERT_EQUALS(p_element->GetNeighbouringElementAtEdgeIndex(edge_local_index).size(), 0u);
                TS_ASSERT_EQUALS(p_element->ContainsEdge(p_edge), true);
            }
        }

        for (unsigned i = 0; i < p_mesh->GetNumEdges(); i++)
        {
            Edge<2>* p_edge = p_mesh->GetEdge(i);
            TS_ASSERT_EQUALS(edge_helper.GetEdge(i), p_edge);
            TS_ASSERT(elements[0]->ContainsEdge(p_edge) || elements[1]->ContainsEdge(p_edge));
        }

        for (unsigned i = 0; i < 2; i++)
        {
            VertexElement<1, 2>* p_element = elements[i];
            p_element->ClearEdges();
            TS_ASSERT_EQUALS(p_element->GetNumEdges(), 0u);
        }
    }

    void TestEdgeProperties()
    {
        unsigned example_edge_index = 1;

        auto p_test_edge = std::make_unique<Edge<2> >(example_edge_index);
        TS_ASSERT_EQUALS(p_test_edge->GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(example_edge_index, p_test_edge->GetIndex());

        // Add nodes to the edge, check validity and centre position
        Node<2> node0(0, ChastePoint<2>(0, 0, 0));
        Node<2> node1(1, ChastePoint<2>(1, 1, 1));

        p_test_edge->SetNodes(&node1, &node0);
        auto edge_node_indices = p_test_edge->GetMapIndex();
        TS_ASSERT_EQUALS(edge_node_indices.first, 0u);
        TS_ASSERT_EQUALS(edge_node_indices.second, 1u);

        p_test_edge->MarkAsDeleted();
        TS_ASSERT(p_test_edge->IsDeleted());

        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.0, -1.0));
        nodes.push_back(new Node<2>(5, false, 1.0, -1.0));

        // Generate two square vertex elements
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[5]);

        std::vector<VertexElement<2, 2>*> elements;
        elements.push_back(new VertexElement<2, 2>(0, nodes_elem_0));
        elements.push_back(new VertexElement<2, 2>(1, nodes_elem_1));

        // Generate a mesh which will automatically build the edges in the constructor
        auto p_mesh = std::make_unique<VertexMesh<2, 2> >(nodes, elements);
        Edge<2>* edge0 = p_mesh->GetElement(0)->GetEdge(0);
        TS_ASSERT_EQUALS(edge0->IsBoundaryEdge(), false);
        Edge<2>* edge1 = p_mesh->GetElement(0)->GetEdge(1);
        TS_ASSERT_EQUALS(edge1->IsBoundaryEdge(), true);

        const std::set<unsigned> elements_contain_edge = edge0->GetNeighbouringElementIndices();
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

            Edge<2>* p_test_edge = new Edge<2>(0);
            p_test_edge->AddElement(0);
            p_test_edge->AddElement(1);

            Node<2> node0(0, ChastePoint<2>(0, 0.0, 0.0));
            Node<2> node1(1, ChastePoint<2>(1, 1.0, 1.0));

            p_test_edge->SetNodes(&node0, &node1);

            // Write the edge to file
            output_arch << p_test_edge;
            delete p_test_edge;
        }

        {
            // Restore the nodes
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            Edge<2>* p_test_edge;
            input_arch >> p_test_edge;

            Node<2>* p_node0 = p_test_edge->GetNode(0);
            Node<2>* p_node1 = p_test_edge->GetNode(1);

            TS_ASSERT_EQUALS(p_node0->GetIndex(), 0u);
            TS_ASSERT_EQUALS(p_node1->GetIndex(), 1u);

            TS_ASSERT_DELTA(p_node0->rGetLocation()[0], 0.0, 1e-9);
            TS_ASSERT_DELTA(p_node0->rGetLocation()[1], 0.0, 1e-9);
            TS_ASSERT_DELTA(p_node1->rGetLocation()[0], 1.0, 1e-9);
            TS_ASSERT_DELTA(p_node1->rGetLocation()[1], 1.0, 1e-9);

            std::set<unsigned> arch_elements = p_test_edge->GetOtherElements(5);
            std::set<unsigned> check_set;
            check_set.insert(0);
            check_set.insert(1);

            TS_ASSERT_EQUALS(arch_elements, check_set);
            TS_ASSERT_EQUALS(p_test_edge->GetIndex(), 0u);

            // Unarchiving has created two new nodes and the edge
            delete p_node0;
            delete p_node1;
            delete p_test_edge;
        }
    }

    void TestEdgeOperations()
    {
        // Testing EdgeOperation and EdgeMapInfo class methods
        EdgeRemapInfo info;
        info.SetSplitProportions(std::vector<double>(2, 0.5));
        EdgeOperation edge_oper(EDGE_OPERATION_ADD, 2, info, false);

        TS_ASSERT_EQUALS(edge_oper.GetOperation(), EDGE_OPERATION_ADD);
        TS_ASSERT_EQUALS(edge_oper.GetElementIndex(), 2u);
        TS_ASSERT_EQUALS(edge_oper.IsElementIndexRemapped(), false);

        edge_oper.SetElementIndex(1);
        edge_oper.SetElementIndex2(5);
        TS_ASSERT_EQUALS(edge_oper.GetElementIndex(), 1u);
        TS_ASSERT_EQUALS(edge_oper.GetElementIndex2(), 5u);

        TS_ASSERT_EQUALS(edge_oper.rGetRemapInfo().GetSplitProportions(), std::vector<double>(2, 0.5));

        const EdgeRemapInfo info_1(std::vector<long>(4, 7), std::vector<unsigned>(4, 8));
        const EdgeRemapInfo info_2(std::vector<long>(3, 5), std::vector<unsigned>(3, 4));
        EdgeOperation edge_oper2(0, 1, info_1, info_2);
        TS_ASSERT_EQUALS(edge_oper2.rGetRemapInfo().GetEdgesMapping(), std::vector<long>(4, 7));
        TS_ASSERT_EQUALS(edge_oper2.rGetRemapInfo().GetEdgesStatus(), std::vector<unsigned>(4, 8));
        TS_ASSERT_EQUALS(edge_oper2.rGetRemapInfo2().GetEdgesMapping(), std::vector<long>(3, 5));
        TS_ASSERT_EQUALS(edge_oper2.rGetRemapInfo2().GetEdgesStatus(), std::vector<unsigned>(3, 4));
        TS_ASSERT_EQUALS(edge_oper2.rGetRemapInfo().GetUnused(), false);
        TS_ASSERT_EQUALS(edge_oper2.rGetRemapInfo2().GetUnused(), false);
    }

    void TestEdgeOperationArchive()
    {
        OutputFileHandler handler("TestEdgeOperation", false);
        const std::string archive_filename = handler.GetOutputDirectoryFullPath() + "edgeOperation.arch";

        {
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            EdgeRemapInfo info;
            info.SetSplitProportions(std::vector<double>(2, 0.5));

            const EdgeRemapInfo info_2(std::vector<long>(3, 5), std::vector<unsigned>(3, 4));
            const EdgeOperation edge_oper(2, 10, info, info_2);

            // Write the edge to file
            output_arch << edge_oper;
        }

        {
            // Restore edge operation
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            EdgeOperation edge_oper;
            input_arch >> edge_oper;

            TS_ASSERT_EQUALS(edge_oper.GetOperation(), EDGE_OPERATION_DIVIDE);
            TS_ASSERT_EQUALS(edge_oper.GetElementIndex(), 2u);
            TS_ASSERT_EQUALS(edge_oper.GetElementIndex2(), 10u);
            TS_ASSERT_EQUALS(edge_oper.IsElementIndexRemapped(), false);

            TS_ASSERT_EQUALS(edge_oper.rGetRemapInfo().GetSplitProportions(), std::vector<double>(2, 0.5));
            TS_ASSERT_EQUALS(edge_oper.rGetRemapInfo2().GetEdgesMapping(), std::vector<long>(3, 5));
            TS_ASSERT_EQUALS(edge_oper.rGetRemapInfo2().GetEdgesStatus(), std::vector<unsigned>(3, 4));
        }
    }

    void TestEdgeGenerateMapIndex()
    {
        TS_ASSERT_EQUALS(std::make_pair(1u, 2u), Edge<2>::GenerateMapIndex(1u, 2u));
        TS_ASSERT_EQUALS(std::make_pair(1u, 2u), Edge<2>::GenerateMapIndex(2u, 1u));
    }
};

#endif //_TESTEDGE_HPP_
