/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTSEMMESH_HPP_
#define TESTSEMMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "SemMesh.hpp"
#include "ArchiveOpener.hpp"

class TestSemMesh : public CxxTest::TestSuite
{
public:

    void Test2DSemMeshConstructor() throw(Exception)
    {
        // Create a two Potts elements containing 100 nodes each.
        std::vector<Node<2>* > nodes_0;
        for (unsigned i=0;i<10; i++)
        {
            for (unsigned j=0; j<10; j++)
            {
                nodes_0.push_back(new Node<2>(j+10*i, double(j/9.0), double(i/9.0)));
            }
        }
        std::vector<Node<2>* > nodes_1;
        for (unsigned i=0;i<10; i++)
        {
            for (unsigned j=0; j<10; j++)
            {
                nodes_1.push_back(new Node<2>(j+10*i+100, 1.05 + double(j/9.0), double(i/9.0)));
            }
        }

        std::vector<PottsElement<2>* > elements;
        elements.push_back(new PottsElement<2>(0, nodes_0));
        elements.push_back(new PottsElement<2>(1, nodes_1));

        std::vector<Node<2>* > nodes;
        nodes.insert(nodes.end(), nodes_0.begin(), nodes_0.end());
        nodes.insert(nodes.end(), nodes_1.begin(), nodes_1.end());

        SemMesh<2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 100u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 100u);

        // Check the nodes know which element they are in.
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetContainingElementIndices().size(), 1u);
        TS_ASSERT_EQUALS((*mesh.GetNode(0)->rGetContainingElementIndices().begin()), 0u);

        TS_ASSERT_EQUALS(mesh.GetNode(100)->rGetContainingElementIndices().size(), 1u);
        TS_ASSERT_EQUALS((*mesh.GetNode(100)->rGetContainingElementIndices().begin()), 1u);

        // Coverage
        TS_ASSERT_EQUALS(mesh.SolveElementMapping(0), 0u);
    }

    void Test3DSemMeshConstructor() throw(Exception)
    {
        // Create a two Potts elements containing 125 nodes each.
        std::vector<Node<3>* > nodes_0;
        for (unsigned i=0;i<5; i++)
        {
            for (unsigned j=0; j<5; j++)
            {
                for (unsigned k=0; k<5; k++)
                {
                    nodes_0.push_back(new Node<3>(j+10*i, double(j/4.0), double(i/4.0), double(k/4.0)));
                }
            }
        }
        std::vector<Node<3>* > nodes_1;
        for (unsigned i=0;i<5; i++)
        {
            for (unsigned j=0; j<5; j++)
            {
                for (unsigned k=0; k<5; k++)
                {
                    nodes_1.push_back(new Node<3>(j+10*i+125, double(j/4.0), double(i/4.0), 1.05 + double(k/4.0)));
                }
            }
        }

        std::vector<PottsElement<3>* > elements;
        elements.push_back(new PottsElement<3>(0, nodes_0));
        elements.push_back(new PottsElement<3>(1, nodes_1));

        std::vector<Node<3>* > nodes;
        nodes.insert(nodes.end(), nodes_0.begin(), nodes_0.end());
        nodes.insert(nodes.end(), nodes_1.begin(), nodes_1.end());

        SemMesh<3> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 125u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 125u);

        // Check the nodes know which element they are in.
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetContainingElementIndices().size(), 1u);
        TS_ASSERT_EQUALS((*mesh.GetNode(0)->rGetContainingElementIndices().begin()), 0u);

        TS_ASSERT_EQUALS(mesh.GetNode(125)->rGetContainingElementIndices().size(), 1u);
        TS_ASSERT_EQUALS((*mesh.GetNode(125)->rGetContainingElementIndices().begin()), 1u);
    }

    void TestSemElementIterator() throw (Exception)
    {
        // Create a two Potts elements containing 100 nodes each.
        std::vector<Node<2>* > nodes_0;
        for (unsigned i=0;i<10; i++)
        {
            for (unsigned j=0; j<10; j++)
            {
                nodes_0.push_back(new Node<2>(j+10*i, double(j/9.0), double(i/9.0)));
            }
        }
        std::vector<Node<2>* > nodes_1;
        for (unsigned i=0;i<10; i++)
        {
            for (unsigned j=0; j<10; j++)
            {
                nodes_1.push_back(new Node<2>(j+10*i+100, 1.05 + double(j/9.0), double(i/9.0)));
            }
        }

        std::vector<PottsElement<2>* > elements;
        elements.push_back(new PottsElement<2>(0, nodes_0));
        elements.push_back(new PottsElement<2>(1, nodes_1));

        std::vector<Node<2>* > nodes;
        nodes.insert(nodes.end(), nodes_0.begin(), nodes_0.end());
        nodes.insert(nodes.end(), nodes_1.begin(), nodes_1.end());

        SemMesh<2> mesh(nodes, elements);

        unsigned counter = 0;
        for (SemMesh<2>::SemElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter, element_index); // assumes the iterator will give elements 0,1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(mesh.GetNumElements(), counter);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), counter);

        // Check that the element iterator correctly handles deleted elements
        mesh.GetElement(0)->MarkAsDeleted();

        counter = 0;
        for (SemMesh<2>::SemElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            // This time check the * operator
            unsigned element_index = (*iter).GetIndex();
            TS_ASSERT_EQUALS(counter+1, element_index); // assumes the iterator will give elements 0,1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(mesh.GetNumElements()-1, counter);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements()-1, counter);

        // For coverage, test with an empty mesh
        SemMesh<2> empty_mesh;

        // Since the mesh is empty, the iterator should be set to mrMesh.mNodes.end() when constructed
        SemMesh<2>::SemElementIterator iter = empty_mesh.GetElementIteratorBegin();

        /*
         * Check that the iterator is now at the end (we need to check this as a double-negative,
         * as we only have a NOT-equals operator defined on the iterator).
         */
        bool iter_is_not_at_end = (iter != empty_mesh.GetElementIteratorEnd());
        TS_ASSERT_EQUALS(iter_is_not_at_end, false);
    }

    void TestAddElement() throw (Exception)
    {
        // Create a two Potts elements containing 100 nodes each.
        unsigned num_elements = 6;
        std::vector<std::vector<Node<2>* > > element_nodes;
        for(unsigned elem_index=0; elem_index < num_elements; elem_index++)
        {
            std::vector<Node<2>* > nodes;
            for (unsigned i=0;i<10; i++)
            {
                for (unsigned j=0; j<10; j++)
                {
                    nodes.push_back(new Node<2>(100*elem_index + 10*i +j, 1.01*elem_index + j/9.0, i/9.0));
                }
            }
            element_nodes.push_back(nodes);
        }

        std::vector<PottsElement<2>* > elements;
        std::vector<Node<2>* > all_nodes;
        for(unsigned i=0; i<num_elements; i++)
        {
            elements.push_back(new PottsElement<2>(i, element_nodes[i]));
            all_nodes.insert(all_nodes.end(), element_nodes[i].begin(), element_nodes[i].end());
        }

        SemMesh<2> mesh(all_nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_elements*100);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6u);

        std::vector<Node<2>* > new_nodes;
        for (unsigned i=0;i<10; i++)
        {
            for (unsigned j=0; j<10; j++)
            {
                new_nodes.push_back(new Node<2>(100*6 + 10*i +j, 1.01*6 + j/9.0, i/9.0));
            }
        }
        // Create a new potts element.
        PottsElement<2>* new_element = new PottsElement<2>(6, new_nodes);

        mesh.AddElement(new_element, new_nodes);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 700u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 7u);

        // Make sure the new nodes are associated with the correct element
        for(unsigned node_index=600; node_index<700; node_index++)
        {
            TS_ASSERT_EQUALS(mesh.GetNode(node_index)->rGetContainingElementIndices().size(), 1u);
            TS_ASSERT_EQUALS((*mesh.GetNode(node_index)->rGetContainingElementIndices().begin()), 6u);
        }
    }

    void TestDeleteElementAndRemesh() throw (Exception)
    {
        // Create a two Potts elements containing 100 nodes each.
       unsigned num_elements = 6;
       std::vector<std::vector<Node<2>* > > element_nodes;
       for(unsigned elem_index=0; elem_index < num_elements; elem_index++)
       {
           std::vector<Node<2>* > nodes;
           for (unsigned i=0;i<10; i++)
           {
               for (unsigned j=0; j<10; j++)
               {
                   nodes.push_back(new Node<2>(100*elem_index + 10*i +j, 1.01*elem_index + j/9.0, i/9.0));
               }
           }
           element_nodes.push_back(nodes);
       }

       std::vector<PottsElement<2>* > elements;
       std::vector<Node<2>* > all_nodes;
       for(unsigned i=0; i<num_elements; i++)
       {
           elements.push_back(new PottsElement<2>(i, element_nodes[i]));
           all_nodes.insert(all_nodes.end(), element_nodes[i].begin(), element_nodes[i].end());
       }

       SemMesh<2> mesh(all_nodes, elements);

       // Delete Element from the end.
       mesh.DeleteElement(2);

       TS_ASSERT_EQUALS(mesh.GetNumElements(), 5u);
       TS_ASSERT_EQUALS(mesh.GetNumNodes(), 500u);

       TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 6u);
       TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 600u);

       mesh.ReMesh();

       TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 5u);
       TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 500u);

       /**
        *   New configuration should look like
        *
        *   Old: 0  1   2   3   4   5
        *   New: 0  1   *   2   3   4
        */

       // Check the elements have been re-shuffled correctly.
       PottsElement<2>* p_element = mesh.GetElement(2);

       for(unsigned i=0; i<p_element->GetNumNodes(); i++)
       {
           TS_ASSERT(p_element->GetNode(i)->GetIndex() > 299);
           TS_ASSERT(p_element->GetNode(i)->GetIndex() < 400);
       }
    }
};

#endif /*TESTSEMMESH_HPP_*/
