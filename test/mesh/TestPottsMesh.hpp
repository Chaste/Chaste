/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTPOTTSMESH_HPP_
#define TESTPOTTSMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "PottsMeshWriter.hpp"
#include "PottsMeshReader.hpp"
#include "PottsMesh.hpp"
#include "PottsMeshGenerator.hpp"
#include "ArchiveOpener.hpp"

class TestPottsMesh : public CxxTest::TestSuite
{
public:
    void TestBasic2dPottsMesh() throw(Exception)
    {
        // Make 6 nodes to assign to two elements
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        basic_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        basic_nodes.push_back(new Node<2>(2, false, 2.0, 0.0));
        basic_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        basic_nodes.push_back(new Node<2>(4, false, 1.0, 1.0));
        basic_nodes.push_back(new Node<2>(5, false, 2.0, 1.0));

        // Make two triangular elements out of these nodes
        std::vector<std::vector<Node<2>*> > nodes_elements(2);
        nodes_elements[0].push_back(basic_nodes[0]);
        nodes_elements[0].push_back(basic_nodes[1]);
        nodes_elements[0].push_back(basic_nodes[3]);

        nodes_elements[1].push_back(basic_nodes[2]);
        nodes_elements[1].push_back(basic_nodes[4]);
        nodes_elements[1].push_back(basic_nodes[5]);

        std::vector<PottsElement<2>*> basic_potts_elements;
        basic_potts_elements.push_back(new PottsElement<2>(0, nodes_elements[0]));
        basic_potts_elements.push_back(new PottsElement<2>(1, nodes_elements[1]));

        // Make a PottsMesh
        PottsMesh<2> basic_potts_mesh(basic_nodes, basic_potts_elements);

        TS_ASSERT_EQUALS(basic_potts_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(basic_potts_mesh.GetNumNodes(), 6u);

        TS_ASSERT_DELTA(basic_potts_mesh.GetNode(2)->rGetLocation()[0], 2.0, 1e-3);
        TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(1)->GetNode(2)->GetIndex(),5u);

        // Check that the nodes know which elements they are in

        // Nodes 0 1 and 3 are only in element 0
        std::set<unsigned> temp_list1;
        temp_list1.insert(0u);
        TS_ASSERT_EQUALS(basic_nodes[0]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(basic_nodes[1]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(basic_nodes[3]->rGetContainingElementIndices(), temp_list1);

        // Node 2 4 and 5 are only in element 1
        std::set<unsigned> temp_list2;
        temp_list2.insert(1u);
        TS_ASSERT_EQUALS(basic_nodes[2]->rGetContainingElementIndices(), temp_list2);
        TS_ASSERT_EQUALS(basic_nodes[4]->rGetContainingElementIndices(), temp_list2);
        TS_ASSERT_EQUALS(basic_nodes[5]->rGetContainingElementIndices(), temp_list2);

        // Test Area and Perimeter of elements
        TS_ASSERT_DELTA(basic_potts_mesh.GetVolumeOfElement(0), 3.0, 1e-12);
        TS_ASSERT_DELTA(basic_potts_mesh.GetSurfaceAreaOfElement(0), 8.0, 1e-12);

        TS_ASSERT_DELTA(basic_potts_mesh.GetVolumeOfElement(1), 3.0, 1e-12);
        TS_ASSERT_DELTA(basic_potts_mesh.GetSurfaceAreaOfElement(1), 8.0, 1e-12);

        // Test GetCentroidOfElements
        TS_ASSERT_DELTA(basic_potts_mesh.GetCentroidOfElement(0)[0], 1.0/3.0, 1e-12);
        TS_ASSERT_DELTA(basic_potts_mesh.GetCentroidOfElement(0)[1], 1.0/3.0, 1e-12);

        TS_ASSERT_DELTA(basic_potts_mesh.GetCentroidOfElement(1)[0], 5.0/3.0, 1e-12);
        TS_ASSERT_DELTA(basic_potts_mesh.GetCentroidOfElement(1)[1], 2.0/3.0, 1e-12);

        // Test GetVectorFromAtoB Method
        c_vector<double, 2> vector=basic_potts_mesh.GetVectorFromAtoB(basic_potts_mesh.GetNode(0)->rGetLocation(), basic_potts_mesh.GetNode(1)->rGetLocation());
        TS_ASSERT_DELTA(vector[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(vector[1], 0.0, 1e-12);

        // Coverage
        TS_ASSERT_EQUALS(basic_potts_mesh.SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(basic_potts_mesh.SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(basic_potts_mesh.SolveBoundaryElementMapping(0), 0u);
        TS_ASSERT_EQUALS(basic_potts_mesh.IsMeshChanging(), true);
    }



    void TestBasic3dPottsMesh() throw(Exception)
	{
		// Make 8 nodes to assign to two elements
		std::vector<Node<3>*> basic_nodes;
		basic_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
		basic_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
		basic_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
		basic_nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 0.0));
		basic_nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 1.0));
		basic_nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 1.0));
		basic_nodes.push_back(new Node<3>(6, false, 0.0, 1.0, 1.0));
	    basic_nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));

		// Make two elements out of these nodes leaving some sites free
		std::vector<std::vector<Node<3>*> > nodes_elements(2);
		nodes_elements[0].push_back(basic_nodes[1]);
		nodes_elements[0].push_back(basic_nodes[2]);
		nodes_elements[0].push_back(basic_nodes[3]);

		nodes_elements[1].push_back(basic_nodes[0]);
		nodes_elements[1].push_back(basic_nodes[4]);
		nodes_elements[1].push_back(basic_nodes[5]);

		std::vector<PottsElement<3>*> basic_potts_elements;
		basic_potts_elements.push_back(new PottsElement<3>(0, nodes_elements[0]));
		basic_potts_elements.push_back(new PottsElement<3>(1, nodes_elements[1]));

		// Make a PottsMesh
		PottsMesh<3> basic_potts_mesh(basic_nodes, basic_potts_elements);

		TS_ASSERT_EQUALS(basic_potts_mesh.GetNumElements(), 2u);
		TS_ASSERT_EQUALS(basic_potts_mesh.GetNumNodes(), 8u);

		TS_ASSERT_DELTA(basic_potts_mesh.GetNode(2)->rGetLocation()[0], 0.0, 1e-3);
		TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(1)->GetNode(2)->GetIndex(), 5u);

		// Check that the nodes know which elements they are in

		// Nodes 1 2 and 3 are only in element 0
		std::set<unsigned> temp_list1;
		temp_list1.insert(0u);
		TS_ASSERT_EQUALS(basic_nodes[1]->rGetContainingElementIndices(), temp_list1);
		TS_ASSERT_EQUALS(basic_nodes[2]->rGetContainingElementIndices(), temp_list1);
		TS_ASSERT_EQUALS(basic_nodes[3]->rGetContainingElementIndices(), temp_list1);

		// Node 0 4 and 5 are only in element 1
		std::set<unsigned> temp_list2;
		temp_list2.insert(1u);
		TS_ASSERT_EQUALS(basic_nodes[0]->rGetContainingElementIndices(), temp_list2);
		TS_ASSERT_EQUALS(basic_nodes[4]->rGetContainingElementIndices(), temp_list2);
		TS_ASSERT_EQUALS(basic_nodes[5]->rGetContainingElementIndices(), temp_list2);

		// Test Area and Volume of elements
		TS_ASSERT_DELTA(basic_potts_mesh.GetVolumeOfElement(0), 3.0, 1e-12);
		TS_ASSERT_DELTA(basic_potts_mesh.GetSurfaceAreaOfElement(0), 14.0, 1e-12);

		TS_ASSERT_DELTA(basic_potts_mesh.GetVolumeOfElement(1), 3.0, 1e-12);
		TS_ASSERT_DELTA(basic_potts_mesh.GetSurfaceAreaOfElement(1), 14.0, 1e-12);

		// Test GetCentroidOfElements
		TS_ASSERT_DELTA(basic_potts_mesh.GetCentroidOfElement(0)[0], 2.0/3.0, 1e-12);
		TS_ASSERT_DELTA(basic_potts_mesh.GetCentroidOfElement(0)[1], 2.0/3.0, 1e-12);
		TS_ASSERT_DELTA(basic_potts_mesh.GetCentroidOfElement(0)[2], 0.0, 1e-12);

		TS_ASSERT_DELTA(basic_potts_mesh.GetCentroidOfElement(1)[0], 1.0/3.0, 1e-12);
		TS_ASSERT_DELTA(basic_potts_mesh.GetCentroidOfElement(1)[1], 0.0, 1e-12);
		TS_ASSERT_DELTA(basic_potts_mesh.GetCentroidOfElement(1)[2], 2.0/3.0, 1e-12);

		// Test GetVectorFromAtoB Method
		c_vector<double, 3> vector=basic_potts_mesh.GetVectorFromAtoB(basic_potts_mesh.GetNode(0)->rGetLocation(), basic_potts_mesh.GetNode(7)->rGetLocation());
		TS_ASSERT_DELTA(vector[0], 1.0, 1e-12);
		TS_ASSERT_DELTA(vector[1], 1.0, 1e-12);
		TS_ASSERT_DELTA(vector[2], 1.0, 1e-12);

		// Coverage
		TS_ASSERT_EQUALS(basic_potts_mesh.SolveNodeMapping(0), 0u);
		TS_ASSERT_EQUALS(basic_potts_mesh.SolveElementMapping(0), 0u);
		TS_ASSERT_EQUALS(basic_potts_mesh.SolveBoundaryElementMapping(0), 0u);
		TS_ASSERT_EQUALS(basic_potts_mesh.IsMeshChanging(), true);
	}



    void TestNodeIterator() throw (Exception)
    {
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);

        // Create mesh
        PottsMesh<2>* p_mesh = generator.GetMesh();

        unsigned counter = 0;
        for (PottsMesh<2>::NodeIterator iter = p_mesh->GetNodeIteratorBegin();
             iter != p_mesh->GetNodeIteratorEnd();
             ++iter)
        {
            unsigned node_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter, node_index); // assumes the iterator will give nodes 0,1..,N in that order
            counter++;
        }
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), counter);

        // Check that the node iterator correctly handles deleted nodes
        p_mesh->GetNode(0)->MarkAsDeleted();

        counter = 0;
        for (PottsMesh<2>::NodeIterator iter = p_mesh->GetNodeIteratorBegin();
             iter != p_mesh->GetNodeIteratorEnd();
             ++iter)
        {
            unsigned node_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter+1, node_index); // assumes the iterator will give nodes 1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(p_mesh->GetNumAllNodes(), counter+1);

        // For coverage, test with an empty mesh
        PottsMesh<2> empty_mesh;

        // Since the mesh is empty, the iterator should be set to mrMesh.mNodes.end() when constructed
        PottsMesh<2>::NodeIterator iter = empty_mesh.GetNodeIteratorBegin();

        // Check that the iterator is now at the end (we need to check this as a double-negative,
        // as we only have a NOT-equals operator defined on the iterator).
        bool iter_is_not_at_end = (iter != empty_mesh.GetNodeIteratorEnd());
        TS_ASSERT_EQUALS(iter_is_not_at_end, false);
    }

    void TestPottsElementIterator() throw (Exception)
    {
        // Create mesh
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        unsigned counter = 0;
        for (PottsMesh<2>::PottsElementIterator iter = p_mesh->GetElementIteratorBegin();
             iter != p_mesh->GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter, element_index); // assumes the iterator will give elements 0,1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), counter);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), counter);

        // Check that the element iterator correctly handles deleted elements
        p_mesh->GetElement(0)->MarkAsDeleted();

        counter = 0;
        for (PottsMesh<2>::PottsElementIterator iter = p_mesh->GetElementIteratorBegin();
             iter != p_mesh->GetElementIteratorEnd();
             ++iter)
        {

        	// This time check the * operator
        	unsigned element_index = (*iter).GetIndex();
            TS_ASSERT_EQUALS(counter+1, element_index); // assumes the iterator will give elements 0,1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(p_mesh->GetNumElements()-1, counter);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements()-1, counter);

        // For coverage, test with an empty mesh
        PottsMesh<2> empty_mesh;

        // Since the mesh is empty, the iterator should be set to mrMesh.mNodes.end() when constructed
        PottsMesh<2>::PottsElementIterator iter = empty_mesh.GetElementIteratorBegin();

        /*
         * Check that the iterator is now at the end (we need to check this as a double-negative,
         * as we only have a NOT-equals operator defined on the iterator).
         */
        bool iter_is_not_at_end = (iter != empty_mesh.GetElementIteratorEnd());
        TS_ASSERT_EQUALS(iter_is_not_at_end, false);
    }

    void TestMeshGetWidthAndBoundingBoxMethod()
    {
        // Create mesh
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Test CalculateBoundingBox() method
        ChasteCuboid<2> bounds=p_mesh->CalculateBoundingBox();
        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], 3.0,   1e-4);
        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], 3.0, 1e-4);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], 0.0,    1e-4);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], 0.0,    1e-4);

        // Test GetWidth() method
        double width = p_mesh->GetWidth(0);
        double height = p_mesh->GetWidth(1);

        TS_ASSERT_DELTA(height, 3.0, 1e-4);
        TS_ASSERT_DELTA(width, 3.0, 1e-4);
    }

    void TestGetMooreNeighbouringNodeIndices2d()
    {
        /* Create a simple potts-based mesh with one element
         * Numbering the nodes as follows:
         *
         *     6----7----8
         *     |    |    |
         *     3----4----5
         *     |    |    |
         *     0----1----2
         */
        PottsMeshGenerator<2> generator(3, 1, 3, 3, 1, 3);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Test bottom left node
        std::set<unsigned> neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 3u);

        std::set<unsigned> expected_neighbouring_sites;
        expected_neighbouring_sites.insert(1);
        expected_neighbouring_sites.insert(3);
        expected_neighbouring_sites.insert(4);
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

        // Test non-corner bottom nodes
        neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(1);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 5u);

        expected_neighbouring_sites.clear();
        expected_neighbouring_sites.insert(0);
        expected_neighbouring_sites.insert(2);
        expected_neighbouring_sites.insert(3);
        expected_neighbouring_sites.insert(4);
        expected_neighbouring_sites.insert(5);
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

        // Test bottom right node
        neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(2);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 3u);

        expected_neighbouring_sites.clear();
        expected_neighbouring_sites.insert(1);
        expected_neighbouring_sites.insert(4);
        expected_neighbouring_sites.insert(5);
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

        // Test non-corner left nodes
        neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(3);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 5u);

        expected_neighbouring_sites.clear();
        expected_neighbouring_sites.insert(0);
        expected_neighbouring_sites.insert(1);
        expected_neighbouring_sites.insert(4);
        expected_neighbouring_sites.insert(6);
        expected_neighbouring_sites.insert(7);
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

        // Test centre node
        neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(4);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 8u);

        expected_neighbouring_sites.clear();
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            if (i != 4)
            {
                expected_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

        // Test non-corner right nodes
        neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(5);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 5u);

        expected_neighbouring_sites.clear();
        expected_neighbouring_sites.insert(1);
        expected_neighbouring_sites.insert(2);
        expected_neighbouring_sites.insert(4);
        expected_neighbouring_sites.insert(7);
        expected_neighbouring_sites.insert(8);
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

        // Test top left node
        neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(6);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 3u);

        expected_neighbouring_sites.clear();
        expected_neighbouring_sites.insert(3);
        expected_neighbouring_sites.insert(4);
        expected_neighbouring_sites.insert(7);
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

        // Test non-corner top nodes
        neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(7);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 5u);

        expected_neighbouring_sites.clear();
        expected_neighbouring_sites.insert(3);
        expected_neighbouring_sites.insert(4);
        expected_neighbouring_sites.insert(5);
        expected_neighbouring_sites.insert(6);
        expected_neighbouring_sites.insert(8);
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

        // Test top right node
        neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(8);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 3u);

        expected_neighbouring_sites.clear();
        expected_neighbouring_sites.insert(4);
        expected_neighbouring_sites.insert(5);
        expected_neighbouring_sites.insert(7);
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);
    }

    void TestGetMooreNeighbouringNodeIndices3d()
    {
        /* Create a simple potts-based mesh with one element
         * Numbering the nodes as follows:
         *
		 *     6------7------8           15-----16-----17            24-----25-----26
		 *     |      |      |           |      |      |             |      |      |
		 *     |      |      |           |      |      |             |      |      |
		 *     3------4------5           12-----13-----14            21-----22-----23
		 *     |      |      |           |      |      |             |      |      |
		 *     |      |      |           |      |      |             |      |      |
		 *     0------1------2           9------10-----11            18-----19-----20
         */

        PottsMeshGenerator<3> generator(3, 1, 3, 3, 1, 3, 3, 1, 3);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        // Test bottom left node
        std::set<unsigned> neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 7u);

        std::set<unsigned> expected_neighbouring_sites;
        expected_neighbouring_sites.insert(1);
        expected_neighbouring_sites.insert(3);
        expected_neighbouring_sites.insert(4);
        expected_neighbouring_sites.insert(9);
        expected_neighbouring_sites.insert(10);
        expected_neighbouring_sites.insert(12);
        expected_neighbouring_sites.insert(13);
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

        // Test middle front
        neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(4);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 17u);

        expected_neighbouring_sites.clear();
        expected_neighbouring_sites.insert(0);
        expected_neighbouring_sites.insert(1);
        expected_neighbouring_sites.insert(2);
        expected_neighbouring_sites.insert(3);
        expected_neighbouring_sites.insert(5);
        expected_neighbouring_sites.insert(6);
        expected_neighbouring_sites.insert(7);
        expected_neighbouring_sites.insert(8);
        expected_neighbouring_sites.insert(9);
        expected_neighbouring_sites.insert(10);
        expected_neighbouring_sites.insert(11);
        expected_neighbouring_sites.insert(12);
        expected_neighbouring_sites.insert(13);
        expected_neighbouring_sites.insert(14);
        expected_neighbouring_sites.insert(15);
        expected_neighbouring_sites.insert(16);
        expected_neighbouring_sites.insert(17);
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

        // Test bottom right node
        neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(13);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 26u);

        expected_neighbouring_sites.clear();
        expected_neighbouring_sites.insert(0);
        expected_neighbouring_sites.insert(1);
        expected_neighbouring_sites.insert(2);
        expected_neighbouring_sites.insert(3);
		expected_neighbouring_sites.insert(4);
        expected_neighbouring_sites.insert(5);
        expected_neighbouring_sites.insert(6);
        expected_neighbouring_sites.insert(7);
        expected_neighbouring_sites.insert(8);
		expected_neighbouring_sites.insert(9);
        expected_neighbouring_sites.insert(10);
        expected_neighbouring_sites.insert(11);
        expected_neighbouring_sites.insert(12);
        expected_neighbouring_sites.insert(14);
		expected_neighbouring_sites.insert(15);
        expected_neighbouring_sites.insert(16);
        expected_neighbouring_sites.insert(17);
        expected_neighbouring_sites.insert(18);
        expected_neighbouring_sites.insert(19);
		expected_neighbouring_sites.insert(20);
        expected_neighbouring_sites.insert(21);
        expected_neighbouring_sites.insert(22);
        expected_neighbouring_sites.insert(23);
        expected_neighbouring_sites.insert(24);
		expected_neighbouring_sites.insert(25);
		expected_neighbouring_sites.insert(26);
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

        // Test non-corner left nodes
        neighbouring_sites = p_mesh->GetMooreNeighbouringNodeIndices(26);
        TS_ASSERT_EQUALS(neighbouring_sites.size(), 7u);

        expected_neighbouring_sites.clear();
        expected_neighbouring_sites.insert(25);
        expected_neighbouring_sites.insert(23);
        expected_neighbouring_sites.insert(22);
        expected_neighbouring_sites.insert(14);
        expected_neighbouring_sites.insert(17);
        expected_neighbouring_sites.insert(16);
        expected_neighbouring_sites.insert(13);
        TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);


    }

    void TestGetVonNeumannNeighbouringNodeIndices2d()
	{
		/* Create a simple potts-based mesh with one element
		 * Numbering the nodes as follows:
		 *
		 *     6----7----8
		 *     |    |    |
		 *     3----4----5
		 *     |    |    |
		 *     0----1----2
		 */
        PottsMeshGenerator<2> generator(3, 1, 3, 3, 1, 3);
		PottsMesh<2>* p_mesh = generator.GetMesh();

		// Test bottom left node
		std::set<unsigned> neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(0);
		TS_ASSERT_EQUALS(neighbouring_sites.size(), 2u);

		std::set<unsigned> expected_neighbouring_sites;
		expected_neighbouring_sites.insert(1);
		expected_neighbouring_sites.insert(3);
		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

		// Test non-corner bottom nodes
		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(1);
		TS_ASSERT_EQUALS(neighbouring_sites.size(), 3u);

		expected_neighbouring_sites.clear();
		expected_neighbouring_sites.insert(0);
		expected_neighbouring_sites.insert(2);
		expected_neighbouring_sites.insert(4);
		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

		// Test bottom right node
		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(2);
		TS_ASSERT_EQUALS(neighbouring_sites.size(), 2u);

		expected_neighbouring_sites.clear();
		expected_neighbouring_sites.insert(1);
		expected_neighbouring_sites.insert(5);
		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

		// Test non-corner left nodes
		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(3);
		TS_ASSERT_EQUALS(neighbouring_sites.size(), 3u);

		expected_neighbouring_sites.clear();
		expected_neighbouring_sites.insert(0);
		expected_neighbouring_sites.insert(4);
		expected_neighbouring_sites.insert(6);
		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

		// Test centre node
		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(4);
		TS_ASSERT_EQUALS(neighbouring_sites.size(), 4u);

		expected_neighbouring_sites.clear();
		expected_neighbouring_sites.insert(1);
		expected_neighbouring_sites.insert(3);
		expected_neighbouring_sites.insert(5);
		expected_neighbouring_sites.insert(7);

		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

		// Test non-corner right nodes
		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(5);
		TS_ASSERT_EQUALS(neighbouring_sites.size(), 3u);

		expected_neighbouring_sites.clear();
		expected_neighbouring_sites.insert(2);
		expected_neighbouring_sites.insert(4);
		expected_neighbouring_sites.insert(8);
		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

		// Test top left node
		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(6);
		TS_ASSERT_EQUALS(neighbouring_sites.size(), 2u);

		expected_neighbouring_sites.clear();
		expected_neighbouring_sites.insert(3);
		expected_neighbouring_sites.insert(7);
		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

		// Test non-corner top nodes
		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(7);
		TS_ASSERT_EQUALS(neighbouring_sites.size(), 3u);

		expected_neighbouring_sites.clear();
		expected_neighbouring_sites.insert(4);
		expected_neighbouring_sites.insert(6);
		expected_neighbouring_sites.insert(8);
		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

		// Test top right node
		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(8);
		TS_ASSERT_EQUALS(neighbouring_sites.size(), 2u);

		expected_neighbouring_sites.clear();
		expected_neighbouring_sites.insert(5);
		expected_neighbouring_sites.insert(7);
		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);
	}

    void TestGetVonNeumannNeighbouringNodeIndices3d()
 	{
		/*
		 * Create a simple Potts-based mesh with one element
		 * Numbering the nodes as follows:
		 *
		 *        0  5
		 *        | /
		 *        |/
		 *  1-----x-----3
		 *      / |
		 *     /  |
		 *    4   2
		 */
        PottsMeshGenerator<3> generator(3, 1, 3, 3, 1, 3, 3, 1, 3);
 		PottsMesh<3>* p_mesh = generator.GetMesh();

 		// Test bottom left node
 		std::set<unsigned> neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(0);
 		TS_ASSERT_EQUALS(neighbouring_sites.size(), 3u);

 		std::set<unsigned> expected_neighbouring_sites;
 		expected_neighbouring_sites.insert(1);
 		expected_neighbouring_sites.insert(3);
 		expected_neighbouring_sites.insert(9);
 		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

 		// Test non-corner bottom nodes
 		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(1);
 		TS_ASSERT_EQUALS(neighbouring_sites.size(), 4u);

 		expected_neighbouring_sites.clear();
 		expected_neighbouring_sites.insert(0);
 		expected_neighbouring_sites.insert(2);
 		expected_neighbouring_sites.insert(4);
 		expected_neighbouring_sites.insert(10);
 		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

 		// Test bottom right node
 		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(2);
 		TS_ASSERT_EQUALS(neighbouring_sites.size(), 3u);

 		expected_neighbouring_sites.clear();
 		expected_neighbouring_sites.insert(1);
 		expected_neighbouring_sites.insert(5);
 		expected_neighbouring_sites.insert(11);
 		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

 		// Test non-corner left nodes
 		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(3);
 		TS_ASSERT_EQUALS(neighbouring_sites.size(), 4u);

 		expected_neighbouring_sites.clear();
 		expected_neighbouring_sites.insert(0);
 		expected_neighbouring_sites.insert(4);
 		expected_neighbouring_sites.insert(6);
 		expected_neighbouring_sites.insert(12);
 		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

 		// Test centre node on the first slice
 		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(4);
 		TS_ASSERT_EQUALS(neighbouring_sites.size(), 5u);

 		expected_neighbouring_sites.clear();
 		expected_neighbouring_sites.insert(1);
 		expected_neighbouring_sites.insert(3);
 		expected_neighbouring_sites.insert(5);
 		expected_neighbouring_sites.insert(7);
 		expected_neighbouring_sites.insert(13);

 		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

 		// Test non-corner right nodes
 		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(5);
 		TS_ASSERT_EQUALS(neighbouring_sites.size(), 4u);

 		expected_neighbouring_sites.clear();
 		expected_neighbouring_sites.insert(2);
 		expected_neighbouring_sites.insert(4);
 		expected_neighbouring_sites.insert(8);
 		expected_neighbouring_sites.insert(14);
 		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

 		// Test top left node
 		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(6);
 		TS_ASSERT_EQUALS(neighbouring_sites.size(), 3u);

 		expected_neighbouring_sites.clear();
 		expected_neighbouring_sites.insert(3);
 		expected_neighbouring_sites.insert(7);
 		expected_neighbouring_sites.insert(15);
 		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

 		// Test non-corner top nodes
 		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(7);
 		TS_ASSERT_EQUALS(neighbouring_sites.size(), 4u);

 		expected_neighbouring_sites.clear();
 		expected_neighbouring_sites.insert(4);
 		expected_neighbouring_sites.insert(6);
 		expected_neighbouring_sites.insert(8);
 		expected_neighbouring_sites.insert(16);
 		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

 		// Test top right node
 		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(8);
 		TS_ASSERT_EQUALS(neighbouring_sites.size(), 3u);

 		expected_neighbouring_sites.clear();
 		expected_neighbouring_sites.insert(5);
 		expected_neighbouring_sites.insert(7);
 		expected_neighbouring_sites.insert(17);
 		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

 		// Test top middle node
 		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(13);
 		TS_ASSERT_EQUALS(neighbouring_sites.size(), 6u);

 		expected_neighbouring_sites.clear();
 		expected_neighbouring_sites.insert(4);
 		expected_neighbouring_sites.insert(22);
 		expected_neighbouring_sites.insert(14);
 		expected_neighbouring_sites.insert(12);
 		expected_neighbouring_sites.insert(10);
 		expected_neighbouring_sites.insert(16);
 		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

		// Test centre node on the back
 		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(22);
 		TS_ASSERT_EQUALS(neighbouring_sites.size(), 5u);

 		expected_neighbouring_sites.clear();
 		expected_neighbouring_sites.insert(21);
 		expected_neighbouring_sites.insert(23);
 		expected_neighbouring_sites.insert(25);
 		expected_neighbouring_sites.insert(19);
 		expected_neighbouring_sites.insert(13);
 		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);

		// Test bottom node on the back
 		neighbouring_sites = p_mesh->GetVonNeumannNeighbouringNodeIndices(19);
 		TS_ASSERT_EQUALS(neighbouring_sites.size(), 4u);

 		expected_neighbouring_sites.clear();
 		expected_neighbouring_sites.insert(18);
 		expected_neighbouring_sites.insert(22);
 		expected_neighbouring_sites.insert(20);
 		expected_neighbouring_sites.insert(10);
 		TS_ASSERT_EQUALS(neighbouring_sites, expected_neighbouring_sites);
 	}

    void Test2dScaleAndTranslate()
    {
        // Create 2D mesh
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        TS_ASSERT_DELTA(p_mesh->GetWidth(0), 3.0, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetWidth(1), 3.0, 1e-4);

        // Squash in the x direction by a factor of 0.5
        p_mesh->Scale(0.5);

        TS_ASSERT_DELTA(p_mesh->GetWidth(0), 1.5, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetWidth(1), 3.0, 1e-4);

        // Stretch in the x and y directions by a factor of 2
        p_mesh->Scale(2.0, 2.0);

        TS_ASSERT_DELTA(p_mesh->GetWidth(0), 3.0, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetWidth(1), 6.0, 1e-4);

        // Test the translate method
        // Pick a certain node and store spatial position
        Node<2>* p_node = p_mesh->GetNode(7);
        ChastePoint<2> original_coordinate = p_node->GetPoint();

        const double x_movement = 1.0;
        const double y_movement = 2.5;

        p_mesh->Translate(x_movement, y_movement);

        ChastePoint<2>  new_coordinate = p_node->GetPoint();

        TS_ASSERT_DELTA(original_coordinate[0], new_coordinate[0] - x_movement, 1e-6);
        TS_ASSERT_DELTA(original_coordinate[1], new_coordinate[1] - y_movement, 1e-6);
    }

    void TestTranslation2DWithUblas()
    {
        // Create 2D mesh
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        c_vector<double, 2> old_location1 = p_mesh->GetNode(4)->rGetLocation();
        c_vector<double, 2> old_location2 = p_mesh->GetNode(9)->rGetLocation();

        // Set translation vector
        c_vector<double, 2> trans_vec;
        trans_vec(0) = 2.0;
        trans_vec(1) = 3.0;

        // Translate
        p_mesh->Translate(trans_vec);
        c_vector<double, 2> new_location1 = p_mesh->GetNode(4)->rGetLocation();
        c_vector<double, 2> new_location2 = p_mesh->GetNode(9)->rGetLocation();

        // Spot check a couple of nodes
        TS_ASSERT_DELTA(new_location1[0], old_location1[0] + 2.0, 1e-6);
        TS_ASSERT_DELTA(new_location1[1], old_location1[1] + 3.0, 1e-6);

        TS_ASSERT_DELTA(new_location2[0], old_location2[0] + 2.0, 1e-6);
        TS_ASSERT_DELTA(new_location2[1], old_location2[1] + 3.0, 1e-6);
    }


    void TestAddElement() throw (Exception)
    {
        // Make 6 nodes to assign to two elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 2.0, 0.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(5, false, 2.0, 1.0));

        // Make two triangular elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[3] = {0, 1, 3};
        unsigned node_indices_elem_1[3] = {2, 4, 5};
        for (unsigned i=0; i<3; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<PottsElement<2>*> elements;
        elements.push_back(new PottsElement<2>(0, nodes_elem_0));

        // Make a vertex mesh
        PottsMesh<2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);

        // Add a new element to the mesh
        mesh.AddElement(new PottsElement<2>(1, nodes_elem_1));

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
    }

    void TestDividePottsElement() throw(Exception)
    {
        {
            // Original Element BELOW new element

            // Make four nodes
            std::vector<Node<2>*> basic_nodes;
            basic_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
            basic_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
            basic_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
            basic_nodes.push_back(new Node<2>(3, false, 1.0, 1.0));

            // Make two rectangular element out of these nodes.
            std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
            nodes_elem_0.push_back(basic_nodes[0]);
            nodes_elem_0.push_back(basic_nodes[2]);
            nodes_elem_1.push_back(basic_nodes[3]);
            nodes_elem_1.push_back(basic_nodes[1]);

            std::vector<PottsElement<2>*> basic_potts_elements;
            basic_potts_elements.push_back(new PottsElement<2>(0, nodes_elem_0));
            basic_potts_elements.push_back(new PottsElement<2>(1, nodes_elem_1));

            // Make a Potts mesh
            PottsMesh<2> basic_potts_mesh(basic_nodes, basic_potts_elements);

            TS_ASSERT_EQUALS(basic_potts_mesh.GetNumElements(), 2u);
            TS_ASSERT_EQUALS(basic_potts_mesh.GetNumNodes(), 4u);

            // Divide elements, putting the original element below the new element
            unsigned new_element_index = basic_potts_mesh.DivideElement(basic_potts_mesh.GetElement(0), true);
            TS_ASSERT_EQUALS(new_element_index, 2u);
            new_element_index = basic_potts_mesh.DivideElement(basic_potts_mesh.GetElement(1), true);
            TS_ASSERT_EQUALS(new_element_index, 3u);

            TS_ASSERT_EQUALS(basic_potts_mesh.GetNumElements(), 4u);

            // Test elements have correct nodes
            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(0)->GetNumNodes(), 1u);
            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);

            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(1)->GetNumNodes(), 1u);
            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(1)->GetNodeGlobalIndex(0), 1u);

            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(2)->GetNumNodes(), 1u);
            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(2)->GetNodeGlobalIndex(0), 2u);

            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(3)->GetNumNodes(), 1u);
            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(3)->GetNodeGlobalIndex(0), 3u);
        }
        {
            // Original element ABOVE new element

            // Make four nodes
            std::vector<Node<2>*> basic_nodes;
            basic_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
            basic_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
            basic_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
            basic_nodes.push_back(new Node<2>(3, false, 1.0, 1.0));

            // Make two rectangular element out of these nodes
            std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
            nodes_elem_0.push_back(basic_nodes[0]);
            nodes_elem_0.push_back(basic_nodes[2]);
            nodes_elem_1.push_back(basic_nodes[3]);
            nodes_elem_1.push_back(basic_nodes[1]);

            std::vector<PottsElement<2>*> basic_potts_elements;
            basic_potts_elements.push_back(new PottsElement<2>(0, nodes_elem_0));
            basic_potts_elements.push_back(new PottsElement<2>(1, nodes_elem_1));

            // Make a Potts mesh
            PottsMesh<2> basic_potts_mesh(basic_nodes, basic_potts_elements);

            TS_ASSERT_EQUALS(basic_potts_mesh.GetNumElements(), 2u);
            TS_ASSERT_EQUALS(basic_potts_mesh.GetNumNodes(), 4u);

            // Divide elements, putting the original element above the new element
            unsigned new_element_index = basic_potts_mesh.DivideElement(basic_potts_mesh.GetElement(0), false);
            TS_ASSERT_EQUALS(new_element_index, 2u);
            new_element_index = basic_potts_mesh.DivideElement(basic_potts_mesh.GetElement(1), false);
            TS_ASSERT_EQUALS(new_element_index, 3u);

            TS_ASSERT_EQUALS(basic_potts_mesh.GetNumElements(), 4u);

            // Test elements have correct nodes
            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(0)->GetNumNodes(), 1u);
            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(0)->GetNodeGlobalIndex(0), 2u);

            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(1)->GetNumNodes(), 1u);
            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(1)->GetNodeGlobalIndex(0), 3u);

            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(2)->GetNumNodes(), 1u);
            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(2)->GetNodeGlobalIndex(0), 0u);

            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(3)->GetNumNodes(), 1u);
            TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(3)->GetNodeGlobalIndex(0), 1u);
        }
        {
            // Testing exceptions

            // Make four nodes
            std::vector<Node<2>*> basic_nodes;
            basic_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
            basic_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
            basic_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
            basic_nodes.push_back(new Node<2>(3, false, 1.0, 1.0));

            // Make one rectangular element out of one of these nodes
            std::vector<Node<2>*> nodes_elem;
            nodes_elem.push_back(basic_nodes[0]);

            std::vector<PottsElement<2>*> basic_potts_elements;
            basic_potts_elements.push_back(new PottsElement<2>(0, nodes_elem));

            // Make a Potts mesh
            PottsMesh<2> basic_potts_mesh(basic_nodes, basic_potts_elements);

            TS_ASSERT_EQUALS(basic_potts_mesh.GetNumElements(), 1u);
            TS_ASSERT_EQUALS(basic_potts_mesh.GetNumNodes(), 4u);

            // Divide element
            TS_ASSERT_THROWS_THIS(basic_potts_mesh.DivideElement(basic_potts_mesh.GetElement(0), false),
                                          "Tried to divide a Potts element with only one node.");
        }
    }

    void TestDeleteAndDividePottsElement() throw(Exception)
    {
        // Make four nodes
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        basic_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        basic_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        basic_nodes.push_back(new Node<2>(3, false, 1.0, 1.0));

        // Make two rectangular element out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        nodes_elem_0.push_back(basic_nodes[0]);
        nodes_elem_0.push_back(basic_nodes[2]);
        nodes_elem_1.push_back(basic_nodes[3]);
        nodes_elem_1.push_back(basic_nodes[1]);

        std::vector<PottsElement<2>*> basic_potts_elements;
        basic_potts_elements.push_back(new PottsElement<2>(0, nodes_elem_0));
        basic_potts_elements.push_back(new PottsElement<2>(1, nodes_elem_1));

        // Make a Potts mesh
        PottsMesh<2> basic_potts_mesh(basic_nodes, basic_potts_elements);

        TS_ASSERT_EQUALS(basic_potts_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(basic_potts_mesh.GetNumNodes(), 4u);

        // Delete element 0
        basic_potts_mesh.DeleteElement(0);

        // Divide remaining element
        unsigned new_element_index = basic_potts_mesh.DivideElement(basic_potts_mesh.GetElement(1));

        // Check index reused
        TS_ASSERT_EQUALS(new_element_index, 0u);

        TS_ASSERT_EQUALS(basic_potts_mesh.GetNumElements(), 2u);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(0)->GetNumNodes(), 1u);
        //TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(0)->GetNodeGlobalIndex(0), 2u);

        TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(1)->GetNumNodes(), 1u);
        //TS_ASSERT_EQUALS(basic_potts_mesh.GetElement(1)->GetNodeGlobalIndex(0), 3u);
    }

    void TestArchive2dPottsMesh()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "potts_mesh_2d.arch";
        ArchiveLocationInfo::SetMeshFilename("potts_mesh");

        // Create mesh
        PottsMeshReader<2> mesh_reader("notforrelease_cell_based/test/data/TestPottsMeshWriter/potts_mesh_2d");
        PottsMesh<2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        AbstractMesh<2,2>* const p_mesh = &mesh;

        /*
         * You need the const above to stop a BOOST_STATIC_ASSERTION failure.
         * This is because the serialization library only allows you to save tracked
         * objects while the compiler considers them const, to prevent the objects
         * changing during the save, and so object tracking leading to wrong results.
         *
         * E.g. A is saved once via pointer, then changed, then saved again. The second
         * save notes that A was saved before, so doesn't write its data again, and the
         * change is lost.
         */

        // Create an output archive
        {
            TS_ASSERT_EQUALS((static_cast<PottsMesh<2>*>(p_mesh))->GetNumNodes(), 6u);
            TS_ASSERT_EQUALS((static_cast<PottsMesh<2>*>(p_mesh))->GetNumElements(), 2u);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // We have to serialize via a pointer here, or the derived class information is lost
            (*p_arch) << p_mesh;
        }

        {
            // De-serialize and compare
            AbstractMesh<2,2>* p_mesh2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_mesh2;

            PottsMesh<2>* p_mesh_original = static_cast<PottsMesh<2>*>(p_mesh);
            PottsMesh<2>* p_mesh_loaded = static_cast<PottsMesh<2>*>(p_mesh2);

            // Compare the loaded mesh against the original

            TS_ASSERT_EQUALS(p_mesh_original->GetNumNodes(), p_mesh_loaded->GetNumNodes());

            for (unsigned node_index=0; node_index<p_mesh_original->GetNumNodes(); node_index++)
            {
                Node<2>* p_node = p_mesh_original->GetNode(node_index);
                Node<2>* p_node2 = p_mesh_loaded->GetNode(node_index);

                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());

                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());

                for (unsigned dimension=0; dimension<2; dimension++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[dimension], p_node2->rGetLocation()[dimension], 1e-4);
                }
            }

            TS_ASSERT_EQUALS(p_mesh_original->GetNumElements(), p_mesh_loaded->GetNumElements());

            for (unsigned elem_index=0; elem_index < p_mesh_original->GetNumElements(); elem_index++)
            {
                TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNumNodes(),
                                 p_mesh_loaded->GetElement(elem_index)->GetNumNodes());

                for (unsigned local_index=0; local_index<p_mesh_original->GetElement(elem_index)->GetNumNodes(); local_index++)
                {
                    TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNodeGlobalIndex(local_index),
                                     p_mesh_loaded->GetElement(elem_index)->GetNodeGlobalIndex(local_index));
                }
            }

            // Tidy up
            delete p_mesh_loaded;
        }
    }
};

#endif /*TESTPOTTSMESH_HPP_*/
