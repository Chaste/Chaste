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

#ifndef TESTMUTABLEMESH_HPP_
#define TESTMUTABLEMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>
#include "MutableMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "RandomNumberGenerator.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include "ArchiveOpener.hpp"
#include <cmath>
#include <vector>

class TestMutableMesh : public CxxTest::TestSuite
{
private:

    template<unsigned DIM>
    void EdgeIteratorTest(std::string meshFilename)
    {
        // Create a simple mesh
        TrianglesMeshReader<DIM,DIM> mesh_reader(meshFilename);
        MutableMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Delete one of the nodes and hence an element
        // to check that iterator skips deleted elements
        // this causes element 0 to be deleted which is a good choice for coverage of the begin method
        mesh.DeleteBoundaryNodeAt(0);

        // Check that we can iterate over the set of edges
        std::set< std::set<unsigned> > edges_visited;

        for (typename MutableMesh<DIM,DIM>::EdgeIterator edge_iterator=mesh.EdgesBegin();
             edge_iterator!=mesh.EdgesEnd();
             ++edge_iterator)
        {
            std::set<unsigned> node_pair;
            node_pair.insert(edge_iterator.GetNodeA()->GetIndex());
            node_pair.insert(edge_iterator.GetNodeB()->GetIndex());

            TS_ASSERT_EQUALS(edges_visited.find(node_pair), edges_visited.end());
            edges_visited.insert(node_pair);
        }

        // Set up expected node pairs
        std::set< std::set<unsigned> > expected_node_pairs;
        for (unsigned i=0; i<mesh.GetNumAllElements(); i++)
        {
            Element<DIM,DIM>* p_element = mesh.GetElement(i);
            if (!p_element->IsDeleted())
            {
                for (unsigned j=0; j<DIM+1; j++)
                {
                    for (unsigned k=0; k<DIM+1; k++)
                    {
                        unsigned node_A = p_element->GetNodeGlobalIndex(j);
                        unsigned node_B = p_element->GetNodeGlobalIndex(k);

                        if (node_A != node_B)
                        {
                            std::set<unsigned> node_pair;
                            node_pair.insert(node_A);
                            node_pair.insert(node_B);

                            expected_node_pairs.insert(node_pair);
                        }
                    }
                }
            }
        }

        TS_ASSERT_EQUALS(edges_visited, expected_node_pairs);
    }

public:

    void TestNodeIterator()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        unsigned counter = 0;

        for (AbstractTetrahedralMesh<2,3>::NodeIterator iter = mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            unsigned node_index = (*iter).GetIndex();
            TS_ASSERT_EQUALS(counter, node_index); // assumes the iterator will give node 0,1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(counter, mesh.GetNumNodes());

        mesh.DeleteNode(0);
        mesh.DeleteNode(10);
        mesh.DeleteNode(70);

        unsigned another_counter = 0;

        for (AbstractTetrahedralMesh<2,3>::NodeIterator iter = mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            another_counter++;
        }

        TS_ASSERT_EQUALS(another_counter+3, counter)
        TS_ASSERT_EQUALS(another_counter+3, mesh.GetNumAllNodes());
        TS_ASSERT_EQUALS(another_counter, mesh.GetNumNodes());
    }

    void TestElementIterator()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        unsigned counter = 0;

        for (AbstractTetrahedralMesh<2,3>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = (*iter).GetIndex();
            TS_ASSERT_EQUALS(counter, element_index); // assumes the iterator will give element 0,1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(counter, mesh.GetNumElements());
    }

    void TestRescaleMeshFromBoundaryNode()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ChastePoint<1> updatedPoint(1.5);
        mesh.RescaleMeshFromBoundaryNode(updatedPoint,10);
        for (int i=0; i<11; i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->GetPoint()[0], 1.5*(i/10.0), 0.001);
        }
    }

    void Test1DSetPoint()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        const int node_index = 3;
        Node<1>* p_node = mesh.GetNode(node_index);

        ChastePoint<1> point=p_node->GetPoint();
        TS_ASSERT_DELTA(point[0], 0.3, 1e-6);

        Element<1,1>* p_element;
        Node<1>::ContainingElementIterator elt_iter = p_node->ContainingElementsBegin();
        p_element = mesh.GetElement(*elt_iter);
        c_matrix<double,1,1> jacobian;
        double jacobian_det;
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.1, 1e-6);

        p_element = mesh.GetElement(*++elt_iter);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.1, 1e-6);

        // Move node 3 from 0.3 (between node 2 at 0.2 and node 4 at 0.4)
        point.SetCoordinate(0, 0.25);
        mesh.SetNode(node_index, point);
        elt_iter = p_node->ContainingElementsBegin();
        p_element = mesh.GetElement(*elt_iter);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.05, 1e-6);

        p_element = mesh.GetElement(*++elt_iter);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.15, 1e-6);

        // Move node 3 from 0.3 (between node 2 at 0.2 and node 4 at 0.4)
        point.SetCoordinate(0, 0.201);
        mesh.SetNode(node_index, point);
        elt_iter = p_node->ContainingElementsBegin();
        p_element = mesh.GetElement(*elt_iter);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.001, 1e-6);

        p_element = mesh.GetElement(*++elt_iter);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.199, 1e-6);

        // Move node 3 so that one element is empty
        point.SetCoordinate(0, 0.200);
        TS_ASSERT_THROWS_THIS(mesh.SetNode(node_index, point),
                "Moving node caused an element to have a non-positive Jacobian determinant");

        // Move node 3 so that one element is negative
        point.SetCoordinate(0, 0.15);
        TS_ASSERT_THROWS_THIS(mesh.SetNode(node_index, point),
                "Moving node caused an element to have a non-positive Jacobian determinant");

        // Move node 3 back (and recover)
        point.SetCoordinate(0, 0.3);
        mesh.SetNode(node_index, point);
        elt_iter = p_node->ContainingElementsBegin();
        p_element = mesh.GetElement(*elt_iter);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.1, 1e-6);

        p_element = mesh.GetElement(*++elt_iter);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.1, 1e-6);

        TS_ASSERT_EQUALS(mesh.IsMeshChanging(), true);
    }

    void Test2DSetPoint()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        const int node_index = 234;
        const int boundary_node_index = 99;
        Node<2>* p_node = mesh.GetNode(node_index);

        // Just focus on one element
        Element<2,2>* p_element;
        Node<2>::ContainingElementIterator elt_iter = p_node->ContainingElementsBegin();
        p_element = mesh.GetElement(*elt_iter);

        ChastePoint<2> point = p_node->GetPoint();
        TS_ASSERT_DELTA(point[0], 0.063497248392600097, 1e-6);
        TS_ASSERT_DELTA(point[1], -0.45483180039309123, 1e-6);

        c_matrix<double,2,2> jacobian;
        double jacobian_det;
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.00907521, 1e-6);

        // Nudge
        point.SetCoordinate(0, 0.06);
        mesh.SetNode(node_index, point);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.00861908, 1e-6);

        // Nudge
        point.SetCoordinate(0, 0.02);
        mesh.SetNode(node_index, point);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.00340215, 1e-6);

        // Nudge
        point.SetCoordinate(0, -0.006);
        mesh.SetNode(node_index, point);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 1.11485e-05, 1e-6);

        // Nudge too far
        point.SetCoordinate(0, -0.0065);
        TS_ASSERT_THROWS_THIS(mesh.SetNode(node_index, point),
                "Moving node caused an element to have a non-positive Jacobian determinant");

        // Put it back
        point.SetCoordinate(0, 0.063497248392600097);
        mesh.SetNode(node_index, point);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.00907521, 1e-6);

        // Now try to move a boundary node
        p_node=mesh.GetNode(boundary_node_index);
        ChastePoint<2> boundary_point = p_node->GetPoint();
        TS_ASSERT_DELTA(boundary_point[0], 0.99211470130000001, 1e-6);
        TS_ASSERT_DELTA(boundary_point[1], -0.12533323360000001, 1e-6);

        Node<2>::ContainingBoundaryElementIterator b_elt_iter = p_node->ContainingBoundaryElementsBegin();
        const BoundaryElement<1,2>* p_boundary_element = mesh.GetBoundaryElement(*b_elt_iter);

        c_vector<double,2> weighted_dir;
        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0628215, 1e-6);

        boundary_point.SetCoordinate(0, 1.0);
        mesh.SetNode(boundary_node_index, boundary_point);
        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0645268, 1e-6);
    }

    void TestMovingNodesIn3D()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double reference_volume = mesh.GetVolume();

        const int interior_node_index = 34;
        Node<3>* p_node = mesh.GetNode(interior_node_index);

        // Just focus on one element
        Element<3,3>* p_element = mesh.GetElement(*p_node->ContainingElementsBegin());
        BoundaryElement<2,3>* p_boundary_element = mesh.GetBoundaryElement(*p_node->ContainingBoundaryElementsBegin());

        ChastePoint<3> point = p_node->GetPoint();
        TS_ASSERT_DELTA(point[0], 1, 1e-6);
        TS_ASSERT_DELTA(point[1], 0.75, 1e-6);
        TS_ASSERT_DELTA(point[2], 0.75, 1e-6);

        c_matrix<double,3,3> jacobian;
        c_vector<double,3> weighted_dir;
        double jacobian_det;

        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.03125, 1e-6);

        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.125, 1e-6);

        // Check the mesh volume hasn't changed
        TS_ASSERT_DELTA(mesh.GetVolume(), reference_volume, 1e-6);

        // Nudge
        point.SetCoordinate(2, 0.9);
        mesh.SetNode(interior_node_index, point);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0125, 1e-6);

        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.05, 1e-6);

        // Nudge
        point.SetCoordinate(2, 0.999);
        mesh.SetNode(interior_node_index, point);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.000125, 1e-6);

        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0005, 1e-6);

        // Nudge
        point.SetCoordinate(2, 0.99999);
        mesh.SetNode(interior_node_index, point);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 1.25e-06, 1e-6);

        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 5.0e-06, 1e-6);

        // Nudge too far
        point.SetCoordinate(2, 1.0);
        TS_ASSERT_THROWS_THIS(mesh.SetNode(interior_node_index, point),
                "Moving node caused a boundary element to have a non-positive Jacobian determinant");

        // Put it back
        point.SetCoordinate(2, 0.75);
        mesh.SetNode(interior_node_index, point);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.03125, 1e-6);

        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.125, 1e-6);

        // Find exterior node
        const int exterior_node_index = 0;
        p_node = mesh.GetNode(exterior_node_index); // this exterior node is at (0,0,0)
        point = p_node->GetPoint();
        TS_ASSERT_DELTA(point[0], 0, 1e-6);
        TS_ASSERT_DELTA(point[1], 0, 1e-6);
        TS_ASSERT_DELTA(point[2], 0, 1e-6);

        // Move exterior node
        point.SetCoordinate(2, -10.0);
        mesh.SetNode(exterior_node_index, point);

        // Check mesh volume has changed
        TS_ASSERT(fabs(mesh.GetVolume() - reference_volume) > 1e-1);
    }

    void Test1DMeshIn2DSetPoint()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/semicircle_outline");
        MutableMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        const int boundary_node_index = 50;
        Node<2>* p_node = mesh.GetNode(boundary_node_index);

        // Just focus on one element
        Element<1,2>* p_element = mesh.GetElement(*p_node->ContainingElementsBegin());
        BoundaryElement<0,2>* p_boundary_element = mesh.GetBoundaryElement(*p_node->ContainingBoundaryElementsBegin());

        ChastePoint<2> point = p_node->GetPoint();
        TS_ASSERT_DELTA(point[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(point[1], 0.0, 1e-6);

        c_vector<double,2> weighted_dir;
        double jacobian_det;

        mesh.GetWeightedDirectionForElement(p_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0628215, 1e-6);

        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 1.0, 1e-6);

        // Nudge left
        point.SetCoordinate(0, -1.5);
        mesh.SetNode(boundary_node_index, point);
        mesh.GetWeightedDirectionForElement(p_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.505885, 1e-6);

        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 1.0, 1e-6);

        // Can't nudge right since an element flips chirality
        point.SetCoordinate(0, -0.5);
        TS_ASSERT_THROWS_THIS(mesh.SetNode(boundary_node_index, point),
                "Moving node caused an subspace element to change direction");

        // Put it back
        point.SetCoordinate(0, -1.0);
        mesh.SetNode(boundary_node_index, point);
        mesh.GetWeightedDirectionForElement(p_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0628215, 1e-6);

        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 1.0, 1e-6);
    }


    void Test2DMeshIn3DSetPoint()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        const int boundary_node_index = 99;
        Node<3>* p_node = mesh.GetNode(boundary_node_index);

        //Just focus on one element
        Element<2,3>* p_element = mesh.GetElement(*p_node->ContainingElementsBegin());
        BoundaryElement<1,3>* p_boundary_element = mesh.GetBoundaryElement(*p_node->ContainingBoundaryElementsBegin());

        ChastePoint<3> point = p_node->GetPoint();
        TS_ASSERT_DELTA(point[0], 0.99211470130000001, 1e-6);
        TS_ASSERT_DELTA(point[1], -0.12533323360000001, 1e-6);
        TS_ASSERT_DELTA(point[2], 0.0, 1e-6);

        c_vector<double,3> weighted_dir;
        double jacobian_det;

        mesh.GetWeightedDirectionForElement(p_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0163772, 1e-6);
        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0628215, 1e-6);

        // Nudge above the plane
        point.SetCoordinate(2, 1e-2);
        mesh.SetNode(boundary_node_index, point);

        mesh.GetWeightedDirectionForElement(p_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0164274, 1e-6);
        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0636124, 1e-6);

        // Nudge it back
        point.SetCoordinate(2, 0.0);
        mesh.SetNode(boundary_node_index, point);

        mesh.GetWeightedDirectionForElement(p_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0163772, 1e-6);
        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0628215, 1e-6);

        // Nudge below the plane
        point.SetCoordinate(2, -1e-2);
        mesh.SetNode(boundary_node_index, point);

        mesh.GetWeightedDirectionForElement(p_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0164274, 1e-6);
        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0636124, 1e-6);

        // Put it back
        point.SetCoordinate(2, 0.0);
        mesh.SetNode(boundary_node_index, point);

        mesh.GetWeightedDirectionForElement(p_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0163772, 1e-6);
        mesh.GetWeightedDirectionForBoundaryElement(p_boundary_element->GetIndex(), weighted_dir, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0628215, 1e-6);

        // Can't nudge to the other side of the circle without changing handedness
        point.SetCoordinate(0, -1.0);
        point.SetCoordinate(2, 0.0);
        TS_ASSERT_THROWS_THIS(mesh.SetNode(boundary_node_index, point),
                "Moving node caused an subspace element to change direction");
    }

    void TestDeletingNodes()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        Node<1>* p_old_rhs_node = mesh.GetNode(10);
        Node<1>* p_old_lhs_node = mesh.GetNode(0);

        MutableMesh<1,1>::BoundaryElementIterator b_elt_iter;
        MutableMesh<1,1>::BoundaryNodeIterator b_node_iter;

        // Delete the right end node
        mesh.DeleteBoundaryNodeAt(10);

        TS_ASSERT(p_old_rhs_node->IsDeleted());

        // Number of *all* nodes & elements should be unchanged, even though we've deleted one,
        // since this is the size of the vector, not the number of active nodes/elements.
        // Yes, it's confusing.
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 11u);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 10u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 9u);

        // Check the boundary lists are correct
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2u);
        b_node_iter = mesh.GetBoundaryNodeIteratorBegin();
        TS_ASSERT_EQUALS((*b_node_iter++)->GetIndex(), 0u);
        TS_ASSERT_EQUALS((*b_node_iter++)->GetIndex(), 9u);

        // NB: New boundary elements not added
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 1u);
        b_elt_iter = mesh.GetBoundaryElementIteratorBegin();
        TS_ASSERT_EQUALS((*b_elt_iter)->GetNumNodes(), 1u);
        TS_ASSERT_EQUALS((*b_elt_iter++)->GetNode(0)->GetIndex(), 0u);

        // Check the new boundary node
        Node<1>* p_new_rhs_node = mesh.GetNode(9);
        TS_ASSERT(p_new_rhs_node->IsBoundaryNode());
        TS_ASSERT_EQUALS(p_new_rhs_node->GetNumContainingElements(), 1u);

        // Only allowed to remove boundary nodes
        TS_ASSERT_THROWS_THIS(mesh.DeleteBoundaryNodeAt(5)," You may only delete a boundary node ");

        // Delete the left end node
        mesh.DeleteBoundaryNodeAt(0);

        TS_ASSERT(p_old_lhs_node->IsDeleted());
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 11u);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 10u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 8u);

        // Check the boundary lists are correct
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2u);
        b_node_iter = mesh.GetBoundaryNodeIteratorBegin();
        TS_ASSERT_EQUALS((*b_node_iter++)->GetIndex(), 9u); // Note that the boundary is now
        TS_ASSERT_EQUALS((*b_node_iter++)->GetIndex(), 1u); // 'reversed'

        // NB: New boundary elements not added
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 0u);

        // Check the new boundary node
        Node<1>* p_new_lhs_node = mesh.GetNode(1);
        TS_ASSERT(p_new_lhs_node->IsBoundaryNode());
        TS_ASSERT_EQUALS(p_new_lhs_node->GetNumContainingElements(), 1u);

        // Check the deleted element/node vectors
    }

    void TestDeleteNodePriorToReMesh()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/circular_fan");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Test it can also delete a boundary node
        mesh.DeleteNodePriorToReMesh(0);
        mesh.DeleteNodePriorToReMesh(11);

        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 100u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 98u);

        NodeMap map(mesh.GetNumNodes());
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 98u);
    }

    void TestAddingAndDeletingNodes()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Add a node at position 0.01
        ChastePoint<1> new_point(0.01);
        Element<1,1>* p_first_element = mesh.GetElement(0);

        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_first_element, new_point));
        TS_ASSERT_EQUALS(p_first_element->GetNode(1)->GetIndex(), 11u);

        // Check the new element is index 10, by comparing nodes
        TS_ASSERT_EQUALS(mesh.GetElement(10)->GetNode(0), p_first_element->GetNode(1));

        // Delete the last node
        mesh.DeleteBoundaryNodeAt(10);

        // Add a node
        ChastePoint<1> new_point2(0.55);
        Element<1,1>* p_sixth_element = mesh.GetElement(5);

        mesh.RefineElement(p_sixth_element, new_point2);

        TS_ASSERT_EQUALS(p_sixth_element->GetNode(1)->GetIndex(), 10u);

        // Check the new element is index 9, by comparing nodes
        TS_ASSERT_EQUALS(mesh.GetElement(9)->GetNode(0), p_sixth_element->GetNode(1));
    }

    void Test1DBoundaryNodeMerger()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_1_element");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetVolume(), 1.0);
        TS_ASSERT_THROWS_THIS(mesh.MoveMergeNode(0, 1),
                "These nodes cannot be merged since they are not neighbours on the boundary");
    }

    void Test1DNodeMerger()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double length = mesh.GetVolume();
        const int node_index = 3;
        const int target_index = 4;
        const int not_neighbour_index = 5;

        // Cannot merge node 3 with node 5 since they are not neighbours
        TS_ASSERT_THROWS_THIS(mesh.MoveMergeNode(node_index, not_neighbour_index),
                "These nodes cannot be merged since they are not neighbours");

        // Merge node 3 with node 4
        mesh.MoveMergeNode(node_index, target_index);

        Element<1,1>* p_element;
        p_element = mesh.GetElement(2);

        c_matrix<double,1,1> jacobian;
        double jacobian_det;
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.2, 1e-6);
        p_element = mesh.GetElement(3);
        mesh.GetJacobianForElement(p_element->GetIndex(), jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0, 1e-6);

        TS_ASSERT_DELTA(length, mesh.GetVolume(), 1e-6);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements() + 1);
    }

    void Test2DNodeMerger()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double area = mesh.GetVolume();

        // Node 432 is supported by a non-convex region
        // Node 206 is the sole reflex vertex in the region
        // Node 172 is not feasible since it is neighbour to the reflex
        const int node_index = 432;
        const int target_index = 206;
        const int not_neighbour_index = 204;
        const int not_feasible_index = 172;

        // Element 309 is shared by the moving node (432), the reflex node (206)
        // the non-feasible node (172) - it will vanish
        // Element 762 is shared by the moving node (432), some other node (205)
        // the non-feasible node (172) - it will increase in size
        c_matrix<double,2,2> jacobian;
        double jacobian_det;
        mesh.GetJacobianForElement(309, jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.00753493, 1e-6);

        mesh.GetJacobianForElement(762, jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.00825652, 1e-6);

        // Cannot merge since they are not neighbours
        TS_ASSERT_THROWS_THIS(mesh.MoveMergeNode(node_index, not_neighbour_index),
                "These nodes cannot be merged since they are not neighbours");

        // Cannot merge since an element goes negative
        // The element 763 shared by moving node (432), reflex node (206) and the
        // other neighbour to the reflex node goes negative
        TS_ASSERT_THROWS_THIS(mesh.MoveMergeNode(node_index, not_feasible_index, false),
                "Moving node caused an element to have a non-positive Jacobian determinant");

        // Added "crossReference=false" to stop elements deregistering
        mesh.MoveMergeNode(node_index, target_index);

        TS_ASSERT_DELTA(area, mesh.GetVolume(), 1e-6);
        mesh.GetJacobianForElement(309, jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0, 1e-6);

        mesh.GetJacobianForElement(762, jacobian, jacobian_det);
        TS_ASSERT_DELTA(jacobian_det, 0.0126728, 1e-6);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements() + 2);
    }

    void Test3DNodeMerger()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double volume = mesh.GetVolume();
        const int node_index = 22; // in the middle
        const int target_index = 310;
        const int not_neighbour_index = 204;
        const int not_feasible_index = 103;

        TS_ASSERT_THROWS_THIS(mesh.MoveMergeNode(node_index, not_neighbour_index),
                "These nodes cannot be merged since they are not neighbours");
        TS_ASSERT_THROWS_THIS(mesh.MoveMergeNode(node_index, not_feasible_index, false),
                "Moving node caused an element to have a non-positive Jacobian determinant");
        // Added "crossReference=false" to stop elements deregistering

        TS_ASSERT_THROWS_NOTHING( mesh.MoveMergeNode(node_index, target_index));
        TS_ASSERT_DELTA(volume, mesh.GetVolume(), 1e-6);

        //Ten elements share 22 and 310.  See:
        /*   Element      N1    N2    N3    N4
                 510     310   348    22   294
                 645      22   328   310   216
                 753      22   329   310   120
                1164     295   310    22   175
                1217     294   310   175    22
                1251     310   336   328    22
                1254     120   310    22   295
                1357     310   336    22   193
                1365      22   329   216   310
                1484     310   348   193    22
        */

        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements() + 10);
    }

    void Test2DBoundaryNodeMergerChangeArea()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double area = mesh.GetVolume();
        double perim = mesh.GetSurfaceArea();
        unsigned num_nodes = mesh.GetNumNodes();
        unsigned num_elements = mesh.GetNumElements();
        unsigned num_boundary_elements = mesh.GetNumBoundaryElements();
        const int node_index = 19;
        const int target_index = 20;
        const int not_boundary_index = 400;

        TS_ASSERT_THROWS_THIS(mesh.MoveMergeNode(node_index, not_boundary_index),"A boundary node can only be moved on to another boundary node");
        mesh.MoveMergeNode(node_index, target_index);

        TS_ASSERT_DELTA(area - mesh.GetVolume(), 1.24e-4, 1e-6);
        TS_ASSERT_DELTA(perim - mesh.GetSurfaceArea(), 6.20e-5, 1e-7);
        TS_ASSERT_EQUALS(num_nodes-mesh.GetNumNodes(), 1u);
        TS_ASSERT_EQUALS(num_elements-mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(num_boundary_elements-mesh.GetNumBoundaryElements(), 1u);
    }

    void Test2DBoundaryNodeMerger()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_800_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double area = mesh.GetVolume();
        double perim = mesh.GetSurfaceArea();
        int num_nodes = mesh.GetNumNodes();
        unsigned num_elements = mesh.GetNumElements();
        unsigned num_boundary_elements = mesh.GetNumBoundaryElements();
        const int node_index = 9;
        const int target_index = 10;
        const int not_boundary_index = 31;

        TS_ASSERT_THROWS_THIS(mesh.MoveMergeNode(node_index, not_boundary_index),"A boundary node can only be moved on to another boundary node");
        mesh.MoveMergeNode(node_index, target_index);

        TS_ASSERT_DELTA(area - mesh.GetVolume(), 0.00, 1e-6);
        TS_ASSERT_DELTA(perim - mesh.GetSurfaceArea(), 0.00, 1e-7);
        TS_ASSERT_EQUALS(num_nodes-mesh.GetNumNodes(), 1u);
        TS_ASSERT_EQUALS(num_elements-mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(num_boundary_elements-mesh.GetNumBoundaryElements(), 1u);

        const int corner_index = 20;
        const int corner_target_index = 19;
        mesh.MoveMergeNode(corner_index, corner_target_index);

        TS_ASSERT_DELTA(area - mesh.GetVolume(), 1.25e-5, 1e-7);
        TS_ASSERT_DELTA(perim - mesh.GetSurfaceArea(), 2.92893e-3, 1e-7);
        TS_ASSERT_EQUALS(num_nodes-mesh.GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(num_elements-mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(num_boundary_elements-mesh.GetNumBoundaryElements(), 2u);
    }

    void Test3DBoundaryNodeMerger()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double volume = mesh.GetVolume();
        double surface = mesh.GetSurfaceArea();
        unsigned num_nodes = mesh.GetNumNodes();
        unsigned num_elements = mesh.GetNumElements();
        unsigned num_boundary_elements = mesh.GetNumBoundaryElements();
        const int node_index = 147;
        const int target_index = 9;
        //const int not_boundary_index=400;

        mesh.MoveMergeNode(node_index, target_index);

        TS_ASSERT_DELTA(volume, mesh.GetVolume(), 1e-7);
        TS_ASSERT_DELTA(surface, mesh.GetSurfaceArea(), 1e-7);
        TS_ASSERT_EQUALS(num_nodes-mesh.GetNumNodes(), 1u);
        TS_ASSERT_EQUALS(num_elements-mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(num_boundary_elements-mesh.GetNumBoundaryElements(), 2u);

        // Can't move corner nodes since this forces some zero volume elements which aren't on the shared list
    }

    void TestCheckVoronoiDisk()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.CheckIsVoronoi(), true);
    }

    void TestCheckVoronoiSquare()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.CheckIsVoronoi(), true);
    }

    void TestCheckCircularFan()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/circular_fan");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.CheckIsVoronoi(5e-3), true);
    }

    void TestCheckMovingMesh()
    {
        MutableMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(1,1,false);

        unsigned top_right=3;
        Node<2>* p_node = mesh.GetNode(top_right);
        ChastePoint<2> point = p_node->GetPoint();
        TS_ASSERT_DELTA(point[0], 1.0, 1e-7);
        TS_ASSERT_DELTA(point[1], 1.0, 1e-7);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 1u);

        for (double x=1.1; x>=0.9; x-=0.01)
        {
            point.SetCoordinate(0,x);
            point.SetCoordinate(1,x);
            mesh.SetNode(top_right, point);

            if (x >= 0.94)
            {
                TS_ASSERT_EQUALS(mesh.CheckIsVoronoi(0.15), true);
            }
            else
            {
                TS_ASSERT_EQUALS(mesh.CheckIsVoronoi(0.15), false);
            }
        }
    }

    void TestDeleteNodes()
    {
        MutableMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2,3);

        TS_ASSERT_EQUALS(mesh.GetVolume(), 6.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 12u);

        // Delete from interior
        mesh.DeleteNode(7);
        TS_ASSERT_EQUALS(mesh.GetVolume(), 6.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 11u);

        // Delete from edge
        mesh.DeleteNode(5);
        TS_ASSERT_EQUALS(mesh.GetVolume(), 6.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 10u);

        // Delete from corner
        mesh.DeleteNode(2);
        TS_ASSERT_EQUALS(mesh.GetVolume(), 5.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9u);

        // Deleting a deleted node should throw an exception
        TS_ASSERT_THROWS_THIS(mesh.DeleteNode(2),"Trying to delete a deleted node");

        // Moving a deleted node should throw an exception
        TS_ASSERT_THROWS_THIS(mesh.MoveMergeNode(2,1),"Trying to move a deleted node");
    }

    void TestDeleteNodeFails()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/HalfSquareWithExtraNode");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_THROWS_THIS(mesh.DeleteNode(0),"Failure to delete node");
    }

    void TestMeshAddNodeAndReMeshMethod()
    {
        MutableMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(1, 1, false);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);

        TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetLocation()[0], 0.0, 1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetLocation()[1], 1.0, 1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(3u)->rGetLocation()[0], 1.0, 1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(3u)->rGetLocation()[1], 1.0, 1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetLocation()[0], 0.0, 1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetLocation()[1], 0.0, 1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetLocation()[0], 1.0, 1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetLocation()[1], 0.0, 1e-7);

        // Test the add node method
        c_vector<double,2> point;
        point[0] = 2.0;
        point[1] = 0.0;
        Node<2>* p_node = new Node<2>(4u, point);
        unsigned new_index = mesh.AddNode(p_node);

        TS_ASSERT_EQUALS(new_index, 4u);
        TS_ASSERT_DELTA(mesh.GetNode(4u)->rGetLocation()[0], 2.0, 1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(4u)->rGetLocation()[1], 0.0, 1e-7);

        // Test the add node and ReMesh method

        point[0] = 2.0;
        point[1] = 1.0;
        Node<2>* p_node2 = new Node<2>(5u, point);

        NodeMap map(mesh.GetNumNodes());
        new_index = mesh.AddNode(p_node2);
        mesh.ReMesh(); // call the version of ReMesh which doesn't use a map
        TS_ASSERT_EQUALS(new_index, 5u);
        TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetLocation()[0], 2.0, 1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetLocation()[1], 1.0, 1e-7);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 4u);
    }

    void TestReindex()
    {
        MutableMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(10, 10);

        unsigned num_old_nodes = mesh.GetNumNodes();

        mesh.DeleteNode(50);
        mesh.DeleteNode(0);

        NodeMap map(num_old_nodes);

        mesh.ReIndex(map);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), (unsigned)(num_old_nodes-2));
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), (unsigned)(num_old_nodes-2));

        TS_ASSERT_EQUALS(map.GetSize(), num_old_nodes);
        TS_ASSERT_EQUALS(map.IsDeleted(50), true);
        TS_ASSERT_EQUALS(map.IsDeleted(0), true);

        for (unsigned i=1; i<num_old_nodes; i++)
        {
            if (i<50)
            {
                TS_ASSERT_EQUALS(map.GetNewIndex(i), (unsigned)(i-1));
            }
            if (i>50)
            {
                TS_ASSERT_EQUALS(map.GetNewIndex(i), (unsigned)(i-2));
            }
        }
    }

    void TestReindex2dMeshIn3dSpace()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        unsigned num_old_nodes = mesh.GetNumNodes();

        mesh.DeleteNode(50);
        mesh.DeleteNode(0);

        NodeMap map(num_old_nodes);

        mesh.ReIndex(map);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), (unsigned)(num_old_nodes-2));
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), (unsigned)(num_old_nodes-2));

        TS_ASSERT_EQUALS(map.GetSize(), num_old_nodes);
        TS_ASSERT_EQUALS(map.IsDeleted(50), true);
        TS_ASSERT_EQUALS(map.IsDeleted(0), true);

        for (unsigned i=1; i<num_old_nodes; i++)
        {
            if (i<50)
            {
                TS_ASSERT_EQUALS(map.GetNewIndex(i), (unsigned)(i-1));
            }
            if (i>50)
            {
                TS_ASSERT_EQUALS(map.GetNewIndex(i), (unsigned)(i-2));
            }
        }
    }

    void TestConstructFromNodes()
    {
        // Create mutable tetrahedral mesh which is Delaunay
        std::vector<Node<3> *> nodes;

        nodes.push_back(new Node<3>(0, true,  0.0,  0.0,  0.0));
        nodes.push_back(new Node<3>(1, true,  1.0,  1.0,  0.0));
        nodes.push_back(new Node<3>(2, true,  1.0,  0.0,  1.0));
        nodes.push_back(new Node<3>(3, true,  0.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(4, false, 0.5,  0.5,  0.5));

        MutableMesh<3,3> mesh(nodes);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 4u);
        TS_ASSERT_DELTA(mesh.GetVolume(), 0.3333, 1e-4);
    }

    void TestEdgeIterator()
    {
        EdgeIteratorTest<3>("mesh/test/data/cube_2mm_12_elements");
        EdgeIteratorTest<2>("mesh/test/data/square_4_elements");
        EdgeIteratorTest<1>("mesh/test/data/1D_0_to_1_10_elements");
    }


    void TestDeleteElement()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        unsigned num_old_nodes = mesh.GetNumNodes();
        unsigned num_old_eles  = mesh.GetNumElements();

        mesh.DeleteElement(18);
        mesh.DeleteElement(0);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), (unsigned)(num_old_eles-2));
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), (unsigned)(num_old_eles));

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), (unsigned)(num_old_nodes));

        //Elements 18, 71, 74, 90, 380, 381 are the only elements using Node 241
        mesh.DeleteElement(71);
        mesh.DeleteElement(74);
        mesh.DeleteElement(90);
        mesh.DeleteElement(380);
        mesh.DeleteElement(381);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), (unsigned)(num_old_nodes - 1));
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), (unsigned)(num_old_nodes));

        NodeMap map(num_old_nodes);

        mesh.ReIndex(map);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), (unsigned)(num_old_nodes - 1));
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), (unsigned)(num_old_nodes - 1));
        TS_ASSERT_EQUALS(mesh.GetNumElements(), (unsigned)(num_old_eles - 7));
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), (unsigned)(num_old_eles - 7));

        TS_ASSERT_EQUALS(map.GetSize(), num_old_nodes);
        TS_ASSERT_EQUALS(map.IsDeleted(241), true);
    }

    void TestDeleteElement1DIn3D()
    {
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/y_branch_3d_mesh");
        MutableMesh<1,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 3u);

        mesh.DeleteElement(1u);

        NodeMap node_map(mesh.GetNumAllNodes());
        mesh.ReIndex(node_map);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 2u);

        //Check that nodes have been shuffled to fill in the indexing
        TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetLocation()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetLocation()[1], -1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetLocation()[2], 0.0, 1e-6);
    }

    void TestArchiving()
    {
        FileFinder archive_dir("archive_mutable_mesh", RelativeTo::ChasteTestOutput);
        std::string archive_file = "mutable_mesh.arch";

        unsigned num_nodes;
        unsigned num_elements;
        unsigned total_num_nodes;
        unsigned total_num_elements;
        std::vector<c_vector<double,2> > node_locations;
        std::vector<c_vector<unsigned,3> > element_nodes;

        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

            AbstractTetrahedralMesh<2,2>* const p_mesh = new MutableMesh<2,2>;
            p_mesh->ConstructFromMeshReader(mesh_reader);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            NodeMap map(p_mesh->GetNumNodes());
            static_cast<MutableMesh<2,2>* >(p_mesh)->ReMesh(map);

            // Record values to test
            for (MutableMesh<2,2>::NodeIterator it = static_cast<MutableMesh<2,2>* >(p_mesh)->GetNodeIteratorBegin();
                it != static_cast<MutableMesh<2,2>* >(p_mesh)->GetNodeIteratorEnd();
                ++it)
            {
                // Mess with the locations of the nodes before archiving.
                // checks that a modified version of the mesh is written to archive.
                (it->rGetModifiableLocation())[0] +=  0.1;
                (it->rGetModifiableLocation())[1] += -0.1;
                node_locations.push_back(it->rGetLocation());
            }

            for (MutableMesh<2,2>::ElementIterator it2 = static_cast<MutableMesh<2,2>* >(p_mesh)->GetElementIteratorBegin();
                it2 != static_cast<MutableMesh<2,2>* >(p_mesh)->GetElementIteratorEnd();
                ++it2)
            {
                c_vector<unsigned, 3> nodes;
                for (unsigned i=0; i<3; i++)
                {
                    nodes[i]=it2->GetNodeGlobalIndex(i);
                }
                element_nodes.push_back(nodes);
            }
            num_nodes = p_mesh->GetNumNodes();
            num_elements = p_mesh->GetNumElements();
            total_num_nodes = p_mesh->GetNumAllNodes();
            total_num_elements = p_mesh->GetNumAllElements();

            (*p_arch) << p_mesh;
            delete p_mesh;
        }

        {
            // Should archive the most abstract class you can to check boost knows what individual classes are.
            // (but here AbstractMesh doesn't have the methods below).
            AbstractTetrahedralMesh<2,2>* p_mesh2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // restore from the archive
            (*p_arch) >> p_mesh2;

            // Check we have the right number of nodes & elements
            TS_ASSERT_EQUALS(num_nodes, p_mesh2->GetNumNodes());
            TS_ASSERT_EQUALS(num_elements, p_mesh2->GetNumElements());
            TS_ASSERT_EQUALS(total_num_nodes, p_mesh2->GetNumAllNodes());
            TS_ASSERT_EQUALS(total_num_elements, p_mesh2->GetNumAllElements());

            // Test recorded node locations
            unsigned counter = 0;
            for (MutableMesh<2,2>::NodeIterator it = static_cast<MutableMesh<2,2>* >(p_mesh2)->GetNodeIteratorBegin();
                it != static_cast<MutableMesh<2,2>* >(p_mesh2)->GetNodeIteratorEnd();
                ++it)
            {
                TS_ASSERT_DELTA(node_locations[counter][0], it->rGetLocation()[0], 1e-9);
                TS_ASSERT_DELTA(node_locations[counter][1], it->rGetLocation()[1], 1e-9);
                counter++;
            }
            counter = 0;
            // Test each element contains the correct nodes.
            for (MutableMesh<2,2>::ElementIterator it2 = static_cast<MutableMesh<2,2>* >(p_mesh2)->GetElementIteratorBegin();
                            it2 != static_cast<MutableMesh<2,2>* >(p_mesh2)->GetElementIteratorEnd();
                            ++it2)
            {
                for (unsigned i=0; i<3; i++)
                {
                    TS_ASSERT_EQUALS(element_nodes[counter][i],it2->GetNodeGlobalIndex(i));
                }
                counter++;
            }

            delete p_mesh2;
        }
    }
};

#endif /*TESTMUTABLEMESH_HPP_*/
