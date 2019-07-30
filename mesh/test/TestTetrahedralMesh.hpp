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

#ifndef _TESTTETRAHEDRALMESH_HPP_
#define _TESTTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include <fstream>
#include <cmath>
#include <vector>
#include <boost/scoped_array.hpp>
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "GenericMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "RandomNumberGenerator.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include "CuboidMeshConstructor.hpp"
#include "ArchiveOpener.hpp"

class TestTetrahedralMesh : public CxxTest::TestSuite
{
public:

    void TestNodeIterator()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        unsigned counter = 0;

        for (AbstractTetrahedralMesh<2,2>::NodeIterator iter = mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            unsigned node_index = (*iter).GetIndex();
            TS_ASSERT_EQUALS(counter, node_index); // assumes the iterator will give node 0,1..,N in that order
            counter++;
        }

        // For coverage, test with an empty mesh
        TetrahedralMesh<2,2> empty_mesh;

        // Since the mesh is empty, the iterator should be set to mrMesh.mNodes.end() when constructed
        AbstractTetrahedralMesh<2,2>::NodeIterator iter = empty_mesh.GetNodeIteratorBegin();

        // We only have a NOT-equals operator defined on the iterator
        TS_ASSERT( !(iter != empty_mesh.GetNodeIteratorEnd()) );
    }

    void TestElementIterator()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        unsigned counter = 0;

        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = (*iter).GetIndex();
            TS_ASSERT_EQUALS(counter, element_index); // assumes the iterator will give element 0,1..,N in that order
            counter++;
        }

        // For coverage, test with an empty mesh
        TetrahedralMesh<2,2> empty_mesh;

        // Since the mesh is empty, the iterator should be set to mrMesh.mNodes.end() when constructed
        AbstractTetrahedralMesh<2,2>::ElementIterator iter = empty_mesh.GetElementIteratorBegin();

        // We only have a NOT-equals operator defined on the iterator
        TS_ASSERT( !(iter != empty_mesh.GetElementIteratorEnd()) );
    }

    void TestMeshConstructionFromMeshReader()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh.GetNumVertices());
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984u);
        TS_ASSERT_EQUALS(mesh.GetNumLocalElements(), 984u);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0],  0.9980267283, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], -0.0627905195, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 0.0, 1e-6);

        TS_ASSERT_DELTA(mesh.GetNodeOrHaloNode(0)->GetPoint()[0],  0.9980267283, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNodeOrHaloNode(0)->GetPoint()[1], -0.0627905195, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNodeOrHaloNode(1)->GetPoint()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNodeOrHaloNode(1)->GetPoint()[1], 0.0, 1e-6);

        // Check first element has the right nodes
        TetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(0), 309u);
        TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(1), 144u);
        TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(2), 310u);
        TS_ASSERT_EQUALS(iter->GetNode(1), mesh.GetNode(144));

        TS_ASSERT_EQUALS(mesh.IsMeshChanging(), false);
        TS_ASSERT_EQUALS(mesh.CalculateMaximumContainingElementsPerProcess(), 8U);
        TS_ASSERT_EQUALS(mesh.CalculateMaximumNodeConnectivityPerProcess(),  9U);

        // Check that there are no halo nodes (coverage)
        std::vector<unsigned> halo_indices;
        mesh.GetHaloNodeIndices(halo_indices);
        unsigned num_halo_nodes = halo_indices.size();
        TS_ASSERT_EQUALS( num_halo_nodes, 0u );
    }

    void TestMeshMoreStatisticsHeart()
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh.CalculateMaximumContainingElementsPerProcess(), 22U);
        TS_ASSERT_EQUALS(mesh.CalculateMaximumNodeConnectivityPerProcess(),  15U);
    }

    void TestMeshMoreStatisticsTumour()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/FromTumourSpheroid");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh.CalculateMaximumContainingElementsPerProcess(), 28U);
        TS_ASSERT_EQUALS(mesh.CalculateMaximumNodeConnectivityPerProcess(),  18U);
    }

    void TestMeshStatisticsSimple()
    {
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(1);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(mesh.CalculateMaximumContainingElementsPerProcess(), 1U);
        TS_ASSERT_EQUALS(mesh.CalculateMaximumNodeConnectivityPerProcess(),  2U);

    }
    void TestMeshConstructionFromMeshReaderIndexedFromOne()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements_indexed_from_1");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984u);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 1.0, 1e-6); // note this mesh is different to disk_984_elements
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 0.9980267283, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 0.0627905195, 1e-6);

        // Check first element has the right nodes
        TetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(0), 309u);
        TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(1), 144u);
        TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(2), 310u);
        TS_ASSERT_EQUALS(iter->GetNode(1), mesh.GetNode(144));
    }

    void Test3dMeshConstructionFromMeshReader()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 51u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 136u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 96u);

        TetrahedralMesh<3,3> mesh;

        try
        {
            mesh.ConstructFromMeshReader(mesh_reader);
        }
        catch (Exception &e)
        {
            std::cout << e.GetMessage() << std::endl;
        }

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 51u);
        //TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 136u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 96u);
        TS_ASSERT_DELTA(mesh.GetVolume(), 1.0, 1e-15);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 6.0, 1e-16);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[0], 0.75, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[1], 0.25, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[2], 0.0, 1e-6);
        TS_ASSERT_EQUALS(mesh.GetNumNodeAttributes(), 0u);
    }

    void TestConstructionFromMeshReaderWithNodeAttributes()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements_with_node_attributes");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 12u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 20u);

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 12u);

        // Check all nodes have 2 attributes
        for (unsigned node_index = 0; node_index < mesh.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(mesh.GetNode(node_index)->rGetNodeAttributes().size(), 2u);
        }
        TS_ASSERT_EQUALS(mesh.GetNumNodeAttributes(), 2u);
        // Check some values
        unsigned probe_node_1 = 0u;
        unsigned probe_node_2 = 8u;

        TS_ASSERT_DELTA(mesh.GetNode(probe_node_1)->rGetNodeAttributes()[0u], 25.2, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(probe_node_1)->rGetNodeAttributes()[1u], 16.3, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(probe_node_2)->rGetNodeAttributes()[0u], 3.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(probe_node_2)->rGetNodeAttributes()[1u], 24.5, 1e-6);
    }

    void Test3dMeshConstructionFromMeshReader2()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_.5mm_1889_elements_irregular");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 425u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 1889u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 436u);

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.GetVolume(), 1.25e-4, 1e-16);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 0.015, 1e-15);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 425u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1889u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 436u);
    }


    void TestMeshWithBoundaryElements()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check for the right number of boundary edges
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 100u);

        // Check all boundary elements have nodes on the boundary
        TetrahedralMesh<2,2>::BoundaryElementIterator it =
            mesh.GetBoundaryElementIteratorBegin();
        while (it != mesh.GetBoundaryElementIteratorEnd())
        {
            for (unsigned i=0; i<(*it)->GetNumNodes(); i++)
            {
                TS_ASSERT((*it)->GetNode(i)->IsBoundaryNode());
            }
            it++;
        }
        TS_ASSERT( mesh.CheckIsConforming() );
    }

    void Test1DClosedMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/circle_outline");
        TetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 100u);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 100u);
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 0u);
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 0u);
        TS_ASSERT_EQUALS( mesh.CalculateMaximumContainingElementsPerProcess(), 2U); //It's 1D

    }

    void Test1DMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/semicircle_outline");
        TetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 51u);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 50u);
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 2u);
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 2u);
    }

    void Test1DBranchedMeshIn3DSpace()
    {
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/branched_1d_in_3d_mesh");
        TetrahedralMesh<1,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 31u);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 30u);
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 3u);
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 3u);
    }

    void Test2DClosedMeshIn3DSpace()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/slab_395_elements");
        TetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 132u);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 224u);
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 0u);
    }

    void Test2DMeshIn3DSpace()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
        TetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 312u);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 522u);

        // Check that the mesh_reader has the culled "faces" (which are edges) (100 instead of 833)
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 100u);
        // These are the 100 edges around the perimeter of the circle
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 100u);

        TS_ASSERT( mesh.CheckIsConforming() );
    }


    void Test1DMeshCrossReference()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        Node<1>* p_node = mesh.GetNode(0);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 1u);
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 1u);
        Node<1>::ContainingElementIterator elt_iter = p_node->ContainingElementsBegin();
        Node<1>::ContainingBoundaryElementIterator b_elt_iter = p_node->ContainingBoundaryElementsBegin();
        TS_ASSERT_EQUALS(*elt_iter, 0u);
        TS_ASSERT_EQUALS(*b_elt_iter, 0u);

        // There is only one boundary element at this end
        TS_ASSERT_EQUALS(++b_elt_iter, p_node->ContainingBoundaryElementsEnd());

        Element<1,1>* p_element = mesh.GetElement(*elt_iter);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 1u);

        c_matrix<double, 1, 1> jacobian;
        double det;
        p_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(det, 0.1, 1e-6);

        c_matrix<double, 1, 1> cached_jacobian;
        double cached_det;
        mesh.GetJacobianForElement(p_element->GetIndex(), cached_jacobian, cached_det);

        TS_ASSERT_EQUALS(cached_det, det);
        TS_ASSERT_EQUALS(jacobian(0,0), cached_jacobian(0,0));

        Node<1>* p_node2 = mesh.GetNode(1);
        TS_ASSERT_EQUALS(p_node2->GetNumContainingElements(), 2u);
        TS_ASSERT_EQUALS(p_node2->GetNumBoundaryElements(), 0u);

        elt_iter = p_node2->ContainingElementsBegin();
        p_element = mesh.GetElement(*elt_iter);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),0u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),1u);

        p_element = mesh.GetElement(*(++elt_iter));
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),1u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),2u);

        p_element->CalculateJacobian(jacobian, det);

        TS_ASSERT_DELTA(det, 0.1, 1e-6);

        mesh.GetJacobianForElement(p_element->GetIndex(), cached_jacobian, cached_det);

        TS_ASSERT_EQUALS(cached_det, det);
        TS_ASSERT_EQUALS(jacobian(0,0), cached_jacobian(0,0));

        // There should be no more containing elements
        TS_ASSERT_EQUALS(++elt_iter, p_node2->ContainingElementsEnd());
    }

    void Test2DMeshCrossReference()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        Node<2>* p_node = mesh.GetNode(234);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 5u);
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 0u);

        Node<2>::ContainingElementIterator elt_iter = p_node->ContainingElementsBegin();
        Element<2,2>* p_element;

        p_element = mesh.GetElement(*elt_iter);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 474u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 290u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2), 234u);

        p_element = mesh.GetElement(*(++elt_iter));
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 234u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 461u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2), 460u);

        p_element = mesh.GetElement(*(++elt_iter));
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 290u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 459u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2), 234u);

        p_element = mesh.GetElement(*(++elt_iter));
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 459u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 461u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2), 234u);

        p_element = mesh.GetElement(*(++elt_iter));
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 460u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 474u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2), 234u);

        // Now look at a boundary node
        p_node = mesh.GetNode(99);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 3u);
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 2u);
        const BoundaryElement<1,2>* p_boundary_element;
        Node<2>::ContainingBoundaryElementIterator b_elt_iter = p_node->ContainingBoundaryElementsBegin();

        p_boundary_element = mesh.GetBoundaryElement(*b_elt_iter);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0), 98u);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1), 99u);
        p_boundary_element = mesh.GetBoundaryElement(*(++b_elt_iter));
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0), 99u);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1), 0u);
    }

    void Test3DMeshCrossReference()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        Node<3>* p_node = mesh.GetNode(34);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 10u);

        Element<3,3>* p_element;
        Node<3>::ContainingElementIterator elt_iter = p_node->ContainingElementsBegin();

        p_element = mesh.GetElement(*elt_iter);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 22u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 34u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2), 33u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3), 10u);

        p_element = mesh.GetElement(*(++elt_iter));
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 22u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 35u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2), 33u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3), 34u);

        // Now look at a boundary node
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 4u);
        const BoundaryElement<2,3>* p_boundary_element;
        Node<3>::ContainingBoundaryElementIterator b_elt_iter = p_node->ContainingBoundaryElementsBegin();
        p_boundary_element = mesh.GetBoundaryElement(*b_elt_iter);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0), 6u);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1), 34u);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2), 24u);
        p_boundary_element = mesh.GetBoundaryElement(*(++b_elt_iter));
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0), 6u);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1), 30u);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2), 34u);
        p_boundary_element = mesh.GetBoundaryElement(*(++b_elt_iter));
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0), 24u);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1), 34u);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2), 10u);
        p_boundary_element = mesh.GetBoundaryElement(*(++b_elt_iter));
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0), 34u);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1), 30u);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2), 10u);
    }

    void TestCalculateDesignatedOwnership()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        unsigned local_faces = 0;
        for (unsigned index = 0; index < mesh.GetNumBoundaryElements(); index++)
        {
            if (mesh.CalculateDesignatedOwnershipOfBoundaryElement(index))
            {
                local_faces++;
            }
        }
        unsigned total_num_faces;
        MPI_Allreduce(&local_faces, &total_num_faces, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        TS_ASSERT_EQUALS(total_num_faces, mesh.GetNumBoundaryElements());

        unsigned local_elements = 0;
        for (unsigned index = 0; index < mesh.GetNumElements(); index++)
        {
            if (mesh.CalculateDesignatedOwnershipOfElement(index))
            {
                local_elements++;
            }
        }
        unsigned total_num_elements;
        MPI_Allreduce(&local_elements, &total_num_elements, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        TS_ASSERT_EQUALS(total_num_elements, mesh.GetNumElements());

    }

    void TestNodePermutation()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        double volume = mesh.GetVolume();
        double surface = mesh.GetSurfaceArea();

        Node<3>* p_node0 = mesh.GetNode(0);
        Node<3>* p_node121 = mesh.GetNode(121);
        Node<3>* p_node125 = mesh.GetNode(125);
        Node<3>* p_node273 = mesh.GetNode(273);
        Node<3>* p_node357 = mesh.GetNode(357);
        Node<3>* p_node35 = mesh.GetNode(35);
        Node<3>* p_node219 = mesh.GetNode(219);
        Node<3>* p_node319 = mesh.GetNode(319);

        mesh.PermuteNodes();

        TS_ASSERT_EQUALS(mesh.GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNode(121)->GetIndex(), 121u);
        TS_ASSERT_EQUALS(mesh.GetNode(125)->GetIndex(), 125u);
        TS_ASSERT_EQUALS(mesh.GetNode(273)->GetIndex(), 273u);

        TS_ASSERT_EQUALS(p_node0->GetIndex(),   344u);
        TS_ASSERT_EQUALS(p_node121->GetIndex(), 193u);
        TS_ASSERT_EQUALS(p_node125->GetIndex(), 203u);
        TS_ASSERT_EQUALS(p_node273->GetIndex(), 266u);
        TS_ASSERT_EQUALS(p_node357->GetIndex(),   6u);
        TS_ASSERT_EQUALS(p_node35->GetIndex(),  188u);
        TS_ASSERT_EQUALS(p_node219->GetIndex(), 199u);
        TS_ASSERT_EQUALS(p_node319->GetIndex(), 160u);
        TS_ASSERT_EQUALS(mesh.GetNode(p_node0->GetIndex()), p_node0);
        TS_ASSERT_EQUALS(mesh.GetNode(p_node121->GetIndex()), p_node121);
        TS_ASSERT_EQUALS(mesh.GetNode(p_node125->GetIndex()), p_node125);
        TS_ASSERT_EQUALS(mesh.GetNode(p_node273->GetIndex()), p_node273);

        TS_ASSERT_DELTA(volume, mesh.GetVolume(), 1e-7);
        TS_ASSERT_DELTA(surface, mesh.GetSurfaceArea(), 1e-7);

    }

    void TestConstructRectangleStagger()
    {
        TetrahedralMesh<2,2> mesh;
        unsigned width = 20;
        unsigned height = 16;

        mesh.ConstructRectangularMesh(width, height, true);

        TS_ASSERT_DELTA(mesh.GetVolume(), width*height, 1e-7);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 2.0*(width+height), 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), ((width+1)*(height+1)));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2*(width + height));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  2*(width+height) );
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2*width*height);

        // Test indexing and connectivity of corners
        Node<2>* p_node;

        // Top-left
        p_node = mesh.GetNode(height*(width+1));
        TS_ASSERT_EQUALS(p_node->rGetLocation()[0], 0.0);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[1], height);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 2u);

        // Top-right
        p_node = mesh.GetNode((height+1)*(width+1)-1);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[0], width);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[1], height);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 2u);

        // Bottom-left
        p_node = mesh.GetNode(0);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[0], 0.0);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 2u);

        // Bottom right
        p_node = mesh.GetNode(width);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[0], width);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 2u);

//      Index with down-loop
//        //Verify some element indices -- top left diagonal goes NW-SE (normal)
//        TS_ASSERT_DELTA(mesh.GetElement(0)->CalculateCentroid()[0],            2.0/3.0, 1e-5);
//        TS_ASSERT_DELTA(mesh.GetElement(0)->CalculateCentroid()[1], (height-1)+2.0/3.0, 1e-5);
//        TS_ASSERT_DELTA(mesh.GetElement(1)->CalculateCentroid()[0],            1.0/3.0, 1e-5);
//        TS_ASSERT_DELTA(mesh.GetElement(1)->CalculateCentroid()[1], (height-1)+1.0/3.0, 1e-5);
//
//        //Verify some element indices -- bottom left diagonal goes SW-NE (stagger)
//        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[0], 1.0/3.0, 1e-5);
//        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[1], 2.0/3.0, 1e-5);
//        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[0], 2.0/3.0, 1e-5);
//        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[1], 1.0/3.0, 1e-5);
        // Verify some element indices -- top left diagonal goes NW-SE (normal)
        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[0],            2.0/3.0, 1e-5);
        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[1], (height-1)+2.0/3.0, 1e-5);
        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[0],            1.0/3.0, 1e-5);
        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[1], (height-1)+1.0/3.0, 1e-5);

        // Verify some element indices -- bottom left diagonal goes SW-NE (stagger)
        TS_ASSERT_EQUALS(height%2, 0u);//If height is even the bottom left is staggered
        TS_ASSERT_DELTA(mesh.GetElement(0)->CalculateCentroid()[0], 1.0/3.0, 1e-5);
        TS_ASSERT_DELTA(mesh.GetElement(0)->CalculateCentroid()[1], 2.0/3.0, 1e-5);
        TS_ASSERT_DELTA(mesh.GetElement(1)->CalculateCentroid()[0], 2.0/3.0, 1e-5);
        TS_ASSERT_DELTA(mesh.GetElement(1)->CalculateCentroid()[1], 1.0/3.0, 1e-5);

        TrianglesMeshWriter<2,2> mesh_writer("","RectangleMeshStagger");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // Test the other version of this method
        TetrahedralMesh<2,2> mesh2;
        mesh2.ConstructRegularSlabMesh(4.0, width, height);
        TS_ASSERT_DELTA(mesh2.GetVolume(), width*height, 1e-7);
        TS_ASSERT_EQUALS(mesh2.GetNumElements(), 4*5*2u);
    }

    void TestConstructRectangleNoStagger()
    {
        TetrahedralMesh<2,2> mesh;
        unsigned width = 38;
        unsigned height = 16;
        mesh.ConstructRectangularMesh(width, height, false);
        TS_ASSERT_DELTA(mesh.GetVolume(), width*height, 1e-7);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 2.0*(width+height), 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), ((width+1)*(height+1)));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  2*(width+height) );
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2*width*height );

        // Test indexing and connectivity of corners
        Node<2>* p_node;

        // Top-left
        p_node = mesh.GetNode(height*(width+1));
        TS_ASSERT_EQUALS(p_node->rGetLocation()[0], 0.0);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[1], height);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 2u);

        // Top-right
        p_node = mesh.GetNode((height+1)*(width+1)-1);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[0], width);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[1], height);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 1u);

        // Bottom-left
        p_node = mesh.GetNode(0);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[0], 0.0);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 1u);

        // Bottom-right
        p_node = mesh.GetNode(width);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[0], width);
        TS_ASSERT_EQUALS(p_node->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 2u);

//      Index with down-loop
//        //Verify some element indices -- top left diagonal goes NW-SE (normal)
//        TS_ASSERT_DELTA(mesh.GetElement(0)->CalculateCentroid()[0],            2.0/3.0, 1e-5);
//        TS_ASSERT_DELTA(mesh.GetElement(0)->CalculateCentroid()[1], (height-1)+2.0/3.0, 1e-5);
//        TS_ASSERT_DELTA(mesh.GetElement(1)->CalculateCentroid()[0],            1.0/3.0, 1e-5);
//        TS_ASSERT_DELTA(mesh.GetElement(1)->CalculateCentroid()[1], (height-1)+1.0/3.0, 1e-5);
//
//        //Verify some element indices -- bottom left diagonal goes NW-SE (normal)
//        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[0], 2.0/3.0, 1e-5);
//        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[1], 2.0/3.0, 1e-5);
//        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[0], 1.0/3.0, 1e-5);
//        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[1], 1.0/3.0, 1e-5);
        //Verify some element indices -- top left diagonal goes NW-SE (normal)
        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[0],            2.0/3.0, 1e-5);
        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[1], (height-1)+2.0/3.0, 1e-5);
        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[0],            1.0/3.0, 1e-5);
        TS_ASSERT_DELTA(mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[1], (height-1)+1.0/3.0, 1e-5);

        // Verify some element indices -- bottom left diagonal goes NW-SE (normal)
        TS_ASSERT_DELTA(mesh.GetElement(0)->CalculateCentroid()[0], 2.0/3.0, 1e-5);
        TS_ASSERT_DELTA(mesh.GetElement(0)->CalculateCentroid()[1], 2.0/3.0, 1e-5);
        TS_ASSERT_DELTA(mesh.GetElement(1)->CalculateCentroid()[0], 1.0/3.0, 1e-5);
        TS_ASSERT_DELTA(mesh.GetElement(1)->CalculateCentroid()[1], 1.0/3.0, 1e-5);

        TrianglesMeshWriter<2,2> mesh_writer("","RectangleMeshNoStagger");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }

    void TestConstructLine()
    {
        TetrahedralMesh<1,1> mesh;
        unsigned width = 39;

        mesh.ConstructLinearMesh(width);

        TS_ASSERT_DELTA(mesh.GetVolume(), width, 1e-7);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 0u, 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), width+1);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  2u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), width);

        TrianglesMeshWriter<1,1> mesh_writer("","LineMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // test other version of the method
        TetrahedralMesh<1,1> mesh2;
        unsigned wrong_step = width/2; // 39/2=19 not 19.5
        TS_ASSERT_THROWS_THIS(mesh2.ConstructRegularSlabMesh(wrong_step, width),
                "Space step does not divide the size of the mesh");

        mesh2.ConstructRegularSlabMesh(width/2.0, width);
        TS_ASSERT_DELTA(mesh2.GetVolume(), width, 1e-7);
        TS_ASSERT_EQUALS(mesh2.GetNumElements(), 2u);
    }

    void TestConstructLineIn3D()
    {
        TetrahedralMesh<1,3> mesh;
        unsigned width = 39;

        mesh.ConstructLinearMesh(width);

        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 0u, 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), width+1);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  2u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), width);

        TrianglesMeshWriter<1,3> mesh_writer("","LineMeshIn3D");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }

    void TestSetOwnerships()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);


        //mesh.SetElementOwnerships();
        mesh.GetDistributedVectorFactory(); //First touch should run SetElementOwnerships

        bool unowned_element = false;
        for (unsigned ele_num=0; ele_num< mesh.GetNumElements(); ele_num++)
        {
            if (mesh.GetElement(ele_num)->GetOwnership() == false)
            {
                unowned_element = true;
                break;
            }
        }
        if (PetscTools::IsSequential())
        {
            //We own everything
            TS_ASSERT_EQUALS(unowned_element, false);
        }
        else
        {
            //We do not own everything
            TS_ASSERT_EQUALS(unowned_element, true);
        }
    }

    void TestOutwardNormal3D()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            BoundaryElement<2,3>* p_b_element = mesh.GetBoundaryElement(i);
            c_vector<double, 3> normal;
            double det;
            p_b_element->CalculateWeightedDirection(normal, det);
            c_vector<double, 3> centroid = p_b_element->CalculateCentroid();
            ChastePoint<3> out(centroid+normal);
            ChastePoint<3> in(centroid-normal);
            TS_ASSERT_THROWS_NOTHING(mesh.GetContainingElementIndex(in));
            TS_ASSERT_THROWS_CONTAINS(mesh.GetContainingElementIndex(out),"is not in mesh"); // full message is "Point (X,Y,Z) is not in mesh - all elements tested"
        }
    }

    void TestCheckOutwardNormals()
    {
        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
            TetrahedralMesh<3,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            mesh.CheckOutwardNormals();
        }
        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements_inward_normal");
            TetrahedralMesh<3,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            TS_ASSERT_THROWS_THIS(mesh.CheckOutwardNormals(), "Inward facing normal in boundary element index 7");
        }
        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            mesh.CheckOutwardNormals();
        }
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRectangularMesh(2,3);
            mesh.CheckOutwardNormals();
        }
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRectangularMesh(2,3,false);
            mesh.CheckOutwardNormals();
        }
        {
            TetrahedralMesh<3,3> mesh;
            mesh.ConstructCuboid(2,3,4);
            mesh.CheckOutwardNormals();
        }
        //These can't be done
        {
            TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
            TetrahedralMesh<2,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            TS_ASSERT_THROWS_THIS(mesh.CheckOutwardNormals(),
                    "Don't have enough information to calculate a normal vector");
        }
        {
            TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            TS_ASSERT_THROWS_THIS(mesh.CheckOutwardNormals(),
                    "1-D mesh has no boundary normals");
        }
    }

    void TestConstructCuboid()
    {
        TetrahedralMesh<3,3> mesh;
        unsigned width = 7;
        unsigned height = 4;
        unsigned depth = 5;

        unsigned num_boundary_nodes =   2*( (width+1)*(height+1) + (width+1)*(depth+1) + (depth+1)*(height+1) )
                                      - 4*(width-1 + height-1 + depth-1)
                                      - 16;

        mesh.ConstructCuboid(width,height,depth);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), ((width+1)*(height+1)*(depth+1)));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), num_boundary_nodes);

        TS_ASSERT_DELTA(mesh.GetVolume(), width*height*depth, 1e-7);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 2.0*(width*height+height*depth+depth*width), 1e-7);

        // Each unit square on the surface is split into 2
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  4*(width*height+height*depth+depth*width));

        // Assuming that each cube is split into 6 tetrahedra
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6*width*height*depth );

        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            BoundaryElement<2,3>* p_b_element = mesh.GetBoundaryElement(i);
            c_vector<double, 3> normal;
            double det;
            p_b_element->CalculateWeightedDirection(normal, det);
            c_vector<double, 3> centroid = p_b_element->CalculateCentroid();
            ChastePoint<3> out(centroid+normal);
            ChastePoint<3> in(centroid-normal);
            normal /= norm_2(normal);
            if (fabs(centroid[0]) < 1e-5)
            {
                TS_ASSERT_DELTA(normal[0], -1.0, 1e-16);
            }
            if (fabs(centroid[0] - width) < 1e-5)
            {
                TS_ASSERT_DELTA(normal[0], 1.0, 1e-16);
            }
            if (fabs(centroid[1]) < 1e-5)
            {
                TS_ASSERT_DELTA(normal[1], -1.0, 1e-16);
            }
            if (fabs(centroid[1] - height) < 1e-5)
            {
                TS_ASSERT_DELTA(normal[1], 1.0, 1e-16);
            }
            if (fabs(centroid[2]) < 1e-5)
            {
                TS_ASSERT_DELTA(normal[2], -1.0, 1e-16);
            }
            if (fabs(centroid[2] - depth) < 1e-5)
            {
                TS_ASSERT_DELTA(normal[2], 1.0, 1e-16);
            }
            TS_ASSERT_THROWS_NOTHING(mesh.GetContainingElementIndex(in))
            TS_ASSERT_THROWS_CONTAINS(mesh.GetContainingElementIndex(out),"is not in mesh"); // full message is "Point (X,Y,Z) is not in mesh - all elements tested"
        }
        TS_ASSERT( mesh.CheckIsConforming() );
        TS_ASSERT_EQUALS(mesh.CalculateMaximumContainingElementsPerProcess(), 24U); //Four surrounding cubes may have all 6 tetrahedra meeting at a node

        TrianglesMeshWriter<3,3> mesh_writer("", "CuboidMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // test the other version of this method
        TetrahedralMesh<3,3> mesh2;
        mesh2.ConstructRegularSlabMesh(1.0, width, height, depth);
        TS_ASSERT_DELTA(mesh2.GetVolume(), width*height*depth, 1e-7);
        TS_ASSERT_EQUALS(mesh2.GetNumElements(), width*height*depth*6u);
    }

    void TestPermute()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[0], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[2], 0.0);

        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[0], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[2], 0.0);

        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[0], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[1], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[2], 0.0);

        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(0)->GetIndex(), 11u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(1)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(2)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(3)->GetIndex(), 0u);

        // Make identity permuation
        std::vector<unsigned> perm;

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            perm.push_back(i);
        }
        // perm is now the identity permuation

        // Rotate first three
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 0;

        // Rotate nodes in the first element
        perm[8] = 11;
        perm[9] = 8;
        perm[11] = 9;

        mesh.PermuteNodes(perm);

        TS_ASSERT_EQUALS(mesh.GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNode(7)->GetIndex(), 7u);

        // 1 was node 0
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[0], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[2], 0.0);

        // 2 was node 1
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[0], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[2], 0.0);

        // 0 was node 2
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[0], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[1], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[2], 0.0);

        // Element 0 new indexes
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(0)->GetIndex(), 9u);  // 9 was 11
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(1)->GetIndex(), 3u);  // 3 is 3
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(2)->GetIndex(), 11u); // 11 was 8
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(3)->GetIndex(), 1u);  // 1 was 0

        // Coverage of GetNodeFromPrePermutationIndex()
        TS_ASSERT_EQUALS(mesh.GetNodeFromPrePermutationIndex(11)->GetIndex(), 9u);
    }

    void TestClear()
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2,3);

        TS_ASSERT_EQUALS(mesh.GetVolume(), 6.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 12u);

        mesh.Clear();

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), 0u);
    }

    void TestGetVectorBetweenPoints()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        c_vector<double, 3> location1 = mesh.GetNode(0)->rGetLocation();
        c_vector<double, 3> location2 = mesh.GetNode(2)->rGetLocation();

        // Test a normal distance calculation
        c_vector<double, 3> vector = mesh.GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 0.2, 1e-7);
        TS_ASSERT_DELTA(vector[1], 0.2, 1e-7);
        TS_ASSERT_DELTA(vector[2], 0.0, 1e-7)
        TS_ASSERT_DELTA(norm_2(vector), sqrt(0.08), 1e-7);
        TS_ASSERT_DELTA(mesh.GetDistanceBetweenNodes(0, 2), sqrt(0.08), 1e-7);

        // And the opposite vector
        vector = mesh.GetVectorFromAtoB(location2, location1);
        TS_ASSERT_DELTA(vector[0], -0.2, 1e-7);
        TS_ASSERT_DELTA(vector[1], -0.2, 1e-7);
        TS_ASSERT_DELTA(vector[2], 0.0, 1e-7);
        TS_ASSERT_DELTA(norm_2(vector), sqrt(0.08), 1e-7);

        // A 3d vector
        location1[0] = 0.5;
        location1[1] = 3.0;
        location1[2] = 1.0;
        location2[0] = 2.5;
        location2[1] = 4.0;
        location2[2] = -3.0;
        vector = mesh.GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], +2.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], +1.0, 1e-7);
        TS_ASSERT_DELTA(vector[2], -4.0, 1e-7);
        TS_ASSERT_DELTA(norm_2(vector), sqrt(21.0), 1e-7);
    }

    void TestMeshGetWidthAndBoundingBoxMethod()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double width = mesh.GetWidth(0u);
        double height = mesh.GetWidth(1u);

        TS_ASSERT_DELTA(width, 2, 1e-6);
        TS_ASSERT_DELTA(height, 2, 1e-6);

        ChasteCuboid<2> bounds=mesh.CalculateBoundingBox();
        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], 1, 1e-6);
        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], 1, 1e-6);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], -1, 1e-6);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], -1, 1e-6);
    }

    void TestPointWeightsInElement1D()
    {
        std::vector<Node<1>*> nodes1d;
        nodes1d.push_back(new Node<1>(0, false, 2.0));
        nodes1d.push_back(new Node<1>(1, false, 2.5));
        Element<1,1> element1d(INDEX_IS_NOT_USED, nodes1d);

        ChastePoint<1> in_point(2.25);
        ChastePoint<1> on_point(2.00);
        ChastePoint<1> out_point(1.25);

        c_vector<double, 2> weights;
        weights = element1d.CalculateInterpolationWeights(on_point);
        TS_ASSERT_EQUALS(weights[0], 1.0);
        TS_ASSERT_EQUALS(weights[1], 0.0);

        weights = element1d.CalculateInterpolationWeights(in_point);
        TS_ASSERT_EQUALS(weights[0], 0.5);
        TS_ASSERT_EQUALS(weights[1], 0.5);

        weights = element1d.CalculateInterpolationWeights(out_point);
        //1.25 = 2.5*2 -1.5 * 2.5
        TS_ASSERT_EQUALS(weights[0], 2.5);
        TS_ASSERT_EQUALS(weights[1], -1.5);

        delete nodes1d[0];
        delete nodes1d[1];
    }

    void TestPointInElement1D()
    {
        std::vector<Node<1>*> nodes1d;
        nodes1d.push_back(new Node<1>(0, false, 2.0));
        nodes1d.push_back(new Node<1>(1, false, 2.5));
        Element<1,1> element1d(INDEX_IS_NOT_USED, nodes1d);

        ChastePoint<1> in_point(2.25);
        ChastePoint<1> on_point(2.00);
        ChastePoint<1> out_point(1.25);
        bool strict = true;
        TS_ASSERT_EQUALS(element1d.IncludesPoint(in_point), true);
        TS_ASSERT_EQUALS(element1d.IncludesPoint(on_point), true);
        TS_ASSERT_EQUALS(element1d.IncludesPoint(on_point, strict), false);
        TS_ASSERT_EQUALS(element1d.IncludesPoint(out_point), false);

        delete nodes1d[0];
        delete nodes1d[1];
    }

    void TestPointinMesh1D()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ChastePoint<1> point1(0.15);
        ChastePoint<1> point2(-0.1);
        ChastePoint<1> point3(0.2);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1),1u);
        TS_ASSERT_THROWS_CONTAINS(mesh.GetContainingElementIndex(point2),"is not in mesh"); // full message is "Point (X,Y,Z) is not in mesh - all elements tested"
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point3),1u);  //in elements 1 and 2

        std::vector<unsigned> indices;
        indices=mesh.GetContainingElementIndices(point1);
        TS_ASSERT_EQUALS(indices.size(), 1u);
        TS_ASSERT_EQUALS(indices[0], 1u);

        indices=mesh.GetContainingElementIndices(point2);
        TS_ASSERT_EQUALS(indices.size(), 0u);

        indices=mesh.GetContainingElementIndices(point3);
        TS_ASSERT_EQUALS(indices.size(), 2u);
        TS_ASSERT_EQUALS(indices[0], 1u);
        TS_ASSERT_EQUALS(indices[1], 2u);

        // Check GetNearestNodeToPoint() methods
        TS_ASSERT_EQUALS(mesh.GetNearestNodeIndex(point1), 1u); // Answer will be the first of nodes 1 and 2 to be checked; 1, in this case
        TS_ASSERT_EQUALS(mesh.GetNearestNodeIndex(point2), 0u); // Point is before start of the fibre, so closest node is the first
        TS_ASSERT_EQUALS(mesh.GetNearestNodeIndex(point3), 2u); // Point is exactly at node 2

    }

    void TestPointWeightsAndInclusion2D()
    {
        std::vector<Node<2>*> nodes2d;
        nodes2d.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes2d.push_back(new Node<2>(1, false, 2.0, 1.0));
        nodes2d.push_back(new Node<2>(2, false, 0.0, 3.0));
        Element<2,2> element2d(INDEX_IS_NOT_USED, nodes2d);

        ChastePoint<2> on_point(0.0, 2.0);
        c_vector<double, 3> weights;
        bool strict = true;
        TS_ASSERT_EQUALS(element2d.IncludesPoint(on_point), true);
        TS_ASSERT_EQUALS(element2d.IncludesPoint(on_point, strict), false);

        weights = element2d.CalculateInterpolationWeights(on_point);
        c_vector<double, 2> xi_on = element2d.CalculateXi(on_point);
        TS_ASSERT_DELTA(weights[0], 1.0/3.0, 1e-5);
        TS_ASSERT_DELTA(weights[1], 0.0, 1e-5);
        TS_ASSERT_DELTA(weights[2], 2.0/3.0, 1e-5);
        TS_ASSERT_DELTA(xi_on[0],  0.0, 1e-5);
        TS_ASSERT_DELTA(xi_on[1],  2.0/3.0, 1e-5);

        ChastePoint<2> in_point(1.0, 1.0);
        TS_ASSERT_EQUALS(element2d.IncludesPoint(in_point), true);
        weights = element2d.CalculateInterpolationWeights(in_point);
        c_vector<double, 2> xi_in = element2d.CalculateXi(in_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(0.0, weights[2]);
        TS_ASSERT_DELTA(xi_in[0], 0.5,1e-12);
        TS_ASSERT_DELTA(xi_in[1], 1.0/6.0,1e-12);

        ChastePoint<2> out_point(1.0, 0.0);
        TS_ASSERT_EQUALS(element2d.IncludesPoint(out_point), false);
        weights = element2d.CalculateInterpolationWeights(out_point);
        c_vector<double, 2> xi_out = element2d.CalculateXi(out_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(weights[2], 0.0);
        TS_ASSERT_DELTA(xi_out[0],0.5,1e-12);
        TS_ASSERT_DELTA(xi_out[1],-1.0/6.0,1e-12);

        delete nodes2d[0];
        delete nodes2d[1];
        delete nodes2d[2];
    }

    void TestPointinMesh2D()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ChastePoint<2> point1(0.051, 0.051);
        ChastePoint<2> point2(0.2, 0.2);
        ChastePoint<2> point3(0.05, 0.05); // node 60 of mesh


        TS_ASSERT_EQUALS(mesh.GetNearestNodeIndex(point1), 60u);
        TS_ASSERT_EQUALS(mesh.GetNearestNodeIndex(point2), 120u);   // Closest node is top right node
        TS_ASSERT_EQUALS(mesh.GetNearestNodeIndex(point3), 60u);

        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1), 110u);
        TS_ASSERT_EQUALS(mesh.GetNearestElementIndex(point1), 110u);
        TS_ASSERT_THROWS_CONTAINS(mesh.GetContainingElementIndex(point2),"is not in mesh"); // full message is "Point (X,Y,Z) is not in mesh - all elements tested"
        TS_ASSERT_EQUALS(mesh.GetNearestElementIndex(point2), 199u); // contains top-right corner
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point3), 89u); // in elements 89,90,91,108,109, 110

        std::vector<unsigned> indices;
        indices = mesh.GetContainingElementIndices(point1);
        TS_ASSERT_EQUALS(indices.size(), 1u);
        TS_ASSERT_EQUALS(indices[0], 110u);

        indices = mesh.GetContainingElementIndices(point2);
        TS_ASSERT_EQUALS(indices.size(), 0u);

        indices = mesh.GetContainingElementIndices(point3);
        TS_ASSERT_EQUALS(indices.size(), 6u);
        TS_ASSERT_EQUALS(indices[0], 89u);
        TS_ASSERT_EQUALS(indices[1], 90u);
        TS_ASSERT_EQUALS(indices[5], 110u);


        std::set<unsigned> test_elements;
        test_elements.insert(0);
        TS_ASSERT_EQUALS(mesh.GetNearestElementIndexFromTestElements(point1, test_elements), 0u);
        test_elements.insert(107);
        test_elements.insert(1);
        test_elements.insert(2);
        TS_ASSERT_EQUALS(mesh.GetNearestElementIndexFromTestElements(point1, test_elements), 107u);
        test_elements.insert(110);
        TS_ASSERT_EQUALS(mesh.GetNearestElementIndexFromTestElements(point1, test_elements), 110u);
    }

    void TestPointInElement3D()
    {
        std::vector<Node<3>*> nodes3d;
        nodes3d.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes3d.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element3d(INDEX_IS_NOT_USED, nodes3d);

        bool strict = true;
        ChastePoint<3> on_point(0., 0.2, 0.);
        c_vector<double, 4> weights;
        TS_ASSERT_EQUALS(element3d.IncludesPoint(on_point), true);
        TS_ASSERT_EQUALS(element3d.IncludesPoint(on_point, strict), false);
        weights = element3d.CalculateInterpolationWeights(on_point);
        c_vector<double, 3> xi_on = element3d.CalculateXi(on_point);

        TS_ASSERT_DELTA(weights[0], 0.8, 1e-5);
        TS_ASSERT_DELTA(weights[1], 0.0, 1e-5);
        TS_ASSERT_DELTA(weights[2], 0.2, 1e-5);
        TS_ASSERT_DELTA(weights[3], 0.0, 1e-5);
        TS_ASSERT_DELTA(xi_on[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(xi_on[1], 0.2, 1e-12);
        TS_ASSERT_DELTA(xi_on[2], 0.0, 1e-12);

        ChastePoint<3> in_point(0.25, 0.25, 0.25);
        TS_ASSERT_EQUALS(element3d.IncludesPoint(in_point), true);

        weights = element3d.CalculateInterpolationWeights(in_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(0.0, weights[2]);
        TS_ASSERT_LESS_THAN(0.0, weights[3]);
        // Weights are non-negative and sum to 1
        TS_ASSERT_DELTA(norm_1(weights), 1.0, 1e-12);

        // If point is in element then projection is redundant - looks the same as above
        weights = element3d.CalculateInterpolationWeightsWithProjection(in_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(0.0, weights[2]);
        TS_ASSERT_LESS_THAN(0.0, weights[3]);
        // Weights are non-negative and sum to 1
        TS_ASSERT_DELTA(norm_1(weights), 1.0, 1e-12);

        ChastePoint<3> out_point(0.1, -10., 0.1);
        TS_ASSERT_EQUALS(element3d.IncludesPoint(out_point), false);

        weights = element3d.CalculateInterpolationWeights(out_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(weights[2], 0.0);
        TS_ASSERT_LESS_THAN(0.0, weights[3]);
        TS_ASSERT_DELTA(weights[0], 10.8, 1e-5);
        TS_ASSERT_DELTA(weights[1], 0.1, 1e-5);
        TS_ASSERT_DELTA(weights[2], -10.0, 1e-5);
        TS_ASSERT_DELTA(weights[3], 0.1, 1e-5);
        // Weights still sum to 1, but one weight is negative
        TS_ASSERT_DELTA(norm_1(weights), 21.0, 1e-12);

        weights = element3d.CalculateInterpolationWeightsWithProjection(out_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN_EQUALS(0.0, weights[2]);
        TS_ASSERT_EQUALS(0.0, weights[2]);
        TS_ASSERT_LESS_THAN(0.0, weights[3]);
        // Weights are non-negative and sum to 1
        TS_ASSERT_DELTA(norm_1(weights), 1.0, 1e-12);

        delete nodes3d[0];
        delete nodes3d[1];
        delete nodes3d[2];
        delete nodes3d[3];
    }

    void TestPointinMesh3D()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_1mm_6000_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ChastePoint<3> point1(0.051, 0.051,0.051);
        ChastePoint<3> point2(0.2, 0.2, 0.2);
        ChastePoint<3> point3(0.050000000000000003, 0.050000000000000003, 0.050000000000000003);
        // Node 665 of mesh
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1), 2992u);
        TS_ASSERT_THROWS_CONTAINS(mesh.GetContainingElementIndex(point2),"is not in mesh"); // full message is "Point (X,Y,Z) is not in mesh - all elements tested"
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point3), 2044u);
        /*in elements 2044, 2047. 2058, 2192, 2268, 2286, 2392, 2414, 2415,
         * 2424, 2426, 2452, 2661, 2704, 2734, 2745, 2846, 2968, 2990, 2992,
         * 3015, 3022, 3024, 3026
         */

        TS_ASSERT_EQUALS(mesh.GetContainingElementIndexWithInitialGuess(point1, 0), 2992u);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndexWithInitialGuess(point1, 2991), 2992u);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndexWithInitialGuess(point1, 2992), 2992u);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndexWithInitialGuess(point1, 2993), 2992u);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndexWithInitialGuess(point1, 5999), 2992u);

        //Note from commemt above that point3 is on the boundary of multiple elements
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndexWithInitialGuess(point3, 2000), 2044u);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndexWithInitialGuess(point3, 2045), 2047u);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndexWithInitialGuess(point3, 3025), 3026u);

        // This should throw because vertex is not contained in any element
        TS_ASSERT_THROWS_CONTAINS(mesh.GetContainingElementIndexWithInitialGuess(point2, 0), "not in mesh - all elements tested");
        // This should throw because vertex is not strictly contained in any element
        TS_ASSERT_THROWS_CONTAINS(mesh.GetContainingElementIndexWithInitialGuess(point3, 0, true), "not in mesh - all elements tested");

        std::vector<unsigned> indices;
        indices = mesh.GetContainingElementIndices(point1);
        TS_ASSERT_EQUALS(indices.size(), 1u);
        TS_ASSERT_EQUALS(indices[0], 2992u);

        indices = mesh.GetContainingElementIndices(point2);
        TS_ASSERT_EQUALS(indices.size(), 0u);

        indices = mesh.GetContainingElementIndices(point3);
        TS_ASSERT_EQUALS(indices.size(), 24u);
        TS_ASSERT_EQUALS(indices[0], 2044u);
        TS_ASSERT_EQUALS(indices[1], 2047u);
        TS_ASSERT_EQUALS(indices[5], 2286u);
        TS_ASSERT_EQUALS(indices[23], 3026u);

        // Test when a suggested set of elements is given
        std::set<unsigned> suggested_elements;
        suggested_elements.insert(2991);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1,false,suggested_elements), 2992u);
        suggested_elements.insert(2992); // should find it quicker, but can't really test that
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1,false,suggested_elements), 2992u);
        suggested_elements.insert(2993);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1,false,suggested_elements), 2992u);

        // Check GetNearestNodeToPoint() methods
        TS_ASSERT_EQUALS(mesh.GetNearestNodeIndex(point1), 665u); // Close to node 665
        TS_ASSERT_EQUALS(mesh.GetNearestNodeIndex(point2), 1330u); // Point is outside mesh, closest node is the one on the outermost corner
        TS_ASSERT_EQUALS(mesh.GetNearestNodeIndex(point3), 665u); // Point is exactly at node 665
    }

    void TestFloatingPointIn3D()
    {
        // There's some weird failing behaviour in the refined mesh test.
        // This test duplicates it.

        TetrahedralMesh<3,3> mesh;

        double third = 1.0L/3.0L;
        mesh.ConstructRegularSlabMesh(third, 1.0, 1.0, 1.0);

        ChastePoint<3> point_on_edge1(5.0/6.0,   0.5,       1.0);
        ChastePoint<3> point_on_edge2(5.0L/6.0L, 0.5,       1.0);
        ChastePoint<3> point_on_edge3(5.0L/6.0L, 0.5L,      1.0L);
        ChastePoint<3> point_on_edge4(5.0L/6.0L, 3.0L/6.0L, 1.0L);
        ChastePoint<3> point_on_edge5(5.0L/6.0L, 0.5L,      6.0L/6.0L);
        ChastePoint<3> point_on_edge6(5.0L/6.0L, 3.0L/6.0L, 6.0L/6.0L);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point_on_edge1), 142u);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point_on_edge2), 142u);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point_on_edge3), 142u);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point_on_edge4), 142u);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point_on_edge5), 142u);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point_on_edge6), 142u);
    }

    void TestGetAngleBetweenNodes()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(0,1),  0.0,      1e-12);
        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(0,2),  M_PI/4,   1e-12);
        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(0,3),  M_PI/2,   1e-12);
        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(1,0),  M_PI,     1e-12);
        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(2,0), -3*M_PI/4, 1e-12);

        TS_ASSERT_THROWS_THIS(mesh.GetAngleBetweenNodes(0,0),"Tried to compute polar angle of (0,0)");
    }

    void TestNodesPerProcessorFile()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Throws because file does not exist
        TS_ASSERT_THROWS_THIS(mesh.ReadNodesPerProcessorFile("dsgund"),"Unable to read nodes per processor file dsgund");

        // Throws because sum of nodes is not equal to the number of nodes in the mesh
        TS_ASSERT_THROWS_THIS(mesh.ReadNodesPerProcessorFile("mesh/test/data/nodes_per_processor_1.txt"),
                "Sum of nodes per processor, 5, not equal to number of nodes in mesh, 4");

        if (PetscTools::GetNumProcs() == 2)
        {
            mesh.ReadNodesPerProcessorFile("mesh/test/data/nodes_per_processor_2.txt");

            if (PetscTools::AmMaster())
            {
                TS_ASSERT_EQUALS(mesh.GetDistributedVectorFactory()->GetLocalOwnership(), 1u);
            }
            else if (PetscTools::AmTopMost())
            {
                TS_ASSERT_EQUALS(mesh.GetDistributedVectorFactory()->GetLocalOwnership(), 3u);
            }
        }
        if (PetscTools::IsSequential())
        {
            // Throws because the file is written for two processes
            TS_ASSERT_THROWS_THIS(mesh.ReadNodesPerProcessorFile("mesh/test/data/nodes_per_processor_2.txt"),
                  "Number of processes doesn't match the size of the nodes-per-processor file")
        }
    }

    void TestReadingMeshesWithRegionsAndGenericReader()
    {
        std::shared_ptr<AbstractMeshReader<1,1> > p_mesh_reader = GenericMeshReader<1,1>("mesh/test/data/1D_0_to_1_10_elements_with_attributes");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(*p_mesh_reader);

        TS_ASSERT_EQUALS(p_mesh_reader->GetNumElementAttributes(), 1u);

        // The first few are unsigned
        for (unsigned i=0; i<5; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetUnsignedAttribute(), i%5+1);
        }
    }

    void TestReadingMeshesWithRegionsElementsAndFaces3D()
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart_nonnegative_flags");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetUnsignedAttribute(), (i+1)%3+1);
        }

       TS_ASSERT_EQUALS(mesh_reader.GetNumFaceAttributes(), 1u);

        bool read_zero_attribute = false;
        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            if (mesh.GetBoundaryElement(i)->GetUnsignedAttribute()==0u)
            {
                read_zero_attribute = true;
            }
            TS_ASSERT_LESS_THAN(mesh.GetBoundaryElement(i)->GetAttribute(), 4.0);
        }
        TS_ASSERT(read_zero_attribute);
    }

    void TestReadingMeshesWithRegionsElementsAndFaces2D()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 0u);

        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetElement(i)->GetAttribute(), 0u, 1e-12);
        }

        TS_ASSERT_EQUALS(mesh_reader.GetNumFaceAttributes(), 1u);

        //In the edge file for this test
        // * All internal edges are marked with 0
        // * All external edges were marked as 1 by triangle
        // * The final edge marker has been edited from 1 to 2
        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            if (i==99)
            {
                TS_ASSERT_EQUALS(mesh.GetBoundaryElement(i)->GetUnsignedAttribute(), 2u);
            }
            else
            {
                TS_ASSERT_EQUALS(mesh.GetBoundaryElement(i)->GetUnsignedAttribute(), 1u);
            }
        }
    }

    void Test1DIn3DBranchWithAttributes()
    {
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/y_branch_3d_mesh");
        TetrahedralMesh<1,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 3u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetAttribute(), 50.0);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetAttribute(), 25.0);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetAttribute(), 25.0);
    }

    void TestCuboidMeshConstructors()
    {
        CuboidMeshConstructor<1> constructor1;
        TetrahedralMesh<1,1> mesh1;
        constructor1.Construct(mesh1, 1, 1.0);
        TS_ASSERT_EQUALS(constructor1.GetWidth(), 1.0);
        TS_ASSERT_EQUALS(mesh1.GetNumNodes(), 9u);
        TS_ASSERT(mesh1.CheckIsConforming());

        CuboidMeshConstructor<2> constructor2;
        TetrahedralMesh<2,2> mesh2;
        constructor2.Construct(mesh2, 1, 1.0);
        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), 9u*9u);

        CuboidMeshConstructor<3> constructor3;
        TetrahedralMesh<3,3> mesh3;
        constructor3.Construct(mesh3, 1, 1.0);
        TS_ASSERT_EQUALS(mesh3.GetNumNodes(), 9u*9u*9u);

        CuboidMeshConstructor<1,3> constructor4;
        TetrahedralMesh<1,3> mesh4;
        constructor4.Construct(mesh4, 1, 1.0);
        TS_ASSERT_EQUALS(mesh4.GetNumNodes(), 9u);

        TS_ASSERT_EQUALS(constructor4.GetNumNodes(), 9u);
    }

    void TestMeshStoresFilename()
    {
        TetrahedralMesh<3,3> mesh;
        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
            mesh.ConstructFromMeshReader(mesh_reader);
        }

        TS_ASSERT_EQUALS(mesh.IsMeshOnDisk(), true);

        std::string mesh_file_base_name = mesh.GetMeshFileBaseName();
        TrianglesMeshReader<3,3> mesh_reader(mesh_file_base_name);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), mesh.GetNumElements());

        TetrahedralMesh<3,3> cuboid_mesh;
        cuboid_mesh.ConstructCuboid(7, 4, 5);

        TS_ASSERT_THROWS_THIS(cuboid_mesh.GetMeshFileBaseName(),"This mesh was not constructed from a file.");
        TS_ASSERT_EQUALS(cuboid_mesh.IsMeshOnDisk(), false);
    }

    void TestCalculateBoundingBox()
    {
        TetrahedralMesh<1,1> mesh1d;
        mesh1d.ConstructLinearMesh(3);
        ChasteCuboid<1> extremes1d = mesh1d.CalculateBoundingBox();
        TS_ASSERT_DELTA(extremes1d.rGetLowerCorner()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(extremes1d.rGetUpperCorner()[0], 3.0, 1e-12);

        mesh1d.GetNode(0)->rGetModifiableLocation()[0] = -0.1;
        mesh1d.Scale(0.5);
        extremes1d = mesh1d.CalculateBoundingBox();
        TS_ASSERT_DELTA(extremes1d.rGetLowerCorner()[0], -0.05, 1e-12);
        TS_ASSERT_DELTA(extremes1d.rGetUpperCorner()[0], 1.5, 1e-12);

        TetrahedralMesh<2,2> mesh2d;
        mesh2d.ConstructRectangularMesh(2,2);
        ChasteCuboid<2> extremes2d = mesh2d.CalculateBoundingBox();
        TS_ASSERT_DELTA(extremes2d.rGetLowerCorner()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(extremes2d.rGetUpperCorner()[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(extremes2d.rGetLowerCorner()[1], 0.0, 1e-12);
        TS_ASSERT_DELTA(extremes2d.rGetUpperCorner()[1], 2.0, 1e-12);

        mesh2d.Rotate(-M_PI/2);
        mesh2d.Scale(2.0,3.0);
        extremes2d = mesh2d.CalculateBoundingBox();
        TS_ASSERT_DELTA(extremes2d.rGetLowerCorner()[0],-4.0, 1e-12);
        TS_ASSERT_DELTA(extremes2d.rGetUpperCorner()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(extremes2d.rGetLowerCorner()[1], 0.0, 1e-12);
        TS_ASSERT_DELTA(extremes2d.rGetUpperCorner()[1], 6.0, 1e-12);

        TetrahedralMesh<3,3> mesh3d;
        mesh3d.ConstructCuboid(4,5,6);
        ChasteCuboid<3> extremes3d = mesh3d.CalculateBoundingBox();
        TS_ASSERT_DELTA(extremes3d.rGetLowerCorner()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(extremes3d.rGetUpperCorner()[0], 4.0, 1e-12);
        TS_ASSERT_DELTA(extremes3d.rGetLowerCorner()[1], 0.0, 1e-12);
        TS_ASSERT_DELTA(extremes3d.rGetUpperCorner()[1], 5.0, 1e-12);
        TS_ASSERT_DELTA(extremes3d.rGetLowerCorner()[2], 0.0, 1e-12);
        TS_ASSERT_DELTA(extremes3d.rGetUpperCorner()[2], 6.0, 1e-12);
    }

    void TestArchiving()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "tetrahedral_mesh.arch";
        ArchiveLocationInfo::SetMeshFilename("tetrahedral_mesh");

        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

            AbstractTetrahedralMesh<2,2>* const p_mesh = new TetrahedralMesh<2,2>;
            p_mesh->ConstructFromMeshReader(mesh_reader);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

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
            TS_ASSERT_EQUALS(p_mesh2->GetNumNodes(), 543u);
            TS_ASSERT_EQUALS(p_mesh2->GetNumElements(), 984u);

            // Check some node co-ordinates
            TS_ASSERT_DELTA(p_mesh2->GetNode(0)->GetPoint()[0],  0.9980267283, 1e-6);
            TS_ASSERT_DELTA(p_mesh2->GetNode(0)->GetPoint()[1], -0.0627905195, 1e-6);
            TS_ASSERT_DELTA(p_mesh2->GetNode(1)->GetPoint()[0], 1.0, 1e-6);
            TS_ASSERT_DELTA(p_mesh2->GetNode(1)->GetPoint()[1], 0.0, 1e-6);

            // Check first element has the right nodes
            TetrahedralMesh<2,2>::ElementIterator iter = p_mesh2->GetElementIteratorBegin();
            TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(0), 309u);
            TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(1), 144u);
            TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(2), 310u);
            TS_ASSERT_EQUALS(iter->GetNode(1), p_mesh2->GetNode(144));

            delete p_mesh2;
        }
    }

    void TestArchiving2din3d()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "surface_mesh.arch";
        ArchiveLocationInfo::SetMeshFilename("surface_mesh");

        {
            TrianglesMeshReader<2,3> mesh_reader("cell_based/test/data/Square2dMeshIn3d/Square2dMeshIn3d");

            AbstractTetrahedralMesh<2,3>* const p_mesh = new TetrahedralMesh<2,3>;
            p_mesh->ConstructFromMeshReader(mesh_reader);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) << p_mesh;
            delete p_mesh;
        }

        {
            // Should archive the most abstract class you can to check boost knows what individual classes are.
            // (but here AbstractMesh doesn't have the methods below).
            AbstractTetrahedralMesh<2,3>* p_mesh2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // restore from the archive
            (*p_arch) >> p_mesh2;

            // Check we have the right number of nodes & elements
            TS_ASSERT_EQUALS(p_mesh2->GetNumNodes(), 4u);
            TS_ASSERT_EQUALS(p_mesh2->GetNumElements(), 2u);

            // Check some node co-ordinates
            TS_ASSERT_DELTA(p_mesh2->GetNode(0)->GetPoint()[0],  0.0, 1e-6);
            TS_ASSERT_DELTA(p_mesh2->GetNode(0)->GetPoint()[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(p_mesh2->GetNode(1)->GetPoint()[0], 1.0, 1e-6);
            TS_ASSERT_DELTA(p_mesh2->GetNode(1)->GetPoint()[1], 0.0, 1e-6);

            // Check first element has the right nodes
            TetrahedralMesh<2,3>::ElementIterator iter = p_mesh2->GetElementIteratorBegin();
            TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(0), 0u);
            TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(1), 1u);
            TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(2), 2u);

            delete p_mesh2;
        }
    }

    void TestDeepCopy()
    {
        TetrahedralMesh<3,3> copy_mesh;

        c_vector <double, 3> node3_location;
        node3_location[0]=0.0;
        node3_location[1]=0.2;
        node3_location[2]=0.0;
        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");

            TetrahedralMesh<3,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), 12u);
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 12u);
            TS_ASSERT_EQUALS(mesh.GetNumElements(), 12u);
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 20u);
            TS_ASSERT_DELTA(mesh.GetVolume(), 0.008, 1e-8);
            TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 0.24, 1e-8);
            TS_ASSERT_DELTA(norm_2(mesh.GetNode(3)->rGetLocation()-node3_location), 0.0, 1e-15)
            TS_ASSERT_EQUALS(mesh.GetNode(11)->GetNumContainingElements(), 6u);
            TS_ASSERT_EQUALS(mesh.GetNode(11)->GetNumBoundaryElements(), 6u);

            copy_mesh.ConstructFromMesh(mesh);
            //Original mesh goes out of scope here, so we'd like the copy not to be a shallow one.
        }
        TS_ASSERT_EQUALS(copy_mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(copy_mesh.GetNumBoundaryNodes(), 12u);
        TS_ASSERT_EQUALS(copy_mesh.GetNumElements(), 12u);
        TS_ASSERT_EQUALS(copy_mesh.GetNumBoundaryElements(), 20u);
        TS_ASSERT_DELTA(copy_mesh.GetVolume(), 0.008, 1e-8);
        TS_ASSERT_DELTA(copy_mesh.GetSurfaceArea(), 0.24, 1e-8);
        TS_ASSERT_DELTA(norm_2(copy_mesh.GetNode(3)->rGetLocation()-node3_location), 0.0, 1e-15)
        TS_ASSERT_EQUALS(copy_mesh.GetNode(11)->GetNumContainingElements(), 6u);
        TS_ASSERT_EQUALS(copy_mesh.GetNode(11)->GetNumBoundaryElements(), 6u);
    }

    void TestNodeExchange()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_800_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        std::vector<std::vector<unsigned> > nodes_to_send_per_process;
        std::vector<std::vector<unsigned> > nodes_to_receive_per_process;
        mesh.CalculateNodeExchange(nodes_to_send_per_process, nodes_to_receive_per_process);

        TS_ASSERT_EQUALS(nodes_to_send_per_process.size(), PetscTools::GetNumProcs());
        TS_ASSERT_EQUALS(nodes_to_receive_per_process.size(), PetscTools::GetNumProcs());
        TS_ASSERT(nodes_to_receive_per_process[PetscTools::GetMyRank()].empty());
        TS_ASSERT(nodes_to_send_per_process[PetscTools::GetMyRank()].empty());

        // Do some communication

        //mesh.rGetDistributedVectorFactory()->rGetGlobalLows();
        for ( unsigned rank_offset = 1; rank_offset < PetscTools::GetNumProcs(); rank_offset++ )
        {
            unsigned send_to      = (PetscTools::GetMyRank() + rank_offset) % (PetscTools::GetNumProcs());
            unsigned receive_from = (PetscTools::GetMyRank() + PetscTools::GetNumProcs()- rank_offset ) % (PetscTools::GetNumProcs());

            MPI_Send( &(nodes_to_send_per_process[send_to][0]),
                      nodes_to_send_per_process[send_to].size(),
                      MPI_UNSIGNED,
                      send_to,
                      0,
                      PETSC_COMM_WORLD );

            boost::scoped_array<unsigned> received(new unsigned[nodes_to_receive_per_process[receive_from].size()]);
            MPI_Status status;

            MPI_Recv( received.get(),
                      nodes_to_receive_per_process[receive_from].size(),
                      MPI_UNSIGNED,
                      receive_from,
                      0,
                      PETSC_COMM_WORLD,
                      &status );

            for ( unsigned i = 0; i < nodes_to_receive_per_process[receive_from].size(); i++ )
            {
                TS_ASSERT_EQUALS( received[i], nodes_to_receive_per_process[receive_from][i] );
            }
        }

//      s  for (unsigned process = 0; process < PetscTools::GetNumProcs(); process++)
//        {
//            PRINT_3_VARIABLES( process,
//                               nodes_to_receive_per_process[process].size(),
//                               nodes_to_send_per_process[process].size() );
//        }
    }


    void TestCalculateEdgeLengths()
    {
        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
            TetrahedralMesh<3,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            c_vector <double, 2> edge_len = mesh.CalculateMinMaxEdgeLengths();
            TS_ASSERT_DELTA(edge_len[0], 0.1, 1e-4);
            TS_ASSERT_DELTA(edge_len[1], 0.2828, 1e-4);
        }
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRectangularMesh(2,3);
            c_vector <double, 2> edge_len = mesh.CalculateMinMaxEdgeLengths();
            TS_ASSERT_DELTA(edge_len[0], 1.0, 1e-5);
            TS_ASSERT_DELTA(edge_len[1], sqrt(2.0), 1e-5);
        }
        {
            TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
            TetrahedralMesh<2,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            c_vector <double, 2> edge_len = mesh.CalculateMinMaxEdgeLengths();
            TS_ASSERT_DELTA(edge_len[0], 0.0628, 1e-4);
            TS_ASSERT_DELTA(edge_len[1], 0.2010, 1e-4);
        }
    }

    void TestConstructSlabMeshWithDimensionSplit()
    {
        double step = 1.0;
        unsigned width = 3;
        unsigned height = 5;
        unsigned depth = 7;

        // In 1D we shouldn't be able to change the split dimension from 0.  (Can only split in x.)
        {
            TetrahedralMesh<1,1> mesh;
            TS_ASSERT_THROWS_THIS(mesh.ConstructRegularSlabMeshWithDimensionSplit(1, step, width), "Cannot split on non-existent dimension");
            mesh.ConstructRegularSlabMeshWithDimensionSplit(0, step, width);
        }

        // In 2D we test that we can split in x (not just y)
        // We also test that a requested y-split gives back the default behaviour (node-for-node)
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(step, width, height);
            TetrahedralMesh<2,2> mesh_with_default_split;
            mesh_with_default_split.ConstructRegularSlabMeshWithDimensionSplit(1, step, width, height);
            TetrahedralMesh<2,2> mesh_with_x_split;
            mesh_with_x_split.ConstructRegularSlabMeshWithDimensionSplit(0, step, width, height);

            // Check that mesh and mesh_with_default_split are identical
            for (AbstractTetrahedralMesh<2,2>::NodeIterator iter = mesh.GetNodeIteratorBegin();
                    iter != mesh.GetNodeIteratorEnd();
                    ++iter)
            {
                unsigned index = iter->GetIndex();
                // Position of this node
                c_vector<double, 2> pos1;
                pos1 = iter->rGetLocation();
                // Position in the other mesh
                c_vector<double, 2> pos2;
                pos2 = mesh_with_default_split.GetNode(index)->rGetLocation();
                TS_ASSERT_DELTA(pos1[0], pos2[0], 1e-5);
                TS_ASSERT_DELTA(pos1[1], pos2[1], 1e-5);
            }

            // Check that the x-split has the same bounding box
            ChasteCuboid<2> bounds = mesh.CalculateBoundingBox();
            ChasteCuboid<2> bounds_with_x_split = mesh_with_x_split.CalculateBoundingBox();
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], bounds_with_x_split.rGetUpperCorner()[0], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], bounds_with_x_split.rGetUpperCorner()[1], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], bounds_with_x_split.rGetLowerCorner()[0], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], bounds_with_x_split.rGetLowerCorner()[1], DBL_EPSILON);

            // Same amount of stuff
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_with_x_split.GetNumNodes());
            TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_with_x_split.GetNumElements());
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_with_x_split.GetNumBoundaryElements());

            // Show that the x-split has a different indexing scheme
            // Normal meshes start at the origin
            c_vector<double, 2> orig1 = mesh_with_default_split.GetNode(0u)->rGetLocation();
            TS_ASSERT_DELTA(orig1[0], 0.0, 1e-5);
            TS_ASSERT_DELTA(orig1[1], 0.0, 1e-5);

            // The new one has the origin at an index height away
            c_vector<double, 2> orig2 = mesh_with_x_split.GetNode(height)->rGetLocation();
            TS_ASSERT_DELTA(orig2[0], 0.0, 1e-5);
            TS_ASSERT_DELTA(orig2[1], 0.0, 1e-5);
        }
        // In 3D we test that we can split in x, y (not just z)
        // We also test that a requested z-split gives back the default behaviour (node-for-node)
        {
            TetrahedralMesh<3,3> mesh;
            mesh.ConstructRegularSlabMesh(step, width, height, depth);
            TetrahedralMesh<3,3> mesh_with_default_split;
            mesh_with_default_split.ConstructRegularSlabMeshWithDimensionSplit(2, step, width, height, depth);
            TetrahedralMesh<3,3> mesh_with_x_split;
            mesh_with_x_split.ConstructRegularSlabMeshWithDimensionSplit(0, step, width, height, depth);
            TetrahedralMesh<3,3> mesh_with_y_split;
            mesh_with_y_split.ConstructRegularSlabMeshWithDimensionSplit(1, step, width, height, depth);

            // Check that mesh and mesh_with_default_split are identical
            for (AbstractTetrahedralMesh<3,3>::NodeIterator iter = mesh.GetNodeIteratorBegin();
                    iter != mesh.GetNodeIteratorEnd();
                    ++iter)
            {
                unsigned index = iter->GetIndex();
                // Position of this node
                c_vector<double, 3> pos1;
                pos1 = iter->rGetLocation();
                // Position in the other mesh
                c_vector<double, 3> pos2;
                pos2 = mesh_with_default_split.GetNode(index)->rGetLocation();
                TS_ASSERT_DELTA(pos1[0], pos2[0], 1e-5);
                TS_ASSERT_DELTA(pos1[1], pos2[1], 1e-5);
                TS_ASSERT_DELTA(pos1[2], pos2[2], 1e-5);
            }


            // Check that the x-split and y-split have the same bounding box
            ChasteCuboid<3> bounds = mesh.CalculateBoundingBox();
            ChasteCuboid<3> bounds_with_x_split = mesh_with_x_split.CalculateBoundingBox();
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], bounds_with_x_split.rGetUpperCorner()[0], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], bounds_with_x_split.rGetUpperCorner()[1], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[2], bounds_with_x_split.rGetUpperCorner()[2], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], bounds_with_x_split.rGetLowerCorner()[0], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], bounds_with_x_split.rGetLowerCorner()[1], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[2], bounds_with_x_split.rGetLowerCorner()[2], DBL_EPSILON);

            ChasteCuboid<3> bounds_with_y_split = mesh_with_y_split.CalculateBoundingBox();
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], bounds_with_y_split.rGetUpperCorner()[0], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], bounds_with_y_split.rGetUpperCorner()[1], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[2], bounds_with_y_split.rGetUpperCorner()[2], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], bounds_with_y_split.rGetLowerCorner()[0], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], bounds_with_y_split.rGetLowerCorner()[1], DBL_EPSILON);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[2], bounds_with_y_split.rGetLowerCorner()[2], DBL_EPSILON);

            // Same amount of stuff
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_with_x_split.GetNumNodes());
            TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_with_x_split.GetNumElements());
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_with_x_split.GetNumBoundaryElements());

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_with_y_split.GetNumNodes());
            TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_with_y_split.GetNumElements());
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_with_y_split.GetNumBoundaryElements());

            // Show that the x-split and y-split have different indexing schemes
            // Normal meshes index in x-direction first
            c_vector<double, 3> xaxis1 = mesh_with_default_split.GetNode(width)->rGetLocation();
            TS_ASSERT_DELTA(xaxis1[0], (double) width, 1e-5);
            TS_ASSERT_DELTA(xaxis1[1], 0.0, 1e-5);
            TS_ASSERT_DELTA(xaxis1[2], 0.0, 1e-5);
            // The x split one indexes in y, then z so x-axis is in top layer
            c_vector<double, 3> xaxis2 = mesh_with_x_split.GetNode(width*(height+1)*(depth+1))->rGetLocation();
            TS_ASSERT_DELTA(xaxis2[0], (double) width, 1e-5);
            TS_ASSERT_DELTA(xaxis2[1], 0.0, 1e-5);
            TS_ASSERT_DELTA(xaxis2[2], 0.0, 1e-5);

            // the y split indexes in z first so x-axis is the end of the first layer
            c_vector<double, 3> xaxis3 = mesh_with_y_split.GetNode(width*(depth+1))->rGetLocation();
            TS_ASSERT_DELTA(xaxis3[0], (double) width, 1e-5);
            TS_ASSERT_DELTA(xaxis3[1], 0.0, 1e-5);
            TS_ASSERT_DELTA(xaxis3[2], 0.0, 1e-5);

        }
    }
};
#endif //_TESTTETRAHEDRALMESH_HPP_
