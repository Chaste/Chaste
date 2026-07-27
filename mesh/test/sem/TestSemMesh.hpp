/*

Copyright (c) 2005-2026, University of Oxford.
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
#include "SemElementGeometry.hpp"
#include "SemMeshReader.hpp"
#include "SemMeshWriter.hpp"
#include "ArchiveOpener.hpp"
#include "OutputFileHandler.hpp"
#include "SemSingleElementMeshGenerator.hpp"
#include "SemMultiElementMeshGenerator.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSemMesh : public CxxTest::TestSuite
{
public:

    void TestSolveNodeMapping()
    {
        ///\todo
    }

    void TestSolveElementMapping()
    {
        ///\todo
    }

    void TestSolveBoundaryElementMapping()
    {
        ///\todo
    }

    void TestDefaultConstructorAndDestructor()
    {
        ///\todo
    }

    void TestGetNumNodes()
    {
        ///\todo
    }

    void TestGetNumElements()
    {
        ///\todo
    }

    void TestGetNumAllElements()
    {
        ///\todo
    }

    void TestGetElementMethods()
    {
        ///\todo
    }

    void TestGetCentroidOfElement()
    {
        // 1D
        {
            SemSingleElementMeshGenerator<1> generator({5}, 2.0);
            auto p_mesh = generator.GetMesh();

            c_vector<double, 1> centroid = p_mesh->GetCentroidOfElement(0);
            TS_ASSERT_DELTA(centroid[0], 0.8, 1e-6);
        }
        // 2D
        {
            SemSingleElementMeshGenerator<2> generator({5, 3}, 4.0);
            auto p_mesh = generator.GetMesh();

            c_vector<double, 2> centroid = p_mesh->GetCentroidOfElement(0);
            TS_ASSERT_DELTA(centroid[0], 1.6, 1e-6);
            TS_ASSERT_DELTA(centroid[1], 0.8, 1e-6);
        }
        // 2D
        {
            SemSingleElementMeshGenerator<3> generator({5, 3, 7}, 6.0);
            auto p_mesh = generator.GetMesh();

            c_vector<double, 3> centroid = p_mesh->GetCentroidOfElement(0);
            TS_ASSERT_DELTA(centroid[0], 2.4, 1e-6);
            TS_ASSERT_DELTA(centroid[1], 1.2, 1e-6);
            TS_ASSERT_DELTA(centroid[2], 3.6, 1e-6);
        }
    }

    void TestConstructFromMeshReader()
    {
        // Build a two-element 2D mesh; 3x3 nodes per element gives a mix of interior
        // (region 0) and boundary (region 1) nodes to exercise region persistence.
        SemMultiElementMeshGenerator<2> generator({3, 3}, {2, 1}, 0.5);
        auto p_original_mesh = generator.GetMesh();

        const unsigned num_nodes = p_original_mesh->GetNumNodes();
        const unsigned num_elements = p_original_mesh->GetNumElements();

        // Sanity check that both regions are actually present in the original mesh
        bool has_interior = false;
        bool has_boundary = false;
        for (unsigned i = 0; i < num_nodes; ++i)
        {
            has_interior = has_interior || (p_original_mesh->GetNode(i)->GetRegion() == 0u);
            has_boundary = has_boundary || (p_original_mesh->GetNode(i)->GetRegion() == 1u);
        }
        TS_ASSERT(has_interior);
        TS_ASSERT(has_boundary);

        // Write the mesh to file, then read it back into a fresh mesh
        SemMeshWriter<2> writer("TestSemMeshRoundTrip", "sem_mesh", true);
        writer.WriteFilesUsingMesh(*p_original_mesh);

        std::string mesh_base = OutputFileHandler::GetChasteTestOutputDirectory() + "TestSemMeshRoundTrip/sem_mesh";
        SemMeshReader<2> reader(mesh_base);
        SemMesh<2> read_mesh;
        read_mesh.ConstructFromMeshReader(reader);

        // Element and node counts are preserved
        TS_ASSERT_EQUALS(read_mesh.GetNumNodes(), num_nodes);
        TS_ASSERT_EQUALS(read_mesh.GetNumElements(), num_elements);

        // Node positions and regions are preserved
        for (unsigned i = 0; i < num_nodes; ++i)
        {
            Node<2>* p_orig = p_original_mesh->GetNode(i);
            Node<2>* p_read = read_mesh.GetNode(i);
            TS_ASSERT_DELTA(p_read->rGetLocation()[0], p_orig->rGetLocation()[0], 1e-6);
            TS_ASSERT_DELTA(p_read->rGetLocation()[1], p_orig->rGetLocation()[1], 1e-6);
            TS_ASSERT_EQUALS(p_read->GetRegion(), p_orig->GetRegion());
        }

        // Per-element node count and node membership are preserved
        for (unsigned e = 0; e < num_elements; ++e)
        {
            SemElement<2>* p_orig = p_original_mesh->GetElement(e);
            SemElement<2>* p_read = read_mesh.GetElement(e);

            TS_ASSERT_EQUALS(p_read->GetNumNodes(), p_orig->GetNumNodes());

            for (unsigned j = 0; j < p_orig->GetNumNodes(); ++j)
            {
                TS_ASSERT_EQUALS(p_read->GetNodeGlobalIndex(j), p_orig->GetNodeGlobalIndex(j));
            }
        }

        // Node -> element membership is rebuilt by RegisterWithNodes() during construction
        for (unsigned e = 0; e < num_elements; ++e)
        {
            SemElement<2>* p_read = read_mesh.GetElement(e);
            for (unsigned j = 0; j < p_read->GetNumNodes(); ++j)
            {
                std::set<unsigned> containing_elements = p_read->GetNode(j)->rGetContainingElementIndices();
                TS_ASSERT_EQUALS(containing_elements.count(e), 1u);
            }
        }
    }

    void TestClear()
    {
        ///\todo
    }

    void TestGetVolumeOfElement()
    {
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0, false, 0.0));
        nodes.push_back(new Node<1>(1, false, 1.0));
        std::vector<Node<1>*> element_nodes = nodes;
        std::vector<SemElement<1>*> elements;
        elements.push_back(new SemElement<1>(0, element_nodes));
        SemMesh<1> mesh(nodes, elements);

        mesh.SetUseExpandedSemSurfaceForVolume(false);
        TS_ASSERT_DELTA(mesh.GetVolumeOfElement(0), 1.0, 1e-6);

        mesh.SetUseExpandedSemSurfaceForVolume(true);
        mesh.SetSemSurfaceExpansionMultiplier(0.5);
        TS_ASSERT_DELTA(mesh.GetVolumeOfElement(0), 2.0, 1e-6);

#ifdef CHASTE_VTK
        std::vector<Node<2>*> nodes_2d;
        nodes_2d.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes_2d.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes_2d.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes_2d.push_back(new Node<2>(3, false, 0.0, 1.0));
        std::vector<Node<2>*> element_nodes_2d = nodes_2d;
        std::vector<SemElement<2>*> elements_2d;
        elements_2d.push_back(new SemElement<2>(0, element_nodes_2d));
        SemMesh<2> mesh_2d(nodes_2d, elements_2d);
        mesh_2d.SetUseExpandedSemSurfaceForVolume(false);
        mesh_2d.SetSemSurfaceExpansionMultiplier(0.0);
        TS_ASSERT_DELTA(mesh_2d.GetVolumeOfElement(0), 1.0, 1e-6);
#else
        std::vector<Node<2>*> nodes_2d;
        nodes_2d.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes_2d.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes_2d.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes_2d.push_back(new Node<2>(3, false, 0.0, 1.0));
        std::vector<Node<2>*> element_nodes_2d = nodes_2d;
        std::vector<SemElement<2>*> elements_2d;
        elements_2d.push_back(new SemElement<2>(0, element_nodes_2d));
        SemMesh<2> mesh_2d(nodes_2d, elements_2d);
        TS_ASSERT_THROWS_CONTAINS(mesh_2d.GetVolumeOfElement(0), "require Chaste to be compiled with VTK");
#endif // CHASTE_VTK
    }

    void TestGetVolumeOfElementThrowsForUnderDefinedElement()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        std::vector<Node<2>*> element_nodes = nodes;
        std::vector<SemElement<2>*> elements;
        elements.push_back(new SemElement<2>(0, element_nodes));
        SemMesh<2> mesh(nodes, elements);

        TS_ASSERT_THROWS_CONTAINS(mesh.GetVolumeOfElement(0), "requires at least three distinct nodes");
    }

    void TestLocalSpacingUsesBoundingBoxEstimate()
    {
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0, false, 0.0));
        nodes.push_back(new Node<1>(1, false, 1.0));
        nodes.push_back(new Node<1>(2, false, 3.0));
        std::vector<Node<1>*> element_nodes = nodes;
        std::vector<SemElement<1>*> elements;
        elements.push_back(new SemElement<1>(0, element_nodes));
        SemMesh<1> mesh(nodes, elements);

        SemElementSurface<1> surface = SemElementGeometry<1>::GenerateSurface(mesh, 0, 2.0, 0.0, false);
        TS_ASSERT_DELTA(surface.LocalSpacing, 1.5, 1e-6);
        TS_ASSERT_DELTA(surface.Alpha, 3.0, 1e-6);
    }

    void TestLocalSpacingUsesBoundingBoxEstimateIn2d()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 4.0, 0.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 0.5));
        nodes.push_back(new Node<2>(4, false, 2.0, 0.5));
        nodes.push_back(new Node<2>(5, false, 4.0, 0.5));
        nodes.push_back(new Node<2>(6, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(7, false, 2.0, 1.0));
        nodes.push_back(new Node<2>(8, false, 4.0, 1.0));

        std::vector<Node<2>*> element_nodes = nodes;
        std::vector<SemElement<2>*> elements;
        elements.push_back(new SemElement<2>(0, element_nodes));

        SemMesh<2> mesh(nodes, elements);

        TS_ASSERT_DELTA(SemElementGeometry<2>::CalculateLocalSpacing(mesh, 0), 1.0, 1e-6);
    }
    void TestGetMeshForVtk()
    {
        ///\todo
    }
};

#endif /*TESTSEMMESH_HPP_*/
