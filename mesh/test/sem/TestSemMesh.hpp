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
#include "ArchiveLocationInfo.hpp"
#include "OutputFileHandler.hpp"
#include "FileFinder.hpp"
#include "SemSingleElementMeshGenerator.hpp"
#include "SemMultiElementMeshGenerator.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSemMesh : public CxxTest::TestSuite
{
public:



    void TestDefaultConstructorAndDestructor()
    {
        SemMesh<2> mesh;

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 0u);

        // Pin the SEM defaults. These are all archived, so a silent change here would alter the
        // behaviour of any restarted simulation.
        TS_ASSERT_DELTA(mesh.GetMaximumInteractionDistance(), 1.0, 1e-12);
        TS_ASSERT_EQUALS(mesh.GetOutputElementSurfacesToVtk(), true);
        TS_ASSERT_DELTA(mesh.GetSemSurfaceAlphaMultiplier(), 2.0, 1e-12);
        TS_ASSERT_DELTA(mesh.GetSemSurfaceExpansionMultiplier(), 0.5, 1e-12);
        TS_ASSERT_EQUALS(mesh.GetUseExpandedSemSurfaceForVolume(), true);

        // An empty mesh must give an immediately-exhausted element iterator
        TS_ASSERT(!(mesh.GetElementIteratorBegin() != mesh.GetElementIteratorEnd()));
    }

    void TestGetNumNodesAndElements()
    {
        // Two elements of 3x3 nodes each
        SemMultiElementMeshGenerator<2> generator({3, 3}, {2, 1}, 0.5);
        auto p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 18u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), 2u);

        // Node indices are the identity mapping, which is what SolveNodeMapping guarantees
        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            TS_ASSERT_EQUALS(p_mesh->GetNode(i)->GetIndex(), i);
        }

        unsigned num_iterated = 0;
        for (auto iter = p_mesh->GetElementIteratorBegin();
             iter != p_mesh->GetElementIteratorEnd();
             ++iter)
        {
            num_iterated++;
        }
        TS_ASSERT_EQUALS(num_iterated, 2u);

        /*
         * Deleting an element does NOT change either count: SemMesh never removes an element or
         * reuses an index, because SemBasedCellPopulation uses the element index as its location
         * index. Only the iterator reflects the deletion. This is the trap that makes
         * GetNumElements() an upper bound on the element index rather than a count of live cells.
         */
        p_mesh->GetElement(0)->MarkAsDeleted();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), 2u);

        num_iterated = 0;
        for (auto iter = p_mesh->GetElementIteratorBegin();
             iter != p_mesh->GetElementIteratorEnd();
             ++iter)
        {
            TS_ASSERT_EQUALS(iter->IsDeleted(), false);
            num_iterated++;
        }
        TS_ASSERT_EQUALS(num_iterated, 1u);

        // The two names are synonyms, which is why neither excludes deleted elements
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), p_mesh->GetNumAllElements());
    }

    void TestGetElementMethods()
    {
        SemMultiElementMeshGenerator<2> generator({3, 3}, {2, 1}, 0.5);
        auto p_mesh = generator.GetMesh();

        for (unsigned e = 0; e < p_mesh->GetNumElements(); ++e)
        {
            SemElement<2>* p_element = p_mesh->GetElement(e);
            TS_ASSERT_EQUALS(p_element->GetIndex(), e);
            TS_ASSERT_EQUALS(p_element->GetNumNodes(), 9u);
        }

        // The iterator visits the elements in index order
        unsigned expected_index = 0;
        for (auto iter = p_mesh->GetElementIteratorBegin();
             iter != p_mesh->GetElementIteratorEnd();
             ++iter)
        {
            TS_ASSERT_EQUALS(iter->GetIndex(), expected_index);
            expected_index++;
        }
        TS_ASSERT_EQUALS(expected_index, 2u);

        // AddElement returns the index of the element it appended
        std::vector<Node<2>*> new_nodes;
        new_nodes.push_back(p_mesh->GetNode(0));
        new_nodes.push_back(p_mesh->GetNode(1));
        SemElement<2>* p_new_element = new SemElement<2>(2, new_nodes);

        TS_ASSERT_EQUALS(p_mesh->AddElement(p_new_element), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 3u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2), p_new_element);
    }

    void TestClear()
    {
        SemMultiElementMeshGenerator<2> generator({3, 3}, {2, 1}, 0.5);
        auto p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 18u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);

        p_mesh->Clear();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), 0u);

        // Clearing an already-empty mesh must be safe: both constructors call Clear(), so this
        // path runs on every mesh that is ever built
        TS_ASSERT_THROWS_NOTHING(p_mesh->Clear());

        // The mesh is reusable afterwards
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        p_mesh->AddNode(nodes[0]);
        p_mesh->AddNode(nodes[1]);
        p_mesh->AddElement(new SemElement<2>(0, nodes));

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 1u);
    }

    void TestGetMeshForVtk()
    {
        SemSingleElementMeshGenerator<2> generator({3, 3}, 1.0);
        auto p_mesh = generator.GetMesh();

        // SemMesh needs no separate VTK representation, so this returns the mesh itself
        TS_ASSERT_EQUALS(p_mesh->GetMeshForVtk(), p_mesh.get());
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
    void TestArchiveSemMesh()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "sem_mesh_2d.arch";
        ArchiveLocationInfo::SetMeshFilename("sem_mesh_2d");

        SemMultiElementMeshGenerator<2> generator({3, 3}, {2, 1}, 0.5);
        auto p_generated_mesh = generator.GetMesh();

        /*
         * Move all five SEM scalars off their default values before archiving. If they were left
         * at the defaults this test would still pass with save()/load() dropping them entirely.
         */
        p_generated_mesh->SetMaximumInteractionDistance(1.75);
        p_generated_mesh->SetOutputElementSurfacesToVtk(false);
        p_generated_mesh->SetSemSurfaceAlphaMultiplier(3.25);
        p_generated_mesh->SetSemSurfaceExpansionMultiplier(0.125);
        p_generated_mesh->SetUseExpandedSemSurfaceForVolume(false);

        AbstractMesh<2, 2>* const p_mesh = p_generated_mesh.get();

        {
            // The const above is needed to stop a BOOST_STATIC_ASSERTION failure
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Serialize via a pointer, or the derived class information is lost
            (*p_arch) << p_mesh;
        }

        {
            AbstractMesh<2, 2>* p_loaded_mesh;

            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) >> p_loaded_mesh;

            SemMesh<2>* p_original = static_cast<SemMesh<2>*>(p_mesh);
            SemMesh<2>* p_loaded = static_cast<SemMesh<2>*>(p_loaded_mesh);

            // The SEM scalar state survives
            TS_ASSERT_DELTA(p_loaded->GetMaximumInteractionDistance(), 1.75, 1e-9);
            TS_ASSERT_EQUALS(p_loaded->GetOutputElementSurfacesToVtk(), false);
            TS_ASSERT_DELTA(p_loaded->GetSemSurfaceAlphaMultiplier(), 3.25, 1e-9);
            TS_ASSERT_DELTA(p_loaded->GetSemSurfaceExpansionMultiplier(), 0.125, 1e-9);
            TS_ASSERT_EQUALS(p_loaded->GetUseExpandedSemSurfaceForVolume(), false);

            // Nodes survive, with their positions and their region labels
            TS_ASSERT_EQUALS(p_loaded->GetNumNodes(), p_original->GetNumNodes());
            for (unsigned i = 0; i < p_original->GetNumNodes(); ++i)
            {
                Node<2>* p_orig_node = p_original->GetNode(i);
                Node<2>* p_loaded_node = p_loaded->GetNode(i);

                TS_ASSERT_EQUALS(p_loaded_node->GetIndex(), p_orig_node->GetIndex());
                TS_ASSERT_EQUALS(p_loaded_node->GetRegion(), p_orig_node->GetRegion());
                for (unsigned dim = 0; dim < 2; ++dim)
                {
                    TS_ASSERT_DELTA(p_loaded_node->rGetLocation()[dim], p_orig_node->rGetLocation()[dim], 1e-4);
                }
            }

            // Elements survive, with their node membership rebuilt by RegisterWithNodes()
            TS_ASSERT_EQUALS(p_loaded->GetNumElements(), p_original->GetNumElements());
            for (unsigned e = 0; e < p_original->GetNumElements(); ++e)
            {
                SemElement<2>* p_orig_element = p_original->GetElement(e);
                SemElement<2>* p_loaded_element = p_loaded->GetElement(e);

                TS_ASSERT_EQUALS(p_loaded_element->GetNumNodes(), p_orig_element->GetNumNodes());
                for (unsigned j = 0; j < p_orig_element->GetNumNodes(); ++j)
                {
                    TS_ASSERT_EQUALS(p_loaded_element->GetNodeGlobalIndex(j),
                                     p_orig_element->GetNodeGlobalIndex(j));
                    TS_ASSERT_EQUALS(p_loaded_element->GetNode(j)->rGetContainingElementIndices().count(e), 1u);
                }
            }

            /*
             * The box collection is deliberately not archived; SemMesh::load() rebuilds it from
             * the loaded node locations using the restored interaction distance. Nothing else in
             * this test would notice if that were dropped, so check it directly -- both calls
             * below assert on a null collection.
             *
             * Note that SetUpBoxCollection() creates the boxes but does not put the nodes in them;
             * UpdateBoxCollection() does that, which is why SemBasedCellPopulation::Update() calls
             * it before every CalculateNodePairs(). Calling CalculateNodePairs() straight after a
             * load, without the update, correctly yields nothing.
             */
            std::vector<std::pair<Node<2>*, Node<2>*> > node_pairs;
            TS_ASSERT_THROWS_NOTHING(p_loaded->UpdateBoxCollection());
            TS_ASSERT_THROWS_NOTHING(p_loaded->CalculateNodePairs(node_pairs));
            TS_ASSERT(!node_pairs.empty());

            delete p_loaded_mesh;
        }
    }
};

#endif /*TESTSEMMESH_HPP_*/
