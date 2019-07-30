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

#ifndef TESTMIXEDDIMENSIONMESH_HPP_
#define TESTMIXEDDIMENSIONMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp" // Must come before other Chaste headers

#include "MixedDimensionMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include "FileComparison.hpp"
#include "Exception.hpp"
#include "ArchiveOpener.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestMixedDimensionMesh : public CxxTest::TestSuite
{
public:
    void TestReadingSquareMesh()
    {
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);
        MixedDimensionMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 121u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 200u);

        /*
         * Cables:
        0       55      56      1.5
        1       56      57      2.5
        2       57      58      3.5
        3       58      59      4.5
        4       59      60      5.5
        5       60      61      6.5
        6       61      62      7.5
        7       62      63      8.5
        8       63      64      9.5
        9       64      65      10.5
         *
         */
        TS_ASSERT_EQUALS(mesh.GetNumCableElements(), 10u);

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(55u))
        {
            // Should be a single cable; number 0
            MixedDimensionMesh<2,2>::CableRangeAtNode cable_range = mesh.GetCablesAtNode(mesh.GetNode(55u));
            MixedDimensionMesh<2,2>::NodeCableIterator iter = cable_range.first;
            iter++;
            TS_ASSERT(iter == cable_range.second);
            TS_ASSERT_EQUALS(cable_range.first->second->GetIndex(), 0u);
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(56u))
        {
            // Should be two cables; numbers 0 & 1
            MixedDimensionMesh<2,2>::CableRangeAtNode cable_range = mesh.GetCablesAtNode(mesh.GetNode(56u));
            MixedDimensionMesh<2,2>::NodeCableIterator iter = cable_range.first;
            iter++;
            iter++;
            TS_ASSERT(iter == cable_range.second);
            std::set<unsigned> indices;
            indices.insert(cable_range.first->second->GetIndex());
            indices.insert((++cable_range.first)->second->GetIndex());
            TS_ASSERT_EQUALS(indices.size(), 2u);
            TS_ASSERT(indices.find(0u) != indices.end());
            TS_ASSERT(indices.find(1u) != indices.end());
        }

        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(mesh.GetNumLocalCableElements(), 10u);

            for (unsigned i=0; i<10u; i++)
            {
                Element<1,2>* p_cable_elt = mesh.GetCableElement(i);
                TS_ASSERT_EQUALS(p_cable_elt->GetNumNodes(), 2u);
                TS_ASSERT_EQUALS(p_cable_elt->GetNodeGlobalIndex(0u), 55u + i);
                TS_ASSERT_EQUALS(p_cable_elt->GetNodeGlobalIndex(1u), 56u + i);
                TS_ASSERT_EQUALS(p_cable_elt->GetNode(0u), mesh.GetNode(55u + i));
                TS_ASSERT_EQUALS(p_cable_elt->GetNode(1u), mesh.GetNode(56u + i));
                TS_ASSERT_DELTA(p_cable_elt->GetAttribute(), i+1.5, 1e-9);
                TS_ASSERT( mesh.CalculateDesignatedOwnershipOfCableElement(i) );
            }

            for (unsigned i=0; i<200u; i++)
            {
                Element<2,2>* p_elt = mesh.GetElement(i);
                TS_ASSERT_EQUALS(p_elt->GetNumNodes(), 3u);
                TS_ASSERT_EQUALS(p_elt->GetNode(0u), mesh.GetNode(p_elt->GetNodeGlobalIndex(0u)));
                TS_ASSERT_EQUALS(p_elt->GetNode(1u), mesh.GetNode(p_elt->GetNodeGlobalIndex(1u)));
                TS_ASSERT_EQUALS(p_elt->GetNode(2u), mesh.GetNode(p_elt->GetNodeGlobalIndex(2u)));
            }
        }
        else if (PetscTools::GetNumProcs() == 2)
        {
            if (PetscTools::GetMyRank() == 0)
            {
                // For a dumb partition, process 0 owns nodes 0 to 59
                // Cables 0, 1, 2, 3, 4, 5

                for (unsigned i=0; i<5u; i++)
                {
                    Element<1,2>* p_cable_elt = mesh.GetCableElement(i);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNumNodes(), 2u);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNodeGlobalIndex(0u), 55u + i);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNodeGlobalIndex(1u), 56u + i);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNode(0u), mesh.GetNode(55u + i));
                    TS_ASSERT_EQUALS(p_cable_elt->GetNode(1u), mesh.GetNodeOrHaloNode(56u + i));
                    TS_ASSERT_DELTA(p_cable_elt->GetAttribute(), i+1.5, 1e-9);
                    TS_ASSERT( mesh.CalculateDesignatedOwnershipOfCableElement(i) ); // Designated owner of all these five, since we own node 0 (lowest index)
                }
                TS_ASSERT_THROWS_THIS(mesh.GetCableElement(6), "Requested cable element 6 does not belong to processor 0");
            }
            else
            {
                // For a dumb partition, process 1 owns nodes 60 to 120
                // Cables  5, 6, 7, 8, 9
                TS_ASSERT_EQUALS(mesh.GetNumLocalCableElements(), 5u);
                TS_ASSERT_THROWS_THIS(mesh.GetCableElement(0), "Requested cable element 0 does not belong to processor 1");

                for (unsigned i=5; i<10u; i++)
                {
                    Element<1,2>* p_cable_elt = mesh.GetCableElement(i);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNumNodes(), 2u);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNodeGlobalIndex(0u), 55u + i);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNodeGlobalIndex(1u), 56u + i);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNode(0u), mesh.GetNodeOrHaloNode(55u + i));
                    TS_ASSERT_EQUALS(p_cable_elt->GetNode(1u), mesh.GetNode(56u + i));
                    TS_ASSERT_DELTA(p_cable_elt->GetAttribute(), i+1.5, 1e-9);

                    // Not designated owner of the first of these as node 0 is owned by process 0
                    if (i==5)
                    {
                        TS_ASSERT( ! mesh.CalculateDesignatedOwnershipOfCableElement(i) );
                    }
                    else
                    {
                        TS_ASSERT( mesh.CalculateDesignatedOwnershipOfCableElement(i) );
                    }

                }
            }
        }
        else if (PetscTools::GetNumProcs() == 3)
        {
            if (PetscTools::GetMyRank() == 0)
            {
                // For a dumb partition, process 0 owns nodes 0 to 39
                TS_ASSERT_EQUALS(mesh.GetNumLocalCableElements(), 0u);
            }
            if (PetscTools::GetMyRank() == 1)
            {
                // For a dumb partition, process 1 owns nodes 40 to 79
                TS_ASSERT_EQUALS(mesh.GetNumLocalCableElements(), 10u);
            }
            if (PetscTools::GetMyRank() == 2)
            {
                // For a dumb partition, process 0 owns nodes 80 to 120
                TS_ASSERT_EQUALS(mesh.GetNumLocalCableElements(), 0u);
            }
        }
        else
        {
            TS_ASSERT_LESS_THAN(mesh.GetNumLocalCableElements(), 11u);
        }
    }


    void TestReadingMeshWithNoCables()
    {
        std::string mesh_base("mesh/test/data/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);
        MixedDimensionMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 121u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 200u);
        TS_ASSERT_EQUALS(mesh.GetNumCableElements(), 0u);
    }

    void TestExceptions()
    {
        // Only TrianglesMeshReader supports cables
        MemfemMeshReader<3,3> memfem_reader("mesh/test/data/Memfem_slab");
        TS_ASSERT_EQUALS(memfem_reader.GetNumCableElements(), 0u);
        TS_ASSERT_EQUALS(memfem_reader.GetNumCableElementAttributes(), 0u);
        TS_ASSERT_THROWS_THIS(memfem_reader.GetNextCableElementData(), "Cable elements are not supported by this mesh format.");
    }


    void TestWritingCableFiles()
    {
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);
        MixedDimensionMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(reader);

        TrianglesMeshWriter<2,2> mesh_writer("TestMixedDimensionMesh", "CableMesh", true);

        mesh_writer.WriteFilesUsingMesh(mesh);

        FileFinder generated("TestMixedDimensionMesh/CableMesh.cable",RelativeTo::ChasteTestOutput);
        FileFinder reference("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements.cable", RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated,reference);
        TS_ASSERT(comparer.CompareFiles());
    }


    void TestWritingCableFilesUsingMeshReader()
    {
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);

        TrianglesMeshWriter<2,2> mesh_writer("TestMixedDimensionMesh", "CableMeshFromReader");
        mesh_writer.WriteFilesUsingMeshReader(reader);
        PetscTools::Barrier("TestWritingCableFilesUsingMeshReader");

        FileFinder generated("TestMixedDimensionMesh/CableMeshFromReader.cable",RelativeTo::ChasteTestOutput);
        FileFinder reference("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements.cable", RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated,reference);
        TS_ASSERT(comparer.CompareFiles());
    }

    void TestWritingBinaryFormat()
    {
        //Read as ascii
        TrianglesMeshReader<2,2> reader("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");

        //Write as binary
        TrianglesMeshWriter<2,2> writer_from_reader("TestMixedDimensionMesh", "CableMeshBinary", false);
        writer_from_reader.SetWriteFilesAsBinary();
        writer_from_reader.WriteFilesUsingMeshReader(reader);

        PetscTools::Barrier();
        //Read created binary file into a mesh
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestMixedDimensionMesh/";
        TrianglesMeshReader<2,2> binary_reader(results_dir + "CableMeshBinary");
        MixedDimensionMesh<2,2> binary_mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        binary_mesh.ConstructFromMeshReader(binary_reader);

        //Read original file into a mesh
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> original_reader(mesh_base);
        MixedDimensionMesh<2,2> original_mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        original_mesh.ConstructFromMeshReader(original_reader);

        //Compare to original
        TS_ASSERT_EQUALS(binary_mesh.GetNumNodes(), original_mesh.GetNumNodes());
        TS_ASSERT_EQUALS(binary_mesh.GetNumElements(), original_mesh.GetNumElements());

        TS_ASSERT_EQUALS(binary_mesh.GetNumCableElements(), original_mesh.GetNumCableElements());

        MixedDimensionMesh<2,2>::CableElementIterator original_iter = original_mesh.GetCableElementIteratorBegin();
        for (MixedDimensionMesh<2,2>::CableElementIterator binary_iter = binary_mesh.GetCableElementIteratorBegin();
             binary_iter != binary_mesh.GetCableElementIteratorEnd();
             ++binary_iter)
        {
            TS_ASSERT_EQUALS((*binary_iter)->GetNumNodes(), (*original_iter)->GetNumNodes());
            TS_ASSERT_EQUALS((*binary_iter)->GetNodeGlobalIndex(0u), (*original_iter)->GetNodeGlobalIndex(0u));
            TS_ASSERT_EQUALS((*binary_iter)->GetNodeGlobalIndex(1u), (*original_iter)->GetNodeGlobalIndex(1u));
            TS_ASSERT_DELTA((*binary_iter)->GetAttribute(), (*original_iter)->GetAttribute(), 1e-12);

            ++original_iter;
        }

        //Write a binary from the original mesh
        TrianglesMeshWriter<2,2> writer_from_mesh("TestMixedDimensionMesh", "CableMeshBinaryFromMesh", false);
        writer_from_mesh.SetWriteFilesAsBinary();
        writer_from_mesh.WriteFilesUsingMesh(original_mesh);

        //Compare the binary written from the reader to the binary written from the mesh
        FileFinder generated(results_dir + "/CableMeshBinary.cable");
        FileFinder reference(results_dir + "/CableMeshBinaryFromMesh.cable");
        FileComparison comparer(generated,reference);
        TS_ASSERT(comparer.CompareFiles());
    }


    void TestGeometryWithMetisPermuation()
    {
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);
        MixedDimensionMesh<2,2> dumb_partition_mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        dumb_partition_mesh.ConstructFromMeshReader(reader);


        for (MixedDimensionMesh<1,2>::CableElementIterator iter = dumb_partition_mesh.GetCableElementIteratorBegin();
             iter != dumb_partition_mesh.GetCableElementIteratorEnd();
             ++iter)
        {
            Element<1,2>& r_element = *(*iter);

            c_matrix<double, 2, 1> jacobian;
            c_matrix<double, 1, 2> inverse_jacobian;
            double jacobian_determinant;

            r_element.CalculateInverseJacobian(jacobian, jacobian_determinant, inverse_jacobian);

            TS_ASSERT_DELTA(jacobian(0,0), 0.01, 1e-6);
            TS_ASSERT_DELTA(jacobian(1,0), 0.0,  1e-6);
            TS_ASSERT_DELTA(inverse_jacobian(0,0), 100, 1e-6);
            TS_ASSERT_DELTA(inverse_jacobian(0,1), 0.0, 1e-6);
            TS_ASSERT_DELTA(jacobian_determinant, 0.01, 1e-6);

            // y value of all the cable nodes should be 0.05
            TS_ASSERT_DELTA(r_element.GetNodeLocation(0,1), 0.05, 1e-6);
            TS_ASSERT_DELTA(r_element.GetNodeLocation(1,1), 0.05, 1e-6);
            // x value at the ends of the cable element can also be anticipated
            TS_ASSERT_DELTA(r_element.GetNodeLocation(0,0), r_element.GetIndex()*0.01, 1e-6);
            TS_ASSERT_DELTA(r_element.GetNodeLocation(1,0), (r_element.GetIndex()+1)*0.01, 1e-6);

        }

        //Test that every cable element has a designated owner
        {
            unsigned local_owned=0u;
            for (unsigned i=0; i<dumb_partition_mesh.GetNumCableElements(); i++)
            {
                if (dumb_partition_mesh.CalculateDesignatedOwnershipOfCableElement(i))
                {
                    local_owned++;
                }
            }
            unsigned total_owned;
            MPI_Allreduce(&local_owned, &total_owned, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
            TS_ASSERT_EQUALS(total_owned, dumb_partition_mesh.GetNumCableElements());
        }

        TrianglesMeshReader<2,2> reader2(mesh_base);
        MixedDimensionMesh<2,2> partitioned_mesh;
        partitioned_mesh.ConstructFromMeshReader(reader2);


        for (MixedDimensionMesh<1,2>::CableElementIterator iter = partitioned_mesh.GetCableElementIteratorBegin();
             iter != partitioned_mesh.GetCableElementIteratorEnd();
             ++iter)
        {
            Element<1,2>& r_element = *(*iter);

            c_matrix<double, 2, 1> jacobian;
            c_matrix<double, 1, 2> inverse_jacobian;
            double jacobian_determinant;

            r_element.CalculateInverseJacobian(jacobian, jacobian_determinant, inverse_jacobian);

            TS_ASSERT_DELTA(jacobian(0,0), 0.01, 1e-6);
            TS_ASSERT_DELTA(jacobian(1,0), 0.0,  1e-6);
            TS_ASSERT_DELTA(inverse_jacobian(0,0), 100, 1e-6);
            TS_ASSERT_DELTA(inverse_jacobian(0,1), 0.0, 1e-6);
            TS_ASSERT_DELTA(jacobian_determinant, 0.01, 1e-6);

            // y value of all the cable nodes should be 0.05
            TS_ASSERT_DELTA(r_element.GetNodeLocation(0,1), 0.05, 1e-6);
            TS_ASSERT_DELTA(r_element.GetNodeLocation(1,1), 0.05, 1e-6);
            // x value at the ends of the cable element can also be anticipated
            TS_ASSERT_DELTA(r_element.GetNodeLocation(0,0), r_element.GetIndex()*0.01, 1e-6);
            TS_ASSERT_DELTA(r_element.GetNodeLocation(1,0), (r_element.GetIndex()+1)*0.01, 1e-6);

        }
        //Test that every cable element has a designated owner
        {
            unsigned local_owned=0u;
            for (unsigned i=0; i<partitioned_mesh.GetNumCableElements(); i++)
            {
                if (partitioned_mesh.CalculateDesignatedOwnershipOfCableElement(i))
                {
                    local_owned++;
                }
            }
            unsigned total_owned;
            MPI_Allreduce(&local_owned, &total_owned, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
            TS_ASSERT_EQUALS(total_owned, partitioned_mesh.GetNumCableElements());
        }
    }

    void TestArchiving()
    {
        OutputFileHandler archive_dir_("mixed_mesh_archive"); // Clear folder
        FileFinder main_archive_dir("mixed_mesh_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "mixed_dimension_mesh.arch";
        ArchiveLocationInfo::SetMeshFilename("mixed_dimension_mesh");

        MixedDimensionMesh<2,2>* p_mesh = new MixedDimensionMesh<2,2>(DistributedTetrahedralMeshPartitionType::DUMB);
        unsigned num_nodes;
        unsigned local_num_nodes;
        unsigned num_elements;
        unsigned num_cable_elements;
        unsigned num_local_cable_elements;
        // archive
        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");

            p_mesh->ConstructFromMeshReader(mesh_reader);
            num_nodes = p_mesh->GetNumNodes();
            local_num_nodes = p_mesh->GetNumLocalNodes();
            num_elements = p_mesh->GetNumElements();
            num_cable_elements = p_mesh->GetNumCableElements();
            num_local_cable_elements = p_mesh->GetNumLocalCableElements();

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(main_archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractTetrahedralMesh<2,2>* const p_mesh_abstract = static_cast<AbstractTetrahedralMesh<2,2>* >(p_mesh);
            (*p_arch) << p_mesh_abstract;
        }

        FileFinder ncl_file("mixed_dimension_mesh.ncl", main_archive_dir);
        TS_ASSERT(ncl_file.Exists());

        // restore
        {
            // Should archive the most abstract class you can to check boost knows what individual classes are.
            // (but here AbstractMesh doesn't have the methods below).
            AbstractTetrahedralMesh<2,2>* p_mesh_abstract2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(main_archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // restore from the archive
            (*p_arch) >> p_mesh_abstract2;
            // Check we have the right number of nodes & elements
            MixedDimensionMesh<2,2>* p_mesh2 = static_cast<MixedDimensionMesh<2,2>*>(p_mesh_abstract2);

            TS_ASSERT_EQUALS(p_mesh2->GetNumNodes(), num_nodes);
            TS_ASSERT_EQUALS(p_mesh2->GetNumLocalNodes(), local_num_nodes);
            TS_ASSERT_EQUALS(p_mesh2->GetNumElements(), num_elements);
            TS_ASSERT_EQUALS(p_mesh2->GetNumCableElements(), num_cable_elements);
            TS_ASSERT_EQUALS(p_mesh2->GetNumLocalCableElements(), num_local_cable_elements);

            // Check elements have the right nodes
            for (unsigned i=0; i<num_cable_elements; i++)
            {
                try
                {
                    Element<1,2>* p_element = p_mesh->GetCableElement(i);
                    Element<1,2>* p_element2 = p_mesh2->GetCableElement(i);
                    TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), p_element2->GetNodeGlobalIndex(0));
                }
                catch(Exception& e)
                {
                    TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
                }
            }
            delete p_mesh2;
        }

        // restore from a single processor archive
        {
            FileFinder archive_dir("mesh/test/data/mixed_mesh_archive", RelativeTo::ChasteSourceRoot);
            if (PetscTools::IsSequential())
            {
                ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
                boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
                AbstractTetrahedralMesh<2,2>* p_mesh_abstract3 = NULL;
                (*p_arch) >> p_mesh_abstract3;

                //Double check that the cables are intact
                MixedDimensionMesh<2,2>* p_mesh3 = static_cast<MixedDimensionMesh<2,2>*>(p_mesh_abstract3);
                Element<1,2>* p_element = p_mesh3->GetCableElement(9);
                TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 64u);
                TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 65u);
                TS_ASSERT_DELTA(p_element->GetAttribute(), 10.5, 1e-8);

                delete p_mesh_abstract3;
            }
            else
            {
                typedef ArchiveOpener<boost::archive::text_iarchive, std::ifstream> InputArchiveOpener;
                if (PetscTools::GetMyRank() > 0)
                {
                    // Should not read this archive because none exists here.
                    TS_ASSERT_THROWS_CONTAINS(InputArchiveOpener arch_opener(archive_dir, archive_file),
                                              "Cannot load secondary archive file:");
                }
                else
                {
                    // Should not read this archive because there are two or more processes and
                    // this archive was written on one process.
                    InputArchiveOpener arch_opener(archive_dir, archive_file);
                    boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
                    AbstractTetrahedralMesh<2,2>* p_mesh3 = NULL;
                    TS_ASSERT_THROWS_THIS((*p_arch) >> p_mesh3,
                                          "This archive was written for a different number of processors");

                }
            }
        }

        delete p_mesh;
    }
};

#endif /*TESTMIXEDDIMENSIONMESH_HPP_*/
