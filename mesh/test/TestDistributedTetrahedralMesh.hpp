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

#ifndef TESTDISTRIBUTEDTETRAHEDRALMESH_HPP_
#define TESTDISTRIBUTEDTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include <sstream>
#include <boost/scoped_array.hpp>

#include "UblasCustomFunctions.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "PetscTools.hpp"
#include "ArchiveOpener.hpp"
#include "FileFinder.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "CmguiMeshWriter.hpp"
#include "FileComparison.hpp"

#include "RandomNumberGenerator.hpp"
#include "Warnings.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestDistributedTetrahedralMesh : public CxxTest::TestSuite
{
private:

    void ReadFileToSetOfLines(const std::string& rFilePath, std::set<std::string>& rSetOfLines)
    {
        std::ifstream filestream(rFilePath.c_str());
        TS_ASSERT(filestream.is_open());
        while (filestream.good())
        {
            std::string line;
            getline(filestream, line);
            if (filestream.fail())
            {
                break;
            }
            if (!(line[0] == '#' || line[0] == '!' || line.substr(0, 10) == "Group name"))
            {
                //Even though both files were created with the same build, they may have slightly different creation
                //times in their provenance line (so we ignore it).
                rSetOfLines.insert(line);
            }
        }
    }

    void ComparePermutedFiles(const std::string& rFilePath1, const std::string& rFilePath2)
    {
        if (!PetscTools::AmMaster())
        {
            //Only the master needs to do this
            return;
        }
        std::set<std::string> lines_from_file_1;
        ReadFileToSetOfLines(rFilePath1, lines_from_file_1);

        std::set<std::string> lines_from_file_2;
        ReadFileToSetOfLines(rFilePath2, lines_from_file_2);

        // Do both files contain the same lines (up to permutations)?
        TS_ASSERT(lines_from_file_1 == lines_from_file_2);
    }

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    void CompareMeshes( DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh1,
                        DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh2 )
    {
        // Check that the same partitioner was used
        TS_ASSERT_EQUALS(rMesh1.mPartitioning, rMesh2.mPartitioning);
        // Check that we have the right number of nodes and elements
        TS_ASSERT_EQUALS(rMesh1.GetNumBoundaryElements(), rMesh2.GetNumBoundaryElements());
        TS_ASSERT_EQUALS(rMesh1.GetNumElements(), rMesh2.GetNumElements());
        TS_ASSERT_EQUALS(rMesh1.GetNumNodes(), rMesh2.GetNumNodes());
        TS_ASSERT_EQUALS(rMesh1.GetNumLocalNodes(), rMesh2.GetNumLocalNodes());
        TS_ASSERT_EQUALS(rMesh1.GetNumLocalElements(), rMesh2.GetNumLocalElements());

        // Check that the nodes and elements of each mesh are identical
        for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator iter = rMesh1.GetElementIteratorBegin();
             iter != rMesh1.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();

            Element<ELEMENT_DIM,SPACE_DIM>* p_element_2 = rMesh2.GetElement(element_index);

            // The elements have the same index and the nodes are located in the same position.
            TS_ASSERT_EQUALS(element_index, p_element_2->GetIndex());
            for (unsigned node_local_index=0; node_local_index < iter->GetNumNodes(); node_local_index++)
            {
                TS_ASSERT_DELTA( norm_2( iter->GetNode(node_local_index)->rGetLocation() -
                                     p_element_2->GetNode(node_local_index)->rGetLocation() ), 0.0, 1e-10 );
            }
        }
    }

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    void CheckEverythingIsAssigned(DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
    {
        /*
         * Check for consistent partitions (i.e. you own or "halo-own" every node in every element you own.
         */
        for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator iter = rMesh.GetElementIteratorBegin();
             iter != rMesh.GetElementIteratorEnd();
             ++iter)
        {
            for (unsigned node_local_index=0; node_local_index<ELEMENT_DIM+1; node_local_index++)
            {
                unsigned node_global_index = iter->GetNodeGlobalIndex(node_local_index);

                TS_ASSERT_THROWS_NOTHING(rMesh.GetNodeOrHaloNode(node_global_index));
            }
        }

        /*
         * Check that nodes are numbered consecutively
         */
        typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator prev_node = rMesh.GetNodeIteratorBegin();
        for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator current_node = ++rMesh.GetNodeIteratorBegin();
             current_node != rMesh.GetNodeIteratorEnd();
             ++prev_node, ++current_node)
        {
            TS_ASSERT_EQUALS(prev_node->GetIndex()+1, current_node->GetIndex())
        }

        /*
         * All the nodes have been assigned
         */
        unsigned total_nodes_this_process = 0;
        {
            const unsigned num_global_nodes = rMesh.GetNumNodes();
            boost::scoped_array<unsigned> nodes_owned(new unsigned[num_global_nodes]);
            for (unsigned index=0; index<num_global_nodes; index++)
            {
                nodes_owned[index]=0u;
            }

            for (unsigned node_id=0; node_id<num_global_nodes;  node_id++)
            {

                try
                {
                     unsigned node_index = rMesh.GetNode(node_id)->GetIndex();
                     TS_ASSERT_EQUALS(node_id, node_index);
                     nodes_owned[node_index] = 1;
                     total_nodes_this_process++;
                }
                catch (Exception &)
                {
                    nodes_owned[node_id] = 0;
                }
            }

            TS_ASSERT_EQUALS(rMesh.GetNumLocalNodes(), total_nodes_this_process);

            // Combine all the local maps by adding them up in the master process
            boost::scoped_array<unsigned> nodes_reduction(new unsigned[num_global_nodes]);
            MPI_Reduce(nodes_owned.get(), nodes_reduction.get(), num_global_nodes, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);

            // Make sure every node is owned at least by one processor
            if (PetscTools::AmMaster())
            {
                for (unsigned node_id=0; node_id<num_global_nodes; node_id++)
                {
                    TS_ASSERT(nodes_reduction[node_id] > 0u);
                }
            }

        }

        /*
         * All elements have been assigned
         */
        unsigned total_elements_this_process = 0;
        {
            const unsigned num_global_elements = rMesh.GetNumElements();
            boost::scoped_array<unsigned> elements_owned(new unsigned[num_global_elements]);

            // Create a local map of the elements this processor owns
            for (unsigned element_id=0; element_id<num_global_elements; element_id++)
            {
                try
                {
                    unsigned element_index = rMesh.GetElement(element_id)->GetIndex();
                    TS_ASSERT_EQUALS(element_id, element_index);

                    elements_owned[element_index] = 1;

                    total_elements_this_process++;
                }
                catch(Exception&)
                {
                    elements_owned[element_id] = 0;
                }
            }

            TS_ASSERT_EQUALS(rMesh.GetNumLocalElements(), total_elements_this_process);

            // Combine all the local maps by adding them up in the master process
            boost::scoped_array<unsigned> elements_reduction(new unsigned[num_global_elements]);
            MPI_Reduce(elements_owned.get(), elements_reduction.get(), num_global_elements, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);

            // Make sure every element is owned at least by one processor
            if (PetscTools::AmMaster())
            {
                for (unsigned element_id=0; element_id<num_global_elements; element_id++)
                {
                    TS_ASSERT(elements_reduction[element_id] > 0);
                }
            }
        }

        /*
         * All boundary elements have been assigned
         */
        unsigned total_b_elements_this_process = 0;
        {
            const unsigned num_global_b_elements = rMesh.GetNumBoundaryElements();
            boost::scoped_array<unsigned> b_elements_owned(new unsigned[num_global_b_elements]);

            // Create a local map of the boundary elements this processor owns
            for (unsigned b_element_id=0; b_element_id<num_global_b_elements; b_element_id++)
            {
                try
                {
                    unsigned b_element_index = rMesh.GetBoundaryElement(b_element_id)->GetIndex();
                    TS_ASSERT_EQUALS(b_element_id, b_element_index);

                    b_elements_owned[b_element_index] = 1;

                    total_b_elements_this_process++;
                }
                catch(Exception&)
                {
                    b_elements_owned[b_element_id] = 0;
                }
            }

            TS_ASSERT_EQUALS(rMesh.GetNumLocalBoundaryElements(), total_b_elements_this_process);

            // Combine all the local maps by adding them up in the master process
            boost::scoped_array<unsigned> b_elements_reduction(new unsigned[num_global_b_elements]);
            MPI_Reduce(b_elements_owned.get(), b_elements_reduction.get(), num_global_b_elements, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);

            // Make sure every boundary element is owned at least by one processor
            if (PetscTools::AmMaster())
            {
                for (unsigned b_element_id=0; b_element_id<num_global_b_elements; b_element_id++)
                {
                    TS_ASSERT(b_elements_reduction[b_element_id] > 0);
                }
            }
        }
        if (total_nodes_this_process != 0)
        {
            TS_ASSERT( 0u != total_nodes_this_process );
            TS_ASSERT( 0u != total_elements_this_process );
            TS_ASSERT( 0u != total_b_elements_this_process );
        }
        else
        {
            // A partitioner may allocate no nodes to a partition if the mesh is small and there are many processes
            // Look out for "You just increased the maxndoms"
            TS_ASSERT(0u == total_nodes_this_process);
            TS_ASSERT(0u == total_elements_this_process);
            TS_ASSERT(0u == total_b_elements_this_process);
        }
    }

public:

    void TestConstructFromMeshReader1D()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_with_attributes");

        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 11u);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 11u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 10u);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            try
            {
                unsigned region = mesh.GetElement(i)->GetUnsignedAttribute();
                TS_ASSERT_EQUALS(region, i%5+1);
                TS_ASSERT_EQUALS(i, mesh.GetElement(i)->GetIndex());
            }
            catch(Exception&)
            {
                // I don't own this element do I?
            }
        }

        //Connectivity is normally 3 (if we own at least 2 local nodes)
        if (mesh.GetNumLocalNodes() >= 2)
        {
            TS_ASSERT_EQUALS(mesh.CalculateMaximumNodeConnectivityPerProcess(), 3U);
        }
    }

    void TestConstructFromMeshReader2DWithoutReordering()
    {
        /*
         * In this test we don't use reordering since we want to check that a TetrahedralMesh and
         * a DistributedTetrahedralMesh create the same geometry from the same file.
         */
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

        DistributedTetrahedralMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB); // No reordering
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543u);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543u);
        TS_ASSERT_EQUALS(mesh.GetDistributedVectorFactory()->GetProblemSize(), 543u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 100u);
        ///\todo There is no way to count the *global* number of boundary nodes
        //TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 100u);

        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            TS_ASSERT(iter->GetOwnership());
        }

        // Check the inverse Jacobian
        c_matrix<double, 2, 2> jacobian;
        double jacobian_determinant;
        c_matrix<double, 2, 2> inverse_jacobian;

        c_matrix<double, 2, 2> element_jacobian;
        double element_jacobian_determinant;
        c_matrix<double, 2, 2> element_inverse_jacobian;

        try
        {
            mesh.GetInverseJacobianForElement(0, jacobian, jacobian_determinant, inverse_jacobian);
            mesh.GetElement(0)->CalculateInverseJacobian(element_jacobian, element_jacobian_determinant, element_inverse_jacobian);

            TS_ASSERT_EQUALS(element_jacobian_determinant, jacobian_determinant);

            for (unsigned row=0; row<2; row++)
            {
                for (unsigned col=0; col<2; col++)
                {
                    TS_ASSERT_EQUALS(element_inverse_jacobian(row,col), inverse_jacobian(row,col));
                }
            }
        }
        catch(Exception&)
        {
            // I don't own this element do I?
        }

        c_vector<double, 2> direction;
        c_vector<double, 2> element_direction;

        try
        {
            mesh.GetWeightedDirectionForBoundaryElement(0, direction, jacobian_determinant);
            mesh.GetBoundaryElement(0)->CalculateWeightedDirection(element_direction, element_jacobian_determinant);

            TS_ASSERT_EQUALS(element_jacobian_determinant, jacobian_determinant);

            for (unsigned row=0; row<2; row++)
            {
                TS_ASSERT_EQUALS(element_direction(row), direction(row));
            }
        }
        catch(Exception&)
        {
            // I don't own this boundary element do I?
        }

        TetrahedralMesh<2,2> seq_mesh;
        seq_mesh.ConstructFromMeshReader(mesh_reader);

        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();

            Element<2,2>* p_sequ_element = seq_mesh.GetElement(element_index);
            TS_ASSERT_EQUALS(element_index, p_sequ_element->GetIndex());

            for (unsigned node_local_index=0; node_local_index < iter->GetNumNodes(); node_local_index++)
            {
                TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(node_local_index),
                                 p_sequ_element->GetNodeGlobalIndex(node_local_index));

                TS_ASSERT_EQUALS(iter->GetNode(node_local_index)->GetPoint()[0],
                                 p_sequ_element->GetNode(node_local_index)->GetPoint()[0]);
            }
        }

        for (DistributedTetrahedralMesh<2,2>::BoundaryElementIterator it=mesh.GetBoundaryElementIteratorBegin();
             it!=mesh.GetBoundaryElementIteratorEnd();
             ++it)
        {
            BoundaryElement<1,2>* p_para_boundary_element = *it;
            unsigned boundary_element_index = p_para_boundary_element->GetIndex();

            BoundaryElement<1,2>* p_sequ_boundary_element = seq_mesh.GetBoundaryElement(boundary_element_index);
            TS_ASSERT_EQUALS(boundary_element_index, p_sequ_boundary_element->GetIndex());

            for (unsigned node_local_index=0; node_local_index < p_para_boundary_element->GetNumNodes(); node_local_index++)
            {
                TS_ASSERT_EQUALS(p_para_boundary_element->GetNodeGlobalIndex(node_local_index),
                                 p_sequ_boundary_element->GetNodeGlobalIndex(node_local_index));

                TS_ASSERT_EQUALS(p_para_boundary_element->GetNode(node_local_index)->GetPoint()[0],
                                 p_sequ_boundary_element->GetNode(node_local_index)->GetPoint()[0]);
            }
        }
    }

    // See #1199
    void TestConstructFromMeshReader2DWithUnevenDistribution()
    {
        unsigned local_nodes = 1u;
        unsigned local_nodes_wrong = 1u;
        unsigned total_nodes = 543u;
        unsigned total_nodes_wrong = 100u;
        if (PetscTools::AmTopMost())
        {
            local_nodes = total_nodes - (PetscTools::GetNumProcs()-1);
            local_nodes_wrong = total_nodes_wrong - (PetscTools::GetNumProcs()-1);
        }

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        DistributedTetrahedralMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB); // No reordering

        // Exceptions
        DistributedVectorFactory* p_wrong_factory1 = new DistributedVectorFactory(PetscTools::GetMyRank(), PetscTools::GetMyRank()+1,
                                                                                  PetscTools::GetNumProcs(), PetscTools::GetNumProcs()+1);
        TS_ASSERT_THROWS_THIS(mesh.SetDistributedVectorFactory(p_wrong_factory1),
                              "The distributed vector factory provided to the mesh is for the wrong number of processes.");
        delete p_wrong_factory1;

        DistributedVectorFactory* p_wrong_factory2 = new DistributedVectorFactory(total_nodes_wrong, local_nodes_wrong);
        mesh.SetDistributedVectorFactory(p_wrong_factory2);
        TS_ASSERT_THROWS_THIS(mesh.ConstructFromMeshReader(mesh_reader),
                              "The distributed vector factory size in the mesh doesn't match the total number of nodes.");
        mesh.mpDistributedVectorFactory=NULL;
        delete p_wrong_factory2;

        // OK call
        DistributedVectorFactory* p_uneven_factory = new DistributedVectorFactory(total_nodes, local_nodes);
        mesh.SetDistributedVectorFactory(p_uneven_factory);
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check the mesh is using the supplied factory
        TS_ASSERT(mesh.GetDistributedVectorFactory() == p_uneven_factory);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), total_nodes);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), total_nodes);
        TS_ASSERT_EQUALS(mesh.GetDistributedVectorFactory()->GetProblemSize(), total_nodes);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 100u);
        TS_ASSERT_EQUALS(mesh.GetNumLocalNodes(), local_nodes);

        // Another exception
        TS_ASSERT_THROWS_THIS(mesh.SetDistributedVectorFactory(p_uneven_factory),
                              "Cannot change the mesh's distributed vector factory once it has been set.");
    }

    void TestConstructFromMeshReader3D()
    {
        /*
         * In this test we let parMETIS reorder the DistributedTetrahedralMesh. We want to check that although
         * the indices of the nodes have changed, the location of the nodes is consistent with a
         * TetrahedralMesh representation of the same mesh.
         */
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 51u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 136u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 96u);

        TetrahedralMesh<3,3> seq_mesh;
        seq_mesh.ConstructFromMeshReader(mesh_reader);

        for (AbstractTetrahedralMesh<3,3>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();

            Element<3,3>* p_sequ_element = seq_mesh.GetElement(element_index);

            // The elements have the same index and the nodes are located in the same position.
            TS_ASSERT_EQUALS(element_index, p_sequ_element->GetIndex());
            for (unsigned node_local_index=0; node_local_index < iter->GetNumNodes(); node_local_index++)
            {
                for (unsigned dim=0; dim<3; dim++)
                {
                    TS_ASSERT_EQUALS(iter->GetNode(node_local_index)->GetPoint()[dim],
                                     p_sequ_element->GetNode(node_local_index)->GetPoint()[dim]);
                }
            }
        }

        for (DistributedTetrahedralMesh<3,3>::BoundaryElementIterator it=mesh.GetBoundaryElementIteratorBegin();
             it!=mesh.GetBoundaryElementIteratorEnd();
             ++it)
        {
            BoundaryElement<2,3>* p_para_boundary_element = *it;
            unsigned boundary_element_index = p_para_boundary_element->GetIndex();

            BoundaryElement<2,3>* p_sequ_boundary_element = seq_mesh.GetBoundaryElement(boundary_element_index);

            // The boundary elements have the same index and the nodes are located in the same position.
            TS_ASSERT_EQUALS(boundary_element_index, p_sequ_boundary_element->GetIndex());
            for (unsigned node_local_index=0; node_local_index < p_para_boundary_element->GetNumNodes(); node_local_index++)
            {
                for (unsigned dim=0; dim<3; dim++)
                {
                    TS_ASSERT_EQUALS(p_para_boundary_element->GetNode(node_local_index)->GetPoint()[dim],
                                     p_sequ_boundary_element->GetNode(node_local_index)->GetPoint()[dim]);
                }
            }
        }

        //Scale it (for coverage)
        mesh.Scale(2.0);

    }

    void TestConstructionFromMeshReaderWithNodeAttributes()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements_with_node_attributes");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 12u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 20u);

        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 12u);

        // Check all nodes have 2 attributes
        for (unsigned node_index = 0; node_index < mesh.GetNumNodes(); node_index++)
        {
            if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(node_index) )
            {
                TS_ASSERT_EQUALS(mesh.GetNode(node_index)->rGetNodeAttributes().size(), 2u);
            }
        }

        // Now check attribute values at two probe nodes
        unsigned probe_node_1 = 0u;
        unsigned probe_node_2 = 8u;

        if (PetscTools::IsParallel())//need to figure out where they end up in permutation
        {
            TS_ASSERT(mesh.rGetNodePermutation().size() > 0);
            probe_node_1 = mesh.rGetNodePermutation()[probe_node_1];
            probe_node_2 = mesh.rGetNodePermutation()[probe_node_2];
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(probe_node_1) )
        {
            TS_ASSERT_DELTA(mesh.GetNode(probe_node_1)->rGetNodeAttributes()[0u], 25.2, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(probe_node_1)->rGetNodeAttributes()[1u], 16.3, 1e-6);
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(probe_node_2) )
        {
            TS_ASSERT_DELTA(mesh.GetNode(probe_node_2)->rGetNodeAttributes()[0u], 3.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(probe_node_2)->rGetNodeAttributes()[1u], 24.5, 1e-6);
        }
    }

    void TestConstructFromMeshReaderWithNclFile()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements_binary");
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<3,3> mesh_writer("WritingNclFile", "cube_136_elements_binary");
        mesh_writer.SetWriteFilesAsBinary();
        mesh_writer.WriteFilesUsingMesh(mesh);

        std::string output_dir = mesh_writer.GetOutputDirectory();

        TrianglesMeshReader<3,3> mesh_reader_ncl(output_dir + "cube_136_elements_binary");
        TS_ASSERT(mesh_reader_ncl.HasNclFile());
        DistributedTetrahedralMesh<3,3> mesh_from_ncl(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh_from_ncl.ConstructFromMeshReader(mesh_reader_ncl);

        CompareMeshes( mesh, mesh_from_ncl );
    }

    void TestRandomShuffle()
    {
        unsigned num_elts = 200;

        std::vector<unsigned> random_order(num_elts);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        p_gen->Reseed(0);
        p_gen->Shuffle(num_elts,random_order);


        if (PetscTools::IsParallel())
        {
            //Check all processes have the same random shuffle
            int num_procs = PetscTools::GetNumProcs();
            int my_rank = PetscTools::GetMyRank();
            int source_rank = (my_rank + num_procs - 1) % num_procs;
            int destination_rank = (my_rank + 1) % num_procs;
            int my_tag;
            int source_tag;

            MPI_Status status;

            for (unsigned element_number = 0; element_number < num_elts; element_number++)
            {
                unsigned my_entry = random_order[element_number];

                my_tag = my_rank + num_elts*element_number;
                source_tag = source_rank + num_elts*element_number;

                //This may not work sequentially on some versions of MPI (MPICH2)
                MPI_Send( &my_entry, 1, MPI_UNSIGNED, destination_rank, my_tag, PETSC_COMM_WORLD );

                unsigned neighbours_entry=0;
                MPI_Recv( &neighbours_entry, 1, MPI_UNSIGNED, source_rank, source_tag, PETSC_COMM_WORLD, &status );
                PetscTools::Barrier();

                TS_ASSERT_EQUALS( my_entry, neighbours_entry );
            }
        }
    }

    /*
     *  If you need to generate a binary mesh from an existing one. Use something like:
     *
     *      TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_1mm_6000_elements");
     *      TrianglesMeshWriter<3,3> mesh_writer("new_binary_mesh", "3D_0_to_1mm_6000_elements_binary");
     *      mesh_writer.SetWriteFilesAsBinary();
     *      mesh_writer.WriteFilesUsingMeshReader(mesh_reader);
     *
     */
    void TestComparePartitionQualities()
    {
        unsigned num_local_nodes_petsc_parmetis, num_local_nodes_parmetis, num_local_nodes_metis_deprecated;

        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_1mm_6000_elements");
            //TrianglesMeshReader<3,3> mesh_reader("heart/test/data/heart");
            DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
            mesh.ConstructFromMeshReader(mesh_reader);
            // There's warning because METIS is deprecated and this is actually a parMETIS partition
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
            Warnings::Instance()->QuietDestroy();

            num_local_nodes_metis_deprecated = mesh.GetNumLocalNodes();
        }

        if (PetscTools::HasParMetis())
        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_1mm_6000_elements_binary");
            //TrianglesMeshReader<3,3> mesh_reader("heart/test/data/heart_binary");
            DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::PETSC_MAT_PARTITION);
            mesh.ConstructFromMeshReader(mesh_reader);

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
            TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());

            CheckEverythingIsAssigned<3,3>(mesh);

            num_local_nodes_petsc_parmetis = mesh.GetNumLocalNodes();
        }
        else
        {
            num_local_nodes_petsc_parmetis = 0u;
        }

        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_1mm_6000_elements_binary");
            //TrianglesMeshReader<3,3> mesh_reader("heart/test/data/heart_binary");
            DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);
            mesh.ConstructFromMeshReader(mesh_reader);

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
            TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());

            CheckEverythingIsAssigned<3,3>(mesh);

            num_local_nodes_parmetis = mesh.GetNumLocalNodes();
            TS_ASSERT_EQUALS(num_local_nodes_metis_deprecated, num_local_nodes_parmetis); // METIS is deprecated so these are the same
        }

        unsigned max_local_nodes_petsc_parmetis;
        unsigned max_local_nodes_parmetis;

        MPI_Allreduce (&num_local_nodes_petsc_parmetis, &max_local_nodes_petsc_parmetis, 1, MPI_UNSIGNED, MPI_MAX, PETSC_COMM_WORLD );
        MPI_Allreduce (&num_local_nodes_parmetis, &max_local_nodes_parmetis, 1, MPI_UNSIGNED, MPI_MAX, PETSC_COMM_WORLD );

        if (PetscTools::AmMaster())
        {
            std::cout << "PETSC PARMETIS\tPARMETIS" << std::endl;
            std::cout << max_local_nodes_petsc_parmetis << "\t\t" << max_local_nodes_parmetis << std::endl;
        }
        PetscTools::Barrier();

        TS_ASSERT(num_local_nodes_petsc_parmetis <= max_local_nodes_parmetis);
        //Watch out for dumb partition and warn about it
        if (PetscTools::IsParallel())
        {
            if (max_local_nodes_petsc_parmetis ==  0u)
            {
                TS_TRACE("Did not use ParMETIS because the interface is not available");
            }
        }
    }

    void TestEverythingIsAssignedParMetisLibraryAsciiFiles()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());

        CheckEverythingIsAssigned<3,3>(mesh);
    }

    void TestEverythingIsAssignedParMetisLibraryBinaryFiles()
    {
        ///\todo This test fails in CheckEverythingIsAssigned with 10 or more processes

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements_binary");
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());

        CheckEverythingIsAssigned<3,3>(mesh);
    }

    void TestEverythingIsAssignedPetscPartition()
    {
        if (PetscTools::HasParMetis())
        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
            DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::PETSC_MAT_PARTITION);
            mesh.ConstructFromMeshReader(mesh_reader);

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
            TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());

            CheckEverythingIsAssigned<3,3>(mesh);
        }
    }

    void TestEverythingIsAssignedPetscPartitionBinaryFiles()
    {
        if (PetscTools::HasParMetis())
        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements_binary");
            DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::PETSC_MAT_PARTITION);
            mesh.ConstructFromMeshReader(mesh_reader);

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
            TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());

            CheckEverythingIsAssigned<3,3>(mesh);
        }
    }


    void TestConstruct3DWithRegions()
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart_nonnegative_flags");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            try
            {
                unsigned region = mesh.GetElement(i)->GetUnsignedAttribute();
                TS_ASSERT_EQUALS(region, (i+1)%3+1);
            }
            catch(Exception&)
            {
                // I don't own this element do I?
            }
        }

        TS_ASSERT_EQUALS(mesh_reader.GetNumFaceAttributes(), 1u);

        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            try
            {
                unsigned region = mesh.GetBoundaryElement(i)->GetUnsignedAttribute();
                TS_ASSERT_LESS_THAN(region, 4u);
            }
            catch(Exception&)
            {
                // I don't own this element do I?
            }
        }
    }


    void TestPartitioningOfEmbeddedDimensionMesh()
    {
        //Shouldn't ever use a partition other than DUMB because it's 1-D
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/branched_1d_in_3d_mesh");
        DistributedTetrahedralMesh<1,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 31u);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 30u);
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 3u);
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 3u);
    }

    /**
     * This test constructs a simple cuboid mesh and divides
     * between two processes
     */
    void TestGeometricPartition()
    {
        unsigned num_procs = PetscTools::GetNumProcs();
        unsigned rank = PetscTools::GetMyRank();

        TetrahedralMesh<3,3> test_mesh;
        test_mesh.ConstructCuboid(num_procs-1,num_procs-1,num_procs-1);

        TrianglesMeshWriter<3,3> mesh_writer("TestGeometricPartition", "TestMesh", false);
        mesh_writer.WriteFilesUsingMesh(test_mesh);


        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<3,3> mesh_reader(output_dir+"TestMesh");

        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::GEOMETRIC);

        TS_ASSERT_THROWS_THIS(mesh.GetProcessRegion(), "Trying to get unset mpSpaceRegion");
        if (PetscTools::IsParallel())   // Won't throw in serial.
        {
            TS_ASSERT_THROWS_THIS(mesh.ConstructFromMeshReader(mesh_reader), "Using GEOMETRIC partition for DistributedTetrahedralMesh with local regions not set. Call SetProcessRegion(ChasteCuboid)");
        }

        // Test throwing an exception if a node doesn't lie in the union of the regions, or regions are not disjoint.
        if (PetscTools::IsParallel())
        {
            ChastePoint<3> lower_wrong(-0.5, -0.5, ((double)rank+0.5));
            ChastePoint<3> upper_wrong((double)(num_procs+1)+0.5, (double)(num_procs+1) + 0.5, ((double)(rank+1)+0.5));

            ChasteCuboid<3> cuboid_wrong(lower_wrong, upper_wrong);

            mesh.SetProcessRegion(&cuboid_wrong);
            ChasteCuboid<3>* test_cuboid = mesh.GetProcessRegion();
            TS_ASSERT_EQUALS(test_cuboid, &cuboid_wrong);

            TS_ASSERT_THROWS_THIS(mesh.ConstructFromMeshReader(mesh_reader), "A node is either not in geometric region, or the regions are not disjoint.");
        }

        ChastePoint<3> lower(-0.5, -0.5, ((double)rank-0.5));
        ChastePoint<3> upper((double)(num_procs+1)+0.5, (double)(num_procs+1) + 0.5, ((double)(rank)+0.5));

        ChasteCuboid<3> cuboid(lower, upper);

        mesh.SetProcessRegion(&cuboid);

        // Make a new mesh reader
        TrianglesMeshReader<3,3> mesh_reader2(output_dir+"TestMesh");

        mesh.ConstructFromMeshReader(mesh_reader2);

        // Check construction is correct.
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_procs*num_procs*num_procs);
        TS_ASSERT_EQUALS(mesh.GetNumLocalNodes(), num_procs*num_procs);

        for (AbstractMesh<3,3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
               node_iter != mesh.GetNodeIteratorEnd();
               ++node_iter)
        {
           TS_ASSERT_LESS_THAN(node_iter->rGetLocation()[2], (double)rank+0.5);
           TS_ASSERT_LESS_THAN((double)rank-0.5, node_iter->rGetLocation()[2]);
        }
    }


    void TestArchiving()
    {
        FileFinder main_archive_dir("distributed_tetrahedral_mesh_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "distributed_tetrahedral_mesh.arch";
        ArchiveLocationInfo::SetMeshFilename("distributed_tetrahedral_mesh");

        DistributedTetrahedralMesh<2,2>* p_mesh = new DistributedTetrahedralMesh<2,2>(DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);
        //std::vector<unsigned> halo_node_indices;
        std::vector<Node<2>*> halo_nodes;
        unsigned num_nodes;
        unsigned local_num_nodes;
        unsigned num_elements;
        // archive
        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

            p_mesh->ConstructFromMeshReader(mesh_reader);
            num_nodes = p_mesh->GetNumNodes();
            local_num_nodes = p_mesh->GetNumLocalNodes();
            num_elements = p_mesh->GetNumElements();

            halo_nodes = p_mesh->mHaloNodes;

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(main_archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractTetrahedralMesh<2,2>* const p_mesh_abstract = static_cast<AbstractTetrahedralMesh<2,2>* >(p_mesh);
            (*p_arch) << p_mesh_abstract;
        }

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
            DistributedTetrahedralMesh<2,2>* p_mesh2 = static_cast<DistributedTetrahedralMesh<2,2>*>(p_mesh_abstract2);

            TS_ASSERT_EQUALS(p_mesh2->GetNumNodes(), num_nodes);
            TS_ASSERT_EQUALS(p_mesh2->GetNumLocalNodes(), local_num_nodes);
            TS_ASSERT_EQUALS(p_mesh2->GetNumElements(), num_elements);

            // Check some node co-ordinates
            try
            {
                Node<2>* p_node1 = p_mesh->GetNode(0);
                Node<2>* p_node2 = p_mesh2->GetNode(0);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], p_node2->GetPoint()[0], 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[1], p_node2->GetPoint()[1], 1e-6);
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            try
            {
                Node<2>* p_node1 = p_mesh->GetNode(500);
                Node<2>* p_node2 = p_mesh2->GetNode(500);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], p_node2->GetPoint()[0], 1e-6);
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            // Check first element has the right nodes
            try
            {
                Element<2,2>* p_element = p_mesh->GetElement(0);
                Element<2,2>* p_element2 = p_mesh2->GetElement(0);
                TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), p_element2->GetNodeGlobalIndex(0));
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            try
            {
                Element<2,2>* p_element = p_mesh->GetElement(500);
                Element<2,2>* p_element2 = p_mesh2->GetElement(500);
                TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), p_element2->GetNodeGlobalIndex(0));
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            // Check the halo nodes are right
            std::vector<Node<2>*> halo_nodes2 = p_mesh2->mHaloNodes;
            TS_ASSERT_EQUALS(halo_nodes2.size(), halo_nodes.size());
            delete p_mesh2;
        }

        // restore from a single processor archive
        {
            FileFinder archive_dir("mesh/test/data/distributed_mesh_archive", RelativeTo::ChasteSourceRoot);
            if (PetscTools::IsSequential())
            {
                ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(
                        archive_dir, "distributed_tetrahedral_mesh.arch");
                boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
                AbstractTetrahedralMesh<2,2>* p_mesh3 = NULL;
                (*p_arch) >> p_mesh3;
                delete p_mesh3;
            }
            else
            {
                typedef ArchiveOpener<boost::archive::text_iarchive, std::ifstream> InputArchiveOpener;
                if (PetscTools::GetMyRank() > 0)
                {
                    // Should not read this archive because none exists here.
                    TS_ASSERT_THROWS_CONTAINS(InputArchiveOpener arch_opener(archive_dir, "distributed_tetrahedral_mesh.arch"),
                                "Cannot load secondary archive file:");
                }
                else
                {
                    // Should not read this archive because there are two or more processes and
                    // this archive was written on one process.
                    InputArchiveOpener arch_opener(archive_dir, "distributed_tetrahedral_mesh.arch");
                    boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
                    AbstractTetrahedralMesh<2,2>* p_mesh3 = NULL;
                    TS_ASSERT_THROWS_THIS((*p_arch) >> p_mesh3,
                                          "This archive was written for a different number of processors");

                }
            }
        }

        delete p_mesh;
    }

    void TestArchivingBinaryMesh()
    {
        FileFinder archive_dir("distributed_tetrahedral_mesh_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "binary_mesh.arch";
        ArchiveLocationInfo::SetMeshFilename("binary_mesh");

        DistributedTetrahedralMesh<3,3>* p_mesh = new DistributedTetrahedralMesh<3,3>(DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);
        //std::vector<unsigned> halo_node_indices;
        std::vector<Node<3>*> halo_nodes;
        unsigned num_nodes;
        unsigned local_num_nodes;
        unsigned num_elements;
        // archive
        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements_binary");
            TS_ASSERT(mesh_reader.HasNclFile());

            p_mesh->ConstructFromMeshReader(mesh_reader);
            num_nodes = p_mesh->GetNumNodes();
            local_num_nodes = p_mesh->GetNumLocalNodes();
            num_elements = p_mesh->GetNumElements();

            halo_nodes = p_mesh->mHaloNodes;

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractTetrahedralMesh<3,3>* const p_mesh_abstract = static_cast<AbstractTetrahedralMesh<3,3>* >(p_mesh);
            (*p_arch) << p_mesh_abstract;
        }

        // restore
        {
            // Should archive the most abstract class you can to check boost knows what individual classes are.
            // (but here AbstractMesh doesn't have the methods below).
            AbstractTetrahedralMesh<3,3>* p_mesh_abstract2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // restore from the archive
            (*p_arch) >> p_mesh_abstract2;
            // Check we have the right number of nodes & elements
            DistributedTetrahedralMesh<3,3>* p_mesh2 = static_cast<DistributedTetrahedralMesh<3,3>*>(p_mesh_abstract2);

            TS_ASSERT_EQUALS(p_mesh2->GetNumNodes(), num_nodes);
            TS_ASSERT_EQUALS(p_mesh2->GetNumLocalNodes(), local_num_nodes);
            TS_ASSERT_EQUALS(p_mesh2->GetNumElements(), num_elements);

            // Check some node co-ordinates
            try
            {
                Node<3>* p_node1 = p_mesh->GetNode(0);
                Node<3>* p_node2 = p_mesh2->GetNode(0);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], p_node2->GetPoint()[0], 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[1], p_node2->GetPoint()[1], 1e-6);
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            try
            {
                Node<3>* p_node1 = p_mesh->GetNode(500);
                Node<3>* p_node2 = p_mesh2->GetNode(500);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], p_node2->GetPoint()[0], 1e-6);
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            // Check first element has the right nodes
            try
            {
                Element<3,3>* p_element = p_mesh->GetElement(0);
                Element<3,3>* p_element2 = p_mesh2->GetElement(0);
                TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), p_element2->GetNodeGlobalIndex(0));
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            try
            {
                Element<3,3>* p_element = p_mesh->GetElement(500);
                Element<3,3>* p_element2 = p_mesh2->GetElement(500);
                TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), p_element2->GetNodeGlobalIndex(0));
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            // Check the halo nodes are right
            std::vector<Node<3>*> halo_nodes2 = p_mesh2->mHaloNodes;
            TS_ASSERT_EQUALS(halo_nodes2.size(), halo_nodes.size());
            delete p_mesh2;
        }

        delete p_mesh;
    }

private:

    template <unsigned DIM>
    void CompareParallelMeshOwnership(DistributedTetrahedralMesh<DIM,DIM> &readMesh, DistributedTetrahedralMesh<DIM,DIM> &constructedMesh)
    {
        // The read mesh has a dumb partition in the test
        TS_ASSERT_EQUALS(readMesh.GetPartitionType(), DistributedTetrahedralMeshPartitionType::DUMB);
        // All constructed meshes have dumb partitioning -- so that they are invariant under archiving
        TS_ASSERT_EQUALS(constructedMesh.GetPartitionType(), DistributedTetrahedralMeshPartitionType::DUMB);
        TS_ASSERT_EQUALS(constructedMesh.GetDistributedVectorFactory()->GetLocalOwnership(),
                         readMesh.GetDistributedVectorFactory()->GetLocalOwnership());
        TS_ASSERT_EQUALS(constructedMesh.GetNumNodes(), readMesh.GetNumNodes());
        TS_ASSERT_EQUALS(constructedMesh.GetNumLocalNodes(), readMesh.GetNumLocalNodes());
        TS_ASSERT_EQUALS(constructedMesh.GetNumBoundaryNodes(), readMesh.GetNumBoundaryNodes());
        TS_ASSERT_EQUALS(constructedMesh.GetNumBoundaryElements(),  readMesh.GetNumBoundaryElements());
        TS_ASSERT_EQUALS(constructedMesh.GetNumElements(), readMesh.GetNumElements());
        TS_ASSERT_EQUALS(constructedMesh.GetNumLocalElements(), readMesh.GetNumLocalElements());

        for (unsigned i=0; i<readMesh.GetNumNodes(); i++)
        {
            try
            {
                unsigned index=constructedMesh.SolveNodeMapping(i);
                // Read mesh didn't throw so owns the node
                TS_ASSERT_THROWS_NOTHING(constructedMesh.GetNode(i));
                TS_ASSERT_EQUALS(index, readMesh.SolveNodeMapping(i));
             }
            catch(Exception&)
            {
                // Read mesh threw so does not own node
                TS_ASSERT_THROWS_CONTAINS(constructedMesh.GetNode(i), "does not belong to processor");
            }
        }

        for (unsigned i=0; i<readMesh.GetNumElements(); i++)
        {
            try
            {
                unsigned index=constructedMesh.SolveElementMapping(i);
                // Read mesh didn't throw so owns the element
                TS_ASSERT_THROWS_NOTHING(constructedMesh.GetElement(i));
                TS_ASSERT_EQUALS(index, readMesh.SolveElementMapping(i));
             }
            catch(Exception&)
            {
                // Read mesh threw so does not own element
                TS_ASSERT_THROWS_CONTAINS(constructedMesh.GetElement(i), "does not belong to processor");
            }
        }

        for (unsigned i=0; i<readMesh.GetNumBoundaryElements(); i++)
        {
            try
            {
                unsigned index = constructedMesh.SolveBoundaryElementMapping(i);
                // Read mesh didn't throw so owns the element
                TS_ASSERT_THROWS_NOTHING(constructedMesh.GetBoundaryElement(i));
                TS_ASSERT_EQUALS(index, readMesh.SolveBoundaryElementMapping(i));
             }
            catch(Exception&)
            {
                // Read mesh threw so does not own element
                TS_ASSERT_THROWS_CONTAINS(constructedMesh.GetBoundaryElement(i), "does not belong to processor");
            }
        }
    }

public:
    void TestConstructLinearMesh()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_with_attributes");
        DistributedTetrahedralMesh<1,1> read_mesh;
        read_mesh.ConstructFromMeshReader(mesh_reader);
        DistributedTetrahedralMesh<1,1> constructed_mesh;
        constructed_mesh.ConstructLinearMesh(10u);

        CompareParallelMeshOwnership(read_mesh, constructed_mesh);

        unsigned owned = constructed_mesh.GetDistributedVectorFactory()->GetLocalOwnership();
        unsigned owned_in_read = read_mesh.GetDistributedVectorFactory()->GetLocalOwnership();
        TS_ASSERT_EQUALS(owned_in_read, owned);
        TS_ASSERT_EQUALS(constructed_mesh.GetNumBoundaryElements(), 2u);
        TS_ASSERT_EQUALS(constructed_mesh.GetNumLocalNodes(), owned);
        // Sequential: Process owns one fewer element than the number of nodes
        // Parallel: End processes own the same as the number of node (since one node is paired with a halo node)
        // Parallel: Middle processes own one more than the number of nodes (since two nodes are paired with a halo nodes)
        unsigned expected_elements = owned+1;
        if (PetscTools::AmMaster())
        {
            expected_elements--;
        }
         if (PetscTools::AmTopMost())
        {
            expected_elements--;
        }
        TS_ASSERT_EQUALS(constructed_mesh.GetNumLocalElements(), expected_elements);

        // Note that boundary nodes are local to the process
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(constructed_mesh.GetNumBoundaryNodes(), 2u);
        }
        else
        {
            TS_ASSERT_LESS_THAN(constructed_mesh.GetNumBoundaryNodes(), 2u);
        }
    }

    void TestConstructLinearMeshVerySmall()
    {
        DistributedTetrahedralMesh<1,1> small_mesh;
        // Coverage hack
        TS_ASSERT_THROWS_THIS(small_mesh.ConstructLinearMesh(0), "There aren't enough nodes to make parallelisation worthwhile");

        // Works with up to 3 processes
        small_mesh.ConstructRegularSlabMesh(10.0, 20.0);
        unsigned owned=small_mesh.GetDistributedVectorFactory()->GetLocalOwnership();
        TS_ASSERT_EQUALS(small_mesh.GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(small_mesh.GetNumLocalNodes(), owned);
        TS_ASSERT_EQUALS(small_mesh.GetNumBoundaryElements(),  2u);
        TS_ASSERT_EQUALS(small_mesh.GetNumElements(), 2u);
        // See logic in earlier test
        unsigned expected_elements=owned+1;

        std::vector<unsigned> halo_indices;
        small_mesh.GetHaloNodeIndices(halo_indices);

        // Check the size
        TS_ASSERT_EQUALS(small_mesh.GetNumHaloNodes(), halo_indices.size());

        // Check no halos in a sequential simulation
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(small_mesh.GetNumHaloNodes(), 0u);
        }

        // Check that iteration does the same thing
        unsigned i = 0;
        for (DistributedTetrahedralMesh<1,1>::HaloNodeIterator it=small_mesh.GetHaloNodeIteratorBegin();
             it != small_mesh.GetHaloNodeIteratorEnd();
             ++it,i++)
        {
            TS_ASSERT_EQUALS(halo_indices[i], (*it)->GetIndex());
        }

        /**
         * 1 Proc:
         * p0:  0 Ow 1 Ow 2 Ow
         * 2 Proc:
         * p0:  0 Ow 1 Ow 2 Ha
         * p1:       1 Ha 2 Ow
         * 3 Proc:
         * p0:  0 Ow 1 Ha
         * p1:  0 Ha 1 Ow 2 Ha
         * p2:       2 Ha 3 Ow
         */
        if (PetscTools::AmMaster())
        {
            expected_elements--;
            //Left processor always owns left node
            TS_ASSERT_EQUALS(small_mesh.GetNode(0),small_mesh.GetNodeOrHaloNode(0));
            TS_ASSERT_DELTA(small_mesh.GetNodeOrHaloNode(0)->rGetLocation()[0], 0.0, 1e-5);
            TS_ASSERT_DELTA(small_mesh.GetNodeOrHaloNode(1)->rGetLocation()[0], 10.0, 1e-5);
            if (PetscTools::IsSequential())
            {
                TS_ASSERT_EQUALS(halo_indices.size(), 0u);
            }
            else
            {
                TS_ASSERT_EQUALS(halo_indices.size(), 1u);
                TS_ASSERT_DELTA(halo_indices[0], 1u, 1u); // Halo is at index 1 (3 procs) or index 2 (2 procs)
            }
        }
        if (PetscTools::AmTopMost() && PetscTools::GetNumProcs() <= 3)
        {
            expected_elements--;
            if (PetscTools::IsSequential())
            {
                TS_ASSERT_EQUALS(halo_indices.size(), 0u);
            }
            else
            {
                TS_ASSERT_EQUALS(halo_indices.size(), 1u);
                TS_ASSERT_EQUALS(halo_indices[0], 1u); //Halo is at index 1 (2 or 3 procs)
                TS_ASSERT_THROWS_CONTAINS(small_mesh.GetNodeOrHaloNode(0), "Requested node/halo");
                // Right processor has  node 1 as halo
                TS_ASSERT_THROWS_CONTAINS(small_mesh.GetNode(1), "does not belong to processor");
                TS_ASSERT_THROWS_NOTHING(small_mesh.GetNodeOrHaloNode(1));//It's a halo
                TS_ASSERT_DELTA(small_mesh.GetNodeOrHaloNode(1)->rGetLocation()[0], 10.0, 1e-5);
                TS_ASSERT_DELTA(small_mesh.GetNodeOrHaloNode(2)->rGetLocation()[0], 20.0, 1e-5);
            }
        }
        if (PetscTools::GetNumProcs() > 3 && PetscTools::GetMyRank()==2)
        {
            // This is the equivalent to "top most"
            expected_elements--;
        }
        if (PetscTools::GetMyRank() < 3)
        {
            TS_ASSERT_EQUALS(small_mesh.GetNumLocalElements(), expected_elements);
        }
        else
        {
            // This process owns nothing
            TS_ASSERT_EQUALS(small_mesh.GetNumLocalNodes(), 0u);
            TS_ASSERT_EQUALS(owned, 0u);
            TS_ASSERT_EQUALS(small_mesh.GetNumLocalElements(), 0u);
            TS_ASSERT_EQUALS(halo_indices.size(), 0u);
        }
    }

    void TestConstructLinearMeshSmall()
    {
        unsigned width = 2;
        // Works well with exactly 3 processors
        if (PetscTools::GetNumProcs() != width + 1)
        {
            TS_TRACE("This test works with exactly 3 processes.");
            return;
        }
        TetrahedralMesh<1,1> base_mesh;
        base_mesh.ConstructLinearMesh(width);
        TrianglesMeshWriter<1,1> mesh_writer("", "linear");
        mesh_writer.WriteFilesUsingMesh(base_mesh);
        PetscTools::Barrier();
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<1,1> mesh_reader(output_dir+"linear");
        DistributedTetrahedralMesh<1,1> read_mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        read_mesh.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<1,1> constructed_mesh;
        constructed_mesh.ConstructLinearMesh(width);

        // Double check
        TS_ASSERT_EQUALS(constructed_mesh.GetNumBoundaryNodes(), read_mesh.GetNumBoundaryNodes());
        CompareParallelMeshOwnership(read_mesh, constructed_mesh);
    }

    void TestConstructRetangularMeshSmall()
    {
        unsigned width = 1;
        unsigned height = 2;
        // Works well with exactly 3 processors
        if (PetscTools::GetNumProcs() != height + 1)
        {
            TS_TRACE("This test works with exactly 3 processes.");
            return;
        }
        TetrahedralMesh<2,2> base_mesh;
        base_mesh.ConstructRectangularMesh(width, height);
        TrianglesMeshWriter<2,2> mesh_writer("", "rectangle");
        mesh_writer.WriteFilesUsingMesh(base_mesh);
        PetscTools::Barrier();
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,2> mesh_reader(output_dir+"rectangle");
        DistributedTetrahedralMesh<2,2> read_mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        read_mesh.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<2,2> constructed_mesh;
        constructed_mesh.ConstructRectangularMesh(width, height, false);

        TS_ASSERT_EQUALS(constructed_mesh.GetNumBoundaryNodes(), read_mesh.GetNumBoundaryNodes());
        CompareParallelMeshOwnership(read_mesh, constructed_mesh);
    }

    void TestConstructRetangularMesh()
    {
        unsigned width = 5;
        unsigned height = 4*PetscTools::GetNumProcs()-1; // 4*NumProcs layers of nodes (ensure dumb partition works in slices)

        TetrahedralMesh<2,2> base_mesh;
        base_mesh.ConstructRectangularMesh(width, height);
        TrianglesMeshWriter<2,2> mesh_writer("", "rectangle");
        mesh_writer.WriteFilesUsingMesh(base_mesh);
        PetscTools::Barrier();
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,2> mesh_reader(output_dir+"rectangle");
        DistributedTetrahedralMesh<2,2> read_mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        read_mesh.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<2,2> constructed_mesh;
        // Coverage
        TS_ASSERT_THROWS_THIS(constructed_mesh.ConstructRectangularMesh(width, 0, false),
                            "There aren't enough nodes to make parallelisation worthwhile");

        // Real mesh construction
        constructed_mesh.ConstructRectangularMesh(width, height, false);

        CompareParallelMeshOwnership(read_mesh, constructed_mesh);

        if (PetscTools::AmTopMost())
        {
            // Verify some element indices -- top left diagonal goes NW-SE (normal)
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[0],            2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[1], (height-1)+2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[0],            1.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[1], (height-1)+1.0/3.0, 1e-5);
        }
        if (PetscTools::AmMaster())
        {
            // Verify some element indices -- bottom left diagonal goes NW-SE (normal)
            TS_ASSERT_DELTA(constructed_mesh.GetElement(0)->CalculateCentroid()[0], 2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(0)->CalculateCentroid()[1], 2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(1)->CalculateCentroid()[0], 1.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(1)->CalculateCentroid()[1], 1.0/3.0, 1e-5);
        }
    }

    void TestConstructRetangularMeshStagger()
    {
        unsigned width = 4;
        unsigned height = 4*PetscTools::GetNumProcs()-1; //4*NumProcs layers of nodes (ensure dumb partition works in slices)

        TetrahedralMesh<2,2> base_mesh;
        base_mesh.ConstructRectangularMesh(width, height, true);
        TrianglesMeshWriter<2,2> mesh_writer("", "rectangle");
        mesh_writer.WriteFilesUsingMesh(base_mesh);
        PetscTools::Barrier();
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,2> mesh_reader(output_dir+"rectangle");
        DistributedTetrahedralMesh<2,2> read_mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        read_mesh.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<2,2> constructed_mesh;
        constructed_mesh.ConstructRectangularMesh(width, height, true);

        CompareParallelMeshOwnership(read_mesh, constructed_mesh);

        if (PetscTools::AmTopMost())
        {
            //Verify some element indices -- top left diagonal goes NW-SE (normal)
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[0],              2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[1],   (height-1)+2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[0],            1.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[1], (height-1)+1.0/3.0, 1e-5);
        }
        if (PetscTools::AmMaster())
        {
            TS_ASSERT_EQUALS(height%2, 1u);//If height is odd the bottom left is not staggered - next one is
            //Verify some element indices -- bottom left diagonal goes SW-NE (stagger)
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2)->CalculateCentroid()[0], 1 + 1.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2)->CalculateCentroid()[1], 2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(3)->CalculateCentroid()[0], 1 + 2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(3)->CalculateCentroid()[1], 1.0/3.0, 1e-5);
        }
    }

    void TestConstructCuboidMesh()
    {
        unsigned width = 2;
        unsigned height = 3;
        unsigned depth = 4*PetscTools::GetNumProcs()-1;

        TetrahedralMesh<3,3> base_mesh;
        base_mesh.ConstructCuboid(width, height, depth);
        TrianglesMeshWriter<3,3> mesh_writer("", "cuboid");
        mesh_writer.WriteFilesUsingMesh(base_mesh);
        PetscTools::Barrier();
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<3,3> mesh_reader(output_dir+"cuboid");
        DistributedTetrahedralMesh<3,3> read_mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        read_mesh.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<3,3> constructed_mesh;
        // Coverage
        TS_ASSERT_THROWS_THIS(constructed_mesh.ConstructCuboid(width, height, 0),
                            "There aren't enough nodes to make parallelisation worthwhile");
        constructed_mesh.ConstructCuboid(width, height, depth);

        CompareParallelMeshOwnership(read_mesh, constructed_mesh);

        // Test the bounding box methods
        ChasteCuboid<3> base_bounding_box=base_mesh.CalculateBoundingBox();
        TS_ASSERT_EQUALS(base_bounding_box.GetWidth(0), (double) width);
        TS_ASSERT_EQUALS(base_bounding_box.GetWidth(1), (double) height);
        TS_ASSERT_EQUALS(base_bounding_box.GetWidth(2), (double) depth);
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(base_bounding_box.GetLongestAxis(), 1U); // Tie between 1 and 2
        }
        else
        {
            TS_ASSERT_EQUALS(base_bounding_box.GetLongestAxis(), 2U); // 2 wins outright
        }
        ChasteCuboid<3> constructed_bounding_box=constructed_mesh.CalculateBoundingBox();
        TS_ASSERT_EQUALS(constructed_bounding_box.GetWidth(0), (double) width);
        TS_ASSERT_EQUALS(constructed_bounding_box.GetWidth(1), (double) height);
        TS_ASSERT_EQUALS(constructed_bounding_box.GetWidth(2), (double) depth);
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(constructed_bounding_box.GetLongestAxis(), 1U); // Tie between 1 and 2
        }
        else
        {
            TS_ASSERT_EQUALS(constructed_bounding_box.GetLongestAxis(), 2U); // 2 wins outright
        }
        TS_ASSERT_EQUALS(constructed_mesh.CalculateMaximumContainingElementsPerProcess(), 24U);  // Four surrounding cubes may have all 6 tetrahedra meeting at a node
        TS_ASSERT_EQUALS(constructed_mesh.CalculateMaximumNodeConnectivityPerProcess(), 15U);  // Four surrounding cubes may have all 6 tetrahedra meeting at a node


        // Test the GetNearestNodeIndex() method
        // This mesh isn't permuted (dumb partition), so we can check against the non-distributed case
        ChastePoint<3> outside_mesh(width*2, depth*2, height*2);  // outside mesh
        ChastePoint<3> within_mesh(width/3, depth/2, height/5);         // somewhere inside the mesh (not necessarily at a node)
        ChastePoint<3> origin_mesh(0.0, 0.0, 0.0);                      // origin
        TS_ASSERT_EQUALS(constructed_mesh.GetNearestNodeIndex(outside_mesh), base_mesh.GetNearestNodeIndex(outside_mesh));
        TS_ASSERT_EQUALS(constructed_mesh.GetNearestNodeIndex(within_mesh), base_mesh.GetNearestNodeIndex(within_mesh));
        TS_ASSERT_EQUALS(constructed_mesh.GetNearestNodeIndex(origin_mesh), base_mesh.GetNearestNodeIndex(origin_mesh));
    }

    void TestNearestNodeIndex2D()
    {

        // Similar to TestPointinMesh2D from TestTetrahedralMesh

        ChastePoint<2> point1(0.051, 0.051);    // close to node 60 in non-dist mesh
        ChastePoint<2> point2(0.2, 0.2);        // outside mesh, closest node is 120 (outermost node)
        ChastePoint<2> point3(0.05, 0.05);      // node 60 of non-dist mesh

        // Read in a mesh and create distributed and non-distributed versions
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        TetrahedralMesh<2,2> non_distributed_mesh;
        DistributedTetrahedralMesh<2,2> distributed_mesh;
        non_distributed_mesh.ConstructFromMeshReader(mesh_reader);
        distributed_mesh.ConstructFromMeshReader(mesh_reader);

        // Check that the closest nodes in the non-dist mesh are the same as those tested in TestPointinMesh2D
        TS_ASSERT_EQUALS(non_distributed_mesh.GetNearestNodeIndex(point1), 60u);
        TS_ASSERT_EQUALS(non_distributed_mesh.GetNearestNodeIndex(point2), 120u);   // Closest node is top right node
        TS_ASSERT_EQUALS(non_distributed_mesh.GetNearestNodeIndex(point3), 60u);


        // If we run on one process, the distributed node ordering is the same as with the
        // no-distributed mesh - we use the same checks
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(distributed_mesh.GetNearestNodeIndex(point1), 60u);
            TS_ASSERT_EQUALS(distributed_mesh.GetNearestNodeIndex(point2), 120u);   // Closest node is top right node
            TS_ASSERT_EQUALS(distributed_mesh.GetNearestNodeIndex(point3), 60u);

        }
        else
        {
            // The order of the nodes is permuted
            // The closest nodes will be permutation[closest node from sequential case]
            std::vector<unsigned> permutation = distributed_mesh.rGetNodePermutation();

            // We'll give some explanatory output - permuted indices of nodes 60 and 120
            if (distributed_mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(permutation[60]))
            {
                double permuted_node_60_x = distributed_mesh.GetNode(permutation[60])->rGetLocation()[0];
                double permuted_node_60_y = distributed_mesh.GetNode(permutation[60])->rGetLocation()[1];
                std::cout << "Original node = 60, permuted node = " << permutation[60]
                          << ", location = " << permuted_node_60_x << ", " << permuted_node_60_y << "\n";
            }
            if (distributed_mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(permutation[120]))
            {
                double permuted_node_120_x = distributed_mesh.GetNode(permutation[120])->rGetLocation()[0];
                double permuted_node_120_y = distributed_mesh.GetNode(permutation[120])->rGetLocation()[1];
                std::cout << "Original node = 120, permuted node = " << permutation[120]
                                          << ", location = " << permuted_node_120_x << ", " << permuted_node_120_y << "\n";
            }

            // Now do the test
            // Get the closest nodes to the same three points using the distributed mesh
            unsigned nearest_node_point_1 = distributed_mesh.GetNearestNodeIndex(point1);
            unsigned nearest_node_point_2 = distributed_mesh.GetNearestNodeIndex(point2);
            unsigned nearest_node_point_3 = distributed_mesh.GetNearestNodeIndex(point3);

            std::cout << "Process " << PetscTools::GetMyRank()
                      << ", node for point 1 = " << nearest_node_point_1
                      << ", node for point 2 = " << nearest_node_point_2
                      << ", node for point 3 = " << nearest_node_point_3 << "\n";

            // Check the results against the permuted node indices
            TS_ASSERT_EQUALS(distributed_mesh.GetNearestNodeIndex(point1), permutation[60]);
            TS_ASSERT_EQUALS(distributed_mesh.GetNearestNodeIndex(point2), permutation[120]);   // Closest node is top right node
            TS_ASSERT_EQUALS(distributed_mesh.GetNearestNodeIndex(point3), permutation[60]);
        }
    }

    void TestNearestNodeIndex1D()
    {
        // Similar to TestPointinMesh2D from TestTetrahedralMesh

        ChastePoint<1> point1(0.012);    // close to node 1 in non-dist mesh
        ChastePoint<1> point2(0.12);     // outside mesh, closest node is 10 (outermost node)
        ChastePoint<1> point3(0.01);     // node 1 of non-dist mesh

        // Read in a mesh and create distributed and non-distributed versions
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1mm_10_elements");
        TetrahedralMesh<1,1> non_distributed_mesh;
        DistributedTetrahedralMesh<1,1> distributed_mesh;
        non_distributed_mesh.ConstructFromMeshReader(mesh_reader);
        distributed_mesh.ConstructFromMeshReader(mesh_reader);

        // Check that the closest nodes in the non-dist mesh are the same as those tested in TestPointinMesh2D
        TS_ASSERT_EQUALS(non_distributed_mesh.GetNearestNodeIndex(point1), 1u);
        TS_ASSERT_EQUALS(non_distributed_mesh.GetNearestNodeIndex(point2), 10u);   // Closest node is top right node
        TS_ASSERT_EQUALS(non_distributed_mesh.GetNearestNodeIndex(point3), 1u);


        assert(distributed_mesh.rGetNodePermutation().size()==0);
        TS_ASSERT_EQUALS(distributed_mesh.GetNearestNodeIndex(point1), 1u);
        TS_ASSERT_EQUALS(distributed_mesh.GetNearestNodeIndex(point2), 10u);   // Closest node is top right node
        TS_ASSERT_EQUALS(distributed_mesh.GetNearestNodeIndex(point3), 1u);
    }


    void TestConstructParallelCuboidMeshGeometricPartition()
    {
        // 1D
        unsigned width = 4*PetscTools::GetNumProcs() - 1;
        ChastePoint<1> lower_1d((double)PetscTools::GetMyRank()*(1+width/(double)PetscTools::GetNumProcs()));
        ChastePoint<1> upper_1d((double)(1+PetscTools::GetMyRank())*(1+width/(double)PetscTools::GetNumProcs()));
        ChasteCuboid<1> cuboid_1d(lower_1d, upper_1d);

        TetrahedralMesh<1,1> base_mesh;
        base_mesh.ConstructLinearMesh(width);
        TrianglesMeshWriter<1,1> mesh_writer("", "linear_geometric");
        mesh_writer.WriteFilesUsingMesh(base_mesh);
        PetscTools::Barrier();
        std::string output_dir_1d = mesh_writer.GetOutputDirectory();

        TrianglesMeshReader<1,1> mesh_reader(output_dir_1d+"linear_geometric");
        DistributedTetrahedralMesh<1,1> read_mesh_1d(DistributedTetrahedralMeshPartitionType::GEOMETRIC);
        read_mesh_1d.SetProcessRegion(&cuboid_1d);
        read_mesh_1d.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<1,1> constructed_mesh_1d(DistributedTetrahedralMeshPartitionType::GEOMETRIC);
        TS_ASSERT_THROWS_THIS(constructed_mesh_1d.ConstructLinearMesh(width), "Space region not set for GEOMETRIC partition of DistributedTetrahedralMesh");
        constructed_mesh_1d.SetProcessRegion(&cuboid_1d);

        constructed_mesh_1d.ConstructLinearMesh(width);

        // Basic check that we have the same number of nodes.
        TS_ASSERT_EQUALS(constructed_mesh_1d.GetNumLocalNodes(),read_mesh_1d.GetNumLocalNodes());

        // 2D
        width = 2;
        unsigned height = 4*PetscTools::GetNumProcs()-1;
        ChastePoint<2> lower_2d(0.0, (double)PetscTools::GetMyRank()*(1+height/(double)PetscTools::GetNumProcs()));
        ChastePoint<2> upper_2d(1+width, (double)(1+PetscTools::GetMyRank())*(1+height/(double)PetscTools::GetNumProcs()));
        ChasteCuboid<2> cuboid_2d(lower_2d, upper_2d);

        TetrahedralMesh<2,2> base_mesh_2d;
        base_mesh_2d.ConstructRectangularMesh(width, height);
        TrianglesMeshWriter<2,2> mesh_writer_2d("", "rectangular_geometric");
        mesh_writer_2d.WriteFilesUsingMesh(base_mesh_2d);
        PetscTools::Barrier();
        std::string output_dir_2d = mesh_writer_2d.GetOutputDirectory();
        TrianglesMeshReader<2,2> mesh_reader_2d(output_dir_2d+"rectangular_geometric");
        DistributedTetrahedralMesh<2,2> read_mesh_2d(DistributedTetrahedralMeshPartitionType::GEOMETRIC);
        read_mesh_2d.SetProcessRegion(&cuboid_2d);
        read_mesh_2d.ConstructFromMeshReader(mesh_reader_2d);

        DistributedTetrahedralMesh<2,2> constructed_mesh_2d(DistributedTetrahedralMeshPartitionType::GEOMETRIC);
    TS_ASSERT_THROWS_THIS(constructed_mesh_2d.ConstructRectangularMesh(width, height), "Space region not set for GEOMETRIC partition of DistributedTetrahedralMesh");
        constructed_mesh_2d.SetProcessRegion(&cuboid_2d);

        constructed_mesh_2d.ConstructRectangularMesh(width, height);

        // Basic check that we have the same number of nodes.
        TS_ASSERT_EQUALS(constructed_mesh_2d.GetNumLocalNodes(),read_mesh_2d.GetNumLocalNodes());


        // 3D
        height = 3;
        unsigned depth = 4*PetscTools::GetNumProcs()-1;
        ChastePoint<3> lower_3d(0.0, 0.0, (double)PetscTools::GetMyRank()*(1+depth/(double)PetscTools::GetNumProcs()));
        ChastePoint<3> upper_3d(1+width, 1+height, (double)(1+PetscTools::GetMyRank())*(1+depth/(double)PetscTools::GetNumProcs()));
        ChasteCuboid<3> cuboid_3d(lower_3d, upper_3d);

        TetrahedralMesh<3,3> base_mesh_3d;
        base_mesh_3d.ConstructCuboid(width, height, depth);
        TrianglesMeshWriter<3,3> mesh_writer_3d("", "cuboid_geometric");
        mesh_writer_3d.WriteFilesUsingMesh(base_mesh_3d);
        PetscTools::Barrier();
        std::string output_dir_3d = mesh_writer_3d.GetOutputDirectory();
        TrianglesMeshReader<3,3> mesh_reader_3d(output_dir_3d+"cuboid_geometric");
        DistributedTetrahedralMesh<3,3> read_mesh_3d(DistributedTetrahedralMeshPartitionType::GEOMETRIC);
        read_mesh_3d.SetProcessRegion(&cuboid_3d);
        read_mesh_3d.ConstructFromMeshReader(mesh_reader_3d);

        DistributedTetrahedralMesh<3,3> constructed_mesh_3d(DistributedTetrahedralMeshPartitionType::GEOMETRIC);
    TS_ASSERT_THROWS_THIS(constructed_mesh_3d.ConstructCuboid(width, height, depth), "Space region not set for GEOMETRIC partition of DistributedTetrahedralMesh");
        constructed_mesh_3d.SetProcessRegion(&cuboid_3d);

        constructed_mesh_3d.ConstructCuboid(width, height, depth);

        // Basic check that we have the same number of nodes.
        TS_ASSERT_EQUALS(constructed_mesh_3d.GetNumLocalNodes(),read_mesh_3d.GetNumLocalNodes());
    }

    void TestConstructLinearMeshSmallest()
    {
        DistributedTetrahedralMesh<1,1> smallest_mesh;
        smallest_mesh.ConstructLinearMesh(1);

        TS_ASSERT_EQUALS(smallest_mesh.GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(smallest_mesh.GetNumBoundaryElements(),  2u);
        TS_ASSERT_EQUALS(smallest_mesh.GetNumElements(), 1u);
        unsigned owned = smallest_mesh.GetDistributedVectorFactory()->GetLocalOwnership();
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalNodes(), 2u);
            TS_ASSERT_EQUALS(owned, 2u);
            TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalElements(), 1u);
        }
        else
        {
            if (PetscTools::GetMyRank() < 2u)
            {
                TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalNodes(), 1u);
                TS_ASSERT_EQUALS(owned, 1u);
                TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalElements(), 1u);
            }
            else
            {
                //Own nothing
                TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalNodes(), 0u);
                TS_ASSERT_EQUALS(owned, 0u);
                TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalElements(), 0u);
            }
        }

        // Note that only processes ranked 0 or 1 will get nodes, so with more than
        // two processes some processes own nothing
        ChastePoint<1> origin(0.0);
        TS_ASSERT_EQUALS(smallest_mesh.GetNearestNodeIndex(origin), 0u);

    }

    void TestConstructRectangularMeshSmallest()
    {
        DistributedTetrahedralMesh<2,2> smallest_mesh;
        smallest_mesh.ConstructRectangularMesh(1,1);

        TS_ASSERT_EQUALS(smallest_mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(smallest_mesh.GetNumBoundaryElements(),  4u);
        TS_ASSERT_EQUALS(smallest_mesh.GetNumElements(), 2u);
        unsigned owned = smallest_mesh.GetDistributedVectorFactory()->GetLocalOwnership();
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalNodes(), 4u);
            TS_ASSERT_EQUALS(owned, 4u);
            TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalElements(), 2u);
        }
        else
        {
            if (PetscTools::GetMyRank() < 2u)
            {
                TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalNodes(), 2u);
                TS_ASSERT_EQUALS(owned, 2u);
                TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalElements(), 2u);
            }
            else
            {
                //Own nothing
                TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalNodes(), 0u);
                TS_ASSERT_EQUALS(owned, 0u);
                TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalElements(), 0u);
            }
        }
    }

    void TestConstructCuboidMeshSmallest()
    {
        DistributedTetrahedralMesh<3,3> smallest_mesh;
        smallest_mesh.ConstructCuboid(1,1,1);

        TS_ASSERT_EQUALS(smallest_mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(smallest_mesh.GetNumBoundaryElements(),  12u);
        TS_ASSERT_EQUALS(smallest_mesh.GetNumElements(), 6u);
        unsigned owned = smallest_mesh.GetDistributedVectorFactory()->GetLocalOwnership();
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalNodes(), 8u);
            TS_ASSERT_EQUALS(owned, 8u);
            TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalElements(), 6u);
        }
        else
        {
            if (PetscTools::GetMyRank() < 2u)
            {
                TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalNodes(), 4u);
                TS_ASSERT_EQUALS(owned, 4u);
                TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalElements(), 6u);
            }
            else
            {
                // Own nothing
                TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalNodes(), 0u);
                TS_ASSERT_EQUALS(owned, 0u);
                TS_ASSERT_EQUALS(smallest_mesh.GetNumLocalElements(), 0u);
            }
        }
    }

    void TestParallelWriting1D()
    {
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_10_elements_with_attributes");
        TetrahedralMesh<1,1> sequential_mesh;
        sequential_mesh.ConstructFromMeshReader(reader);
        TrianglesMeshWriter<1,1> mesh_writer1("TestDistributedMeshWriter", "seq_line_10_elements");
        mesh_writer1.WriteFilesUsingMesh(sequential_mesh);

        DistributedTetrahedralMesh<1,1> distributed_mesh(DistributedTetrahedralMeshPartitionType::DUMB); //Makes sure that there is no permutation
        AbstractTetrahedralMesh<1,1> *p_distributed_mesh = &distributed_mesh; //Hide the fact that it's distributed from the compiler
        distributed_mesh.ConstructFromMeshReader(reader);
        TrianglesMeshWriter<1,1> mesh_writer2("TestDistributedMeshWriter", "par_line_10_elements", false);
        mesh_writer2.WriteFilesUsingMesh(*p_distributed_mesh);

        std::string output_dir = mesh_writer1.GetOutputDirectory();
        FileComparison(output_dir + "par_line_10_elements.node",
                       output_dir + "seq_line_10_elements.node").CompareFiles();
        FileComparison(output_dir + "par_line_10_elements.ele",
                       output_dir + "seq_line_10_elements.ele").CompareFiles();
    }

    void TestNodeExchange()
    {
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        std::vector<std::vector<unsigned> > nodes_to_send_per_process;
        std::vector<std::vector<unsigned> > nodes_to_receive_per_process;
        mesh.CalculateNodeExchange(nodes_to_send_per_process, nodes_to_receive_per_process);

        TS_ASSERT_EQUALS(nodes_to_send_per_process.size(), PetscTools::GetNumProcs());
        TS_ASSERT_EQUALS(nodes_to_receive_per_process.size(), PetscTools::GetNumProcs());
        TS_ASSERT(nodes_to_receive_per_process[PetscTools::GetMyRank()].empty());
        TS_ASSERT(nodes_to_send_per_process[PetscTools::GetMyRank()].empty());

        // Do some communication

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
    }

    void TestParallelWriting3D()
    {
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> sequential_mesh;
        sequential_mesh.ConstructFromMeshReader(reader);
        TrianglesMeshWriter<3,3> mesh_writer1("TestDistributedMeshWriter", "seq_cube_2mm_12_elements", false);
        mesh_writer1.WriteFilesUsingMesh(sequential_mesh);

        DistributedTetrahedralMesh<3,3> distributed_mesh(DistributedTetrahedralMeshPartitionType::DUMB); //Makes sure that there is no permutation
        AbstractTetrahedralMesh<3,3> *p_distributed_mesh = &distributed_mesh; //Hide the fact that it's distributed from the compiler

        distributed_mesh.ConstructFromMeshReader(reader);
        TrianglesMeshWriter<3,3> mesh_writer2("TestDistributedMeshWriter", "par_cube_2mm_12_elements", false);
        mesh_writer2.WriteFilesUsingMesh(*p_distributed_mesh);

        std::string output_dir = mesh_writer1.GetOutputDirectory();

        std::vector<std::string> files_to_compare;
        files_to_compare.push_back("node");
        files_to_compare.push_back("ele");
        files_to_compare.push_back("face");

        for (unsigned i=0; i<files_to_compare.size(); i++)
        {
            std::cout << "Comparing ." << files_to_compare[i] << std::endl;
            FileFinder generated_parallel(output_dir + "/par_cube_2mm_12_elements." + files_to_compare[i]);
            FileFinder generated_sequential(output_dir + "/seq_cube_2mm_12_elements." + files_to_compare[i]);
            FileComparison comparer(generated_parallel,generated_sequential);
            TS_ASSERT(comparer.CompareFiles());
        }
    }

    void TestEfficientParallelWriting3D()
    {
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> sequential_mesh;
        sequential_mesh.ConstructFromMeshReader(reader);

        DistributedTetrahedralMesh<3,3> distributed_mesh(DistributedTetrahedralMeshPartitionType::DUMB); //Makes sure that there is no permutation
        AbstractTetrahedralMesh<3,3> *p_distributed_mesh = &distributed_mesh; //Hide the fact that it's distributed from the compiler

        distributed_mesh.ConstructFromMeshReader(reader);

        // Meshalyzer
        MeshalyzerMeshWriter<3,3> mesh_writer("TestDistributedMeshWriter", "seq_cube_2mm_12_elements", false);
        mesh_writer.WriteFilesUsingMesh(sequential_mesh);
        std::string output_dir = mesh_writer.GetOutputDirectory();
        MeshalyzerMeshWriter<3,3> mesh_writer_par("TestDistributedMeshWriter", "par_efficient_cube_2mm_12_elements", false);
        mesh_writer_par.WriteFilesUsingMesh(*p_distributed_mesh, false);  //"false == Don't preserve element ordering
        {
            FileFinder parallel(output_dir + "/par_efficient_cube_2mm_12_elements.pts");
            FileFinder sequential(output_dir + "/seq_cube_2mm_12_elements.pts");
            FileComparison comparer(parallel,sequential);
            TS_ASSERT(comparer.CompareFiles());

            ComparePermutedFiles(output_dir + "/seq_cube_2mm_12_elements.tetras", output_dir + "/par_efficient_cube_2mm_12_elements.tetras");
            ComparePermutedFiles(output_dir + "/seq_cube_2mm_12_elements.tri", output_dir + "/par_efficient_cube_2mm_12_elements.tri");
        }

        // Cool graphics
        MeshalyzerMeshWriter<3,3> mesh_writer_cg("TestDistributedMeshWriter", "seq_cube_2mm_12_elements_cg", false, true); //Don't clean, Do use CG
        mesh_writer_cg.WriteFilesUsingMesh(sequential_mesh);
        MeshalyzerMeshWriter<3,3> mesh_writer_par_cg("TestDistributedMeshWriter", "par_efficient_cube_2mm_12_elements_cg", false, true); //Don't clean, Do use CG
        mesh_writer_par_cg.WriteFilesUsingMesh(*p_distributed_mesh, false);  //"false == Don't preserve element ordering
        {
            // cg output is indexed from 1, but the pts file doesn't have indices
            FileFinder parallel(output_dir + "/par_efficient_cube_2mm_12_elements_cg.pts");
            FileFinder sequential(output_dir + "/seq_cube_2mm_12_elements_cg.pts");
            FileComparison comparer(parallel,sequential);
            TS_ASSERT(comparer.CompareFiles());
        }

        //cmgui
        CmguiMeshWriter<3,3> cmgui_writer("TestDistributedMeshWriter", "seq_cube_2mm_12_elements_cmgui", false); //Don't clean
        cmgui_writer.WriteFilesUsingMesh(sequential_mesh);
        CmguiMeshWriter<3,3> cmgui_writer_par("TestDistributedMeshWriter", "par_efficient_cube_2mm_12_elements_cmgui", false);
        cmgui_writer_par.WriteFilesUsingMesh(*p_distributed_mesh, false);  //"false == Don't preserve element ordering
        {
            FileFinder parallel(output_dir + "/par_efficient_cube_2mm_12_elements_cmgui.exnode");
            FileFinder sequential(output_dir + "/seq_cube_2mm_12_elements_cmgui.exnode");
            FileComparison comparer(parallel,sequential);
            comparer.SetIgnoreLinesBeginningWith("Group name:");
            TS_ASSERT(comparer.CompareFiles());
            ComparePermutedFiles(output_dir + "/seq_cube_2mm_12_elements_cmgui.exelem", output_dir + "/par_efficient_cube_2mm_12_elements_cmgui.exelem");
        }
    }

    void TestArchiveOfConstructedMesh()
    {
        FileFinder archive_dir("distributed_tetrahedral_mesh_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "distributed_rectangle.arch";
        ArchiveLocationInfo::SetMeshFilename("distributed_rectangle");

        DistributedTetrahedralMesh<2,2>* p_mesh = new DistributedTetrahedralMesh<2,2>;
        //std::vector<unsigned> halo_node_indices;
        std::vector<Node<2>*> halo_nodes;
        unsigned num_nodes;
        unsigned local_num_nodes;
        unsigned num_elements;

        unsigned width = 4;
        unsigned height = 4*PetscTools::GetNumProcs()-1; // 4*NumProcs layers of nodes (ensure dumb partition works in slices)
        // Archive
        {
            p_mesh->ConstructRectangularMesh(width, height);
            num_nodes = p_mesh->GetNumNodes();
            local_num_nodes = p_mesh->GetNumLocalNodes();
            num_elements = p_mesh->GetNumElements();

            halo_nodes = p_mesh->mHaloNodes;

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractTetrahedralMesh<2,2>* const p_mesh_abstract = static_cast<AbstractTetrahedralMesh<2,2>* >(p_mesh);
            (*p_arch) << p_mesh_abstract;
        }

        // Restore
        {
            // Should archive the most abstract class you can to check boost knows what individual classes are.
            // (but here AbstractMesh doesn't have the methods below).
            AbstractTetrahedralMesh<2,2>* p_mesh_abstract2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_mesh_abstract2;

            // Check we have the right number of nodes & elements
            DistributedTetrahedralMesh<2,2>* p_mesh2 = static_cast<DistributedTetrahedralMesh<2,2>*>(p_mesh_abstract2);

            TS_ASSERT_EQUALS(p_mesh2->GetNumNodes(), num_nodes);
            TS_ASSERT_EQUALS(p_mesh2->GetNumLocalNodes(), local_num_nodes);
            TS_ASSERT_EQUALS(p_mesh2->GetNumElements(), num_elements);

            CompareParallelMeshOwnership(*p_mesh, *p_mesh2);

            // Check some node co-ordinates
            try
            {
                Node<2>* p_node1 = p_mesh->GetNode(0);
                Node<2>* p_node2 = p_mesh2->GetNode(0);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], p_node2->GetPoint()[0], 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[1], p_node2->GetPoint()[1], 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], 0.0, 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[1], 0.0, 1e-6);
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            try
            {
                Node<2>* p_node1 = p_mesh->GetNode((width+1)*(height+1)-1);
                Node<2>* p_node2 = p_mesh2->GetNode((width+1)*(height+1)-1);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], p_node2->GetPoint()[0], 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], p_node2->GetPoint()[0], 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], (double) width, 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[1], (double) height, 1e-6);
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            // Check first element has the right nodes
            try
            {
                Element<2,2>* p_element = p_mesh->GetElement(0);
                Element<2,2>* p_element2 = p_mesh2->GetElement(0);
                TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), p_element2->GetNodeGlobalIndex(0));
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            try
            {
                Element<2,2>* p_element = p_mesh->GetElement(width*height);
                Element<2,2>* p_element2 = p_mesh2->GetElement(500);
                TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), p_element2->GetNodeGlobalIndex(0));
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            // Check the halo nodes are right
            std::vector<Node<2>*> halo_nodes2 = p_mesh2->mHaloNodes;
            TS_ASSERT_EQUALS(halo_nodes2.size(), halo_nodes.size());
            delete p_mesh2;
        }
        delete p_mesh;
    }

    void TestLoadBadFacesException()
    {
        DistributedTetrahedralMesh<3,3> distributed_mesh_bad;
        TrianglesMeshReader<3,3> mesh_reader_bad("mesh/test/data/cube_21_nodes_side/Cube21_bad_faces"); // 5x5x5mm cube (internode distance = 0.25mm)
        if (PetscTools::IsSequential())
        {
            distributed_mesh_bad.ConstructFromMeshReader(mesh_reader_bad);
            TS_ASSERT_EQUALS(distributed_mesh_bad.GetNumNodes(), 9261u); // 21x21x21 nodes
            TS_ASSERT_EQUALS(distributed_mesh_bad.GetNumElements(), 48000u);
            TS_ASSERT_EQUALS(distributed_mesh_bad.GetNumBoundaryElements(), 4800u);
        }
        else
        {
            TS_ASSERT_THROWS_CONTAINS(distributed_mesh_bad.ConstructFromMeshReader(mesh_reader_bad), "Face does not appear in element file (Face ");
        }

        DistributedTetrahedralMesh<3,3> distributed_mesh_good;
        TrianglesMeshReader<3,3> mesh_reader_good("mesh/test/data/cube_21_nodes_side/Cube21"); // 5x5x5mm cube (internode distance = 0.25mm)
        distributed_mesh_good.ConstructFromMeshReader(mesh_reader_good);
        TS_ASSERT_EQUALS(distributed_mesh_good.GetNumNodes(), 9261u); // 21x21x21 nodes
        TS_ASSERT_EQUALS(distributed_mesh_good.GetNumElements(), 48000u);
        TS_ASSERT_EQUALS(distributed_mesh_good.GetNumBoundaryElements(), 4800u);
    }

    void TestWritingDistributedMeshBinary()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_152_elements_v2");

        // This test will only pass if the node and element orderings are preserved (i.e. dumb partition)
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<3,3> mesh_writer("WritingDistributedMeshBinary", "3dDistributedMesh");
        mesh_writer.SetWriteFilesAsBinary();

        mesh_writer.WriteFilesUsingMesh(mesh);

        std::string output_dir = mesh_writer.GetOutputDirectory();
//        TrianglesMeshReader<3,3> mesh_reader2(output_dir + "3dDistributedMesh");

//        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 312u);
//        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 522u);
//        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), mesh_reader.GetNumNodes());
//        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), mesh_reader.GetNumElements());

        // Test for connectivity
        ///\todo #1621 use the mesh reader when it's written
        FileFinder generated(output_dir + "3dDistributedMesh.ncl");
        FileFinder reference("mesh/test/data/cube_2mm_152_elements_binary_v2.ncl",RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated,reference);
        TS_ASSERT(comparer.CompareFiles());

    }

    void TestCheckOutwardNormals()
    {
        {
            DistributedTetrahedralMesh<2,2> mesh;
            mesh.ConstructRectangularMesh(2, 4);
            mesh.CheckOutwardNormals();
        }
        {
            DistributedTetrahedralMesh<3,3> mesh;
            mesh.ConstructCuboid(2, 2, 4);
            mesh.CheckOutwardNormals();
        }
    }

    void TestCalculateEdgeLengths()
    {
        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
            DistributedTetrahedralMesh<3,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            c_vector <double, 2> edge_len = mesh.CalculateMinMaxEdgeLengths();
            TS_ASSERT_DELTA(edge_len[0], 0.1, 1e-4);
            TS_ASSERT_DELTA(edge_len[1], 0.2828, 1e-4);
        }
        {
            DistributedTetrahedralMesh<2,2> mesh;
            mesh.ConstructRectangularMesh(2,3);
            c_vector <double, 2> edge_len = mesh.CalculateMinMaxEdgeLengths();
            TS_ASSERT_DELTA(edge_len[0], 1.0, 1e-5);
            TS_ASSERT_DELTA(edge_len[1], sqrt(2.0), 1e-5);
        }
        {
            TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
            DistributedTetrahedralMesh<2,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            c_vector <double, 2> edge_len = mesh.CalculateMinMaxEdgeLengths();
            TS_ASSERT_DELTA(edge_len[0], 0.0628, 1e-4);
            TS_ASSERT_DELTA(edge_len[1], 0.2010, 1e-4);
        }
    }

    void TestPartitionHasNoElementsWithAllHaloNodes()
    {
        TrianglesMeshReader<3,3> reader("heart/test/data/box_shaped_heart/box_heart");
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);
        mesh.ConstructFromMeshReader(reader);

        // Loop over all nodes in all elements. If any element is made of all-halo-nodes, hit assertion.
        for (AbstractTetrahedralMesh<3,3>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            bool any_local = false;
            for (unsigned node_local_index=0; node_local_index < 4u; node_local_index++)
            {
                unsigned node_global_index = iter->GetNodeGlobalIndex(node_local_index);
                if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(node_global_index))
                {
                    any_local = true;
                }
            }
            TS_ASSERT(any_local);
        }
    }
    /* HOW_TO_TAG Mesh
     * Construct a distributed regular mesh (rectangle in 2D or cuboid in 3D) which does not have a default
     * split plane.  The default is for parallel code to split 2-D meshes into slices in the y-dimension
     * and 3-D meshes in the z-dimension.
     */
    void TestConstructSlabMeshWithDimensionSplit()
    {
        double step = 1.0;
        unsigned width = 3;
        unsigned height = 5;
        unsigned depth = 7;

        // In 1D we shouldn't be able to change the split dimension from 0.  (Can only split in x.)
        {
            DistributedTetrahedralMesh<1,1> mesh;
            TS_ASSERT_THROWS_THIS(mesh.ConstructRegularSlabMeshWithDimensionSplit(1, step, width), "Cannot split on non-existent dimension");
            mesh.ConstructRegularSlabMeshWithDimensionSplit(0, step, width);
         }

        {
            // Splitting dimension is not specified
            DistributedTetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(step, width, height);
            // Splitting dimension is y-dimension (the default)
            DistributedTetrahedralMesh<2,2> mesh_with_default_split;
            mesh_with_default_split.ConstructRegularSlabMeshWithDimensionSplit(1, step, width, height);
            // Splitting dimension is x-dimension
            DistributedTetrahedralMesh<2,2> mesh_with_x_split;
            mesh_with_x_split.ConstructRegularSlabMeshWithDimensionSplit(0, step, width, height);

            // Check that mesh and mesh_with_default_split are identical
            CompareMeshes(mesh, mesh_with_default_split);

            // Check that the y-split has the same bounding box
            ChasteCuboid<2> bounds = mesh.CalculateBoundingBox();
            ChasteCuboid<2> bounds_with_x_split = mesh_with_x_split.CalculateBoundingBox();
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], bounds_with_x_split.rGetUpperCorner()[0], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], bounds_with_x_split.rGetUpperCorner()[1], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], bounds_with_x_split.rGetLowerCorner()[0], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], bounds_with_x_split.rGetLowerCorner()[1], 1e-6);

            // Same amount of stuff
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_with_x_split.GetNumNodes());
            TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_with_x_split.GetNumElements());
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_with_x_split.GetNumBoundaryElements());

            // Show that the y-split has a different indexing scheme
            // Normal meshes start at the origin
            if (mesh_with_default_split.GetDistributedVectorFactory()->IsGlobalIndexLocal(0u))
            {
                c_vector<double, 2> orig1 = mesh_with_default_split.GetNode(0u)->rGetLocation();
                TS_ASSERT_DELTA(orig1[0], 0.0, 1e-5);
                TS_ASSERT_DELTA(orig1[1], 0.0, 1e-5);
            }

            // The new one has the origin at an index height away
            if (mesh_with_x_split.GetDistributedVectorFactory()->IsGlobalIndexLocal(height))
            {
                c_vector<double, 2> orig2 = mesh_with_x_split.GetNode(height)->rGetLocation();
                TS_ASSERT_DELTA(orig2[0], 0.0, 1e-5);
                TS_ASSERT_DELTA(orig2[1], 0.0, 1e-5);
            }

            // Check that we are really splitting on x, by checking every process' split reaches the top of y.
            double max_x_with_x_split = DBL_MIN;
            double max_y_with_x_split = DBL_MIN;
            for (AbstractTetrahedralMesh<2,2>::NodeIterator iter = mesh_with_x_split.GetNodeIteratorBegin();
                 iter != mesh_with_x_split.GetNodeIteratorEnd();
                 ++iter)
            {
                c_vector<double, 2> pos;
                pos = iter->rGetLocation();
                max_x_with_x_split = std::max(max_x_with_x_split, pos[0]);
                max_y_with_x_split = std::max(max_y_with_x_split, pos[1]);
            }
            TS_ASSERT_DELTA(max_y_with_x_split, height, 1e-6);
            if (PetscTools::AmTopMost())
            {
                TS_ASSERT_DELTA(max_x_with_x_split, width, 1e-6);
            }
            else
            {
                TS_ASSERT_LESS_THAN(max_x_with_x_split, width - step/2.0);
            }

            // We also sanity check that the default split does split on y.
            double max_x_with_default_split = DBL_MIN;
            double max_y_with_default_split = DBL_MIN;
            for (AbstractTetrahedralMesh<2,2>::NodeIterator iter = mesh_with_default_split.GetNodeIteratorBegin();
                 iter != mesh_with_default_split.GetNodeIteratorEnd();
                 ++iter)
            {
                c_vector<double, 2> pos;
                pos = iter->rGetLocation();
                max_x_with_default_split = std::max(max_x_with_default_split, pos[0]);
                max_y_with_default_split = std::max(max_y_with_default_split, pos[1]);
            }
            TS_ASSERT_DELTA(max_x_with_default_split, width, 1e-6);
            if (PetscTools::AmTopMost())
            {
                TS_ASSERT_DELTA(max_y_with_default_split, height, 1e-6);
            }
            else
            {
                TS_ASSERT_LESS_THAN(max_y_with_default_split, height - step/2.0);
            }

            c_matrix<double, 2, 2> dummy;
            // Check that all the elements (including halos) are the same size -- TetrahedralMesh has an active RefreshMesh but DistributedTetrahedralMesh does not
            for (AbstractTetrahedralMesh<2, 2>::ElementIterator iter = mesh_with_x_split.GetElementIteratorBegin();
                 iter != mesh_with_x_split.GetElementIteratorEnd();
                 ++iter)
            {
                double jacobian_det;
                iter->CalculateJacobian(dummy, jacobian_det);
                TS_ASSERT_DELTA(jacobian_det, 1.0, 1e-6)
            }

        }

        {
            DistributedTetrahedralMesh<3,3> mesh;
            mesh.ConstructRegularSlabMesh(step, width, height, depth);
            DistributedTetrahedralMesh<3,3> mesh_with_default_split;
            mesh_with_default_split.ConstructRegularSlabMeshWithDimensionSplit(2, step, width, height, depth);
            DistributedTetrahedralMesh<3,3> mesh_with_x_split;
            mesh_with_x_split.ConstructRegularSlabMeshWithDimensionSplit(0, step, width, height, depth);
            DistributedTetrahedralMesh<3,3> mesh_with_y_split;
            mesh_with_y_split.ConstructRegularSlabMeshWithDimensionSplit(1, step, width, height, depth);

            // Check that mesh and mesh_with_default_split are identical
            CompareMeshes(mesh, mesh_with_default_split);

            // Check that the x-split and y-split have the same bounding box
            ChasteCuboid<3> bounds = mesh.CalculateBoundingBox();
            ChasteCuboid<3> bounds_with_x_split = mesh_with_x_split.CalculateBoundingBox();
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], bounds_with_x_split.rGetUpperCorner()[0], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], bounds_with_x_split.rGetUpperCorner()[1], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[2], bounds_with_x_split.rGetUpperCorner()[2], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], bounds_with_x_split.rGetLowerCorner()[0], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], bounds_with_x_split.rGetLowerCorner()[1], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[2], bounds_with_x_split.rGetLowerCorner()[2], 1e-6);
            ChasteCuboid<3> bounds_with_y_split = mesh_with_y_split.CalculateBoundingBox();
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], bounds_with_y_split.rGetUpperCorner()[0], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], bounds_with_y_split.rGetUpperCorner()[1], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetUpperCorner()[2], bounds_with_y_split.rGetUpperCorner()[2], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], bounds_with_y_split.rGetLowerCorner()[0], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], bounds_with_y_split.rGetLowerCorner()[1], 1e-6);
            TS_ASSERT_DELTA(bounds.rGetLowerCorner()[2], bounds_with_y_split.rGetLowerCorner()[2], 1e-6);

            // Show that the splits have different indexing schemes
            // Normal meshes start at the origin
            if (mesh_with_default_split.GetDistributedVectorFactory()->IsGlobalIndexLocal(width))
            {
                c_vector<double, 3> xaxis1 = mesh_with_default_split.GetNode(width)->rGetLocation();
                TS_ASSERT_DELTA(xaxis1[0], (double) width, 1e-5);
                TS_ASSERT_DELTA(xaxis1[1], 0.0, 1e-5);
                TS_ASSERT_DELTA(xaxis1[2], 0.0, 1e-5);
            }

            // The x split one indexes in y, then z so x-axis is in top layer
            unsigned expected_xaxis2_index = width*(depth+1)*(height+1);
            if (mesh_with_x_split.GetDistributedVectorFactory()->IsGlobalIndexLocal(expected_xaxis2_index))
            {
                c_vector<double, 3> xaxis2 = mesh_with_x_split.GetNode(expected_xaxis2_index)->rGetLocation();
                TS_ASSERT_DELTA(xaxis2[0], (double) width, 1e-5);
                TS_ASSERT_DELTA(xaxis2[1], 0.0, 1e-5);
                TS_ASSERT_DELTA(xaxis2[2], 0.0, 1e-5);
            }
            // the y split indexes in z first so x-axis is the end of the first layer
            unsigned expected_xaxis3_index = width*(depth+1);
            if (mesh_with_y_split.GetDistributedVectorFactory()->IsGlobalIndexLocal(expected_xaxis3_index))
            {
                c_vector<double, 3> xaxis3 = mesh_with_y_split.GetNode(expected_xaxis3_index)->rGetLocation();
                TS_ASSERT_DELTA(xaxis3[0], (double) width, 1e-5);
                TS_ASSERT_DELTA(xaxis3[1], 0.0, 1e-5);
                TS_ASSERT_DELTA(xaxis3[2], 0.0, 1e-5);
            }

            // Check that we are really splitting on x, by checking every process' split reaches the top of y,z.
            double max_x_with_x_split = DBL_MIN;
            double max_y_with_x_split = DBL_MIN;
            double max_z_with_x_split = DBL_MIN;
            for (AbstractTetrahedralMesh<3,3>::NodeIterator iter = mesh_with_x_split.GetNodeIteratorBegin();
                 iter != mesh_with_x_split.GetNodeIteratorEnd();
                 ++iter)
            {
                c_vector<double, 3> pos;
                pos = iter->rGetLocation();
                max_x_with_x_split = std::max(max_x_with_x_split, pos[0]);
                max_y_with_x_split = std::max(max_y_with_x_split, pos[1]);
                max_z_with_x_split = std::max(max_z_with_x_split, pos[2]);
            }
            TS_ASSERT_DELTA(max_y_with_x_split, height, 1e-6);
            TS_ASSERT_DELTA(max_z_with_x_split, depth, 1e-6);
            if (PetscTools::AmTopMost())
            {
                TS_ASSERT_DELTA(max_x_with_x_split, width, 1e-6);
            }
            else
            {
                TS_ASSERT_LESS_THAN(max_x_with_x_split, width - step/2.0);
            }

            // Check that we are really splitting on y, by checking every process' split reaches the top of x,z.
            double max_x_with_y_split = DBL_MIN;
            double max_y_with_y_split = DBL_MIN;
            double max_z_with_y_split = DBL_MIN;
            for (AbstractTetrahedralMesh<3,3>::NodeIterator iter = mesh_with_y_split.GetNodeIteratorBegin();
                 iter != mesh_with_y_split.GetNodeIteratorEnd();
                 ++iter)
            {
                c_vector<double, 3> pos;
                pos = iter->rGetLocation();
                max_x_with_y_split = std::max(max_x_with_y_split, pos[0]);
                max_y_with_y_split = std::max(max_y_with_y_split, pos[1]);
                max_z_with_y_split = std::max(max_z_with_y_split, pos[2]);
            }
            TS_ASSERT_DELTA(max_x_with_y_split, width, 1e-6);
            TS_ASSERT_DELTA(max_z_with_y_split, depth, 1e-6);
            if (PetscTools::AmTopMost())
            {
                TS_ASSERT_DELTA(max_y_with_y_split, height, 1e-6);
            }
            else
            {
                TS_ASSERT_LESS_THAN(max_y_with_y_split, height - step/2.0);
            }

            // We also sanity check that the default split does split on z.
            double max_x_with_default_split = DBL_MIN;
            double max_y_with_default_split = DBL_MIN;
            double max_z_with_default_split = DBL_MIN;
            for (AbstractTetrahedralMesh<3,3>::NodeIterator iter = mesh_with_default_split.GetNodeIteratorBegin();
                 iter != mesh_with_default_split.GetNodeIteratorEnd();
                 ++iter)
            {
                c_vector<double, 3> pos;
                pos = iter->rGetLocation();
                max_x_with_default_split = std::max(max_x_with_default_split, pos[0]);
                max_y_with_default_split = std::max(max_y_with_default_split, pos[1]);
                max_z_with_default_split = std::max(max_z_with_default_split, pos[2]);
            }
            TS_ASSERT_DELTA(max_x_with_default_split, width, 1e-6);
            TS_ASSERT_DELTA(max_y_with_default_split, height, 1e-6);
            if (PetscTools::AmTopMost())
            {
                TS_ASSERT_DELTA(max_z_with_default_split, depth, 1e-6);
            }
            else
            {
                TS_ASSERT_LESS_THAN(max_z_with_default_split, depth - step/2.0);
            }

            c_matrix<double, 3, 3> dummy;
            // Check that all the elements (including halos) are the same size -- TetrahedralMesh has an active RefreshMesh but DistributedTetrahedralMesh does not
            for (AbstractTetrahedralMesh<3, 3>::ElementIterator iter = mesh_with_x_split.GetElementIteratorBegin();
                 iter != mesh_with_x_split.GetElementIteratorEnd();
                 ++iter)
            {
                double jacobian_det;
                iter->CalculateJacobian(dummy, jacobian_det);
                TS_ASSERT_DELTA(jacobian_det, 1.0, 1e-6)
            }
            for (AbstractTetrahedralMesh<3, 3>::ElementIterator iter = mesh_with_y_split.GetElementIteratorBegin();
                 iter != mesh_with_y_split.GetElementIteratorEnd();
                 ++iter)
            {
                double jacobian_det;
                iter->CalculateJacobian(dummy, jacobian_det);
                TS_ASSERT_DELTA(jacobian_det, 1.0, 1e-6)
            }
        }
    }
};

#endif /*TESTDISTRIBUTEDTETRAHEDRALMESH_HPP_*/
