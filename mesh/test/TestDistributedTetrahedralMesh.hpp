/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef TESTDISTRIBUTEDTETRAHEDRALMESH_HPP_
#define TESTDISTRIBUTEDTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <sstream>

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

#include "RandomNumberGenerator.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestDistributedTetrahedralMesh : public CxxTest::TestSuite
{
private:

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    void CompareMeshes( DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh1,
                        DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh2 )
    {
        // Check that we have the right number of nodes and elements
        TS_ASSERT_EQUALS(rMesh1.GetNumBoundaryElements(), rMesh2.GetNumBoundaryElements());
        TS_ASSERT_EQUALS(rMesh1.GetNumElements(), rMesh2.GetNumElements());
        TS_ASSERT_EQUALS(rMesh1.GetNumNodes(), rMesh2.GetNumNodes());

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
            unsigned nodes_owned[num_global_nodes];
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
                catch (Exception &e)
                {
                    nodes_owned[node_id] = 0;
                }
            }

            TS_ASSERT_EQUALS(rMesh.GetNumLocalNodes(), total_nodes_this_process);

            // Combine all the local maps by adding them up in the master process
            unsigned nodes_reduction[num_global_nodes];
            MPI_Reduce(&nodes_owned, &nodes_reduction, num_global_nodes, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);

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
            unsigned elements_owned[num_global_elements];

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
                catch(Exception& e)
                {
                    elements_owned[element_id] = 0;
                }
            }

            TS_ASSERT_EQUALS(rMesh.GetNumLocalElements(), total_elements_this_process);

            // Combine all the local maps by adding them up in the master process
            unsigned elements_reduction[num_global_elements];
            MPI_Reduce(&elements_owned, &elements_reduction, num_global_elements, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);

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
            unsigned b_elements_owned[num_global_b_elements];

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
                catch(Exception& e)
                {
                    b_elements_owned[b_element_id] = 0;
                }
            }

            TS_ASSERT_EQUALS(rMesh.GetNumLocalBoundaryElements(), total_b_elements_this_process);

            // Combine all the local maps by adding them up in the master process
            unsigned b_elements_reduction[num_global_b_elements];
            MPI_Reduce(&b_elements_owned, &b_elements_reduction, num_global_b_elements, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);

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
            //Metis may allocate no nodes to a partition if the mesh is small and there are many processes
            // Look out for "You just increased the maxndoms"
            TS_ASSERT( 0u == total_nodes_this_process );
            TS_ASSERT( 0u == total_elements_this_process );
            TS_ASSERT( 0u == total_b_elements_this_process );
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
                unsigned region = mesh.GetElement(i)->GetRegion();
                TS_ASSERT_EQUALS(region, i%5+1);
                TS_ASSERT_EQUALS(i, mesh.GetElement(i)->GetIndex());
            }
            catch(Exception& e)
            {
                // I don't own this element do I?
            }
        }
        TS_ASSERT_EQUALS(mesh.CalculateMaximumNodeConnectivityPerProcess(), 3U);
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
        catch(Exception& e)
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
        catch(Exception& e)
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
         * In this test we let METIS reorder the DistributedTetrahedralMesh. We want to check that although
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

    void TestConstructionFromMeshReaderWithNodeAttributes() throw(Exception)
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

        if (mesh.rGetNodePermutation().size() > 0)//need to figure out where they end up in permutation
        {
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

    void TestConstructFromMeshReaderWithBinaryFiles()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements_binary");
        TrianglesMeshReader<3,3> mesh_reader_ascii("mesh/test/data/cube_136_elements");

        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<3,3> mesh_from_ascii;
        mesh_from_ascii.ConstructFromMeshReader(mesh_reader_ascii);

        CompareMeshes(mesh, mesh_from_ascii);
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

    void TestEverythingIsAssignedMetisLibrary()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());

        CheckEverythingIsAssigned<3,3>(mesh);
    }

    void TestRandomShuffle() throw (Exception)
    {
        unsigned num_elts = 200;

        std::vector<unsigned> random_order(num_elts);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        p_gen->Reseed(0);
        p_gen->Shuffle(num_elts,random_order);

        unsigned my_entry;
        unsigned neighbours_entry;

        int num_procs = PetscTools::GetNumProcs();
        int my_rank = PetscTools::GetMyRank();
        int source_rank = (my_rank + num_procs - 1) % num_procs;
        int destination_rank = (my_rank + 1) % num_procs;
        int my_tag;
        int source_tag;

        MPI_Status status;

        for (unsigned element_number = 0; element_number < num_elts; element_number++)
        {
            my_entry = random_order[element_number];

            my_tag = my_rank + num_elts*element_number;
            source_tag = source_rank + num_elts*element_number;

            MPI_Send( &my_entry, 1, MPI_UNSIGNED, destination_rank, my_tag, PETSC_COMM_WORLD );
            MPI_Recv( &neighbours_entry, 1, MPI_UNSIGNED, source_rank, source_tag, PETSC_COMM_WORLD, &status );
            PetscTools::Barrier();

            TS_ASSERT_EQUALS( my_entry, neighbours_entry );
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
        unsigned num_local_nodes_petsc_parmetis, num_local_nodes_binary, num_local_nodes_metis, num_total_nodes;

        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_1mm_6000_elements");
            //TrianglesMeshReader<3,3> mesh_reader("heart/test/data/heart");
            DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);
            mesh.ConstructFromMeshReader(mesh_reader);

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
            TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());

            CheckEverythingIsAssigned<3,3>(mesh);

            num_local_nodes_metis = mesh.GetNumLocalNodes();
            num_total_nodes=mesh.GetNumNodes();
        }

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

        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_1mm_6000_elements_binary");
            //TrianglesMeshReader<3,3> mesh_reader("heart/test/data/heart_binary");
            DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);
            mesh.ConstructFromMeshReader(mesh_reader);

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
            TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());

            CheckEverythingIsAssigned<3,3>(mesh);

            num_local_nodes_binary = mesh.GetNumLocalNodes();
        }

        unsigned max_local_nodes_metis;
        unsigned max_local_nodes_petsc_parmetis;
        unsigned max_local_nodes_binary;

        MPI_Allreduce (&num_local_nodes_metis, &max_local_nodes_metis, 1, MPI_UNSIGNED, MPI_MAX, PETSC_COMM_WORLD );
        MPI_Allreduce (&num_local_nodes_petsc_parmetis, &max_local_nodes_petsc_parmetis, 1, MPI_UNSIGNED, MPI_MAX, PETSC_COMM_WORLD );
        MPI_Allreduce (&num_local_nodes_binary, &max_local_nodes_binary, 1, MPI_UNSIGNED, MPI_MAX, PETSC_COMM_WORLD );

        if (PetscTools::AmMaster())
        {
            std::cout << "METIS\tPETSC PARMETIS\tPARMETIS BINARY" << std::endl;
            std::cout << max_local_nodes_metis << "\t" << max_local_nodes_petsc_parmetis << "\t\t" << max_local_nodes_binary << std::endl;
        }
        PetscTools::Barrier();

        TS_ASSERT(num_local_nodes_petsc_parmetis <= max_local_nodes_binary);
        //Watch out for dumb partition and warn about it
        if (PetscTools::IsParallel())
        {
            //Dumb partition is ceil(n/p), ceil(n/p), .... [ceil(n/p) + n - p*ceil(n/p)]
            //i.e. most processes get ceil(n/p)  = floor((n+p-1)/p)
            unsigned max_in_dumb_partition = (num_total_nodes + PetscTools::GetNumProcs() - 1)/PetscTools::GetNumProcs();
            if (max_local_nodes_petsc_parmetis ==  max_in_dumb_partition)
            {
                TS_TRACE("That was dumb partition -- it did not use ParMETIS");
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
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::PETSC_MAT_PARTITION);
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());

        CheckEverythingIsAssigned<3,3>(mesh);
    }

    void TestEverythingIsAssignedPetscPartitionBinaryFiles()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements_binary");
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::PETSC_MAT_PARTITION);
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());

        CheckEverythingIsAssigned<3,3>(mesh);
    }


    void TestConstruct3DWithRegions() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart_nonnegative_flags");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            try
            {
                unsigned region = mesh.GetElement(i)->GetRegion();
                TS_ASSERT_EQUALS(region, (i+1)%3+1);
            }
            catch(Exception& e)
            {
                // I don't own this element do I?
            }
        }

        TS_ASSERT_EQUALS(mesh_reader.GetNumFaceAttributes(), 1u);

        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            try
            {
                unsigned region = mesh.GetBoundaryElement(i)->GetRegion();
                TS_ASSERT_LESS_THAN(region, 4u);
            }
            catch(Exception& e)
            {
                // I don't own this element do I?
            }
        }
    }

    void TestMetisPartitioning()
    {
        EXIT_IF_SEQUENTIAL;

        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
            DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);
            mesh.ConstructFromMeshReader(mesh_reader);

            // Check that each processor owns the number of nodes corresponding to its METIS partition
            unsigned local_nodes = mesh.GetDistributedVectorFactory()->GetLocalOwnership();
            TS_ASSERT_EQUALS(local_nodes, mesh.GetNumLocalNodes());

            TS_ASSERT_EQUALS(mesh.GetPartitionType(), DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);
        }

        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
            DistributedTetrahedralMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);
            mesh.ConstructFromMeshReader(mesh_reader);

            // Check that each processor owns the number of nodes corresponding to its METIS partition
            unsigned local_nodes = mesh.GetDistributedVectorFactory()->GetLocalOwnership();
            TS_ASSERT_EQUALS(local_nodes, mesh.GetNumLocalNodes());
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


    void TestArchiving() throw(Exception)
    {
        FileFinder archive_dir("distributed_tetrahedral_mesh_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "distributed_tetrahedral_mesh.arch";
        ArchiveLocationInfo::SetMeshFilename("distributed_tetrahedral_mesh");

        DistributedTetrahedralMesh<2,2>* p_mesh = new DistributedTetrahedralMesh<2,2>(DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);
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

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
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
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
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
            if ( PetscTools::IsSequential() )
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

    void TestArchivingBinaryMesh() throw(Exception)
    {
        FileFinder archive_dir("distributed_tetrahedral_mesh_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "binary_mesh.arch";
        ArchiveLocationInfo::SetMeshFilename("binary_mesh");

        DistributedTetrahedralMesh<3,3>* p_mesh = new DistributedTetrahedralMesh<3,3>(DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);
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
            catch(Exception& e)
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
            catch(Exception& e)
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
            catch(Exception& e)
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
        /* Compare
          grep -ab "#" /tmp/chaste/testoutput/TestDistributedMeshWriter/seq_line_10_elements.node
         219:# Created by Chaste version 1.1.8525 on Wed, 31 Mar 2010 14:13:58 +0000.  Chaste was built on Wed, 31 Mar 2010 14:10:35 +0000 by machine (uname) 'Linux userpc59.comlab.ox.ac.uk 2.6.24-27-generic #1 SMP Fri Mar 12 00:52:19 UTC 2010 x86_64' using settings: default, shared libraries.
         grep -ab "#" /tmp/$USER/testoutput/TestDistributedMeshWriter/seq_line_10_elements.node
         i.e. Bytes after 219 are provenance data and will change -- one second clock difference will be spotted
         grep -ab "#" /tmp/$USER/testoutput/TestDistributedMeshWriter/seq_line_10_elements.ele
         */
        TS_ASSERT_EQUALS(system(("cmp -n 219 " + output_dir + "/par_line_10_elements.node "+ output_dir + "/seq_line_10_elements.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp -n 88 " + output_dir + "/par_line_10_elements.ele "+ output_dir + "/seq_line_10_elements.ele").c_str()), 0);
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

            unsigned received[nodes_to_receive_per_process[receive_from].size()];
            MPI_Status status;

            MPI_Recv( received,
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

//        for (unsigned process = 0; process < PetscTools::GetNumProcs(); process++)
//        {
//            PRINT_3_VARIABLES( process,
//                               nodes_to_receive_per_process[process].size(),
//                               nodes_to_send_per_process[process].size() );
//        }
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

        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + output_dir + "/par_cube_2mm_12_elements.node "+ output_dir + "/seq_cube_2mm_12_elements.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + output_dir + "/par_cube_2mm_12_elements.ele "+ output_dir + "/seq_cube_2mm_12_elements.ele").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + output_dir + "/par_cube_2mm_12_elements.face "+ output_dir + "/seq_cube_2mm_12_elements.face").c_str()), 0);
    }

    void TestEfficientParallelWriting3D()
    {
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> sequential_mesh;
        sequential_mesh.ConstructFromMeshReader(reader);
        MeshalyzerMeshWriter<3,3> mesh_writer("TestDistributedMeshWriter", "seq_cube_2mm_12_elements", false);
        mesh_writer.WriteFilesUsingMesh(sequential_mesh);
        MeshalyzerMeshWriter<3,3> mesh_writer_cg("TestDistributedMeshWriter", "seq_cube_2mm_12_elements_cg", false, true); //Don't clean, Do use CG
        mesh_writer_cg.WriteFilesUsingMesh(sequential_mesh);
        CmguiMeshWriter<3,3> cmgui_writer("TestDistributedMeshWriter", "seq_cube_2mm_12_elements_cmgui", false); //Don't clean
        cmgui_writer.WriteFilesUsingMesh(sequential_mesh);

        DistributedTetrahedralMesh<3,3> distributed_mesh(DistributedTetrahedralMeshPartitionType::DUMB); //Makes sure that there is no permutation
        AbstractTetrahedralMesh<3,3> *p_distributed_mesh = &distributed_mesh; //Hide the fact that it's distributed from the compiler

        distributed_mesh.ConstructFromMeshReader(reader);
        MeshalyzerMeshWriter<3,3> mesh_writer_par("TestDistributedMeshWriter", "par_efficient_cube_2mm_12_elements", false);
        mesh_writer_par.WriteFilesUsingMesh(*p_distributed_mesh, false);  //"false == Don't preserve element ordering

        MeshalyzerMeshWriter<3,3> mesh_writer_par_cg("TestDistributedMeshWriter", "par_efficient_cube_2mm_12_elements_cg", false, true); //Don't clean, Do use CG
        mesh_writer_par_cg.WriteFilesUsingMesh(*p_distributed_mesh, false);  //"false == Don't preserve element ordering

        CmguiMeshWriter<3,3> cmgui_writer_par("TestDistributedMeshWriter", "par_efficient_cube_2mm_12_elements_cmgui", false);
        cmgui_writer_par.WriteFilesUsingMesh(*p_distributed_mesh, false);  //"false == Don't preserve element ordering

        std::string output_dir = mesh_writer.GetOutputDirectory();

        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + output_dir + "/par_efficient_cube_2mm_12_elements.pts "+ output_dir + "/seq_cube_2mm_12_elements.pts").c_str()), 0);

        // cg output is indexed from 1, but the pts file doesn't have indices
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + output_dir + "/par_efficient_cube_2mm_12_elements_cg.pts "+ output_dir + "/seq_cube_2mm_12_elements_cg.pts").c_str()), 0);

        //cmgui
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" -I \"Group name:\" " + output_dir + "/par_efficient_cube_2mm_12_elements_cmgui.exnode "+ output_dir + "/seq_cube_2mm_12_elements_cmgui.exnode").c_str()), 0);

        // Master process sorts element and face file and the rest wait before comparing.
        if (PetscTools::AmMaster())
        {
            system(("sort " + output_dir + "/seq_cube_2mm_12_elements.tetras > " + output_dir + "seq_sorted.tetras").c_str());
            system(("sort " + output_dir + "/par_efficient_cube_2mm_12_elements.tetras > " + output_dir + "par_eff_sorted.tetras").c_str());

            system(("sort " + output_dir + "/seq_cube_2mm_12_elements.tri > " + output_dir + "seq_sorted.tri").c_str());
            system(("sort " + output_dir + "/par_efficient_cube_2mm_12_elements.tri > " + output_dir + "par_eff_sorted.tri").c_str());

            //for the cmgui, we employ some grep trickery to sort the files. We create one file per element containing the element number and the nodes that make it.
            for (unsigned elem_index = 1; elem_index<=sequential_mesh.GetNumAllElements(); elem_index++)
            {
                std::stringstream ss;
                ss << elem_index;
                std::string elem_string(ss.str());
                system(("grep -i '" + elem_string + " 0 0' " + output_dir + "par_efficient_cube_2mm_12_elements_cmgui.exelem > " + output_dir + "element_" + elem_string + "_efficient ").c_str());
                system(("grep -i '" + elem_string + " 0 0' " + output_dir + "seq_cube_2mm_12_elements_cmgui.exelem > " + output_dir + "element_" + elem_string + "_sequential ").c_str());
            }
        }
        PetscTools::Barrier();

        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + output_dir + "seq_sorted.tetras " + output_dir + "par_eff_sorted.tetras").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + output_dir + "seq_sorted.tri " + output_dir + "par_eff_sorted.tri").c_str()), 0);

        //compare teh cmgui element, one element at a time.
        for (unsigned elem_index = 1; elem_index<=sequential_mesh.GetNumAllElements(); elem_index++)
        {
            std::stringstream ss;
            ss << elem_index;
            std::string elem_string(ss.str());
            TS_ASSERT_EQUALS(system(("diff " + output_dir + "element_" + elem_string + "_efficient " + output_dir + "element_" + elem_string + "_sequential").c_str()), 0);
        }
    }

    void TestArchiveOfConstructedMesh() throw(Exception)
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

    void TestLoadBadFacesException() throw (Exception)
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
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + output_dir + "3dDistributedMesh.ncl mesh/test/data/cube_2mm_152_elements_binary_v2.ncl").c_str()), 0);
    }

    void TestCheckOutwardNormals() throw (Exception)
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
};

#endif /*TESTDISTRIBUTEDTETRAHEDRALMESH_HPP_*/
