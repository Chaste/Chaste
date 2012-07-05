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

#ifndef TESTDISTRIBUTEDQUADRATICMESH_HPP_
#define TESTDISTRIBUTEDQUADRATICMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <set>
#include <vector>

#include "DistributedQuadraticMesh.hpp"
#include "NodePartitioner.hpp"
#include "QuadraticMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscTools.hpp"

#include "PetscSetupAndFinalize.hpp"


class TestDistributedQuadraticMesh : public CxxTest::TestSuite
{
public:
    void TestDumbMeshPartitioning() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1, false);
        QuadraticMesh<2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        std::set<unsigned> nodes_owned;

        NodePartitioner<2, 2>::DumbPartitioning(mesh, nodes_owned);

        if(PetscTools::GetNumProcs() == 1)
        {
            TS_ASSERT_EQUALS(nodes_owned.size(), 289u);
        }
        else if(PetscTools::GetNumProcs() == 2)
        {
            if(PetscTools::GetMyRank() == 0 )
            {
                TS_ASSERT_EQUALS(nodes_owned.size(), 145u);
            }
            else //PetscTools::GetMyRank() == 1
            {
                TS_ASSERT_EQUALS(nodes_owned.size(), 144u);
            }
        }
        else if(PetscTools::GetNumProcs() == 3)
        {
            if(PetscTools::GetMyRank() == 0 )
            {
                TS_ASSERT_EQUALS(nodes_owned.size(), 97u);
            }
            else if(PetscTools::GetMyRank() == 1 )
            {
                TS_ASSERT_EQUALS(nodes_owned.size(), 96u);
            }
            else //PetscTools::GetMyRank() == 2
            {
                TS_ASSERT_EQUALS(nodes_owned.size(), 96u);
            }
        }
    }

    void TestMetisMeshPartitioning() throw (Exception)
    {
        EXIT_IF_SEQUENTIAL //Doesn't make sense to try and partition in sequential

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1, false);
        std::vector<unsigned> nodes_permutation;
        std::set<unsigned> nodes_owned;
        std::vector<unsigned> processor_offset;

        typedef NodePartitioner<2, 2> Partitioner2D;

        TS_ASSERT_THROWS_THIS(Partitioner2D::MetisLibraryPartitioning(mesh_reader, nodes_permutation, nodes_owned, processor_offset),
                              "Metis cannot partition a quadratic mesh.");
    }

    void TestPetscMatrixPartitioning() throw (Exception)
    {
        EXIT_IF_SEQUENTIAL //Doesn't make sense to try and partition in sequential

        if(!PetscTools::HasParMetis())
        {
            std::cout << "\n\nWarning: PETSc support for ParMetis is not installed. Mesh partitioning not tested." << std::endl;
            return;
        }

        //// Note: if this test is run with square_128_elements_quadratic and num_procs > 2, the
        //// partition returned is junk. Have verified that this is just due to the mesh being
        //// very coarse - if the following is used to create a finer mesh, the partition looks
        //// fine for num_procs up to 6 (at least).
        //QuadraticMesh<2> mesh0(0.01,1.0,1.0);
        //TrianglesMeshWriter<2,2> writer("", "quad_mesh");
        //writer.WriteFilesUsingMesh(mesh0);
        //TrianglesMeshReader<2,2> mesh_reader("../../../../tmp/rafb/testoutput/quad_mesh",2,1,false);

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1, false);

        std::vector<unsigned> nodes_permutation;
        std::set<unsigned> nodes_owned;
        std::vector<unsigned> processor_offset;

        typedef NodePartitioner<2, 2> Partitioner2D;
        Partitioner2D::PetscMatrixPartitioning(mesh_reader, nodes_permutation, nodes_owned, processor_offset);

        if(PetscTools::GetNumProcs() != 2)
        {
            return;
        }

        QuadraticMesh<2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Compute the centre of mass (average position of all the nodes) for each processor,
        // and check that all the nodes are within a certain distance of the centre of mass.
        c_vector<double,2> centre_of_mass = zero_vector<double>(2);
        unsigned counter = 0;
        for(std::set<unsigned>::iterator iter = nodes_owned.begin();
            iter != nodes_owned.end();
            ++iter)
        {
            double x = mesh.GetNode(*iter)->rGetLocation()[0];
            double y = mesh.GetNode(*iter)->rGetLocation()[1];

            centre_of_mass(0) += x;
            centre_of_mass(1) += y;

            counter++;
        }
        centre_of_mass(0) /= counter;
        centre_of_mass(1) /= counter;

        for(std::set<unsigned>::iterator iter = nodes_owned.begin();
                    iter != nodes_owned.end();
                    ++iter)
        {
            double dx = mesh.GetNode(*iter)->rGetLocation()[0] - centre_of_mass(0);
            double dy = mesh.GetNode(*iter)->rGetLocation()[1] - centre_of_mass(1);

            double dist_to_centre_of_mass = sqrt(dx*dx + dy*dy);
            TS_ASSERT_LESS_THAN(dist_to_centre_of_mass, 0.66);
        }

        //// For visualising the partition - see the matlab / octave code below
        //OutputFileHandler handler("TestDistributedQuadMeshPartitioning");
        //std::stringstream ss;
        //ss << "res_" << PetscTools::GetNumProcs() << "_" << PetscTools::GetMyRank() << ".txt";
        //out_stream p_file = handler.OpenOutputFile(ss.str());
        //
        //for(std::set<unsigned>::iterator iter = nodes_owned.begin();
        //    iter != nodes_owned.end();
        //    ++iter)
        //{
        //    *p_file << mesh.GetNode(*iter)->rGetLocation()[0] << " "
        //              << mesh.GetNode(*iter)->rGetLocation()[1] << std::endl;
        //}
        //p_file->close();

        /*
        function viz_partition(N)
        col = {'*', 'r*', 'k*', 'm*', 'g*', 'y*'};
        figure; hold on;
        for i=0:N-1
          file = ['/tmp/rafb/testoutput/TestDistributedQuadMeshPartitioning/res_',num2str(N),'_',num2str(i),'.txt'];
          d = load(file);
          plot(d(:,1),d(:,2),col{i+1});
        end;
        */
    }


    void TestConstructFromMeshReader2D() throw (Exception)
    {
        /*
         * Note that these mesh files have
         *  - quadratic elements
         *  - linear edges
         *  - edge file doesn't say which element the edge belongs too
         */
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic", 2, 1, false);
        DistributedQuadraticMesh<2> mesh(DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 289u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 128u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 32u);

        /**
         * Currently we get
         *         NumLocalElements [min,max]
         *         r16019(nnodes=3) r16020(nnodes=2)
         *  Procs  PARMETIS_LIBRARY PARMETIS_LIBRARY   PETSC_MAT_PARTITION DUMB
         * -----   ---------------- ----------          ------------------- ----
         * 1       [128, 128]                                               (one process gets 128 up until 7 procs)
         * 2       [106, 109]       [66, 83]            [72, 72]            [93, 128]
         * 3       [86,  92]        [49, 65]            [77, 95]            [63, 128]
         * 4       [60,  94]        [39, 47]            [46, 74]            [49, 128]
         * 5       [50,  68]        [30, 52]            [35, 64]            [39, 128]
         * 
         */
        // Check that it is not a dumb partition.
        // (Dumb partitions with few processes require ownership of all the mesh by at least one process
        if (PetscTools::IsParallel())
        {
            TS_ASSERT_LESS_THAN(mesh.GetNumLocalElements(), mesh.GetNumElements());
        }


        TS_ASSERT_EQUALS(mesh.GetDistributedVectorFactory()->GetProblemSize(), 289u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 32u);

        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            TS_ASSERT(iter->GetOwnership());
        }

        //Compare the DistributedQuadraticMesh to a normal QuadraticMesh
        const std::vector<unsigned>& r_node_perm = mesh.rGetNodePermutation();
        QuadraticMesh<2> seq_mesh;
        seq_mesh.ConstructFromMeshReader(mesh_reader);

        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
            iter != mesh.GetElementIteratorEnd();
            ++iter)
        {
            unsigned element_index = iter->GetIndex();

            Element<2,2>* p_sequ_element = seq_mesh.GetElement(element_index);
            TS_ASSERT_EQUALS(element_index, p_sequ_element->GetIndex());
            TS_ASSERT_EQUALS(iter->GetNumNodes(), p_sequ_element->GetNumNodes());

            for (unsigned node_local_index=0; node_local_index < iter->GetNumNodes(); node_local_index++)
            {
                unsigned node_global_index = p_sequ_element->GetNodeGlobalIndex(node_local_index);
                if (!r_node_perm.empty())
                {
                    node_global_index = r_node_perm[node_global_index];
                }
                TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(node_local_index), node_global_index);

                TS_ASSERT_EQUALS(iter->GetNode(node_local_index)->GetPoint()[0],
                                 p_sequ_element->GetNode(node_local_index)->GetPoint()[0]);
            }
        }

        for (AbstractTetrahedralMesh<2,2>::BoundaryElementIterator it=mesh.GetBoundaryElementIteratorBegin();
            it!=mesh.GetBoundaryElementIteratorEnd();
            ++it)
        {
            BoundaryElement<1,2>* p_para_boundary_element = *it;
            unsigned boundary_element_index = p_para_boundary_element->GetIndex();

            BoundaryElement<1,2>* p_sequ_boundary_element = seq_mesh.GetBoundaryElement(boundary_element_index);
            TS_ASSERT_EQUALS(boundary_element_index, p_sequ_boundary_element->GetIndex());
            TS_ASSERT_EQUALS(p_para_boundary_element->GetNumNodes(), p_sequ_boundary_element->GetNumNodes());
            TS_ASSERT_EQUALS(p_para_boundary_element->GetNumNodes(), 3u); //Quadratic edge

            for (unsigned node_local_index=0; node_local_index < p_para_boundary_element->GetNumNodes(); node_local_index++)
            {
                unsigned node_global_index = p_sequ_boundary_element->GetNodeGlobalIndex(node_local_index);
                if (!r_node_perm.empty())
                {
                    node_global_index = r_node_perm[node_global_index];
                }
                TS_ASSERT_EQUALS(p_para_boundary_element->GetNodeGlobalIndex(node_local_index), node_global_index);

                TS_ASSERT_EQUALS(p_para_boundary_element->GetNode(node_local_index)->GetPoint()[0],
                                 p_sequ_boundary_element->GetNode(node_local_index)->GetPoint()[0]);
            }
        }
    }
    void TestConstructFromMeshReader2DWithPetscSupport() throw (Exception)
    {
        if (!PetscTools::HasParMetis())
        {
            std::cout << "\n\nWarning: PETSc support for ParMetis is not installed. Mesh partitioning not tested." << std::endl;
            return;
        }

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic", 2, 1, false);
        DistributedQuadraticMesh<2> mesh(DistributedTetrahedralMeshPartitionType::PETSC_MAT_PARTITION);

        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 289u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 128u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 152u);

        // Check that it is not a dumb partition.
        // (Dumb partitions with few processes require ownership of all the mesh by at least one process
        if (PetscTools::IsParallel())
        {
            TS_ASSERT_LESS_THAN(mesh.GetNumLocalElements(), mesh.GetNumElements());
        }
    }

    void TestConstructFromMeshReader3DContainingElementInFile() throw (Exception)
    {
        /*
         * Note that these mesh files have
         *  - quadratic elements
         *  - linear edges
         *  - face file does say which element the face belongs too
         */
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_152_elements_v3", 2, 1, true);
        DistributedQuadraticMesh<3> mesh(DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 335u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 152u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 116u);

        // Check that it is not a dumb partition.
        // (Dumb partitions with few processes require ownership of all the mesh by at least one process
        if (PetscTools::IsParallel())
        {
            TS_ASSERT_LESS_THAN(mesh.GetNumLocalElements(), mesh.GetNumElements());
        }

        BoundaryElement<2,3>* p_face = *(mesh.GetBoundaryElementIteratorBegin());
        TS_ASSERT_EQUALS(p_face->GetNumNodes(), 6u);
    }

};


#endif // TESTDISTRIBUTEDQUADRATICMESH_HPP_
