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
#include <cassert>
#include <algorithm>

#include "Exception.hpp"
#include "NodePartitioner.hpp"
#include "PetscMatTools.hpp"
#include "PetscTools.hpp"
#include "Timer.hpp"
#include "TrianglesMeshReader.hpp"
#include "Warnings.hpp"
#include "petscao.h"
#include "petscmat.h"

#include <parmetis.h>

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodePartitioner<ELEMENT_DIM, SPACE_DIM>::DumbPartitioning(AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                                               std::set<unsigned>& rNodesOwned)
{
    //Note that if there is not DistributedVectorFactory in the mesh, then a dumb partition based
    //on rMesh.GetNumNodes() is applied automatically.
    //If there is a DistributedVectorFactory then that one will be returned
    if (rMesh.GetDistributedVectorFactory()->GetProblemSize() != rMesh.GetNumNodes())
    {
        EXCEPTION("The distributed vector factory size in the mesh doesn't match the total number of nodes.");
    }

    for (unsigned node_index = rMesh.GetDistributedVectorFactory()->GetLow();
         node_index < rMesh.GetDistributedVectorFactory()->GetHigh();
         node_index++)
    {
         rNodesOwned.insert(node_index);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodePartitioner<ELEMENT_DIM, SPACE_DIM>::PetscMatrixPartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                              std::vector<unsigned>& rNodePermutation,
                                              std::set<unsigned>& rNodesOwned,
                                              std::vector<unsigned>& rProcessorsOffset)
{
    assert(PetscTools::IsParallel());
    assert(ELEMENT_DIM==2 || ELEMENT_DIM==3);      // LCOV_EXCL_LINE // Metis works with triangles and tetras

    if (!PetscTools::HasParMetis()) //We must have ParMetis support compiled into Petsc
    {
        WARN_ONCE_ONLY("PETSc support for ParMetis is not installed.");      // LCOV_EXCL_LINE
    }

    unsigned num_nodes = rMeshReader.GetNumNodes();
    PetscInt num_procs = PetscTools::GetNumProcs();
    unsigned local_proc_index = PetscTools::GetMyRank();

    unsigned num_elements = rMeshReader.GetNumElements();
    unsigned num_local_elements = num_elements / num_procs;
    unsigned first_local_element = num_local_elements * local_proc_index;
    if (PetscTools::AmTopMost())
    {
        // Take the excess elements
        num_local_elements += num_elements - (num_local_elements * num_procs);
    }

    PetscTools::Barrier();
    Timer::Reset();

    /*
     * Create PETSc matrix which will have 1 for adjacent nodes.
     */
    Mat connectivity_matrix;
    ///\todo #1216 change the number 54 below (row nonzero allocation) to be nonmagic
    PetscTools::SetupMat(connectivity_matrix, num_nodes, num_nodes, 54, PETSC_DECIDE, PETSC_DECIDE, false);

    if (!rMeshReader.IsFileFormatBinary())
    {
        // Advance the file pointer to the first element I own
        for (unsigned element_index = 0; element_index < first_local_element; element_index++)
        {
            rMeshReader.GetNextElementData();
        }
    }

    // In the loop below we assume that there exist edges between any pair of nodes in an element. This is
    // a reasonable assumption for triangles and tetrahedra. This won't be the case for quadratic simplices,
    // squares or hexahedra (or higher order elements), leading to slightly suboptimal partitions in these
    // cases.
    // We allow ELEMENT_DIM smaller than SPACE_DIM in case this is a 2D mesh in
    // a 3D space.
    assert(SPACE_DIM >= ELEMENT_DIM);

    for (unsigned element_index = 0; element_index < num_local_elements; element_index++)
    {
        ElementData element_data;

        if (rMeshReader.IsFileFormatBinary())
        {
            element_data = rMeshReader.GetElementData(first_local_element + element_index);
        }
        else
        {
            element_data = rMeshReader.GetNextElementData();
        }

        for (unsigned i=0; i<element_data.NodeIndices.size(); i++)
        {
            unsigned row = element_data.NodeIndices[i];
            for (unsigned j=0; j<i; j++)
            {
                unsigned col = element_data.NodeIndices[j];
                MatSetValue(connectivity_matrix, row, col, 1.0, INSERT_VALUES);
                MatSetValue(connectivity_matrix, col, row, 1.0, INSERT_VALUES);
            }
        }
    }
    /// \todo: This assembly is likely to generate many communications. Try to interleave other operations by executing them between Begin() and End().
    PetscMatTools::Finalise(connectivity_matrix);

    PetscTools::Barrier();
    if (PetscTools::AmMaster())
    {
        Timer::PrintAndReset("Connectivity matrix assembly");
    }

    rMeshReader.Reset();

    PetscInt connectivity_matrix_lo;
    PetscInt connectivity_matrix_hi;
    MatGetOwnershipRange(connectivity_matrix, &connectivity_matrix_lo, &connectivity_matrix_hi);

    unsigned num_local_nodes = connectivity_matrix_hi - connectivity_matrix_lo;

    /* PETSc MatCreateMPIAdj and parMETIS rely on adjacency arrays
     * Named here:      xadj            adjncy
     * parMETIS name:   xadj            adjncy    (see S4.2 in parMETIS manual)
     * PETSc name:      i               j         (see MatCreateMPIAdj in PETSc manual)
     *
     * The adjacency information of all nodes is listed in the main array, adjncy.  Since each node
     * has a variable number of adjacent nodes, the array xadj is used to store the index (in adjncy) where
     * this information starts.  Since xadj[i] is the start of node i's information, xadj[i+1] marks the end.
     *
     *
     */
    MatInfo matrix_info;
    MatGetInfo(connectivity_matrix, MAT_LOCAL, &matrix_info);
    unsigned local_num_nz = (unsigned) matrix_info.nz_used;

    size_t size = (num_local_nodes+1)*sizeof(PetscInt);
    void* ptr;
    PetscMalloc(size, &ptr);
    PetscInt* xadj = (PetscInt*) ptr;
    size = local_num_nz*sizeof(PetscInt);
    PetscMalloc(size, &ptr);
    PetscInt* adjncy = (PetscInt*) ptr;

    PetscInt row_num_nz;
    const PetscInt* column_indices;

    xadj[0]=0;
    for (PetscInt row_global_index=connectivity_matrix_lo; row_global_index<connectivity_matrix_hi; row_global_index++)
    {
        MatGetRow(connectivity_matrix, row_global_index, &row_num_nz, &column_indices, PETSC_NULL);

        unsigned row_local_index = row_global_index - connectivity_matrix_lo;
        xadj[row_local_index+1] = xadj[row_local_index] + row_num_nz;
        for (PetscInt col_index=0; col_index<row_num_nz; col_index++)
        {
           adjncy[xadj[row_local_index] + col_index] =  column_indices[col_index];
        }

        MatRestoreRow(connectivity_matrix, row_global_index, &row_num_nz,&column_indices, PETSC_NULL);
    }

    PetscTools::Destroy(connectivity_matrix);

    // Convert to an adjacency matrix
    Mat adj_matrix;
    MatCreateMPIAdj(PETSC_COMM_WORLD, num_local_nodes, num_nodes, xadj, adjncy, PETSC_NULL, &adj_matrix);

    PetscTools::Barrier();
    if (PetscTools::AmMaster())
    {
        Timer::PrintAndReset("Adjacency matrix creation");
    }

    // Get PETSc to call ParMETIS
    MatPartitioning part;
    MatPartitioningCreate(PETSC_COMM_WORLD, &part);
    MatPartitioningSetAdjacency(part, adj_matrix);
    MatPartitioningSetFromOptions(part);
    IS new_process_numbers;

    MatPartitioningApply(part, &new_process_numbers);
    MatPartitioningDestroy(PETSC_DESTROY_PARAM(part));

    /// It seems to be free-ing xadj and adjncy as a side effect
    PetscTools::Destroy(adj_matrix);

    PetscTools::Barrier();
    if (PetscTools::AmMaster())
    {
        Timer::PrintAndReset("PETSc-ParMETIS call");
    }

    // Figure out who owns what - processor offsets
    PetscInt* num_nodes_per_process = new PetscInt[num_procs];
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
    ISPartitioningCount(new_process_numbers, num_procs, num_nodes_per_process);
#else
    ISPartitioningCount(new_process_numbers, num_nodes_per_process);
#endif

    rProcessorsOffset.resize(num_procs);
    rProcessorsOffset[0] = 0;
    for (PetscInt i=1; i<num_procs; i++)
    {
        rProcessorsOffset[i] = rProcessorsOffset[i-1] + num_nodes_per_process[i-1];
    }
    unsigned my_num_nodes = num_nodes_per_process[local_proc_index];
    delete[] num_nodes_per_process;

    // Figure out who owns what - new node numbering
    IS new_global_node_indices;
    ISPartitioningToNumbering(new_process_numbers, &new_global_node_indices);

    // Index sets only give local information, we want global
    AO ordering;
    AOCreateBasicIS(new_global_node_indices, PETSC_NULL /* natural ordering */, &ordering);

    // The locally owned range under the new numbering
    PetscInt* local_range = new PetscInt[my_num_nodes];
    for (unsigned i=0; i<my_num_nodes; i++)
    {
        local_range[i] = rProcessorsOffset[local_proc_index] + i;
    }
    AOApplicationToPetsc(ordering, my_num_nodes, local_range);
    //AOView(ordering, PETSC_VIEWER_STDOUT_WORLD);

    // Fill in rNodesOwned
    rNodesOwned.insert(local_range, local_range + my_num_nodes);
    delete[] local_range;

    // Once we know the offsets we can compute the permutation vector
    PetscInt* global_range = new PetscInt[num_nodes];
    for (unsigned i=0; i<num_nodes; i++)
    {
        global_range[i] = i;
    }
    AOPetscToApplication(ordering, num_nodes, global_range);

    rNodePermutation.resize(num_nodes);
    std::copy(global_range, global_range+num_nodes, rNodePermutation.begin());
    delete[] global_range;

    AODestroy(PETSC_DESTROY_PARAM(ordering));
    ISDestroy(PETSC_DESTROY_PARAM(new_process_numbers));
    ISDestroy(PETSC_DESTROY_PARAM(new_global_node_indices));

    PetscTools::Barrier();
    if (PetscTools::AmMaster())
    {
        Timer::PrintAndReset("PETSc-ParMETIS output manipulation");
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodePartitioner<ELEMENT_DIM, SPACE_DIM>::GeometricPartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                    std::vector<unsigned>& rNodePermutation,
                                    std::set<unsigned>& rNodesOwned,
                                    std::vector<unsigned>& rProcessorsOffset,
                                    ChasteCuboid<SPACE_DIM>* pRegion)
{
    PetscTools::Barrier();

    unsigned num_nodes =  rMeshReader.GetNumNodes();

    // Work out where each node should lie.
    std::vector<unsigned> node_ownership(num_nodes, 0);
    std::vector<unsigned> node_index_ownership(num_nodes, UNSIGNED_UNSET);

    for (unsigned node=0; node < num_nodes; node++)
    {
        std::vector<double> location = rMeshReader.GetNextNode();

        // Make sure it is the correct size
        assert(location.size() == SPACE_DIM);

        // Establish whether it lies in the domain. ChasteCuboid::DoesContain is
        // insufficient for this as it treats all boundaries as open.
        const ChastePoint<SPACE_DIM>& lower = pRegion->rGetLowerCorner();
        const ChastePoint<SPACE_DIM>& upper = pRegion->rGetUpperCorner();

        bool does_contain = true;

        for (unsigned d=0; d<SPACE_DIM; d++)
        {
            bool boundary_check;
            boundary_check = ((location[d] > lower[d]) || sqrt((location[d]-lower[d])*(location[d]-lower[d])) < DBL_EPSILON);
            does_contain = (does_contain && boundary_check && (location[d] < upper[d]));
        }

        if (does_contain)
        {
            node_ownership[node] = 1;
            rNodesOwned.insert(node);
            node_index_ownership[node] = PetscTools::GetMyRank();
        }
    }

    // Make sure each node will be owned by exactly on process
    std::vector<unsigned> global_ownership(num_nodes, 0);
    std::vector<unsigned> global_index_ownership(num_nodes, UNSIGNED_UNSET);

    MPI_Allreduce(&node_ownership[0], &global_ownership[0], num_nodes, MPI_UNSIGNED, MPI_LXOR, PETSC_COMM_WORLD);
    for (unsigned i=0; i<num_nodes; i++)
    {
        if (global_ownership[i] == 0)
        {
            EXCEPTION("A node is either not in geometric region, or the regions are not disjoint.");
        }
    }

    // Create the permutation and offset vectors.
    MPI_Allreduce(&node_index_ownership[0], &global_index_ownership[0], num_nodes, MPI_UNSIGNED, MPI_MIN, PETSC_COMM_WORLD);

    for (unsigned proc=0; proc<PetscTools::GetNumProcs(); proc++)
    {
        rProcessorsOffset.push_back(rNodePermutation.size());
        for (unsigned node=0; node<num_nodes; node++)
        {
            if (global_index_ownership[node] == proc)
            {
                rNodePermutation.push_back(node);
            }
        }
    }
    assert(rNodePermutation.size() == num_nodes);
}

// Explicit instantiation
template class NodePartitioner<1,1>;
template class NodePartitioner<1,2>;
template class NodePartitioner<1,3>;
template class NodePartitioner<2,2>;
template class NodePartitioner<2,3>;
template class NodePartitioner<3,3>;
