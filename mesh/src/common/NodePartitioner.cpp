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
#include <cassert>

#include "petscao.h"

#include "NodePartitioner.hpp"
#include "PetscMatTools.hpp"
#include "PetscTools.hpp"
#include "Timer.hpp"





template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodePartitioner<ELEMENT_DIM, SPACE_DIM>::PetscMatrixPartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                              std::vector<unsigned>& rNodesPermutation,
                                              std::set<unsigned>& rNodesOwned,
                                              std::vector<unsigned>& rProcessorsOffset)
{
    assert(PetscTools::IsParallel());
    assert(ELEMENT_DIM==2 || ELEMENT_DIM==3); // Metis works with triangles and tetras

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

    if ( ! rMeshReader.IsFileFormatBinary() )
    {
        // Advance the file pointer to the first element I own
        for (unsigned element_index = 0; element_index < first_local_element; element_index++)
        {
            ElementData element_data = rMeshReader.GetNextElementData();
        }
    }

    // In the loop below we assume that there exist edges between any pair of nodes in an element. This is
    // a reasonable assumption for triangles and tetrahedra. This won't be the case for squares or hexahedra
    // (or higher order elements). We allow ELEMENT_DIM smaller than SPACE_DIM in case this is a 2D mesh in
    // a 3D space.
    assert(SPACE_DIM >= ELEMENT_DIM);

    for (unsigned element_index = 0; element_index < num_local_elements; element_index++)
    {
        ElementData element_data;

        if ( rMeshReader.IsFileFormatBinary() )
        {
            element_data = rMeshReader.GetElementData(first_local_element + element_index);
        }
        else
        {
            element_data = rMeshReader.GetNextElementData();
        }

        for (unsigned i=0; i<ELEMENT_DIM+1; i++)
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

    MatInfo matrix_info;
    MatGetInfo(connectivity_matrix, MAT_LOCAL, &matrix_info);
    unsigned local_num_nz = (unsigned) matrix_info.nz_used;

    size_t size = (num_local_nodes+1)*sizeof(PetscInt);
    void* ptr;
    PetscMalloc(size, &ptr);
    PetscInt* local_ia = (PetscInt*) ptr;
    size = local_num_nz*sizeof(PetscInt);
    PetscMalloc(size, &ptr);
    PetscInt* local_ja = (PetscInt*) ptr;

    PetscInt row_num_nz;
    const PetscInt* column_indices;

    local_ia[0]=0;
    for (PetscInt row_global_index=connectivity_matrix_lo; row_global_index<connectivity_matrix_hi; row_global_index++)
    {
        MatGetRow(connectivity_matrix, row_global_index, &row_num_nz, &column_indices, PETSC_NULL);

        unsigned row_local_index = row_global_index - connectivity_matrix_lo;
        local_ia[row_local_index+1] = local_ia[row_local_index] + row_num_nz;
        for (PetscInt col_index=0; col_index<row_num_nz; col_index++)
        {
           local_ja[local_ia[row_local_index] + col_index] =  column_indices[col_index];
        }

        MatRestoreRow(connectivity_matrix, row_global_index, &row_num_nz,&column_indices, PETSC_NULL);
    }

    PetscTools::Destroy(connectivity_matrix);

    // Convert to an adjacency matrix
    Mat adj_matrix;
    MatCreateMPIAdj(PETSC_COMM_WORLD, num_local_nodes, num_nodes, local_ia, local_ja, PETSC_NULL, &adj_matrix);

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

    /// It seems to be free-ing local_ia and local_ja as a side effect
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
    ///\todo do something smarter with iterators...
    for (unsigned i=0; i<my_num_nodes; i++)
    {
        rNodesOwned.insert(local_range[i]);
    }
    delete[] local_range;

    // Once we know the offsets we can compute the permutation vector
    PetscInt* global_range = new PetscInt[num_nodes];
    for (unsigned i=0; i<num_nodes; i++)
    {
        global_range[i] = i;
    }
    AOPetscToApplication(ordering, num_nodes, global_range);

    rNodesPermutation.resize(num_nodes);

    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        rNodesPermutation[node_index] = global_range[node_index];
    }
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
////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class NodePartitioner<1,1>;
template class NodePartitioner<1,2>;
template class NodePartitioner<1,3>;
template class NodePartitioner<2,2>;
template class NodePartitioner<2,3>;
template class NodePartitioner<3,3>;




