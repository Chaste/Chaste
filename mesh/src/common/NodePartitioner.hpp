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

#ifndef NODEPARTITIONER_HPP_
#define NODEPARTITIONER_HPP_

#include <set>

#include "AbstractMesh.hpp"
#include "AbstractMeshReader.hpp"

/**
 * Static methods to allow node-wise partitioning of meshes.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class NodePartitioner
{
public:

    /**
      * Specialised method to compute a parallel partitioning of a given mesh
      * (called by ComputeMeshPartitioning, based on the value of mMetisPartitioning)
      *
      * @param rMesh is the original mesh (so that we can set the DistributedVectorFactory up
      * @param rNodesOwned is an empty set to be filled with the indices of nodes owned by this process
      */
     static void DumbPartitioning(AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                  std::set<unsigned>& rNodesOwned);

    /**
     * Method to compute a parallel partitioning of a given mesh.
     *
     * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
     * @param rNodePermutation is the vector to be filled with node permutation information.
     * @param rNodesOwned is an empty set to be filled with the indices of nodes owned by this process
     * @param rProcessorsOffset a vector of length NumProcs to be filled with the index of the lowest indexed node owned by each process
     *
     */
    static void PetscMatrixPartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                        std::vector<unsigned>& rNodePermutation,
                                        std::set<unsigned>& rNodesOwned,
                                        std::vector<unsigned>& rProcessorsOffset);
    /**
     * Specialised method to compute the partition of a mesh based on geometric partitioning
     *
     * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
     * @param rNodePermutation is the vector to be filled with node permutation information.
     * @param rNodesOwned is an empty set to be filled with the indices of nodes owned by this process
     * @param rProcessorsOffset a vector of length NumProcs to be filled with the index of the lowest indexed node owned by each process
     * @param pRegion the local region owned by this process.
     *
     */
    static void GeometricPartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                        std::vector<unsigned>& rNodePermutation,
                                        std::set<unsigned>& rNodesOwned,
                                        std::vector<unsigned>& rProcessorsOffset,
                                        ChasteCuboid<SPACE_DIM>* pRegion);


private:
};

#endif // NODEPARTITIONER_HPP_
