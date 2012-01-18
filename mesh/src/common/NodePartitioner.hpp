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

#ifndef NODEPARTITIONER_HPP_
#define NODEPARTITIONER_HPP_

#include <set>

#include "AbstractMeshReader.hpp"

/**
 * Static methods to allow node-wise partitioning of meshes.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class NodePartitioner
{
public:
    /**
     * Method to compute a parallel partitioning of a given mesh.
     *
     * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
     * @param rNodesPermutation is the vector to be filled with node permutation information.
     * @param rNodesOwned is an empty set to be filled with the indices of nodes owned by this process
     * @param rProcessorsOffset a vector of length NumProcs to be filled with the index of the lowest indexed node owned by each process
     *
     */
    static void PetscMatrixPartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                        std::vector<unsigned>& rNodesPermutation,
                                        std::set<unsigned>& rNodesOwned,
                                        std::vector<unsigned>& rProcessorsOffset);

private:
};

#endif // NODEPARTITIONER_HPP_
