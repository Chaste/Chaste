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

#ifndef DISTRIBUTEDQUADRATICMESH_HPP_
#define DISTRIBUTEDQUADRATICMESH_HPP_

#include <map>
#include <vector>
#include <set>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "DistributedTetrahedralMesh.hpp"
#include "Node.hpp"
#include "AbstractMeshReader.hpp"
#include "DistributedTetrahedralMeshPartitionType.hpp"
#include "QuadraticMeshHelper.hpp"

#define UNASSIGNED_NODE UINT_MAX

#include <parmetis.h>

/**
 * Parallel implementation of a quadratic mesh
 * Nodes are distributed such that each process has
 * - A set of nodes (possibly reordered) with contiguous global indices
 * - A local copy of all the elements supporting those nodes
 * - A local copy of ghost/halo nodes which are all the nodes used in the supporting elements, but not owned outright.
 */
template<unsigned DIM>
class DistributedQuadraticMesh : public DistributedTetrahedralMesh<DIM, DIM>
{
    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the mesh.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<DistributedTetrahedralMesh<DIM, DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param partitioningMethod  defaults to PARMETIS_LIBRARY, Nb: This is should have the same default behaviour as DistributedTetrahedralMesh
     */
    DistributedQuadraticMesh(DistributedTetrahedralMeshPartitionType::type partitioningMethod=DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);

    /**
     * Destructor.
     */
    virtual ~DistributedQuadraticMesh();

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<DIM,DIM>& rMeshReader);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DistributedQuadraticMesh)

namespace boost
{
namespace serialization
{
/**
 * Record number of processors when saving...
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const DistributedQuadraticMesh<DIM> * t, const unsigned int file_version)
{
    unsigned num_procs = PetscTools::GetNumProcs();
    const DistributedTetrahedralMeshPartitionType::type partition_type = t->GetPartitionType();
    ar << num_procs;
    ar << partition_type;
}

/**
 * De-serialize constructor parameters and initialise a DistributedQuadraticMesh,
 * checking the number of processors is the same.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, DistributedQuadraticMesh<DIM> * t, const unsigned int file_version)
{
    unsigned num_procs;
    DistributedTetrahedralMeshPartitionType::type partition_type;

    ar >> num_procs;
    ar >> partition_type;

    // Invoke inplace constructor to initialise instance
    /// \todo #1199  Lots of stuff can't cope if we re-partition
    //::new(t)DistributedTetrahedralMesh<DIM, DIM>(partition_type);
    ::new(t)DistributedQuadraticMesh<DIM>(DistributedTetrahedralMeshPartitionType::DUMB);

    /*
     * The exception needs to be thrown after the call to ::new(t), or Boost will try
     * to free non-allocated memory when the exception is thrown.
     */
    if (DistributedVectorFactory::CheckNumberOfProcessesOnLoad() &&
        num_procs != PetscTools::GetNumProcs())
    {
        NEVER_REACHED;
    }

}
}
} // namespace ...

#endif /*DISTRIBUTEDQUADRATICMESH_HPP_*/
