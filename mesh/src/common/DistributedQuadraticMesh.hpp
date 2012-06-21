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
        archive & boost::serialization::base_object<AbstractTetrahedralMesh<DIM, DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * \todo Change default partitioningMethod to PETSC_MAT_PARTITION when this class is complete
     *
     * @param partitioningMethod  defaults to DUMB, Nb: This is in contrast to the default behaviour of DistributedTetrahedralMesh, which defaults to METIS_LIBRARY
     */
    DistributedQuadraticMesh(DistributedTetrahedralMeshPartitionType::type partitioningMethod=DistributedTetrahedralMeshPartitionType::DUMB);

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

private:

    /**
     * Top level method for making 2D edges have 3 nodes not 2 and making 3D faces have 6 nodes not 3  (ie linear to quadratic).
     * @param pMeshReader Pointer to the reader. Only used if boundaryElemFileHasContainElementInfo==true (can be null if not).
     */
    void AddNodesToBoundaryElements(TrianglesMeshReader<DIM,DIM>* pMeshReader);

    /**
     * This method adds the given node (defined by an element and a node index)
     * to the given boundary element, and also sets the node as a boundary
     * element and adds it to the std::vector of boundary elements.
     *
     * @param pBoundaryElement  pointer to a boundary element in the mesh
     * @param pElement  pointer to an element in the mesh
     * @param internalNode  index of a node in the mesh
     */
    void AddNodeToBoundaryElement(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                  Element<DIM,DIM>* pElement,
                                  unsigned internalNode);

    /**
     * Given a face in an element (defined by giving an element and the opposite
     * node number to the face) that corresponds to a given boundary element,
     * this method adds in the face's internal nodes to the boundary element
     * (in the correct order).
     *
     * @param pBoundaryElement  pointer to a boundary element in the mesh
     * @param pElement  pointer to an element in the mesh
     * @param nodeIndexOppositeToFace  index of a node in the mesh
     */
    void AddExtraBoundaryNodes(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                               Element<DIM,DIM>* pElement,
                               unsigned nodeIndexOppositeToFace);

    /**
     * Nasty helper method for AddNodeToBoundaryElement() in 3D.
     *
     * This method takes in the three vertices of a face which match the given boundary
     * element, and figure out if the order of the nodes in the face is reversed in
     * the boundary element (returned in the bool 'rReverse'). Also, the offset between
     * the first node in the face (as given to this method) and the first node in
     * the boundary element is computed (returned in the variable 'rOffset'). Offset
     * should then be applied before reverse to match the face nodes to the boundary
     * element nodes.
     *
     * \todo document these parameters
     *
     * @param boundaryElemNode0
     * @param boundaryElemNode1
     * @param pElement
     * @param node0
     * @param node1
     * @param node2
     * @param rOffset
     * @param rReverse
     */
    void HelperMethod1(unsigned boundaryElemNode0, unsigned boundaryElemNode1,
                       Element<DIM,DIM>* pElement,
                       unsigned node0, unsigned node1, unsigned node2,
                       unsigned& rOffset,
                       bool& rReverse);

    /**
     * Nasty helper method for AddNodeToBoundaryElement() in 3D.
     *
     * This method takes the three internal nodes for some face in some element,
     * applies the given offset and reverse (see HelperMethod1) to them, to get
     * the ordered internal nodes which should given to the boundary element.
     * It then calls AddNodeToBoundaryElement with each of the three internal nodes.
     *
     * \todo document these parameters
     *
     * @param pBoundaryElement
     * @param pElement
     * @param internalNode0
     * @param internalNode1
     * @param internalNode2
     * @param offset
     * @param reverse
     */
    void HelperMethod2(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                       Element<DIM,DIM>* pElement,
                       unsigned internalNode0, unsigned internalNode1, unsigned internalNode2,
                       unsigned offset,
                       bool reverse);



protected:

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
    Archive & ar, const DistributedQuadraticMesh<DIM> * t, const BOOST_PFTO unsigned int file_version)
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
        EXCEPTION("This archive was written for a different number of processors");
    }

}
}
} // namespace ...

#endif /*DISTRIBUTEDQUADRATICMESH_HPP_*/
