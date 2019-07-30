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

#ifndef MIXEDDIMENSIONMESH_HPP_
#define MIXEDDIMENSIONMESH_HPP_

#include "DistributedTetrahedralMesh.hpp"
#include "AbstractMeshReader.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A tetrahedral mesh that also supports embedded 1D cable elements.
 *
 * Could be used for Purkinje or blood vessels, etc.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MixedDimensionMesh : public DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * Constructor.
     *
     * @param partitioningMethod  defaults to METIS_LIBRARY, but in 1-D is always overridden in this constructor to be the DUMB partition
     */
    MixedDimensionMesh(DistributedTetrahedralMeshPartitionType::type partitioningMethod=DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);

    /**
     * Destructor - cleans up the cables
     *
     */
    ~MixedDimensionMesh();

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader);

   /**
     * Add the most recently constructed cable element to the global->local cable element mapping
     *
     * @param index is the global index of cable element to be registered
     */
    void RegisterCableElement(unsigned index);

    /**
     * @return the number of cable elements.
     */
    unsigned GetNumCableElements() const;

    /**
     * @return the number of cable elements on this process.
     */
    unsigned GetNumLocalCableElements() const;

    /**
     * Get the cable element with a given index in the mesh.
     *
     * @param globalElementIndex the global index of the cable element
     * @return a pointer to the cable element.
     */
    Element<1u, SPACE_DIM>* GetCableElement(unsigned globalElementIndex) const;

    /**
     * Determine whether or not the current process owns node 0 of this cable element (tie breaker to determine which process writes
     * to file for when two or more share ownership of a cable element).
     * @return true if designated owner
     * @param globalElementIndex is the global index of the cable element
     */
     bool CalculateDesignatedOwnershipOfCableElement( unsigned globalElementIndex );

     /** Iterator type over #mNodeToCablesMapping. */
     typedef typename std::multimap<const Node<SPACE_DIM>*, Element<1u, SPACE_DIM>*>::iterator NodeCableIterator;

     /** The type returned by GetCablesAtNode. */
     typedef std::pair<NodeCableIterator, NodeCableIterator> CableRangeAtNode;

     /**
      * @return the cables that are attached to the given node.
      *
      * @param pNode a node to find the adjoining cables of
      * @return the adjoining cables.
      */
     CableRangeAtNode GetCablesAtNode(const Node<SPACE_DIM>* pNode);

private:
    /** The elements making up the 1D cables */
    std::vector<Element<1u, SPACE_DIM>*> mCableElements;

    /** The global number of cables over all processes*/
    unsigned mNumCableElements;

    /** A map from global cable index to local index used by this process. */
    std::map<unsigned, unsigned> mCableElementsMapping;

    /** Records which cables are attached to each node. */
    std::multimap<const Node<SPACE_DIM>*, Element<1u, SPACE_DIM>*> mNodeToCablesMapping;

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
        archive & boost::serialization::base_object<DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

    /** Definition of cable element Iterator type. */
    typedef typename std::vector<Element<1, SPACE_DIM> *>::const_iterator CableElementIterator;

    /**
     * @return a pointer to the first boundary element in the mesh.
     */
    CableElementIterator GetCableElementIteratorBegin() const;

    /**
     * @return a pointer to *one past* the last boundary element in the mesh
     * (for consistency with STL iterators).
     */
    CableElementIterator GetCableElementIteratorEnd() const;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MixedDimensionMesh)

namespace boost
{
namespace serialization
{
/**
 * Record number of processors when saving...
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    unsigned num_procs = PetscTools::GetNumProcs();
    const DistributedTetrahedralMeshPartitionType::type partition_type = t->GetPartitionType();
    ar << num_procs;
    ar << partition_type;
}

/**
 * De-serialize constructor parameters and initialise a MixedDimensionMesh,
 * checking the number of processors is the same.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    unsigned num_procs;
    DistributedTetrahedralMeshPartitionType::type partition_type;

    ar >> num_procs;
    ar >> partition_type;

    // Invoke inplace constructor to initialise instance
    /// \todo #1199  Lots of stuff can't cope if we re-partition
    //::new(t)MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>(partition_type);
    ::new(t)MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>(DistributedTetrahedralMeshPartitionType::DUMB);

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

#endif /*MIXEDDIMENSIONMESH_HPP_*/
