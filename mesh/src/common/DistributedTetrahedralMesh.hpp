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

#ifndef DISTRIBUTEDTETRAHEDRALMESH_HPP_
#define DISTRIBUTEDTETRAHEDRALMESH_HPP_

#include <map>
#include <vector>
#include <set>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractTetrahedralMesh.hpp"
#include "Node.hpp"
#include "AbstractMeshReader.hpp"
#include "DistributedTetrahedralMeshPartitionType.hpp"

#define UNASSIGNED_NODE UINT_MAX


/**
 * Parallel implementation of a mesh
 * Nodes are distributed such that each process has
 * A set of nodes (possibly reordered) with contiguous global indices
 * A local copy of all the elements supporting those nodes
 * A local copy of ghost/halo nodes which are all the nodes used in the supporting elements, but not owned outright.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class DistributedTetrahedralMesh : public AbstractTetrahedralMesh< ELEMENT_DIM, SPACE_DIM>
{
    friend class TestDistributedTetrahedralMesh;
    friend class TestDistributedQuadraticMesh;
private:

    /** The total number of elements in the mesh. */
    unsigned mTotalNumElements;

    /** The total number of boundary elements in the mesh. */
    unsigned mTotalNumBoundaryElements;

    /** The total number of nodes in the mesh. */
    unsigned mTotalNumNodes;

    /** Vector of pointer to halo nodes used by this process. */
    std::vector<Node<SPACE_DIM>* > mHaloNodes;

    /** A map from node global index to local index used by this process. */
    std::map<unsigned, unsigned> mNodesMapping;

    /** A map from halo node global index to local index used by this process. */
    std::map<unsigned, unsigned> mHaloNodesMapping;

    /** A map from element global index to local index used by this process. */
    std::map<unsigned, unsigned> mElementsMapping;

    /** A map from boundary element global index to local index used by this process. */
    std::map<unsigned, unsigned> mBoundaryElementsMapping;

    /** The region of space owned by this process, if using geometric partition. */
    ChasteCuboid<SPACE_DIM>* mpSpaceRegion;

    /** Partitioning method. */
    DistributedTetrahedralMeshPartitionType::type mPartitioning;

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
        archive & boost::serialization::base_object<AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> >(*this);
    }


    /**
     * Sets the ownership of each element according to which nodes are owned by the
     * process.
     *
     * Information on node ownership comes from the distributed vector factory and
     * an element is "owned" if one or more of its nodes are owned
     */
    void SetElementOwnerships();


public:

    /**
     * Constructor.
     *
     * @param partitioningMethod  defaults to PARMETIS_LIBRARY, but in 1-D is always overridden in this constructor to be the DUMB partition
     */
    DistributedTetrahedralMesh(DistributedTetrahedralMeshPartitionType::type partitioningMethod=DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);

    /**
     * Destructor.
     */
    virtual ~DistributedTetrahedralMesh();

    /**
     * Specify the node distribution across processes.
     * This also makes sure we don't try to use METIS to partition the mesh.
     *
     * @param pFactory a factory to use for this mesh
     */
    void SetDistributedVectorFactory(DistributedVectorFactory* pFactory);

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    virtual void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader);

    /**
     * @return the number of nodes that are entirely owned by the local process.
     * (Does not include halo nodes).
     */
    unsigned GetNumLocalNodes() const;

    /**
     * @return the number of nodes that are halo owned by the local process.
     */
    unsigned GetNumHaloNodes() const;

    /**
     * @return the number of Elements which are owned by this process (have at least one entirely
     * locally-owned node).
     */
    unsigned GetNumLocalElements() const;

    /**
     * @return the number of Boundary Elements which are owned by this process (have at least one entirely
     * locally-owned node).
     */
    unsigned GetNumLocalBoundaryElements() const;

    /**
     * @return the total number of nodes that are actually in use (globally).
     */
    unsigned GetNumNodes() const;

    /**
     * @return the total number of nodes that are actually in use (globally).
     */
    unsigned GetNumAllNodes() const;

    /**
     * @return the total number of elements that are actually in use (globally).
     */
    unsigned GetNumElements() const;

    /**
     * @return the type of mesh partitioning that is being used...
     *
     * serialization uses this method.
     */
    DistributedTetrahedralMeshPartitionType::type GetPartitionType() const;

    /**
     * @return the total number of boundary elements that are actually in use (globally).
     */
    unsigned GetNumBoundaryElements() const;

    /**
     * Utility method to give the functionality of iterating through the halo nodes of a process
     * @param rHaloIndices  A vector to fill with the global indices of the nodes which are locally halos
     */
    void GetHaloNodeIndices(std::vector<unsigned>& rHaloIndices) const;

    /**
     * Set the local region of space to be owned by this process
     *
     * @param pRegion The region, defined by a ChasteCuboid.
     */
    void SetProcessRegion(ChasteCuboid<SPACE_DIM>* pRegion);

    /**
     * Get the local region of space owned by this process
     *
     * @return mSpaceRegion
     */
    ChasteCuboid<SPACE_DIM>* GetProcessRegion();

    /**
     * Determine whether or not the current process owns node 0 of this element (tie breaker to determine which process writes
     * to file for when two or more share ownership of an element).
     *
     * @return true if process is designated owner
     * @param elementIndex is the global index of the element
     */
    bool CalculateDesignatedOwnershipOfElement( unsigned elementIndex );

    /**
     * Determine whether or not the current process owns node 0 of this boundary element (tie breaker to determine which process writes
     * to file for when two or more share ownership of a face).
     * @return true if process is designated owner
     * @param faceIndex is the global index of the face
     */
    bool CalculateDesignatedOwnershipOfBoundaryElement( unsigned faceIndex );

    /**
     * Construct a 1D linear grid on [0,width]
     *
     * Throws if there are more processes than the number of nodes (width+1)
     *
     * @param width  width of the mesh (in the x-direction)
     */
     void ConstructLinearMesh(unsigned width);

    /**
     * Construct a 2D rectangular grid on [0,width]x[0,height].
     *
     * Diagonals can be staggered so that there is no preferred
     * diffusion propagation direction.
     *
     * Distributed version splits the mesh in layers in the y-direction.
     * That is, the zeroth process will own from y=0 to about y=height/num_procs etc.
     *
     * @param width  width of the mesh (in the x-direction)
     * @param height  height of the mesh (in the y-direction)
     * @param stagger  whether the mesh should 'jumble' up the elements (defaults to true)
     */
    void ConstructRectangularMesh(unsigned width, unsigned height, bool stagger=true);

    /**
     * Construct a 3D cuboid grid on [0,width]x[0,height]x[0,depth].
     *
     * Distributed version splits the mesh in layers in the z-direction.
     * That is, the zeroth process will own from z=0 to about z=depth/num_procs etc.
     *
     * @param width  width of the mesh (in the x-direction)
     * @param height  height of the mesh (in the y-direction)
     * @param depth  depth of the mesh (in the z-direction)
     */
    void ConstructCuboid(unsigned width, unsigned height, unsigned depth);
    /**
     * Scale the mesh - uses the parent class for scaling the nodes.  This derived specialisation
     * is for scaling halo nodes.
     *
     * @param xFactor is the scale in the x-direction (defaults to 1.0)
     * @param yFactor is the scale in the y-direction (defaults to 1.0)
     * @param zFactor is the scale in the z-direction (defaults to 1.0)
     */
    virtual void Scale(const double xFactor=1.0, const double yFactor=1.0, const double zFactor=1.0);

    /**
     * @return the local pointer to a node which is
     * either owned or in the halo of this process.
     *
     * We first search halo node (as there are fewer),
     * then search totally owned nodes.  Otherwise throw.
     *
     * @param index the global index of the node
     */
    Node<SPACE_DIM>* GetNodeOrHaloNode(unsigned index) const;

    /**
     * Calculate the bounding box (width extremes for all dimensions of the mesh).
     * Override for distributed case
     *
     * @return The minimum and maximum co-ordinates of any node in each dimension
     *
     */
    virtual ChasteCuboid<SPACE_DIM> CalculateBoundingBox() const;

    /** GetNearestNodeIndex iterates through all nodes in the mesh and returns the global index
      * with the smallest distance to the provided point.
      *
      * This method is overridden in the distributed case to return the global node index.
      *
      * @param rTestPoint reference to the point
      * @return node index
      */
    virtual unsigned GetNearestNodeIndex(const ChastePoint<SPACE_DIM>& rTestPoint);

    /**
     * Computes the minimum and maximum lengths of the edges in the mesh.
     * This method overrides the default implementation in the parent class
     * \todo Should be const
     *
     * @return The minimum and maximum edge lengths in the mesh
     *
     */
    virtual c_vector<double, 2> CalculateMinMaxEdgeLengths();

    /**
     * Do a general mesh rotation with a positive determinant orthonormal rotation matrix.
     * This is the rotation method that actually does the work.
     * This override is because class has halo nodes.
     *
     * @param rotationMatrix is a Ublas rotation matrix of the correct form
     */
    void Rotate(c_matrix<double , SPACE_DIM, SPACE_DIM> rotationMatrix);

    /**
     * Translate the mesh given the displacement vector.
     * This is the translation method that actually does the work.
     * This override is because class has halo nodes.
     *
     * @param rDisplacement is a translation vector of the correct size
     */
    void Translate(const c_vector<double, SPACE_DIM>& rDisplacement);


protected:
    /**
     * Overridden solve node mapping method.
     *
     * @param index the global index of the node
     * @return local index
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Overridden solve element mapping method.
     *
     * @param index the global index of the element
     * @return local index
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Overridden solve boundary element mapping method.
     *
     * @param index the global index of the boundary element
     * @return local index
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;
private:

    /**
     * Add the most recently constructed node to the global->local node mapping
     *
     * @param index is the global index of node to be registered
     */
    void RegisterNode(unsigned index);

    /**
     * Add the most recently constructed halo node to the global->local halo node mapping
     *
     * @param index is the global index of halo node to be registered
     */
    void RegisterHaloNode(unsigned index);

    /**
     * Add the most recently constructed element to the global->local element mapping
     *
     * @param index is the global index of element to be registered
     */
    void RegisterElement(unsigned index);

     /**
     * Add the most recently constructed boundary element to the global->local boundary element mapping
     *
     * @param index is the global index of boundary element to be registered
     */
    void RegisterBoundaryElement(unsigned index);


    /**
     * Compute a parallel partitioning of a given mesh
     * using specialised methods below based on the value
     * of mPartitioning
     *
     * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
     * @param rNodesOwned is a set to be filled with the indices of nodes owned by this process
     * @param rHaloNodesOwned is a set to be filled with the indices of halo nodes owned by this process
     * @param rElementsOwned is a set to be filled with the indices of elements owned by this process
     * @param rProcessorsOffset a vector of length NumProcs to be filled with the index of the lowest indexed node owned by each process
     */
    void ComputeMeshPartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                 std::set<unsigned>& rNodesOwned,
                                 std::set<unsigned>& rHaloNodesOwned,
                                 std::set<unsigned>& rElementsOwned,
                                 std::vector<unsigned>& rProcessorsOffset);

    /**
      * Specialised method to compute a parallel partitioning of a given mesh with the ParMetis library
      * (called by ComputeMeshPartitioning, based on the value of mPartitioning)
      *
      * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
      * @param rElementsOwned is an empty set to be filled with the indices of elements owned by this process
      * @param rNodesOwned is an empty set to be filled with the indices of nodes owned by this process
      * @param rHaloNodesOwned is an empty set to be filled with the indices of halo nodes owned by this process
      * @param rProcessorsOffset a vector of length NumProcs to be filled with the index of the lowest indexed node owned by each process
      *
      */
     void ParMetisLibraryNodeAndElementPartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                          std::set<unsigned>& rElementsOwned,
                                          std::set<unsigned>& rNodesOwned,
                                          std::set<unsigned>& rHaloNodesOwned,
                                          std::vector<unsigned>& rProcessorsOffset);

    /**
     * Reorder the node indices in this mesh by applying the permutation
     * give in mNodePermutation.
     *
     * The node indexed with "i" will be re-assigned with the new index mNodePermutation[i]
     */
    void ReorderNodes();

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

public:
    ///\todo #1494, this iterator needs to be dereferenced twice because it is an STL iterator to a pointer.
    // The other iterators aren't (so only need to be dereferenced once). Consistency would be good...

    /** Definition of halo node Iterator type. */
    typedef typename std::vector<Node<SPACE_DIM> *>::const_iterator HaloNodeIterator;

    /**
     * @return an iterator to the first halo node in the mesh.
     */
    HaloNodeIterator GetHaloNodeIteratorBegin() const;

    /**
     * @return an iterator to one past the last halo node in the mesh.
     */
    HaloNodeIterator GetHaloNodeIteratorEnd() const;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(DistributedTetrahedralMesh)

namespace boost
{
namespace serialization
{
/**
 * Record number of processors when saving...
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    unsigned num_procs = PetscTools::GetNumProcs();
    const DistributedTetrahedralMeshPartitionType::type partition_type = t->GetPartitionType();
    ar << num_procs;
    ar << partition_type;
}

/**
 * De-serialize constructor parameters and initialise a DistributedTetrahedralMesh,
 * checking the number of processors is the same.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    unsigned num_procs;
    DistributedTetrahedralMeshPartitionType::type partition_type;

    ar >> num_procs;
    ar >> partition_type;

    // Invoke inplace constructor to initialise instance
    /// \todo #1199  Lots of stuff can't cope if we re-partition
    //::new(t)DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>(partition_type);
    ::new(t)DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>(DistributedTetrahedralMeshPartitionType::DUMB);

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

#endif /*DISTRIBUTEDTETRAHEDRALMESH_HPP_*/
