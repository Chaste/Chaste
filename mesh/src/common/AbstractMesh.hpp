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

#ifndef ABSTRACTMESH_HPP_
#define ABSTRACTMESH_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include "UblasVectorInclude.hpp"
#include "UblasMatrixInclude.hpp"

#include <vector>
#include <string>
#include <cassert>

#include "Node.hpp"
#include "DistributedVectorFactory.hpp"
#include "ProcessSpecificArchive.hpp"
#include "ChasteCuboid.hpp"

#include <boost/utility.hpp>

/**
 * Abstract base class for all meshes.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractMesh : private boost::noncopyable
{
    friend class TestDistributedTetrahedralMesh;
    template <unsigned A_DIMENSION> friend class NodesOnlyMesh; //NodesOnlyMesh is able to grab the node information in order to copy
    template <unsigned A_DIMENSION> friend class QuadraticMeshHelper;

private:
    /**
     * Pure virtual solve node mapping method. For a node with a given global
     * index, get the local index used by this process.
     *
     * Overridden in TetrahedralMesh DistributedTetrahedralMesh and Vertex Mesh classes.
     *
     * @param index the global index of the node
     * @return local index
     */
    virtual unsigned SolveNodeMapping(unsigned index) const = 0;

    /** Needed for serialization. */
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
        archive & mMeshChangesDuringSimulation;
        (*ProcessSpecificArchive<Archive>::Get()) & mpDistributedVectorFactory;
    }

protected:  // Give access of these variables to subclasses

    /** Vector of pointers to nodes in the mesh. */
    std::vector<Node<SPACE_DIM> *> mNodes;

    /** Vector of pointers to boundary nodes in the mesh. */
    std::vector<Node<SPACE_DIM> *> mBoundaryNodes;

    /**
     * DistributedVectorFactory capable of reproducing the
     * given number of nodes owned by each processor.
     */
    DistributedVectorFactory* mpDistributedVectorFactory;

    /**
     * Vector containing node permutation information.
     *
     *  When empty (most meshes) there is no node permutation
     *  When non-empty (parallel distributed meshes) then for a given original_index
     *  #mNodePermutation[original_index] holds the new assigned index of that node in memory
     */
    std::vector<unsigned> mNodePermutation;

    /**
     * If the mesh is constructed from file using a MeshReader, this member
     * variable stores the base name of these files.
     */
    std::string mMeshFileBaseName;

    /**
     * Whether this mesh changes during simulation (used to know whether to write a new one to file)
     */
    bool mMeshChangesDuringSimulation;

    /**
     * Does nothing.  Used in derived classes which have elements
     */
    virtual void SetElementOwnerships();

    /**
     * Calculate a bounding box from a set of nodes. A generalised version of the public
     * CalculateBoundingBox method.
     *
     * @param rNodes the list of nodes to calculate the bounding box for.
     * @return the bounding box in the form of a ChasteCuboid.
     */
    ChasteCuboid<SPACE_DIM> CalculateBoundingBox(const std::vector<Node<SPACE_DIM>* >& rNodes) const;

public:

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

    /** Definition of boundary node Iterator type. */
    typedef typename std::vector<Node<SPACE_DIM> *>::const_iterator BoundaryNodeIterator;

    /** Forward declaration */
    class NodeIterator;

    /**
     * @return an iterator to the first node in the mesh.
     *
     * @param skipDeletedNodes whether to include deleted nodes
     */
    inline NodeIterator GetNodeIteratorBegin(bool skipDeletedNodes=true);

    /**
     * @return an iterator to one past the last node in the mesh.
     */
    inline NodeIterator GetNodeIteratorEnd();

    //////////////////////////////////////////////////////////////////////
    //                             Methods                              //
    //////////////////////////////////////////////////////////////////////

    /**
     * Constructor.
     */
    AbstractMesh();

    /**
     * Virtual destructor, since this class has virtual methods.
     */
    virtual ~AbstractMesh();

    /**
     * @return the number of nodes that are actually in use.
     *
     * Overridden in MutableMesh and DistributedTetrahedralMesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * @return the number of boundary nodes in the mesh.
     */
    unsigned GetNumBoundaryNodes() const;

    /**
     * @return the total number of nodes (including those marked as deleted).
     */
    virtual unsigned GetNumAllNodes() const;

    /**
     * @return the  number of attributes on each node.
     * Note, this implementation assumes that all nodes have the same number of attributes
     * so that the first node in the container is indicative of the others.
     */
    unsigned GetNumNodeAttributes() const;

    /**
     * Get the node with a given index in the mesh.
     *
     * @param index the global index of the node
     * @return a pointer to the node.
     */
    Node<SPACE_DIM>* GetNode(unsigned index) const;

    /**
     * Get the node with a given index in the mesh (synonym of GetNode() unless overridden in a distributed mesh).
     *
     * @param index the global index of the node
     * @return a pointer to the node.
     */
    virtual Node<SPACE_DIM>* GetNodeOrHaloNode(unsigned index) const;

    /**
     * Get the node with a given index in the mesh, prior to any node permutation
     * being applied.  For non-permuted meshes, this will have the same effect
     * as GetNode.
     *
     * This method is intended for use by the archiving code, to enable checkpoint
     * migration, so that we can load the correct cells and boundary conditions
     * after the mesh has been re-partitioned.
     *
     * If unsure, use GetNode in preference to this method!
     *
     * @param index the global index of the node prior to a permutation being applied
     * @return a pointer to the node
     */
    Node<SPACE_DIM>* GetNodeFromPrePermutationIndex(unsigned index) const;

    /**
     * Read in the number of nodes per processor from file.
     *
     * @param rNodesPerProcessorFile the name of the file
     */
    virtual void ReadNodesPerProcessorFile(const std::string& rNodesPerProcessorFile);

    /**
     * @return the DistributedVectorFactory.
     */
    virtual DistributedVectorFactory * GetDistributedVectorFactory();

    /**
     * Set method for #mpDistributedVectorFactory.
     * Must be called before the mesh is used for anything.  This only actually
     * impacts the DistributedTetrahedralMesh subclass, in which the supplied factory
     * is then used to specify the node distribution among the processes.
     *
     * @param pFactory a factory to use for this mesh
     */
    virtual void SetDistributedVectorFactory(DistributedVectorFactory* pFactory);

    /**
     * Permute the nodes so that they appear in a different order in mNodes
     * (and their mIndex's are altered accordingly).
     */
    virtual void PermuteNodes();

    /**
     * @return a pointer to the first boundary node in the mesh.
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorBegin() const;

    /**
     * @return a pointer to *one past* the last boundary node in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorEnd() const;

    /**
     * @return mMeshFileBaseName.
     */
    std::string GetMeshFileBaseName() const;

    /**
     * Get whether this mesh was read from file.
     *
     * @return whether this mesh was read from file
     */
    bool IsMeshOnDisk() const;

    /**
     * @return #mNodePermutation.
     *
     *  When empty (most meshes) there is no node permutation
     *  When non-empty (parallel distributed meshes) then for a given original_index
     *  #mNodePermutation[original_index] holds the new assigned index of that node in memory
     *
     */
    const std::vector<unsigned>& rGetNodePermutation() const;

    /**
     * @return a vector between two points in space.
     *
     * This method is overridden in some daughter classes (e.g. Cylindrical2dMesh).
     *
     * @param rLocationA a c_vector of coordinates
     * @param rLocationB a c_vector of coordinates
     *
     * @return c_vector from location A to location B.
     */
    virtual c_vector<double, SPACE_DIM> GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocationA,
                                                          const c_vector<double, SPACE_DIM>& rLocationB);

    /**
     * Return the distance between two nodes.
     *
     * This method calls GetVectorFromAtoB(), which is overridden in some
     * daughter classes (e.g. Cylindrical2dMesh).
     *
     * @param indexA a node index
     * @param indexB a node index
     *
     * @return distance between two nodes.
     */
    double GetDistanceBetweenNodes(unsigned indexA, unsigned indexB);

    /**
     * Calculate the 'width' of any dimension of the mesh.
     *
     * This method is overridden in some daughter classes (e.g. Cylindrical2dMesh).
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     */
    virtual double GetWidth(const unsigned& rDimension) const;

    /**
     * Calculate the bounding box (width extremes for all dimensions of the mesh).
     * Overridden in distributed case
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
      * This method uses GetVectorFromAtoB distance and hence may return a correct solution
      * in non-Euclidean space, but only if this method is overridden in a subclass
      * (see e.g. Cylindrical2dMesh for an example of this).
      *
      * @param rTestPoint reference to the point
      * @return node index
      */
    virtual unsigned GetNearestNodeIndex(const ChastePoint<SPACE_DIM>& rTestPoint);

    /**
     * Scale the mesh.
     *
     * @param xFactor is the scale in the x-direction (defaults to 1.0)
     * @param yFactor is the scale in the y-direction (defaults to 1.0)
     * @param zFactor is the scale in the z-direction (defaults to 1.0)
     */
    virtual void Scale(const double xFactor=1.0, const double yFactor=1.0, const double zFactor=1.0);

    /**
     * Translate the mesh given the displacement vector.
     * This is the translation method that actually does the work.
     * Should be overridden when the child class has halo nodes.
     *
     * @param rDisplacement is a translation vector of the correct size
     */
    virtual void Translate(const c_vector<double, SPACE_DIM>& rDisplacement);

    /**
     * Translate the mesh given the coordinate displacements separately.
     *
     * @param xMovement is the x-displacement (defaults to 0.0)
     * @param yMovement is the y-displacement (defaults to 0.0)
     * @param zMovement is the z-displacement (defaults to 0.0)
     */
    void Translate(const double xMovement=0.0, const double yMovement=0.0, const double zMovement=0.0);

    /**
     * Do a general mesh rotation with a positive determinant orthonormal rotation matrix.
     * This is the rotation method that actually does the work.
     * Should be overridden when the child class has halo nodes.
     *
     * @param rotationMatrix is a Ublas rotation matrix of the correct form
     */
    virtual void Rotate(c_matrix<double , SPACE_DIM, SPACE_DIM> rotationMatrix);

    /**
     * Do an angle axis rotation.
     *
     * @param axis is the axis of rotation (does not need to be normalised)
     * @param angle is the angle of rotation in radians
     */
    void Rotate(c_vector<double,3> axis, double angle);

    /**
     * Rotate the mesh about the x-axis.
     *
     * @param theta is the angle of rotation in radians
     */
    void RotateX(const double theta);

    /**
     * Rotate the mesh about the y-axis.
     *
     * @param theta is the angle of rotation in radians
     */
    void RotateY(const double theta);

    /**
     * Rotate the mesh about the z-axis.
     *
     * @param theta is the angle of rotation in radians
     */
    void RotateZ(const double theta);

    /**
     * Rotating a 2D mesh equates that rotation around the z-axis.
     *
     * @param theta is the angle of rotation in radians
     */
    void Rotate(double theta);

    /**
     * This method allows the mesh properties to be re-calculated after
     * one or more nodes have been moved.
     */
    virtual void RefreshMesh();

    /**
     * @return Whether the mesh changes (used in archiving).
     */
    bool IsMeshChanging() const;

    /**
     * @return Iterates through local nodes and finds the maximum number of containing elements for all locally owned nodes
     * Useful for determining FEM matrix fill.
     */
    unsigned CalculateMaximumContainingElementsPerProcess() const;

    /**
     * Set whether the mesh has been modified since it was read from file.
     * This prevents the archiving code just blithely storing the original,
     * unmodified, mesh.
     */
    void SetMeshHasChangedSinceLoading();

    //////////////////////////////////////////////////////////////////////
    //                         Nested classes                           //
    //////////////////////////////////////////////////////////////////////

    /**
     * A smart iterator over the nodes in the mesh.
     */
    class NodeIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current node.
         * @return reference
         * Make sure to use a reference for the result to avoid copying nodes unnecessarily.
         */
        inline Node<SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline Node<SPACE_DIM>* operator->();

        /**
         * @return Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator& rOther);

        /**
         * Prefix increment operator.
         * @return incremented object
         */
        inline NodeIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * AbstractMesh::GetNodeIteratorBegin and AbstractMesh::GetNodeIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param nodeIter where to start iterating
         * @param skipDeletedNodes whether to include deleted nodes
         */
        NodeIterator(AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                     typename std::vector<Node<SPACE_DIM> *>::iterator nodeIter,
                     bool skipDeletedNodes=true);
    private:
        /** The mesh we're iterating over. */
        AbstractMesh& mrMesh;

        /** The actual node iterator. */
        typename std::vector<Node<SPACE_DIM> *>::iterator mNodeIter;

        /** Whether to skip deleted nodes. */
        bool mSkipDeletedNodes;

        /**
         * Helper method to say when we're at the end.
         * @return true if at end
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this node.
         * @return true if allowed
         */
        inline bool IsAllowedNode();
    };
};

TEMPLATED_CLASS_IS_ABSTRACT_2_UNSIGNED(AbstractMesh)

//////////////////////////////////////////////////////////////////////////////
//      NodeIterator class implementation - most methods are inlined        //
//////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNodeIteratorBegin(
        bool skipDeletedNodes)
{
    return NodeIterator(*this, mNodes.begin(), skipDeletedNodes);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNodeIteratorEnd()
{
    return NodeIterator(*this, mNodes.end());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>& AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::operator*()
{
    assert(!IsAtEnd());
    return **mNodeIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::operator->()
{
    assert(!IsAtEnd());
    return *mNodeIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::operator!=(const typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator& rOther)
{
    return mNodeIter != rOther.mNodeIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator& AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::operator++()
{
    do
    {
        ++mNodeIter;
    }
    while (!IsAtEnd() && !IsAllowedNode());

    return (*this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::NodeIterator(
        AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
        typename std::vector<Node<SPACE_DIM> *>::iterator nodeIter,
        bool skipDeletedNodes)
    : mrMesh(rMesh),
      mNodeIter(nodeIter),
      mSkipDeletedNodes(skipDeletedNodes)
{
    if (mrMesh.mNodes.size() == 0)
    {
        // Cope with empty meshes
        mNodeIter = mrMesh.mNodes.end();
    }
    else
    {
        // Make sure we start at an allowed node
        if (mNodeIter == mrMesh.mNodes.begin() && !IsAllowedNode())
        {
            ++(*this);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::IsAtEnd()
{
    return mNodeIter == mrMesh.mNodes.end();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::IsAllowedNode()
{
    return !(mSkipDeletedNodes && (*this)->IsDeleted());
}


#endif /*ABSTRACTMESH_HPP_*/
