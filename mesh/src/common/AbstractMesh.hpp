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
class AbstractMesh : boost::noncopyable
{
    friend class TestDistributedTetrahedralMesh;
    template <unsigned A_DIMENSION> friend class NodesOnlyMesh; //NodesOnlyMesh is able to grab the node information in order to copy
private:
    /**
     * Pure virtual solve node mapping method. For a node with a given global
     * index, get the local index used by this process.
     *
     * Overridden in TetrahedralMesh DistributedTetrahedralMesh and Vertex Mesh classes.
     *
     * @param index the global index of the node
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

    /** DistributedVectorFactory capable of reproducing the given number of nodes owned by each processor. */
    DistributedVectorFactory* mpDistributedVectorFactory;

    /** Vector containing node permutation information.
     *  When empty (most meshes) there is no node permutation
     *  When non-empty (parallel distributed meshes) then for a given original_index
     *  #mNodesPermutation[original_index] holds the new assigned index of that node in memory
     */
    std::vector<unsigned> mNodesPermutation;

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
public:

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

    /** Definition of boundary node Iterator type. */
    typedef typename std::vector<Node<SPACE_DIM> *>::const_iterator BoundaryNodeIterator;

    /** Forward declaration */
    class NodeIterator;

    /**
     * Get an iterator to the first node in the mesh.
     *
     * @param skipDeletedNodes whether to include deleted nodes
     */
    inline NodeIterator GetNodeIteratorBegin(bool skipDeletedNodes=true);

    /**
     * Get an iterator to one past the last node in the mesh.
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
     * Get the number of nodes that are actually in use.
     *
     * Overridden in MutableMesh and DistributedTetrahedralMesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * Get the number of boundary nodes in the mesh.
     */
    unsigned GetNumBoundaryNodes() const;

    /**
     * Get the total number of nodes (including those marked as deleted).
     */
    virtual unsigned GetNumAllNodes() const;

    /**
     * Get the node with a given index in the mesh.
     *
     * @param index the global index of the node
     * @return a pointer to the node.
     */
    Node<SPACE_DIM>* GetNode(unsigned index) const;

    /**
     * Get the node with a given index in the mesh (synonym of GetNode()).
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
     * Get method for DistributedVectorFactory.
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
     * Return a pointer to the first boundary node in the mesh.
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorBegin() const;

    /**
     * Return a pointer to *one past* the last boundary node in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorEnd() const;

    /**
     * Get method for mMeshFileBaseName.
     */
    std::string GetMeshFileBaseName() const;

    /**
     * Get whether this mesh was read from file.
     *
     * @return whether this mesh was read from file
     */
    bool IsMeshOnDisk() const;

    /**
     * Get method for #mNodesPermutation.
     *
     *  When empty (most meshes) there is no node permutation
     *  When non-empty (parallel distributed meshes) then for a given original_index
     *  #mNodesPermutation[original_index] holds the new assigned index of that node in memory
     *
     */
    const std::vector<unsigned>& rGetNodePermutation() const;

    /**
     * Return a vector between two points in space.
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
     * Calculate the bounding box (width extremes for all dimensions of the mesh.
     * Overridden in Distribute case
     *
     * @return The minimum and maximum co-ordinates of any node in each dimension
     *
     */
    virtual ChasteCuboid<SPACE_DIM> CalculateBoundingBox() const;

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
      *
      * @param rDisplacement is a translation vector of the correct size
      */
     void Translate(const c_vector<double, SPACE_DIM>& rDisplacement);

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
      *
      * @param rotationMatrix is a Ublas rotation matrix of the correct form
      */
     void Rotate(c_matrix<double , SPACE_DIM, SPACE_DIM> rotationMatrix);

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
         *
         * Make sure to use a reference for the result to avoid copying nodes unnecessarily.
         */
        inline Node<SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         */
        inline Node<SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator& rOther);

        /**
         * Prefix increment operator.
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
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this node.
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
bool AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::operator!=(const AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator& rOther)
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
