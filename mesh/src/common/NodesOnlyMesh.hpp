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

#ifndef NODESONLYMESH_HPP_
#define NODESONLYMESH_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>

#include "SmartPointers.hpp"
#include "PetscTools.hpp"
#include "DistributedBoxCollection.hpp"
#include "MutableMesh.hpp"
/**
 * Mesh class for storing lists of nodes (no elements). This inherits from MutableMesh
 * because we want to be able to add and delete nodes.
 */
template<unsigned SPACE_DIM>
class NodesOnlyMesh: public MutableMesh<SPACE_DIM, SPACE_DIM>
{
private:

    friend class TestNodesOnlyMesh;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archives the member variables of the object which have to be preserved
     * during its lifetime.
     *
     * Note that we must archive any member variables FIRST so that this
     * method can call a ReMesh (to convert from TrianglesMeshReader input
     * format into our native format).
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & mMaximumInteractionDistance;
        archive & mMinimumNodeDomainBoundarySeparation;
        std::vector<unsigned> indices = GetAllNodeIndices();
        archive & indices;
        archive & boost::serialization::base_object<MutableMesh<SPACE_DIM, SPACE_DIM> >(*this);
    }

    /**
     * Load member variables of the object which have to be preserved
     * during its lifetime.
     *
     * Note that we must archive any member variables FIRST so that this
     * method can call a ReMesh (to convert from TrianglesMeshReader input
     * format into our native format).
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & mMaximumInteractionDistance;
        archive & mMinimumNodeDomainBoundarySeparation;
        std::vector<unsigned> indices;
        archive & indices;
        archive & boost::serialization::base_object<MutableMesh<SPACE_DIM, SPACE_DIM> >(*this);
        // Re-index the nodes according to what we've just read
        assert(GetNumNodes() == indices.size());
        this->mNodesMapping.clear();
        for (unsigned i=0; i<this->mNodes.size(); i++)
        {
            unsigned new_index = indices[i];
            this->mNodes[i]->SetIndex(new_index);
            this->mNodesMapping[new_index] = i;
        }
        mMaxAddedNodeIndex = *(std::max_element(indices.begin(), indices.end()));
        mIndexCounter = mMaxAddedNodeIndex + 1; // Next available fresh index
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    /** Vector of shared-pointers to halo nodes used by this process. */
    std::vector<boost::shared_ptr<Node<SPACE_DIM> > > mHaloNodes;

    /** Nodes separated by a distance less than mMaximumInteractionDistance are neighbours. */
    double mMaximumInteractionDistance;

    /** A map from node global index to local index in mNodes. */
    std::map<unsigned, unsigned> mNodesMapping;

    /** A map from halo node global index to local index in mHaloNodes. */
    std::map<unsigned, unsigned> mHaloNodesMapping;

    /** A counter of the number of fresh node indices used on this process. Ensures unique indices in parallel. */
    unsigned mIndexCounter;

    /** A minimum separation to maintain between nodes and the boundary of mpBoxCollection, which grows with the mesh. */
    double mMinimumNodeDomainBoundarySeparation;

    /** A list of the global indices of nodes that have been deleted from this process and can be reused. */
    std::vector<unsigned> mDeletedGlobalNodeIndices;

    /** A list of global indices of nodes that need to be moved to the right hand process. */
    std::vector<unsigned> mNodesToSendRight;

    /** A list of global indices of nodes that need to be moved to the left hand process. */
    std::vector<unsigned> mNodesToSendLeft;

    /**A list of flags showing which initial nodes passed to ConstructNodesWithoutMesh
     * were created on this process. */
    std::vector<bool> mLocalInitialNodes;

    /** A variable to keep track of added node indices so we know the largest on this process. Used to
     * create a large enough NodeMap when remeshing. */
    unsigned mMaxAddedNodeIndex;

    /** A pointer to the DistributedBoxCollection. */
    DistributedBoxCollection<SPACE_DIM>* mpBoxCollection;

    /** Whether to calculate node neighbours in the box collection. Switch off for efficiency */
    bool mCalculateNodeNeighbours;

    /**
     * Calculate the next unique global index available on this
     * process. Uses a hashing function to ensure that a unique
     * index is given to every node.
     *
     * For example for 3 process they will have access to the following
     * integers:
     *
     * Proc 0:  0   3   6   9   12  ...
     *
     * Proc 1:  1   4   7   10   13  ...
     *
     * Proc 2:  2   5   8   11   14  ...
     *
     * Deleted node indices can be locally re-used.
     *
     * Deleted node indices of nodes that have *moved* process cannot be re-used.
     *
     * @return the next available index.
     */
    unsigned GetNextAvailableIndex();

    /** Increase box collection size by adding an extra row/face to each edge. */
    void EnlargeBoxCollection();

    /**
     * @return Whether any node (deleted or otherwise) is close to the domain boundary.
     */
    bool IsANodeCloseToDomainBoundary();

     /**
      * Set up a DistributedBoxCollection by calculating the correct domain size from the node locations.
      *
      * @param rNodes the nodes that will be contained in the box collection.
      */
     void SetUpBoxCollection(const std::vector<Node<SPACE_DIM>* >& rNodes);

     /**
      * Remove all nodes that return mIsDeleted as true.
      *
      * @param map the NodeMap to record which nodes have been removed.
      */
     void RemoveDeletedNodes(NodeMap& map);

     /** Make sure that node indices match their location, and update mNodesMapping. */
     void UpdateNodeIndices();

     /**
      * Add pNewNode to the mesh, maintaining its current global index. Called by AddNode and AddMovedNode.
      *
      * @param pNewNode the new node to add to this mesh.
      */
     void AddNodeWithFixedIndex(Node<SPACE_DIM>* pNewNode);

protected:

    /**  Clear the BoxCollection  */
    void ClearBoxCollection();

    /**
     * Set up the box collection. Overridden in subclasses to implement periodicity.
     *
     * @param cutOffLength the cut off length for node neighbours.
     * @param domainSize the size of the domain containing the nodes.
     * @param numLocalRows the number of rows that should be owned by this process.
     * @param isPeriodic whether the DistributedBoxCollection should be periodic.
     */
     virtual void SetUpBoxCollection(double cutOffLength, c_vector<double, 2*SPACE_DIM> domainSize, int numLocalRows = PETSC_DECIDE, bool isPeriodic = false);

     /** @return mpBoxCollection */
     DistributedBoxCollection<SPACE_DIM>* GetBoxCollection();

public:

    /**  Default constructor to initialise BoxCollection to NULL.  */
    NodesOnlyMesh();

    /**  Over-written destructor to delete pointer to BoxCollection. */
    virtual ~NodesOnlyMesh();

    /**
     * Construct the mesh using only nodes. No mesh is created, but the nodes are stored.
     * The original vector of nodes is deep-copied: new node objects are made with are
     * independent of the pointers in the input so that they can be safely deleted.
     *
     * If this is the only way of constructing a mesh of this type, then we can be certain that
     * elements and boundary elements are always unused.
     *
     * @param rNodes a vector of pointers to nodes.
     * @param maxInteractionDistance the distance that defines node neighbours in CalculateNodePairs.
     */
    void ConstructNodesWithoutMesh(const std::vector<Node<SPACE_DIM>*>& rNodes, double maxInteractionDistance);

    /**
     * Construct the mesh using only nodes. No mesh is created, but the nodes are stored.
     * The original vector of nodes is deep-copied: new node objects are made with are
     * independent of the pointers in the input so that they can be safely deleted.
     *
     * If this is the only way of constructing a mesh of this type, then we can be certain that
     * elements and boundary elements are always unused.
     *
     * @param rNodes a vector of shared pointers to nodes.
     * @param maxInteractionDistance the distance that defines node neighbours in CalculateNodePairs.
     */
    void ConstructNodesWithoutMesh(const std::vector<boost::shared_ptr<Node<SPACE_DIM> > >& rNodes, double maxInteractionDistance);

    /**
     * A Helper method to enable you to construct a nodes-only mesh by stripping the nodes
     * TetrahedralMesh, this calls the ConstructNodesWithoutMesh method with the nodes
     *
     * @param rGeneratingMesh any mesh with nodes, used to generate the NodesOnlyMesh.
     * @param maxInteractionDistance the distance that defines node neighbours in CalculateNodePairs.
     */
    void ConstructNodesWithoutMesh(const AbstractMesh<SPACE_DIM,SPACE_DIM>& rGeneratingMesh, double maxInteractionDistance);

    /** @return whether each initial node given to ConstructNodesWithoutMesh is owned by this process. */
    std::vector<bool>& rGetInitiallyOwnedNodes();

    /**  Overridden Clear() method for NodesOnlyMesh.  */
    void Clear();

    /**
     * Overridden solve node mapping method
     *
     * @param index the global index of the node.
     *
     * @return the local index of the node in mNodes.
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Overridden method to get node or halo node from the mesh.
     *
     * @param index the global index of the node.
     *
     * @return a pointer to the node.
     */
    Node<SPACE_DIM>* GetNodeOrHaloNode(unsigned index) const;

    /**
     * A method to identify whether a location is owned in the parallel space decomposition.
     * @param location the location to test
     * @return whether the point is owned.
     */
    bool IsOwned(c_vector<double, SPACE_DIM>& location);

    /**
     * @return the local number of nodes that are actually in use.
     * Does not include halo nodes.
     */
    unsigned GetNumNodes() const;

    /**
     * Get the largest node global index on this process.
     *
     * @return the maximum node index.
     */
    virtual unsigned GetMaximumNodeIndex();

    /**
     * Set the maximum node interaction distance.
     *
     * @param maxDistance the new maximum distance.
     */
    void SetMaximumInteractionDistance(double maxDistance);

    /**
     * @return mMaxInteractionDistance.
     */
    double GetMaximumInteractionDistance();

    /**
     * Overridden GetWidth method to work in parallel.
     *
     * @param rDimension the dimension along which to get the width.
     * @return the width.
     */
    double GetWidth(const unsigned& rDimension) const;

    /**
     * Set whether to calculate node neighbours for the rNodeNeigbours set in CalculateNodePairs. Switch off for efficiency
     *
     * @param calculateNodeNeighbours whether to store the neighbours.
     */
    void SetCalculateNodeNeighbours(bool calculateNodeNeighbours);

    /**
     * Calculate pairs of nodes from interior boxes using the BoxCollection.
     *
     * @param rNodePairs reference to the set of node pairs to populate.
     */
    void CalculateInteriorNodePairs(std::vector<std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> >& rNodePairs);

    /**
     * Calculate pairs of nodes from boxes on the process boundary using the BoxCollection
     *
     * @param rNodePairs reference to the set of node pairs to populate.
     */
    void CalculateBoundaryNodePairs(std::vector<std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> >& rNodePairs);

    /**
     * Overridden ReMesh() method. Since only Nodes are stored, this method cleans up mNodes by
     * removing nodes marked as deleted and reallocating mNodes to 'fill the gaps'.
     *
     * @param rMap a reference to a NodeMap which associates the indices of the old mesh
     * with the new mesh. It should be large enough to contain all node indices.
     */
    void ReMesh(NodeMap& rMap);

    /**
     * Set the initial box collection without passing nodes to the mesh. Used for memory efficient construction in parallel.
     * @param domainSize the initial domain size of the mesh.
     * @param maxInteractionDistance the max interaction distance between nodes.
     */
    void SetInitialBoxCollection(const c_vector<double, 2*SPACE_DIM> domainSize, double maxInteractionDistance);

    /**
     * Clear the old box collection and set up a new one if necessary.
     */
    void UpdateBoxCollection();

    /**
     * Check whether any nodes are close to the edge of the box collection and increase
     * the size of it if necessary.
     */
    void ResizeBoxCollection();

    /**
     * Iterate through each node and add it to its appropriate box.
     */
    void AddNodesToBoxes();

    /**
     * Iterate through each halo node and add it to its appropriate box.
     */
    void AddHaloNodesToBoxes();

    /**
     * Work out which nodes lie outside the local domain and add their indices to the vectors #mNodesToSendLeft and #mNodesToSendRight.
     */
    void CalculateNodesOutsideLocalDomain();

    /**
     * @return #mNodesToSendLeft.
     */
    std::vector<unsigned>& rGetNodesToSendLeft();

    /**
     * @return #mNodesToSendRight.
     */
    std::vector<unsigned>& rGetNodesToSendRight();

    /**
     * @return the indices of halo nodes, owned by this process, on the right hand boundary.
     */
    std::vector<unsigned>& rGetHaloNodesToSendRight();

    /**
     * @return the indices of halo nodes, owned by this process, on the left hand boundary.
     */
    std::vector<unsigned>& rGetHaloNodesToSendLeft();

    /**
     * Add a temporary halo node on this process.
     * @param pNewNode a shared pointer to the new node to add.
     */
    void AddHaloNode(boost::shared_ptr<Node<SPACE_DIM> > pNewNode);

    /**
     * Delete all the halo nodes on this process.
     */
    void ClearHaloNodes();

    /**
     * Overridden SetNode() method.
     *
     * @param nodeIndex is the index of the node to be moved
     * @param point is the new target location of the node
     * @param concreteMove is set to false if we want to skip the signed area tests in the parent Class Note this should always be false here
     */
    void SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point, bool concreteMove = false);

    /**
     * Overridden AddNode() method.
     *
     * @param pNewNode  pointer to the new node.
     * @return index of new node.
     */
    unsigned AddNode(Node<SPACE_DIM>* pNewNode);

    /**
     * Add a node to this process that has moved from another process.
     * @param pMovedNode the node to add to this mesh.
     */
    void AddMovedNode(boost::shared_ptr<Node<SPACE_DIM> > pMovedNode);

    /**
     * Overridden DeleteNode() method.
     *
     * @param index of the node to be deleted
     */
    void DeleteNode(unsigned index);

    /**
     * Make a clean delete of a node that has moved off this process.
     *
     * @param index the global index of the node moving off this process.
     */
    void DeleteMovedNode(unsigned index);

    /**
     * Set the value of mMinimumNodeDomainBoundarySeparation.
     *
     * @param separation the new value for the separation.
     */
    void SetMinimumNodeDomainBoundarySeparation(double separation);

    /**
     * Re-allocate the underlaying BoxCollection rows based on the load-balance algorithm implemented
     * in the box collection.
     */
    void LoadBalanceMesh();

    /**
     * Overridden ConstructFromMeshReader to correctly assign global node indices on load.
     *
     * @param rMeshReader the mesh reader for input.
     */
    void ConstructFromMeshReader(AbstractMeshReader<SPACE_DIM, SPACE_DIM>& rMeshReader);

    /**
     * Get all node indices in order of appearance.
     * @return Node vector of node indices for this process ignoring all delete nodes
     */
    std::vector<unsigned> GetAllNodeIndices() const;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodesOnlyMesh)

#endif /*NODESONLYMESH_HPP_*/
