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

#include <map>
#include "NodesOnlyMesh.hpp"
#include "ChasteCuboid.hpp"

template<unsigned SPACE_DIM>
NodesOnlyMesh<SPACE_DIM>::NodesOnlyMesh()
        : MutableMesh<SPACE_DIM, SPACE_DIM>(),
          mMaximumInteractionDistance(1.0),
          mIndexCounter(0u),
          mMinimumNodeDomainBoundarySeparation(1.0),
          mMaxAddedNodeIndex(0u),
          mpBoxCollection(nullptr),
          mCalculateNodeNeighbours(true)
{
}

template<unsigned SPACE_DIM>
NodesOnlyMesh<SPACE_DIM>::~NodesOnlyMesh()
{
    Clear();
    ClearBoxCollection();
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ConstructNodesWithoutMesh(const std::vector<Node<SPACE_DIM>*>& rNodes, double maxInteractionDistance)
{
    assert(maxInteractionDistance > 0.0 && maxInteractionDistance < DBL_MAX);
    mMaximumInteractionDistance = maxInteractionDistance;

    mMinimumNodeDomainBoundarySeparation = mMaximumInteractionDistance;

    Clear();

    SetUpBoxCollection(rNodes);

    mLocalInitialNodes.resize(rNodes.size(), false);

    for (unsigned i=0; i<rNodes.size(); i++)
    {
        if (mpBoxCollection->IsOwned(rNodes[i]))
        {
            assert(!rNodes[i]->IsDeleted());

            mLocalInitialNodes[i] = true;

            // Create a copy of the node, sharing its location
            c_vector<double, SPACE_DIM> location = rNodes[i]->rGetLocation();
            Node<SPACE_DIM>* p_node_copy = new Node<SPACE_DIM>(GetNextAvailableIndex(), location);

            p_node_copy->SetRadius(0.5);

            // If the original node has attributes, then copy these
            if (rNodes[i]->HasNodeAttributes())
            {
                p_node_copy->rGetNodeAttributes() = rNodes[i]->rGetNodeAttributes();
            }

            this->mNodes.push_back(p_node_copy);

            // Update the node map
            mNodesMapping[p_node_copy->GetIndex()] = this->mNodes.size()-1;
        }
    }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ConstructNodesWithoutMesh(const std::vector<boost::shared_ptr<Node<SPACE_DIM> > >& rNodes, double maxInteractionDistance)
{
    // This is not efficient. It should replace the corresponding raw ptr method if SetUpBoxCollection and Chaste Cuboid methods are changed to take shared ptrs.
    std::vector<Node<SPACE_DIM>*> temp_nodes(rNodes.size());
    for(unsigned idx=0; idx<rNodes.size(); idx++)
    {
        temp_nodes[idx] = rNodes[idx].get();
    }

    ConstructNodesWithoutMesh(temp_nodes, maxInteractionDistance);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ConstructNodesWithoutMesh(const AbstractMesh<SPACE_DIM,SPACE_DIM>& rGeneratingMesh, double maxInteractionDistance)
{
    ConstructNodesWithoutMesh(rGeneratingMesh.mNodes, maxInteractionDistance);
}

template<unsigned SPACE_DIM>
std::vector<bool>& NodesOnlyMesh<SPACE_DIM>::rGetInitiallyOwnedNodes()
{
    return mLocalInitialNodes;
}

template<unsigned SPACE_DIM>
unsigned NodesOnlyMesh<SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    std::map<unsigned, unsigned>::const_iterator node_position = mNodesMapping.find(index);

    if (node_position == mNodesMapping.end())
    {
        EXCEPTION("Requested node " << index << " does not belong to process " << PetscTools::GetMyRank());
    }

    return node_position->second;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::Clear()
{
    // Call Clear() on the parent class
    MutableMesh<SPACE_DIM,SPACE_DIM>::Clear();

    // Clear the nodes mapping
    mNodesMapping.clear();

    mIndexCounter = 0;
}

template<unsigned SPACE_DIM>
DistributedBoxCollection<SPACE_DIM>* NodesOnlyMesh<SPACE_DIM>::GetBoxCollection()
{
    return mpBoxCollection;
}

template <unsigned SPACE_DIM>
Node<SPACE_DIM>* NodesOnlyMesh<SPACE_DIM>::GetNodeOrHaloNode(unsigned index) const
{
    Node<SPACE_DIM>* p_node;

    std::map<unsigned, unsigned>::const_iterator node_position = mHaloNodesMapping.find(index);

    if (node_position != mHaloNodesMapping.end())
    {
        p_node = mHaloNodes[node_position->second].get();
    }
    else
    {
        p_node = this->GetNode(index);
    }

    assert(p_node != nullptr);

    return p_node;
}

template<unsigned SPACE_DIM>
bool NodesOnlyMesh<SPACE_DIM>::IsOwned(c_vector<double, SPACE_DIM>& location)
{
    return mpBoxCollection->IsOwned(location);
}

template<unsigned SPACE_DIM>
unsigned NodesOnlyMesh<SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size() - this->mDeletedNodeIndices.size();
}

template<unsigned SPACE_DIM>
unsigned NodesOnlyMesh<SPACE_DIM>::GetMaximumNodeIndex()
{
    return std::max(mIndexCounter* PetscTools::GetNumProcs() + PetscTools::GetMyRank(), mMaxAddedNodeIndex);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::SetMaximumInteractionDistance(double maxDistance)
{
    mMaximumInteractionDistance = maxDistance;
}

template<unsigned SPACE_DIM>
double NodesOnlyMesh<SPACE_DIM>::GetMaximumInteractionDistance()
{
    return mMaximumInteractionDistance;
}

template<unsigned SPACE_DIM>
double NodesOnlyMesh<SPACE_DIM>::GetWidth(const unsigned& rDimension) const
{
    double local_width = AbstractMesh<SPACE_DIM, SPACE_DIM>::GetWidth(rDimension);
    double global_width;

    MPI_Allreduce(&local_width, &global_width, 1, MPI_DOUBLE, MPI_MAX, PetscTools::GetWorld());

    return global_width;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::SetCalculateNodeNeighbours(bool calculateNodeNeighbours)
{
    mCalculateNodeNeighbours = calculateNodeNeighbours;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::CalculateInteriorNodePairs(std::vector<std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> >& rNodePairs)
{
    assert(mpBoxCollection);

    mpBoxCollection->CalculateInteriorNodePairs(this->mNodes, rNodePairs);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::CalculateBoundaryNodePairs(std::vector<std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> >& rNodePairs)
{
    assert(mpBoxCollection);

    mpBoxCollection->CalculateBoundaryNodePairs(this->mNodes, rNodePairs);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ReMesh(NodeMap& map)
{
    map.ResetToIdentity();

    RemoveDeletedNodes(map);

    this->mDeletedNodeIndices.clear();
    this->mAddedNodes = false;

    UpdateNodeIndices();

    this->SetMeshHasChangedSinceLoading();
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::RemoveDeletedNodes(NodeMap& map)
{
    typename std::vector<Node<SPACE_DIM>* >::iterator node_iter = this->mNodes.begin();
    while (node_iter != this->mNodes.end())
    {
        if ((*node_iter)->IsDeleted())
        {
            map.SetDeleted((*node_iter)->GetIndex());

            mNodesMapping.erase((*node_iter)->GetIndex());

            // Free memory before erasing the pointer from the list of nodes.
            delete (*node_iter);
            node_iter = this->mNodes.erase(node_iter);
        }
        else
        {
            ++node_iter;
        }
    }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::UpdateNodeIndices()
{
    mNodesMapping.clear();
    for (unsigned location_in_vector=0; location_in_vector < this->mNodes.size(); location_in_vector++)
    {
        unsigned global_index = this->mNodes[location_in_vector]->GetIndex();
        mNodesMapping[global_index] = location_in_vector;
    }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::CalculateNodesOutsideLocalDomain()
{
    mNodesToSendRight.clear();
    mNodesToSendLeft.clear();

    for (typename AbstractMesh<SPACE_DIM, SPACE_DIM>::NodeIterator node_iter = this->GetNodeIteratorBegin();
            node_iter != this->GetNodeIteratorEnd();
            ++node_iter)
    {
        unsigned owning_process = mpBoxCollection->GetProcessOwningNode(&(*node_iter));
        if (owning_process == PetscTools::GetMyRank())
        {
            // Do nothing.
        }
        else if (owning_process == PetscTools::GetMyRank() + 1)
        {
            mNodesToSendRight.push_back(node_iter->GetIndex());
        }
        else if (owning_process == PetscTools::GetMyRank() - 1)
        {
            mNodesToSendLeft.push_back(node_iter->GetIndex());
        }
    }
}

template<unsigned SPACE_DIM>
std::vector<unsigned>& NodesOnlyMesh<SPACE_DIM>::rGetNodesToSendLeft()
{
    return mNodesToSendLeft;
}

template<unsigned SPACE_DIM>
std::vector<unsigned>& NodesOnlyMesh<SPACE_DIM>::rGetNodesToSendRight()
{
    return mNodesToSendRight;
}

template<unsigned SPACE_DIM>
std::vector<unsigned>& NodesOnlyMesh<SPACE_DIM>::rGetHaloNodesToSendRight()
{
    return mpBoxCollection->rGetHaloNodesRight();
}

template<unsigned SPACE_DIM>
std::vector<unsigned>& NodesOnlyMesh<SPACE_DIM>::rGetHaloNodesToSendLeft()
{
    return mpBoxCollection->rGetHaloNodesLeft();
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::AddNodeWithFixedIndex(Node<SPACE_DIM>* pNewNode)
{
    unsigned location_in_nodes_vector = 0;

    if (this->mDeletedNodeIndices.empty())
    {
        this->mNodes.push_back(pNewNode);
        location_in_nodes_vector = this->mNodes.size() - 1;
    }
    else
    {
        location_in_nodes_vector = this->mDeletedNodeIndices.back();
        this->mDeletedNodeIndices.pop_back();
        delete this->mNodes[location_in_nodes_vector];
        this->mNodes[location_in_nodes_vector] = pNewNode;
    }

    this->mAddedNodes = true;

    mMaxAddedNodeIndex = (pNewNode->GetIndex() > mMaxAddedNodeIndex) ? pNewNode->GetIndex() : mMaxAddedNodeIndex;

    // Update mNodesMapping
    mNodesMapping[pNewNode->GetIndex()] = location_in_nodes_vector;

    // Then update cell radius to default.
    pNewNode->SetRadius(0.5);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::AddHaloNode(boost::shared_ptr<Node<SPACE_DIM> > pNewNode)
{
    mHaloNodes.push_back(pNewNode);
    mHaloNodesMapping[pNewNode->GetIndex()] = mHaloNodes.size() - 1;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ClearHaloNodes()
{
    mHaloNodes.clear();

    mHaloNodesMapping.clear();
}

template<unsigned SPACE_DIM>
unsigned NodesOnlyMesh<SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNewNode)
{
    unsigned fresh_global_index = GetNextAvailableIndex();
    pNewNode->SetIndex(fresh_global_index);

    AddNodeWithFixedIndex(pNewNode);

    return fresh_global_index;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point, bool concreteMove)
{
    // concreteMove should always be false for a NodesOnlyMesh as there are no elements to check
    assert(!concreteMove);

    // Update the node's location
    this->GetNode(nodeIndex)->SetPoint(point);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::AddMovedNode(boost::shared_ptr<Node<SPACE_DIM> > pMovedNode)
{
    // Make a deep copy of this node pointer so that it isn't accidentally deleted.
    unsigned index = pMovedNode->GetIndex();
    c_vector<double, SPACE_DIM> location = pMovedNode->rGetLocation();

    Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(index, location);

    if (pMovedNode->HasNodeAttributes())
    {
        double radius = pMovedNode->GetRadius();
        p_node->SetRadius(radius);

        unsigned region = pMovedNode->GetRegion();
        p_node->SetRegion(region);

        bool is_particle = pMovedNode->IsParticle();
        p_node->SetIsParticle(is_particle);

        for (unsigned i=0; i<pMovedNode->GetNumNodeAttributes(); i++)
        {
            double attribute = pMovedNode->rGetNodeAttributes()[i];
            p_node->AddNodeAttribute(attribute);
        }
    }

    AddNodeWithFixedIndex(p_node);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::DeleteNode(unsigned index)
{
    if (this->GetNode(index)->IsDeleted())
    {
        EXCEPTION("Trying to delete a deleted node");
    }

    unsigned local_index = SolveNodeMapping(index);

    this->mNodes[local_index]->MarkAsDeleted();
    this->mDeletedNodeIndices.push_back(local_index);
    mDeletedGlobalNodeIndices.push_back(index);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::DeleteMovedNode(unsigned index)
{
    DeleteNode(index);

    // Remove index from deleted indices, as moved indices must not be re-used.
    mDeletedGlobalNodeIndices.pop_back();
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::SetMinimumNodeDomainBoundarySeparation(double separation)
{
    assert(!(separation < 0.0));

    mMinimumNodeDomainBoundarySeparation = separation;
}

template<unsigned SPACE_DIM>
unsigned NodesOnlyMesh<SPACE_DIM>::GetNextAvailableIndex()
{
    unsigned index;

    if (!this->mDeletedGlobalNodeIndices.empty())
    {
        index = this->mDeletedGlobalNodeIndices.back();
        this->mDeletedGlobalNodeIndices.pop_back();
    }
    else
    {
        unsigned counter = mIndexCounter;
        mIndexCounter++;
        index = counter * PetscTools::GetNumProcs() + PetscTools::GetMyRank();
    }

    return index;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::EnlargeBoxCollection()
{
    assert(mpBoxCollection);

    int num_local_rows = mpBoxCollection->GetNumLocalRows();
    int new_local_rows = num_local_rows + (int)(PetscTools::AmTopMost()) + (int)(PetscTools::AmMaster());

    c_vector<double, 2*SPACE_DIM> current_domain_size = mpBoxCollection->rGetDomainSize();
    c_vector<double, 2*SPACE_DIM> new_domain_size = current_domain_size;

    double fudge = 1e-14;
    // We don't enlarge the x direction if periodic
    unsigned d0 = ( mpBoxCollection->GetIsPeriodicInX() ) ? 1 : 0;
    for (unsigned d=d0; d < SPACE_DIM; d++)
    {
        new_domain_size[2*d] = current_domain_size[2*d] - (mMaximumInteractionDistance - fudge);
        new_domain_size[2*d+1] = current_domain_size[2*d+1] + (mMaximumInteractionDistance - fudge);
    }
    SetUpBoxCollection(mMaximumInteractionDistance, new_domain_size, new_local_rows);
}

template<unsigned SPACE_DIM>
bool NodesOnlyMesh<SPACE_DIM>::IsANodeCloseToDomainBoundary()
{
    assert(mpBoxCollection);

    int is_local_node_close = 0;
    c_vector<double, 2*SPACE_DIM> domain_boundary = mpBoxCollection->rGetDomainSize();

    // We ignore the x direction if the domain is periodic in x
    unsigned d0 = ( mpBoxCollection->GetIsPeriodicInX() ) ? 1 : 0;

    for (typename AbstractMesh<SPACE_DIM, SPACE_DIM>::NodeIterator node_iter = this->GetNodeIteratorBegin();
         node_iter != this->GetNodeIteratorEnd();
         ++node_iter)
    {
        // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
        c_vector<double, SPACE_DIM> location;
        location = node_iter->rGetLocation();

        for (unsigned d=d0; d<SPACE_DIM; d++)
        {
            if (location[d] < (domain_boundary[2*d] + mMinimumNodeDomainBoundarySeparation) ||  location[d] > (domain_boundary[2*d+1] - mMinimumNodeDomainBoundarySeparation))
            {
                is_local_node_close = 1;
                break;
            }
        }
        if (is_local_node_close)
        {
            break;  // Saves checking every node if we find one close to the boundary
        }
    }

    // Synchronise between processes
    int is_any_node_close = 0;
    MPI_Allreduce(&is_local_node_close, &is_any_node_close, 1, MPI_INT, MPI_SUM, PetscTools::GetWorld());

    return (is_any_node_close > 0);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ClearBoxCollection()
{
    if (mpBoxCollection)
    {
        delete mpBoxCollection;
    }
    mpBoxCollection = nullptr;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::SetInitialBoxCollection(const c_vector<double, 2*SPACE_DIM> domainSize, double maxInteractionDistance)
{
    this->SetUpBoxCollection(maxInteractionDistance, domainSize);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::SetUpBoxCollection(const std::vector<Node<SPACE_DIM>* >& rNodes)
{
    ClearBoxCollection();

    ChasteCuboid<SPACE_DIM> bounding_box = this->CalculateBoundingBox(rNodes);

    c_vector<double, 2*SPACE_DIM> domain_size;
    for (unsigned i=0; i < SPACE_DIM; i++)
    {
        domain_size[2*i] = bounding_box.rGetLowerCorner()[i] - 1e-14;
        domain_size[2*i+1] = bounding_box.rGetUpperCorner()[i] + 1e-14;
    }

    SetUpBoxCollection(mMaximumInteractionDistance, domain_size);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::SetUpBoxCollection(double cutOffLength, c_vector<double, 2*SPACE_DIM> domainSize, int numLocalRows, bool isPeriodic)
{
     ClearBoxCollection();

     mpBoxCollection = new DistributedBoxCollection<SPACE_DIM>(cutOffLength, domainSize, isPeriodic, numLocalRows);
     mpBoxCollection->SetupLocalBoxesHalfOnly();
     mpBoxCollection->SetCalculateNodeNeighbours(mCalculateNodeNeighbours);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::AddNodesToBoxes()
{
     // Put the nodes in the boxes.
     for (typename AbstractMesh<SPACE_DIM, SPACE_DIM>::NodeIterator node_iter = this->GetNodeIteratorBegin();
               node_iter != this->GetNodeIteratorEnd();
               ++node_iter)
     {
          unsigned box_index = mpBoxCollection->CalculateContainingBox(&(*node_iter));
          mpBoxCollection->rGetBox(box_index).AddNode(&(*node_iter));
     }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::AddHaloNodesToBoxes()
{
    // Add halo nodes
    for (typename std::vector<boost::shared_ptr<Node<SPACE_DIM> > >::iterator halo_node_iter = mHaloNodes.begin();
            halo_node_iter != mHaloNodes.end();
            ++halo_node_iter)
    {
        unsigned box_index = mpBoxCollection->CalculateContainingBox((*halo_node_iter).get());
        mpBoxCollection->rGetHaloBox(box_index).AddNode((*halo_node_iter).get());
    }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::UpdateBoxCollection()
{
    assert(mpBoxCollection);

    // Remove node pointers from boxes in BoxCollection.
    mpBoxCollection->EmptyBoxes();

    AddNodesToBoxes();

    mpBoxCollection->UpdateHaloBoxes();
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ResizeBoxCollection()
{
    if (!mpBoxCollection)
    {
        SetUpBoxCollection(this->mNodes);
    }

    while (IsANodeCloseToDomainBoundary())
    {
        EnlargeBoxCollection();
    }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::LoadBalanceMesh()
{
    std::vector<int> local_node_distribution = mpBoxCollection->CalculateNumberOfNodesInEachStrip();

    unsigned new_rows = mpBoxCollection->LoadBalance(local_node_distribution);

    c_vector<double, 2*SPACE_DIM> current_domain_size = mpBoxCollection->rGetDomainSize();

    // This ensures the domain will stay the same size.
    double fudge = 1e-14;
    for (unsigned d=0; d < SPACE_DIM; d++)
    {
        current_domain_size[2*d] = current_domain_size[2*d] + fudge;
        current_domain_size[2*d+1] = current_domain_size[2*d+1] - fudge;
    }
    SetUpBoxCollection(mMaximumInteractionDistance, current_domain_size, new_rows);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ConstructFromMeshReader(AbstractMeshReader<SPACE_DIM, SPACE_DIM>& rMeshReader)
{
    TetrahedralMesh<SPACE_DIM, SPACE_DIM>::ConstructFromMeshReader(rMeshReader);

    // Set the correct global node indices
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->SetIndex(GetNextAvailableIndex());
    }
}

template<unsigned SPACE_DIM>
std::vector<unsigned> NodesOnlyMesh<SPACE_DIM>::GetAllNodeIndices() const
{
    std::vector<unsigned> indices(GetNumNodes()); // GetNumNodes = mNodes - mDeletedNodes
    unsigned live_index=0;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        // Only use nodes which are not deleted
        if (!this->mNodes[i]->IsDeleted())
        {
            indices[live_index] = this->mNodes[i]->GetIndex();
            live_index++;
        }
    }
    return indices;
}

// Explicit instantiation
template class NodesOnlyMesh<1>;
template class NodesOnlyMesh<2>;
template class NodesOnlyMesh<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodesOnlyMesh)
