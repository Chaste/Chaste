/*

Copyright (c) 2005-2013, University of Oxford.
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
          mpBoxCollection(NULL),
          mTotalNumNodes(0u),
          mIndexCounter(0u)
{
}

template<unsigned SPACE_DIM>
NodesOnlyMesh<SPACE_DIM>::~NodesOnlyMesh()
{
    Clear();
    ClearBoxCollection();
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ConstructNodesWithoutMesh(const std::vector<Node<SPACE_DIM>*>& rNodes, double maxInteractionDistance, double domainPadding)
{
    assert(maxInteractionDistance > 0.0 && maxInteractionDistance < DBL_MAX);
    mMaximumInteractionDistance = maxInteractionDistance;

    mMinimumNodeDomainBoundarySeparation = mMaximumInteractionDistance;

    this->Clear();
    mpBoxCollection = NULL;

    ChasteCuboid<SPACE_DIM> bounding_box = this->CalculateBoundingBox(rNodes);

    ///\todo #2323 make a constructor for BoxCollection that accepts a ChasteCuboid.
    c_vector<double, 2*SPACE_DIM> domain_size;
    for (unsigned i=0; i < SPACE_DIM; i++)
    {
        domain_size[2*i] = bounding_box.rGetLowerCorner()[i];
        domain_size[2*i+1] = bounding_box.rGetUpperCorner()[i];
    }

    mpBoxCollection = new BoxCollection<SPACE_DIM>(mMaximumInteractionDistance, domain_size);

    for (unsigned i=0; i<rNodes.size(); i++)
    {
        assert(!rNodes[i]->IsDeleted());
        c_vector<double, SPACE_DIM> location = rNodes[i]->rGetLocation();

        Node<SPACE_DIM>* p_node_copy = new Node<SPACE_DIM>(i, location);
        p_node_copy->SetRadius(0.5);    // Default value.

        this->mNodes.push_back(p_node_copy);

        // Update the node map
        mNodesMapping[p_node_copy->GetIndex()] = this->mNodes.size()-1;

        mIndexCounter++;
    }

    mTotalNumNodes = rNodes.size();
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ConstructNodesWithoutMesh(const AbstractMesh<SPACE_DIM,SPACE_DIM>& rGeneratingMesh, double maxInteractionDistance, double domainPadding)
{
    ConstructNodesWithoutMesh(rGeneratingMesh.mNodes, maxInteractionDistance, domainPadding);
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
unsigned NodesOnlyMesh<SPACE_DIM>::GetNextAvailableIndex()
{
//    if(!this->mDeletedGlobalNodeIndices.empty())
//    {
//        unsigned index = this->mDeletedGlobalNodeIndices.back();
//        this->mDeletedGlobalNodeIndices.pop_back();
//        return index;
//    }
//    else
    {
        unsigned counter = mIndexCounter;
        mIndexCounter++;
        return counter * PetscTools::GetNumProcs() + PetscTools::GetMyRank();
    }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::Clear()
{
    // Call Clear() on the parent class
    MutableMesh<SPACE_DIM,SPACE_DIM>::Clear();

    // Clear the nodes mapping
    mNodesMapping.clear();

    // Set the global number of nodes to zero.
    mTotalNumNodes = 0u;

    // Clear the box collection.
    ClearBoxCollection();
}

template<unsigned SPACE_DIM>
unsigned NodesOnlyMesh<SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size() - this->mDeletedNodeIndices.size();
}

template<unsigned SPACE_DIM>
unsigned NodesOnlyMesh<SPACE_DIM>::GetGlobalNumNodes() const
{
    return mTotalNumNodes;
}

template<unsigned SPACE_DIM>
double NodesOnlyMesh<SPACE_DIM>::GetMaximumInteractionDistance()
{
    return mMaximumInteractionDistance;
}

template<unsigned SPACE_DIM>
BoxCollection<SPACE_DIM>* NodesOnlyMesh<SPACE_DIM>::GetBoxCollection()
{
    return mpBoxCollection;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ClearBoxCollection()
{
    if (mpBoxCollection != NULL)
    {
        delete mpBoxCollection;
    }
    mpBoxCollection = NULL;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::SetUpBoxCollection(double cutOffLength, c_vector<double, 2*SPACE_DIM> domainSize)
{
    if (mpBoxCollection)
    {
        delete mpBoxCollection;
    }

    mpBoxCollection = new BoxCollection<SPACE_DIM>(cutOffLength, domainSize);
    mpBoxCollection->SetupLocalBoxesHalfOnly();

    //Put the nodes in the boxes.
    for (typename AbstractMesh<SPACE_DIM, SPACE_DIM>::NodeIterator node_iter = this->GetNodeIteratorBegin();
            node_iter != this->GetNodeIteratorEnd();
            ++node_iter)
    {
        unsigned box_index = mpBoxCollection->CalculateContainingBox(&(*node_iter));
        mpBoxCollection->rGetBox(box_index).AddNode(&(*node_iter));
    }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::CalculateNodePairs(std::set<std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> >& rNodePairs, std::map<unsigned, std::set<unsigned> >& rNodeNeighbours)
{
    assert(mpBoxCollection != NULL);
    mpBoxCollection->CalculateNodePairs(this->mNodes, rNodePairs, rNodeNeighbours);
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ReMesh(NodeMap& map)
{
    // Store the node locations
    std::vector<c_vector<double, SPACE_DIM> > old_node_locations;
    std::vector<double> old_cell_radii;
    std::vector<bool> old_is_particle;

    unsigned new_index = 0;
    for (unsigned i=0; i<this->GetNumAllNodes(); i++)
    {
        if (this->mNodes[i]->IsDeleted())
        {
            map.SetDeleted(i);
        }
        else
        {
            map.SetNewIndex(i, new_index);
            old_node_locations.push_back(this->mNodes[i]->rGetLocation());
            old_cell_radii.push_back(this->mNodes[i]->GetRadius());
            old_is_particle.push_back(this->mNodes[i]->IsParticle());

            new_index++;
        }
    }
    // Remove current data
    this->Clear();

    // Construct the nodes
    for (unsigned node_index=0; node_index<old_node_locations.size(); node_index++)
    {
        Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index, old_node_locations[node_index], false);
        p_node->SetRadius(old_cell_radii[node_index]);
        if (old_is_particle[node_index])
        {
            p_node->SetIsParticle(true);
        }
        this->mNodes.push_back(p_node);

        mNodesMapping[p_node->GetIndex()] = this->mNodes.size()-1;
    }

    // Update the global number of nodes
    mTotalNumNodes = this->mNodes.size();
}

template<unsigned SPACE_DIM>
unsigned NodesOnlyMesh<SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNewNode)
{
    // Call method on parent class
    unsigned new_node_index = MutableMesh<SPACE_DIM, SPACE_DIM>::AddNode(pNewNode);

    // update mNodesMapping
    mNodesMapping[pNewNode->GetIndex()] = new_node_index;

    // Then update cell radius to default.
    pNewNode->SetRadius(0.5);

    // Increase the counter for the number of nodes.
    mTotalNumNodes += 1;

    return new_node_index;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::DeleteNode(unsigned index)
{
    if (this->mNodes[index]->IsDeleted())
    {
        EXCEPTION("Trying to delete a deleted node");
    }

    this->mNodes[index]->MarkAsDeleted();
    this->mDeletedNodeIndices.push_back(index);
//    mDeletedGlobalNodeIndices.push_back(index);
    mNodesMapping.erase(index);

    // Decrease the counter for the number of nodes.
    mTotalNumNodes -= 1;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::EnlargeBoxCollection()
{
    assert(mpBoxCollection);

    c_vector<double, 2*SPACE_DIM> current_domain_size = mpBoxCollection->rGetDomainSize();
    c_vector<double, 2*SPACE_DIM> new_domain_size;

    for (unsigned d=0; d < SPACE_DIM; d++)
    {
        new_domain_size[2*d] = current_domain_size[2*d] - mMaximumInteractionDistance;
        new_domain_size[2*d+1] = current_domain_size[2*d+1] + mMaximumInteractionDistance;
    }

    SetUpBoxCollection(mMaximumInteractionDistance, new_domain_size);
}

template<unsigned SPACE_DIM>
bool NodesOnlyMesh<SPACE_DIM>::IsANodeCloseToDomainBoundary()
{
    assert (mpBoxCollection);

    bool is_any_node_close = false;
    c_vector<double, 2*SPACE_DIM> domain_boundary = mpBoxCollection->rGetDomainSize();

    for (typename AbstractMesh<SPACE_DIM, SPACE_DIM>::NodeIterator node_iter = this->GetNodeIteratorBegin();
            node_iter != this->GetNodeIteratorEnd();
            ++node_iter)
    {
        c_vector<double, SPACE_DIM> location = node_iter->rGetLocation();

        for (unsigned d=0; d < SPACE_DIM; d++)
        {
            if (location[d] < (domain_boundary[2*d] + mMinimumNodeDomainBoundarySeparation) ||  location[d] > (domain_boundary[2*d+1] - mMinimumNodeDomainBoundarySeparation))
            {
                is_any_node_close = true;
                break;
            }
        }
        if (is_any_node_close)
        {
            break;  // Saves checking every node if we find one close to the boundary.
        }
    }

    return is_any_node_close;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::SetMinimumNodeDomainBoundarySeparation(double separation)
{
    assert(!(separation < 0.0));

    mMinimumNodeDomainBoundarySeparation = separation;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class NodesOnlyMesh<1>;
template class NodesOnlyMesh<2>;
template class NodesOnlyMesh<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodesOnlyMesh)
