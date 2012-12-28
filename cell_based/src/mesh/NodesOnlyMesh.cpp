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

#include <map>
#include "NodesOnlyMesh.hpp"

template<unsigned SPACE_DIM>
NodesOnlyMesh<SPACE_DIM>::NodesOnlyMesh()
        : MutableMesh<SPACE_DIM, SPACE_DIM>(),
          mpBoxCollection(NULL),
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
void NodesOnlyMesh<SPACE_DIM>::ConstructNodesWithoutMesh(const std::vector<Node<SPACE_DIM>*>& rNodes)
{
    this->Clear();
    mpBoxCollection = NULL;

    for (unsigned i=0; i<rNodes.size(); i++)
    {
        assert(!rNodes[i]->IsDeleted());
        c_vector<double, SPACE_DIM> location = rNodes[i]->rGetLocation();

        Node<SPACE_DIM>* p_node_copy = new Node<SPACE_DIM>(i, location);
        this->mNodes.push_back(p_node_copy);

        mCellRadii[i] = 0.5;

        mIndexCounter++;
    }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ConstructNodesWithoutMesh(const AbstractMesh<SPACE_DIM,SPACE_DIM>& rGeneratingMesh)
{
    ConstructNodesWithoutMesh(rGeneratingMesh.mNodes);
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

    // Clear the cell radii
    mCellRadii.clear();
}

template<unsigned SPACE_DIM>
double NodesOnlyMesh<SPACE_DIM>::GetCellRadius(unsigned index)
{
    std::map<unsigned, double>::iterator local =  mCellRadii.find(index);
    std::map<unsigned, double>::const_iterator halo = mHaloCellRadii.find(index);

    if(local != mCellRadii.end())
    {
        return local->second;
    }
    else if(halo != mHaloCellRadii.end())
    {
        return halo->second;
    }
    else
    {
        EXCEPTION( "Requested radius of a node which is not set. Either does not lie on this process as a node or halo node, or has not been set.");
    }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::SetCellRadius(unsigned index, double radius)
{
//    // Make sure that we own the node
//    if(mNodesMapping.find(index) == mNodesMapping.end())
//    {
//        #define COVERAGE_IGNORE
//        EXCEPTION("Trying to set the radius of node which does not lie on this process");
//        #undef COVERAGE_IGNORE
//    }

    // Set the radius
    mCellRadii[index] = radius;
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
void NodesOnlyMesh<SPACE_DIM>::SetMaximumInteractionDistance(double maximumInteractionDistance)
{
    mMaximumInteractionDistance = maximumInteractionDistance;
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
    bool copy_radii = !mCellRadii.empty();

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
            if (copy_radii)
            {
                old_cell_radii.push_back(mCellRadii[i]);
            }

            new_index++;
        }
    }
    // Remove current data
    this->Clear();

    // Replace radius data
    for (unsigned i=0; i<old_cell_radii.size(); i++)
    {
        mCellRadii[i] = old_cell_radii[i];
    }

    // Construct the nodes
    for (unsigned node_index=0; node_index<old_node_locations.size(); node_index++)
    {
        Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index, old_node_locations[node_index], false);
        this->mNodes.push_back(p_node);
    }
}

template<unsigned SPACE_DIM>
unsigned NodesOnlyMesh<SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNewNode)
{
    // Call method on parent class
    unsigned new_node_index = MutableMesh<SPACE_DIM, SPACE_DIM>::AddNode(pNewNode);

    // Then update mCellRadii
    SetCellRadius(new_node_index, 0.5);

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

    /**
     * Note: we may not need to update mCellRadii here, since if the
     * node index is ever re-used when a new node is added, mCellRadii
     * will be updated correctly.
     */
    mCellRadii.erase(index);
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
