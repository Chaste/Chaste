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
		  mpBoxCollection(NULL)
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

        mCellRadii.push_back(0.5);
    }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ConstructNodesWithoutMesh(const AbstractMesh<SPACE_DIM,SPACE_DIM>& rGeneratingMesh)
{
    ConstructNodesWithoutMesh(rGeneratingMesh.mNodes);
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
    assert(index < mCellRadii.size());
    return mCellRadii[index];
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::SetCellRadius(unsigned index, double radius)
{
    assert(index < mCellRadii.size());
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
	if(mpBoxCollection != NULL)
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
    for (unsigned i=0; i< this->GetNumNodes(); i++)
    {
        unsigned box_index = mpBoxCollection->CalculateContainingBox(this->GetNode(i));
        mpBoxCollection->rGetBox(box_index).AddNode(this->GetNode(i));
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
    mCellRadii = old_cell_radii;

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
    if (new_node_index >= mCellRadii.size())
    {
        mCellRadii.resize(new_node_index+1);
    }
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
    mCellRadii[index] = DOUBLE_UNSET;
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
