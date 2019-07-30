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

#include "PottsMesh.hpp"
#include "RandomNumberGenerator.hpp"


template<unsigned DIM>
PottsMesh<DIM>::PottsMesh(std::vector<Node<DIM>*> nodes,
                          std::vector<PottsElement<DIM>*> pottsElements,
                          std::vector<std::set<unsigned> > vonNeumannNeighbouringNodeIndices,
                          std::vector<std::set<unsigned> > mooreNeighbouringNodeIndices)
{
    // Reset member variables and clear mNodes, mElements.
    Clear();

    // Verify the same size of nodes and neighbour information.
    if ((vonNeumannNeighbouringNodeIndices.size() != nodes.size()) || (mooreNeighbouringNodeIndices.size() != nodes.size()))
    {
        EXCEPTION("Nodes and neighbour information for a Potts mesh need to be the same length.");
    }
    mVonNeumannNeighbouringNodeIndices = vonNeumannNeighbouringNodeIndices;
    mMooreNeighbouringNodeIndices = mooreNeighbouringNodeIndices;

    // Populate mNodes and mElements
    for (unsigned node_index=0; node_index<nodes.size(); node_index++)
    {
        Node<DIM>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }
    for (unsigned elem_index=0; elem_index<pottsElements.size(); elem_index++)
    {
        PottsElement<DIM>* p_temp_element = pottsElements[elem_index];
        mElements.push_back(p_temp_element);
    }

    // Register elements with nodes
    for (unsigned index=0; index<mElements.size(); index++)
    {
        PottsElement<DIM>* p_element = mElements[index];

        unsigned element_index = p_element->GetIndex();
        unsigned num_nodes_in_element = p_element->GetNumNodes();

        for (unsigned node_index=0; node_index<num_nodes_in_element; node_index++)
        {
            p_element->GetNode(node_index)->AddElement(element_index);
        }
    }

    this->mMeshChangesDuringSimulation = true;
}

template<unsigned DIM>
PottsMesh<DIM>::PottsMesh()
{
    this->mMeshChangesDuringSimulation = true;
    Clear();
}

template<unsigned DIM>
PottsMesh<DIM>::~PottsMesh()
{
    Clear();
}

template<unsigned DIM>
unsigned PottsMesh<DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size());
    return index;
}

template<unsigned DIM>
unsigned PottsMesh<DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size());
    return index;
}

template<unsigned DIM>
unsigned PottsMesh<DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    return index;
}

template<unsigned DIM>
void PottsMesh<DIM>::Clear()
{
    // Delete elements
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    mElements.clear();

    // Delete nodes
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }
    this->mNodes.clear();

    mDeletedElementIndices.clear();

    // Delete neighbour info
    //mVonNeumannNeighbouringNodeIndices.clear();
    //mMooreNeighbouringNodeIndices.clear();
}

template<unsigned DIM>
unsigned PottsMesh<DIM>::GetNumNodes() const
{
    return this->mNodes.size();
}

template<unsigned DIM>
unsigned PottsMesh<DIM>::GetNumElements() const
{
    return mElements.size() - mDeletedElementIndices.size();
}

template<unsigned DIM>
unsigned PottsMesh<DIM>::GetNumAllElements() const
{
    return mElements.size();
}

template<unsigned DIM>
PottsElement<DIM>* PottsMesh<DIM>::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}

template<unsigned DIM>
c_vector<double, DIM> PottsMesh<DIM>::GetCentroidOfElement(unsigned index)
{
    PottsElement<DIM>* p_element = GetElement(index);
    unsigned num_nodes_in_element = p_element->GetNumNodes();

    ///\todo This should probably be returning the nearest node
    c_vector<double, DIM> centroid = zero_vector<double>(DIM);

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        // Find location of current node and add it to the centroid
        centroid += p_element->GetNodeLocation(local_index);
    }

    centroid /= num_nodes_in_element;

    return centroid;
}

template<unsigned DIM>
double PottsMesh<DIM>::GetVolumeOfElement(unsigned index)
{
    PottsElement<DIM>* p_element = GetElement(index);
    double element_volume = (double) p_element->GetNumNodes();

    return element_volume;
}

template<unsigned DIM>
double PottsMesh<DIM>::GetSurfaceAreaOfElement(unsigned index)
{
    ///\todo not implemented in 3d yet
    assert(DIM==2 || DIM==3); // LCOV_EXCL_LINE

    // Helper variables
    PottsElement<DIM>* p_element = GetElement(index);
    unsigned num_nodes = p_element->GetNumNodes();

    double surface_area = 0.0;
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        std::set<unsigned> neighbouring_node_indices = GetVonNeumannNeighbouringNodeIndices(p_element->GetNode(node_index)->GetIndex());
        unsigned local_edges = 2*DIM;
        for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
             iter != neighbouring_node_indices.end();
             iter++)
        {
            std::set<unsigned> neighbouring_node_element_indices = this->mNodes[*iter]->rGetContainingElementIndices();

            if (!(neighbouring_node_element_indices.empty()) && (local_edges!=0))
            {
                unsigned neighbouring_node_element_index = *(neighbouring_node_element_indices.begin());
                if (neighbouring_node_element_index == index)
                {
                    local_edges--;
                }
            }
        }
        surface_area += local_edges;
    }
    return surface_area;
}

template<unsigned DIM>
std::set<unsigned> PottsMesh<DIM>::GetMooreNeighbouringNodeIndices(unsigned nodeIndex)
{
    return mMooreNeighbouringNodeIndices[nodeIndex];
}

template<unsigned DIM>
std::set<unsigned> PottsMesh<DIM>::GetVonNeumannNeighbouringNodeIndices(unsigned nodeIndex)
{
    return mVonNeumannNeighbouringNodeIndices[nodeIndex];
}

template<unsigned DIM>
void PottsMesh<DIM>::DeleteElement(unsigned index)
{
    // Mark this element as deleted; this also updates the nodes containing element indices
    this->mElements[index]->MarkAsDeleted();
    mDeletedElementIndices.push_back(index);
}

template<unsigned DIM>
void PottsMesh<DIM>::RemoveDeletedElements()
{
    // Remove any elements that have been removed and re-order the remaining ones
    unsigned num_deleted_elements = mDeletedElementIndices.size();

    for (unsigned index = num_deleted_elements; index>0; index--)
    {
        unsigned deleted_elem_index = mDeletedElementIndices[index-1];
        delete mElements[deleted_elem_index];
        mElements.erase(mElements.begin()+deleted_elem_index);
        for (unsigned elem_index=deleted_elem_index; elem_index<mElements.size(); elem_index++)
        {
            mElements[elem_index]->ResetIndex(elem_index);
        }
    }
    mDeletedElementIndices.clear();
}

template<unsigned DIM>
void PottsMesh<DIM>::DeleteNode(unsigned index)
{
    //Mark node as deleted so we don't consider it when iterating over nodes
    this->mNodes[index]->MarkAsDeleted();

    //Remove from Elements
    std::set<unsigned> containing_element_indices = this->mNodes[index]->rGetContainingElementIndices();

    for (std::set<unsigned>::iterator iter = containing_element_indices.begin();
         iter != containing_element_indices.end();
         iter++)
    {
        assert(mElements[*iter]->GetNumNodes() > 0);
        if (mElements[*iter]->GetNumNodes() == 1)
        {
            DeleteElement(*iter);
        }
        else
        {
            this->mElements[*iter]->DeleteNode(this->mElements[*iter]->GetNodeLocalIndex(index));
        }
    }

    // Remove from connectivity
    mVonNeumannNeighbouringNodeIndices[index].clear();
    mMooreNeighbouringNodeIndices[index].clear();

    assert(mVonNeumannNeighbouringNodeIndices.size()==mMooreNeighbouringNodeIndices.size());
    for (unsigned node_index = 0;
         node_index < mVonNeumannNeighbouringNodeIndices.size();
         node_index++)
    {
        // Remove node "index" from the Von Neuman neighbourhood of node "node_index".
        mVonNeumannNeighbouringNodeIndices[node_index].erase(index);
        mMooreNeighbouringNodeIndices[node_index].erase(index);

        // Check there's still connectivity for the other non-deleted nodes
        if (!this->mNodes[node_index]->IsDeleted())
        {
            assert(!mVonNeumannNeighbouringNodeIndices[node_index].empty());
            assert(!mMooreNeighbouringNodeIndices[node_index].empty());
        }
    }

    // Remove node from mNodes and renumber all the elements and nodes
    delete this->mNodes[index];
    this->mNodes.erase(this->mNodes.begin()+index);
    unsigned num_nodes = GetNumNodes();
    mVonNeumannNeighbouringNodeIndices.erase(mVonNeumannNeighbouringNodeIndices.begin()+index);
    mMooreNeighbouringNodeIndices.erase(mMooreNeighbouringNodeIndices.begin()+index);

    assert(mVonNeumannNeighbouringNodeIndices.size()==num_nodes);
    assert(mMooreNeighbouringNodeIndices.size()==num_nodes);

    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        // Reduce the index of all nodes greater than  node "index"
        if (node_index >= index)
        {
            assert(this->mNodes[node_index]->GetIndex() == node_index+1);
            this->mNodes[node_index]->SetIndex(node_index);
        }
        assert(this->mNodes[node_index]->GetIndex() == node_index);

        // Reduce the index of all nodes greater than  node "index"
        // in the Moore and Von Neuman neighbourhoods.
        std::set<unsigned> von_neuman = mVonNeumannNeighbouringNodeIndices[node_index];
        mVonNeumannNeighbouringNodeIndices[node_index].clear();
        for (std::set<unsigned>::iterator iter = von_neuman.begin();
             iter != von_neuman.end();
             iter++)
        {
            if (*iter >= index)
            {
                mVonNeumannNeighbouringNodeIndices[node_index].insert(*iter-1);
            }
            else
            {
                mVonNeumannNeighbouringNodeIndices[node_index].insert(*iter);
            }
        }
        std::set<unsigned> moore = mMooreNeighbouringNodeIndices[node_index];
        mMooreNeighbouringNodeIndices[node_index].clear();
        for (std::set<unsigned>::iterator iter = moore.begin();
             iter != moore.end();
             iter++)
        {
            if (*iter >= index)
            {
                mMooreNeighbouringNodeIndices[node_index].insert(*iter-1);
            }
            else
            {
                mMooreNeighbouringNodeIndices[node_index].insert(*iter);
            }
        }
    }
    // Finally remove any elememts that have been removed
    assert(mDeletedElementIndices.size() <= 1); // Should have at most one element to remove
    if (mDeletedElementIndices.size() == 1)
    {
        unsigned deleted_elem_index = mDeletedElementIndices[0];
        delete mElements[deleted_elem_index];
        mElements.erase(mElements.begin()+deleted_elem_index);
        mDeletedElementIndices.clear();

        for (unsigned elem_index=deleted_elem_index; elem_index<GetNumElements(); elem_index++)
        {
            mElements[elem_index]->ResetIndex(elem_index);
        }
    }
}

template<unsigned DIM>
unsigned PottsMesh<DIM>::DivideElement(PottsElement<DIM>* pElement,
                                       bool placeOriginalElementBelow)
{
    /// Not implemented in 1d
    assert(DIM==2 || DIM==3); // LCOV_EXCL_LINE

    // Store the number of nodes in the element (this changes when nodes are deleted from the element)
    unsigned num_nodes = pElement->GetNumNodes();

    if (num_nodes < 2)
    {
        EXCEPTION("Tried to divide a Potts element with only one node. Cell dividing too often given dynamic parameters.");
    }

    // Copy the nodes in this element
    std::vector<Node<DIM>*> nodes_elem;
    for (unsigned i=0; i<num_nodes; i++)
    {
        nodes_elem.push_back(pElement->GetNode(i));
    }

    // Get the index of the new element
    unsigned new_element_index;
    if (mDeletedElementIndices.empty())
    {
        new_element_index = this->mElements.size();
    }
    else
    {
        new_element_index = mDeletedElementIndices.back();
        mDeletedElementIndices.pop_back();
        delete this->mElements[new_element_index];
    }

    // Add the new element to the mesh
    AddElement(new PottsElement<DIM>(new_element_index, nodes_elem));

    /**
     * Remove the correct nodes from each element. If placeOriginalElementBelow is true,
     * place the original element below (in the y direction or z in 3d) the new element; otherwise,
     * place it above.
     */
    unsigned half_num_nodes = num_nodes/2; // This will round down
    assert(half_num_nodes > 0);
    assert(half_num_nodes < num_nodes);

    // Find lowest element
    ///\todo this could be more efficient
    double height_midpoint_1 = 0.0;
    double height_midpoint_2 = 0.0;
    unsigned counter_1 = 0;
    unsigned counter_2 = 0;

    for (unsigned i=0; i<num_nodes; i++)
    {
        if (i<half_num_nodes)
        {
            height_midpoint_1 += pElement->GetNode(i)->rGetLocation()[DIM - 1];
            counter_1++;
        }
        else
        {
            height_midpoint_2 += pElement->GetNode(i)->rGetLocation()[DIM -1];
            counter_2++;
        }
    }
    height_midpoint_1 /= (double)counter_1;
    height_midpoint_2 /= (double)counter_2;

    for (unsigned i=num_nodes; i>0; i--)
    {
        if (i-1 >= half_num_nodes)
        {
            if (height_midpoint_1 < height_midpoint_2)
            {
                if (placeOriginalElementBelow)
                {
                    pElement->DeleteNode(i-1);
                }
                else
                {
                    this->mElements[new_element_index]->DeleteNode(i-1);
                }
            }
            else
            {
                if (placeOriginalElementBelow)
                {
                    this->mElements[new_element_index]->DeleteNode(i-1);
                }
                else
                {
                    pElement->DeleteNode(i-1);
                }
            }
        }
        else // i-1 < half_num_nodes
        {
            if (height_midpoint_1 < height_midpoint_2)
            {
                if (placeOriginalElementBelow)
                {
                    this->mElements[new_element_index]->DeleteNode(i-1);
                }
                else
                {
                    pElement->DeleteNode(i-1);
                }
            }
            else
            {
                if (placeOriginalElementBelow)
                {
                    pElement->DeleteNode(i-1);
                }
                else
                {
                    this->mElements[new_element_index]->DeleteNode(i-1);
                }
            }
        }
    }

    return new_element_index;
}

template<unsigned DIM>
unsigned PottsMesh<DIM>::AddElement(PottsElement<DIM>* pNewElement)
{
    unsigned new_element_index = pNewElement->GetIndex();

    if (new_element_index == this->mElements.size())
    {
        this->mElements.push_back(pNewElement);
    }
    else
    {
        this->mElements[new_element_index] = pNewElement;
    }
    pNewElement->RegisterWithNodes();
    return pNewElement->GetIndex();
}

template<unsigned DIM>
std::set<unsigned> PottsMesh<DIM>::GetNeighbouringElementIndices(unsigned elementIndex)
{
    // Helper variables
    PottsElement<DIM>* p_element = this->GetElement(elementIndex);
    unsigned num_nodes = p_element->GetNumNodes();

    // Create a set of neighbouring element indices
    std::set<unsigned> neighbouring_element_indices;

    // Loop over nodes owned by this element
    for (unsigned local_index=0; local_index<num_nodes; local_index++)
    {
        // Get a pointer to this node
        Node<DIM>* p_node = p_element->GetNode(local_index);

        // Find the indices of the elements owned by neighbours of this node

        // Loop over neighbouring nodes. Only want Von Neuman neighbours (i.e N,S,E,W) as need to share an edge
        std::set<unsigned> neighbouring_node_indices = GetVonNeumannNeighbouringNodeIndices(p_node->GetIndex());

         // Iterate over these neighbouring nodes
         for (std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
              neighbour_iter != neighbouring_node_indices.end();
              ++neighbour_iter)
         {
             std::set<unsigned> neighbouring_node_containing_elem_indices = this->GetNode(*neighbour_iter)->rGetContainingElementIndices();

             assert(neighbouring_node_containing_elem_indices.size()<2); // Either in element or in medium

             if (neighbouring_node_containing_elem_indices.size()==1) // Node is in an element
             {
                 // Add this element to the neighbouring elements set
                 neighbouring_element_indices.insert(*(neighbouring_node_containing_elem_indices.begin()));
             }
         }
    }

    // Lastly remove this element's index from the set of neighbouring element indices
    neighbouring_element_indices.erase(elementIndex);

    return neighbouring_element_indices;
}

template<unsigned DIM>
void PottsMesh<DIM>::ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rMeshReader)
{
    assert(rMeshReader.HasNodePermutation() == false);

    // Store numbers of nodes and elements
    unsigned num_nodes = rMeshReader.GetNumNodes();
    unsigned num_elements = rMeshReader.GetNumElements();

    // Reserve memory for nodes
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    // Add nodes
    std::vector<double> node_data;
    for (unsigned i=0; i<num_nodes; i++)
    {
        node_data = rMeshReader.GetNextNode();
        unsigned is_boundary_node = (bool) node_data[DIM];
        node_data.pop_back();
        this->mNodes.push_back(new Node<DIM>(i, node_data, is_boundary_node));
    }

    rMeshReader.Reset();

    // Reserve memory for nodes
    mElements.reserve(rMeshReader.GetNumElements());

    // Add elements
    for (unsigned elem_index=0; elem_index<num_elements; elem_index++)
    {
        // Get the data for this element
        ElementData element_data = rMeshReader.GetNextElementData();

        // Get the nodes owned by this element
        std::vector<Node<DIM>*> nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j=0; j<num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        // Use nodes and index to construct this element
        PottsElement<DIM>* p_element = new PottsElement<DIM>(elem_index, nodes);
        mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            double attribute_value = element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }

    // If we are just using a mesh reader, then there is no neighbour information (see #1932)
    if (mVonNeumannNeighbouringNodeIndices.empty())
    {
        mVonNeumannNeighbouringNodeIndices.resize(num_nodes);
    }
    if (mMooreNeighbouringNodeIndices.empty())
    {
        mMooreNeighbouringNodeIndices.resize(num_nodes);
    }
}

// Explicit instantiation
template class PottsMesh<1>;
template class PottsMesh<2>;
template class PottsMesh<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PottsMesh)

