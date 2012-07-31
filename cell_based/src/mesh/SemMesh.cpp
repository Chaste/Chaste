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

#include "SemMesh.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"
#include <list>

template<unsigned DIM>
SemMesh<DIM>::SemMesh( std::vector<Node<DIM>*> nodes,
                        std::vector<PottsElement<DIM>*> pottsElements)
{
    // SEM model only defined in 2 and 3D
    assert(DIM > 1);

    // Reset member variables and clear mNodes, mElements.
    Clear();

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
SemMesh<DIM>::SemMesh()
{
    assert(DIM > 1);
    this->mMeshChangesDuringSimulation = true;
    Clear();
}

template<unsigned DIM>
SemMesh<DIM>::~SemMesh()
{
    Clear();
}

template<unsigned DIM>
unsigned SemMesh<DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size());
    return index;
}

template<unsigned DIM>
unsigned SemMesh<DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size());
    return index;
}

template<unsigned DIM>
unsigned SemMesh<DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    return index;
}

template<unsigned DIM>
void SemMesh<DIM>::Clear()
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
    mDeletedNodeIndices.clear();
}

template<unsigned DIM>
void SemMesh<DIM>::ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rMeshReader)
{
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
        unsigned is_boundary_node = (unsigned) node_data[DIM];
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
         unsigned attribute_value = element_data.AttributeValue;
         p_element->SetAttribute(attribute_value);
     }
    }
}

template<unsigned DIM>
unsigned SemMesh<DIM>::GetNumNodes() const
{
    return this->mNodes.size() - mDeletedNodeIndices.size();
}

template<unsigned DIM>
unsigned SemMesh<DIM>::GetNumAllNodes() const
{
    return this->mNodes.size();
}

template<unsigned DIM>
unsigned SemMesh<DIM>::GetNumElements() const
{
    return mElements.size() - mDeletedElementIndices.size();
}

template<unsigned DIM>
unsigned SemMesh<DIM>::GetNumAllElements() const
{
    return mElements.size();
}

template<unsigned DIM>
PottsElement<DIM>* SemMesh<DIM>::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}

template<unsigned DIM>
void SemMesh<DIM>::ReMesh()
{
    std::vector<unsigned> old_element_indices;
    std::vector< std::vector<unsigned> > old_element_node_indices;
    std::vector< std::vector<c_vector<double, DIM> > > old_element_node_location;

    // For each element either note that it was deleted, or store the information needed to re-construct it, i.e. index and nodes.
    for(typename std::vector<PottsElement<DIM>* >::iterator element_iterator = this->mElements.begin();
            element_iterator != this->mElements.end();
            ++element_iterator)
    {
        if((*element_iterator)->IsDeleted())
        {
            // Marked as deleted - don't store the nodes.
        }
        else
        {
            // Store the index and node location.
            old_element_indices.push_back( (*element_iterator)->GetIndex() );

            std::vector<unsigned> local_node_indices;
            std::vector< c_vector<double, DIM> > local_node_locations;
            for (unsigned local_index=0; local_index < (*element_iterator)->GetNumNodes(); local_index++)
            {
                Node<DIM>* p_node = (*element_iterator)->GetNode(local_index);
                local_node_indices.push_back(p_node->GetIndex());
                local_node_locations.push_back(p_node->rGetLocation());
            }
            old_element_node_indices.push_back(local_node_indices);
            old_element_node_location.push_back(local_node_locations);
        }
    }

    // Clear all the local information.
    Clear();

    // Reconstruct the elements
    for (unsigned index=0; index < old_element_indices.size(); index++)
    {
        std::vector<Node<DIM>* > new_element_nodes;
        for(unsigned node=0; node < old_element_node_indices[index].size(); node++)
        {
            c_vector<double, DIM> new_location = old_element_node_location[index][node];
            unsigned new_index = old_element_node_indices[index][node];
            Node<DIM>* p_new_node = new Node<DIM>(new_index, new_location);
            new_element_nodes.push_back(p_new_node);
            this->mNodes.push_back(p_new_node);
        }

        this->mElements.push_back(new PottsElement<DIM>(old_element_indices[index], new_element_nodes));
    }

    // Register elements with nodes
    for (unsigned index=0; index < this->mElements.size(); index++)
    {
        PottsElement<DIM>* p_element = this->mElements[index];

        unsigned element_index = p_element->GetIndex();
        unsigned num_nodes_in_element = p_element->GetNumNodes();

        for (unsigned node_index=0; node_index<num_nodes_in_element; node_index++)
        {
            p_element->GetNode(node_index)->AddElement(element_index);
        }
    }
}

template<unsigned DIM>
void SemMesh<DIM>::DeleteElement(unsigned index)
{
    // Delete the nodes.
    for(unsigned node_index = 0; node_index < this->mElements[index]->GetNumNodes(); node_index++)
    {
        Node<DIM>* p_node = this->mElements[index]->GetNode(node_index);
        p_node->MarkAsDeleted();
        mDeletedNodeIndices.push_back(node_index);
    }

    // Mark this element as deleted; this also updates the nodes containing element indices
    this->mElements[index]->MarkAsDeleted();
    mDeletedElementIndices.push_back(index);
}

template<unsigned DIM>
unsigned SemMesh<DIM>::AddElement(PottsElement<DIM>* pNewElement, std::vector<Node<DIM>* > newNodes)
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

    // Add the nodes
    for(unsigned index=0; index<newNodes.size(); index++)
    {
        // Set index of the node correctly
        newNodes[index]->SetIndex(this->mNodes.size());
        this->mNodes.push_back(newNodes[index]);
    }

    pNewElement->RegisterWithNodes();
    return pNewElement->GetIndex();
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class SemMesh<1>;
template class SemMesh<2>;
template class SemMesh<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemMesh)
