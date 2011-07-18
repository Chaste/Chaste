/*

Copyright (C) University of Oxford, 2005-2011

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

#include "PottsMesh.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"
#include <list>

PottsMesh::PottsMesh(std::vector<Node<2>*> nodes, std::vector<PottsElement<2>*> pottsElements)
{
    // Reset member variables and clear mNodes and mElements
    Clear();

    // Populate mNodes and mElements
    for (unsigned node_index=0; node_index<nodes.size(); node_index++)
    {
        Node<2>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }
    for (unsigned elem_index=0; elem_index<pottsElements.size(); elem_index++)
    {
        PottsElement<2>* p_temp_element = pottsElements[elem_index];
        mElements.push_back(p_temp_element);
    }

    // Register elements with nodes
    for (unsigned index=0; index<mElements.size(); index++)
    {
        PottsElement<2>* p_element = mElements[index];

        unsigned element_index = p_element->GetIndex();
        unsigned num_nodes_in_element = p_element->GetNumNodes();

        for (unsigned node_index=0; node_index<num_nodes_in_element; node_index++)
        {
            p_element->GetNode(node_index)->AddElement(element_index);
        }
    }

    this->mMeshChangesDuringSimulation = true;
}

PottsMesh::PottsMesh()
{
    this->mMeshChangesDuringSimulation = true;
    Clear();
}

PottsMesh::~PottsMesh()
{
    Clear();
}

unsigned PottsMesh::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size());
    return index;
}

unsigned PottsMesh::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size());
    return index;
}

unsigned PottsMesh::SolveBoundaryElementMapping(unsigned index) const
{
    return index;
}

void PottsMesh::Clear()
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
}

unsigned PottsMesh::GetNumNodes() const
{
    return this->mNodes.size();
}

unsigned PottsMesh::GetNumElements() const
{
    return mElements.size() - mDeletedElementIndices.size();
}

unsigned PottsMesh::GetNumAllElements() const
{
    return mElements.size();
}

PottsElement<2>* PottsMesh::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}

c_vector<double, 2> PottsMesh::GetCentroidOfElement(unsigned index)
{
    PottsElement<2>* p_element = GetElement(index);
    unsigned num_nodes_in_element = p_element->GetNumNodes();

    ///\todo This should probably be returning the nearest node
    c_vector<double, 2> centroid = zero_vector<double>(2);

    double temp_centroid_x = 0;
    double temp_centroid_y = 0;

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        // Find location of current node
        c_vector<double, 2> current_node = p_element->GetNodeLocation(local_index);

        temp_centroid_x += current_node[0];
        temp_centroid_y += current_node[1];
    }

    centroid(0) = temp_centroid_x/num_nodes_in_element;
    centroid(1) = temp_centroid_y/num_nodes_in_element;

    return centroid;
}

///// \cond Get Doxygen to ignore, since it's confused by these templates
//template<>
//void PottsMesh<2,2>::ConstructFromMeshReader(AbstractMeshReader<2,2>& rMeshReader)
///// \endcond Get Doxygen to ignore, since it's confused by these templates
//{
//    // Store numbers of nodes and elements
//    unsigned num_nodes = rMeshReader.GetNumNodes();
//    unsigned num_elements = rMeshReader.GetNumElements();
//
//    // Reserve memory for nodes
//    this->mNodes.reserve(num_nodes);
//
//    rMeshReader.Reset();
//
//    // Add nodes
//    std::vector<double> node_data;
//    for (unsigned i=0; i<num_nodes; i++)
//    {
//        node_data = rMeshReader.GetNextNode();
//        unsigned is_boundary_node = (unsigned) node_data[2];
//        node_data.pop_back();
//        this->mNodes.push_back(new Node<2>(i, node_data, is_boundary_node));
//    }
//
//    rMeshReader.Reset();
//
//    // Reserve memory for nodes
//    mElements.reserve(rMeshReader.GetNumElements());
//
//    // Add elements
//    for (unsigned elem_index=0; elem_index<num_elements; elem_index++)
//    {
//        // Get the data for this element
//        ElementData element_data = rMeshReader.GetNextElementData();
//
//        // Get the nodes owned by this element
//        std::vector<Node<2>*> nodes;
//        unsigned num_nodes_in_element = element_data.NodeIndices.size();
//        for (unsigned j=0; j<num_nodes_in_element; j++)
//        {
//            assert(element_data.NodeIndices[j] < this->mNodes.size());
//            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
//        }
//
//        // Use nodes and index to construct this element
//        PottsElement<2><2,2>* p_element = new PottsElement<2><2,2>(elem_index, nodes);
//        mElements.push_back(p_element);
//
//        if (rMeshReader.GetNumElementAttributes() > 0)
//        {
//            assert(rMeshReader.GetNumElementAttributes() == 1);
//            unsigned attribute_value = element_data.AttributeValue;
//            p_element->SetRegion(attribute_value);
//        }
//    }
//}

c_vector<double, 2> PottsMesh::GetVectorFromAtoB(const c_vector<double, 2>& rLocationA, const c_vector<double, 2>& rLocationB)
{
    c_vector<double, 2> vector;

    vector = AbstractMesh<2, 2>::GetVectorFromAtoB(rLocationA, rLocationB);

    return vector;
}

double PottsMesh::GetVolumeOfElement(unsigned index)
{
    // Get pointer to this element
    PottsElement<2>* p_element = GetElement(index);

    double element_volume = (double) p_element->GetNumNodes();

    return element_volume;
}

double PottsMesh::GetSurfaceAreaOfElement(unsigned index)
{
    //\todo this is not correct need to work this out from the number of free boundaries. See #1683.
    double surface_area = 0.0;

    return surface_area;
}

std::set<unsigned> PottsMesh::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
    // Create a set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices;

    double width = this->GetWidth(0);
    unsigned nodes_across = (unsigned)width + 1; // Getting the number of nodes along x axis of mesh

    double height = this->GetWidth(1);
    unsigned nodes_up = (unsigned)height + 1;

    /*
     * This stores the available neighbours using the following numbering:
     *
     *  1-----0-----7
     *  |     |     |
     *  |     |     |
     *  2-----x-----6
     *  |     |     |
     *  |     |     |
     *  3-----4-----5
     */

    // Create a vector of possible neighbouring node indices
    std::vector<unsigned> neighbour_indices_vector(8, nodeIndex);
    neighbour_indices_vector[0] += nodes_across;
    neighbour_indices_vector[1] += nodes_across - 1;
    neighbour_indices_vector[2] -= 1;
    neighbour_indices_vector[3] -= nodes_across + 1;
    neighbour_indices_vector[4] -= nodes_across;
    neighbour_indices_vector[5] -= nodes_across - 1;
    neighbour_indices_vector[6] += 1;
    neighbour_indices_vector[7] += nodes_across + 1;

    // Work out whether this node lies on any edge of the mesh
    bool on_south_edge = (nodeIndex < nodes_across);
    bool on_north_edge = (nodeIndex > nodes_up*(nodes_across - 1)-1);
    bool on_west_edge = (nodeIndex%nodes_across == 0);
    bool on_east_edge = (nodeIndex%nodes_across == nodes_across - 1);

    // Create a vector of booleans for which neighbours are available
    // Use the order N, NW, W, SW, S, SE, E, NE
    std::vector<bool> available_neighbours = std::vector<bool>(8, true);
    available_neighbours[0] = !on_north_edge;
    available_neighbours[1] = !(on_north_edge || on_west_edge);
    available_neighbours[2] = !on_west_edge;
    available_neighbours[3] = !(on_south_edge || on_west_edge);
    available_neighbours[4] = !on_south_edge;
    available_neighbours[5] = !(on_south_edge || on_east_edge);
    available_neighbours[6] = !on_east_edge;
    available_neighbours[7] = !(on_north_edge || on_east_edge);

    // Using neighbour_indices_vector and available_neighbours, store the indices of all available neighbours to the set all_neighbours
    for (unsigned i=0; i<8; i++)
    {
        if (available_neighbours[i])
        {
            neighbouring_node_indices.insert(neighbour_indices_vector[i]);
        }
    }
    return neighbouring_node_indices;
}

void PottsMesh::DeleteElement(unsigned index)
{
    // Mark this element as deleted; this also updates the nodes containing element indices
    this->mElements[index]->MarkAsDeleted();
    mDeletedElementIndices.push_back(index);
}

unsigned PottsMesh::DivideElement(PottsElement<2>* pElement,
                                  bool placeOriginalElementBelow)
{
    // Store the number of nodes in the element (this changes when nodes are deleted from the element)
    unsigned num_nodes = pElement->GetNumNodes();

    if (num_nodes < 2)
    {
        EXCEPTION("Tried to divide a Potts element with only one node.");
    }

    // Copy the nodes in this element
    std::vector<Node<2>*> nodes_elem;
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
    AddElement(new PottsElement<2>(new_element_index, nodes_elem));

    /**
     * Remove the correct nodes from each element. If placeOriginalElementBelow is true,
     * place the original element below (in the y direction) the new element; otherwise,
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
            height_midpoint_1 += pElement->GetNode(i)->rGetLocation()[1];
            counter_1++;
        }
        else
        {
            height_midpoint_2 += pElement->GetNode(i)->rGetLocation()[1];
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

unsigned PottsMesh::AddElement(PottsElement<2>* pNewElement)
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

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

//template class PottsMesh<1,1>;
//template class PottsMesh<1,2>;
//template class PottsMesh<1,3>;
//template class PottsMesh<2,2>;
//template class PottsMesh<2,3>;
//template class PottsMesh<3,3>;


//// Serialization for Boost >= 1.36
//#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_ALL_DIMS(PottsMesh)
