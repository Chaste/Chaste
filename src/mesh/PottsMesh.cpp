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

PottsMesh::PottsMesh(std::vector<Node<2>*> nodes, std::vector<PottsElement*> pottsElements)
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
        PottsElement* p_temp_vertex_element = pottsElements[elem_index];
        mElements.push_back(p_temp_vertex_element);
    }

    // Register elements with nodes
    for (unsigned index=0; index<mElements.size(); index++)
    {
        PottsElement* p_element = mElements[index];

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
}


unsigned PottsMesh::GetNumNodes() const
{
    return this->mNodes.size();
}


unsigned PottsMesh::GetNumElements() const
{
    return mElements.size();
}


unsigned PottsMesh::GetNumAllElements() const
{
    return mElements.size();
}


PottsElement* PottsMesh::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}


c_vector<double, 2> PottsMesh::GetCentroidOfElement(unsigned index)
{
    PottsElement* p_element = GetElement(index);
    unsigned num_nodes_in_element = p_element->GetNumNodes();


    // This should probably be returning the nearest node
    c_vector<double, 2> centroid = zero_vector<double>(2);

    double temp_centroid_x = 0;
    double temp_centroid_y = 0;

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        // Find locations of current node
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
//        PottsElement<2,2>* p_element = new PottsElement<2,2>(elem_index, nodes);
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
    PottsElement* p_element = GetElement(index);

    double element_volume = (double) p_element->GetNumNodes();

    return element_volume;
}


double PottsMesh::GetSurfaceAreaOfElement(unsigned index)
{
    // TODO this is not correct need to work this out from the number of free boudaries
    double surface_area = 0.0;

    return surface_area;
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
