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
#include "Debug.hpp"

template<unsigned DIM>
PottsMesh<DIM>::PottsMesh(std::vector<Node<DIM>*> nodes, std::vector<PottsElement<DIM>*> pottsElements)
{
    // Reset member variables and clear mNodes and mElements
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

    // Calculate and store the neeighbourhoods.
    CaclulateNeighbouringNodeIndices();

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

template<unsigned DIM>
c_vector<double, DIM> PottsMesh<DIM>::GetVectorFromAtoB(const c_vector<double, DIM>& rLocationA, const c_vector<double, DIM>& rLocationB)
{
    c_vector<double, DIM> vector = AbstractMesh<DIM, DIM>::GetVectorFromAtoB(rLocationA, rLocationB);

    return vector;
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
	assert(DIM==2 || DIM==3);

	// Get pointer to this element
	PottsElement<DIM>* p_element = GetElement(index);

    double surface_area = 0.0;
    for (unsigned node_index=0; node_index< p_element->GetNumNodes(); node_index++)
    {
    	std::set<unsigned> neighbouring_node_indices = GetVonNeumannNeighbouringNodeIndices(p_element->GetNode(node_index)->GetIndex());
    	unsigned local_edges = 2*DIM;
    	for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
			 iter != neighbouring_node_indices.end();
			 iter++)
    	{
    		std::set<unsigned> neighbouring_node_element_indices = this->mNodes[*iter]->rGetContainingElementIndices();

    		if (neighbouring_node_element_indices.size()>0 && local_edges>0)
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
void PottsMesh<DIM>::CaclulateNeighbouringNodeIndices()
{
    assert(DIM==2 || DIM==3);

    mMooreNeighbouringNodeIndices.resize(GetNumNodes());
    mVonNeumannNeighbouringNodeIndices.resize(GetNumNodes());

    c_vector<unsigned, DIM> num_nodes_across;

    for(unsigned i=0;i<DIM;i++)
    {
        num_nodes_across(i)=(unsigned)this->GetWidth(i) + 1;
    }

    for (unsigned node_index=0; node_index<GetNumNodes(); node_index++)
    {
        // Clear the set of neighbouring node indices

        mMooreNeighbouringNodeIndices[node_index].clear();

        switch (DIM)
        {
            case 2:
            {
                assert(DIM == 2);
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
                std::vector<unsigned> moore_neighbour_indices_vector(8, node_index);
                moore_neighbour_indices_vector[0] += num_nodes_across(0);
                moore_neighbour_indices_vector[1] += num_nodes_across(0) - 1;
                moore_neighbour_indices_vector[2] -= 1;
                moore_neighbour_indices_vector[3] -= num_nodes_across(0) + 1;
                moore_neighbour_indices_vector[4] -= num_nodes_across(0);
                moore_neighbour_indices_vector[5] -= num_nodes_across(0) - 1;
                moore_neighbour_indices_vector[6] += 1;
                moore_neighbour_indices_vector[7] += num_nodes_across(0) + 1;

                // Work out whether this node lies on any edge of the mesh
                bool on_south_edge = (node_index < num_nodes_across(0));
                bool on_north_edge = (node_index > num_nodes_across(0)*(num_nodes_across(1) - 1) - 1);
                bool on_west_edge = (node_index%num_nodes_across(0) == 0);
                bool on_east_edge = (node_index%num_nodes_across(0) == num_nodes_across(0) - 1);

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
                        assert(moore_neighbour_indices_vector[i] < this->mNodes.size());
                        mMooreNeighbouringNodeIndices[node_index].insert(moore_neighbour_indices_vector[i]);
                    }
                }
                break;
            }
            case 3:
            {
                assert(DIM ==3 );
                /*
                 * This stores the available neighbours using the following numbering:
                 *                      FRONT           BACK
                 *  1-----0-----7   10-----9---16   19----18---25
                 *  |     |     |   |     |     |   |     |     |
                 *  |     |     |   |     |     |   |     |     |
                 *  2-----x-----6   11----8-----15  20----17---24
                 *  |     |     |   |     |     |   |     |     |
                 *  |     |     |   |     |     |   |     |     |
                 *  3-----4-----5   12----13----14  21---22----23
                 */

                // Create a vector of possible neighbouring node indices
                std::vector<unsigned> moore_neighbour_indices_vector(26, node_index);
                moore_neighbour_indices_vector[0] += num_nodes_across(0);
                moore_neighbour_indices_vector[1] += num_nodes_across(0) - 1;
                moore_neighbour_indices_vector[2] -= 1;
                moore_neighbour_indices_vector[3] -= num_nodes_across(0) + 1;
                moore_neighbour_indices_vector[4] -= num_nodes_across(0);
                moore_neighbour_indices_vector[5] -= num_nodes_across(0) - 1;
                moore_neighbour_indices_vector[6] += 1;
                moore_neighbour_indices_vector[7] += num_nodes_across(0) + 1;
                moore_neighbour_indices_vector[8] -= num_nodes_across(0)*num_nodes_across(1);
                for( unsigned i=9; i<17; i++)
                {
                    moore_neighbour_indices_vector[i] = moore_neighbour_indices_vector[i-9]-num_nodes_across(0)*num_nodes_across(1);
                }
                moore_neighbour_indices_vector[17] += num_nodes_across(0)*num_nodes_across(1);
                for (unsigned i=18; i<26; i++)
                {
                    moore_neighbour_indices_vector[i]=moore_neighbour_indices_vector[i-18]+num_nodes_across(0)*num_nodes_across(1);
                }

                // Work out whether this node lies on any edge of the mesh
                bool on_south_edge = (node_index%(num_nodes_across(0)*num_nodes_across(1))<num_nodes_across(0));
                bool on_north_edge = (node_index%(num_nodes_across(0)*num_nodes_across(1))>(num_nodes_across(0)*num_nodes_across(1)-num_nodes_across(0)-1));
                bool on_west_edge = (node_index%num_nodes_across(0) == 0);
                bool on_east_edge = (node_index%num_nodes_across(0) == num_nodes_across(0) - 1);
                bool on_front_edge = (node_index < num_nodes_across(0)*num_nodes_across(1)-1);
                bool on_back_edge = (node_index > num_nodes_across(0)*num_nodes_across(1)*num_nodes_across(2)-num_nodes_across(0)*num_nodes_across(1)-1);

                // Create a vector of booleans for which neighbours are available
                // Use the order N, NW, W, SW, S, SE, E, NE
                std::vector<bool> available_neighbours = std::vector<bool>(26, true);
                available_neighbours[0] = !on_north_edge;
                available_neighbours[1] = !(on_north_edge || on_west_edge);
                available_neighbours[2] = !on_west_edge;
                available_neighbours[3] = !(on_south_edge || on_west_edge);
                available_neighbours[4] = !on_south_edge;
                available_neighbours[5] = !(on_south_edge || on_east_edge);
                available_neighbours[6] = !on_east_edge;
                available_neighbours[7] = !(on_north_edge || on_east_edge);
                available_neighbours[8] = !(on_front_edge);
                for (unsigned i=9; i<17; i++)
                {
                    available_neighbours[i] = (available_neighbours[i-9] && !(on_front_edge));
                }
                available_neighbours[17] = !(on_back_edge);
                for (unsigned i=18; i<26; i++)
                {
                    available_neighbours[i] = (available_neighbours[i-18] && !(on_back_edge));
                }

                // Using neighbour_indices_vector and available_neighbours, store the indices of all available neighbours to the set all_neighbours
                for (unsigned i=0; i<26; i++)
                {
                    if (available_neighbours[i] && moore_neighbour_indices_vector[i] < num_nodes_across(0)*num_nodes_across(1)*num_nodes_across(2))
                    {
                        assert(moore_neighbour_indices_vector[i] < this->mNodes.size());
                        mMooreNeighbouringNodeIndices[node_index].insert(moore_neighbour_indices_vector[i]);
                    }
                }
                break;
            }
            default:
                NEVER_REACHED;
        }

        // Clear the set of neighbouring node indices
        mVonNeumannNeighbouringNodeIndices[node_index].clear();

        switch (DIM)
        {
            case 2:
            {
                assert(DIM == 2);
                /*
                 * This stores the available neighbours using the following numbering:
                 *
                 *        0
                 *        |
                 *        |
                 *  1-----x-----3
                 *        |
                 *        |
                 *        2
                 */
                // Create a vector of possible neighbouring node indices
                std::vector<unsigned> von_neumann_neighbour_indices_vector(4, node_index);
                von_neumann_neighbour_indices_vector[0] += num_nodes_across(0);
                von_neumann_neighbour_indices_vector[1] -= 1;
                von_neumann_neighbour_indices_vector[2] -= num_nodes_across(0);
                von_neumann_neighbour_indices_vector[3] += 1;

                // Work out whether this node lies on any edge of the mesh
                bool on_south_edge = (node_index < num_nodes_across(0));
                bool on_north_edge = (node_index > num_nodes_across(0)*(num_nodes_across(1) - 1) - 1);
                bool on_west_edge = (node_index%num_nodes_across(0) == 0);
                bool on_east_edge = (node_index%num_nodes_across(0) == num_nodes_across(0) - 1);

                // Create a vector of booleans for which neighbours are available
                // Use the order N, W, S, E
                std::vector<bool> available_neighbours = std::vector<bool>(2*DIM, true);
                available_neighbours[0] = !on_north_edge;
                available_neighbours[1] = !on_west_edge;
                available_neighbours[2] = !on_south_edge;
                available_neighbours[3] = !on_east_edge;

                // Using von_neumann_neighbour_indices_vector and available_neighbours, store the indices of all available neighbours to the set all_neighbours
                for (unsigned i=0; i<4; i++)
                {
                    if (available_neighbours[i])
                    {
                        assert(von_neumann_neighbour_indices_vector[i] < this->mNodes.size());
                        mVonNeumannNeighbouringNodeIndices[node_index].insert(von_neumann_neighbour_indices_vector[i]);
                    }
                }
                break;
            }
            case 3:
            {
                assert(DIM == 3);

                /*
                 * This stores the available neighbours using the following numbering:
                 *
                 *        0  5
                 *        | /
                 *        |/
                 *  1-----x-----3
                 *      / |
                 *     /  |
                 *    4   2
                 */
                // Create a vector of possible neighbouring node indices
                std::vector<unsigned> von_neumann_neighbour_indices_vector(6, node_index);
                von_neumann_neighbour_indices_vector[0] += num_nodes_across(0);
                von_neumann_neighbour_indices_vector[1] -= 1;
                von_neumann_neighbour_indices_vector[2] -= num_nodes_across(0);
                von_neumann_neighbour_indices_vector[3] += 1;
                von_neumann_neighbour_indices_vector[4] -= num_nodes_across(0)*num_nodes_across(1);
                von_neumann_neighbour_indices_vector[5] += num_nodes_across(0)*num_nodes_across(1);

                // Work out whether this node lies on any edge of the mesh
                bool on_south_edge = (node_index%(num_nodes_across(0)*num_nodes_across(1))<num_nodes_across(0));
                bool on_north_edge = (node_index%(num_nodes_across(0)*num_nodes_across(1))>(num_nodes_across(0)*num_nodes_across(1)-num_nodes_across(0)-1));
                bool on_west_edge = (node_index%num_nodes_across(0)== 0);
                bool on_east_edge = (node_index%num_nodes_across(0) == num_nodes_across(0) - 1);
                bool on_front_edge = (node_index < num_nodes_across(0)*num_nodes_across(1)-1);
                bool on_back_edge = (node_index > num_nodes_across(0)*num_nodes_across(1)*num_nodes_across(2)-num_nodes_across(0)*num_nodes_across(1)-1);

                // Create a vector of booleans for which neighbours are available
                // Use the order N, W, S, E, F, B
                std::vector<bool> available_neighbours = std::vector<bool>(2*DIM, true);
                available_neighbours[0] = !on_north_edge;
                available_neighbours[1] = !on_west_edge;
                available_neighbours[2] = !on_south_edge;
                available_neighbours[3] = !on_east_edge;
                available_neighbours[4] = !on_front_edge;
                available_neighbours[5] = !on_back_edge;

                // Using von_neumann_neighbour_indices_vector and available_neighbours, store the indices of all available neighbours to the set all_neighbours
                for (unsigned i=0; i<6; i++)
                {
                    if (available_neighbours[i] && von_neumann_neighbour_indices_vector[i]<num_nodes_across(0)*num_nodes_across(1)*num_nodes_across(2))
                    {
                        assert(von_neumann_neighbour_indices_vector[i] < this->mNodes.size());
                        mVonNeumannNeighbouringNodeIndices[node_index].insert(von_neumann_neighbour_indices_vector[i]);
                    }
                }
                break;
            }
            default:
                NEVER_REACHED;
        }
    }
}

template<unsigned DIM>
void PottsMesh<DIM>::DeleteElement(unsigned index)
{
    // Mark this element as deleted; this also updates the nodes containing element indices
    this->mElements[index]->MarkAsDeleted();
    mDeletedElementIndices.push_back(index);
}

template<unsigned DIM>
unsigned PottsMesh<DIM>::DivideElement(PottsElement<DIM>* pElement,
                                  bool placeOriginalElementBelow)
{
	///\todo not implemented in 3d yet
	assert(DIM==2);

    // Store the number of nodes in the element (this changes when nodes are deleted from the element)
    unsigned num_nodes = pElement->GetNumNodes();

    if (num_nodes < 2)
    {
        EXCEPTION("Tried to divide a Potts element with only one node.");
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
void PottsMesh<DIM>::ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rMeshReader)
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
            p_element->SetRegion(attribute_value);
        }
    }

    // Calculate and store the neeighbourhoods.
    CaclulateNeighbouringNodeIndices();
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class PottsMesh<1>;
template class PottsMesh<2>;
template class PottsMesh<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PottsMesh)
