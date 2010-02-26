/*

Copyright (C) University of Oxford, 2005-2010

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

#include "VertexMesh.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                               std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertexElements,
                                               double cellRearrangementThreshold,
                                               double edgeDivisionThreshold,
                                               double t2Threshold)
    : mCellRearrangementThreshold(cellRearrangementThreshold),
      mEdgeDivisionThreshold(edgeDivisionThreshold),
      mT2Threshold(t2Threshold)
{
    assert(cellRearrangementThreshold > 0.0);
    assert(edgeDivisionThreshold > 0.0);
    assert(t2Threshold > 0.0);

    // Reset member variables and clear mNodes and mElements
    Clear();

    // Populate mNodes and mElements
    for (unsigned node_index=0; node_index<nodes.size(); node_index++)
    {
        Node<SPACE_DIM>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }
    for (unsigned elem_index=0; elem_index<vertexElements.size(); elem_index++)
    {
        VertexElement<ELEMENT_DIM,SPACE_DIM>* p_temp_vertex_element = vertexElements[elem_index];
        mElements.push_back(p_temp_vertex_element);
    }

    // Register elements with nodes
    for (unsigned index=0; index<mElements.size(); index++)
    {
        VertexElement<ELEMENT_DIM,SPACE_DIM>* p_temp_vertex_element = mElements[index];
        for (unsigned node_index=0; node_index<p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<SPACE_DIM>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }

    this->mMeshChangesDuringSimulation = true;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh()
    : mCellRearrangementThreshold(0.01), // Overwritten as soon as archiving is complete
      mEdgeDivisionThreshold(DBL_MAX), // Overwritten as soon as archiving is complete
      mT2Threshold(0.001) // Overwritten as soon as archiving is complete
{
    this->mMeshChangesDuringSimulation = true;
    Clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::~VertexMesh()
{
    Clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size());
    return index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size());
    return index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    /// \todo sort out boundary elements in a vertex mesh
//    assert(index < this->mBoundaryElements.size() );
    return index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCellRearrangementThreshold() const
{
    return mCellRearrangementThreshold;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetEdgeDivisionThreshold() const
{
    return mEdgeDivisionThreshold;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetT2Threshold() const
{
    return mT2Threshold;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetCellRearrangementThreshold(double cellRearrangementThreshold)
{
    mCellRearrangementThreshold = cellRearrangementThreshold;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetEdgeDivisionThreshold(double edgeDivisionThreshold)
{
    mEdgeDivisionThreshold = edgeDivisionThreshold;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetT2Threshold(double t2Threshold)
{
    mT2Threshold = t2Threshold;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    mDeletedNodeIndices.clear();
    mDeletedElementIndices.clear();

    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }

    this->mNodes.clear();
    mElements.clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size() - mDeletedNodeIndices.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size() - mDeletedElementIndices.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements() const
{
    return mElements.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM,SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetAreaOfElement(unsigned index)
{
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2);
    #undef COVERAGE_IGNORE

    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    c_vector<double, SPACE_DIM> first_node;
    c_vector<double, SPACE_DIM> current_node;
    c_vector<double, SPACE_DIM> anticlockwise_node;

    unsigned num_nodes_in_element = p_element->GetNumNodes();

    double element_area = 0;

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        // Find locations of current node and anticlockwise node
        current_node = p_element->GetNodeLocation(local_index);
        anticlockwise_node = p_element->GetNodeLocation((local_index+1)%num_nodes_in_element);

        element_area += 0.5*(current_node[0]*anticlockwise_node[1] - anticlockwise_node[0]*current_node[1]);
    }

    return element_area;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetPerimeterOfElement(unsigned index)
{
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2);
    #undef COVERAGE_IGNORE

    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    unsigned current_node_index;
    unsigned anticlockwise_node_index;
    unsigned num_nodes_in_element = p_element->GetNumNodes();

    double element_perimeter = 0;

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        // Find locations of current node and anticlockwise node
        current_node_index = p_element->GetNodeGlobalIndex(local_index);
        anticlockwise_node_index = p_element->GetNodeGlobalIndex((local_index+1)%num_nodes_in_element);

        element_perimeter += this->GetDistanceBetweenNodes(current_node_index, anticlockwise_node_index);
    }

    return element_perimeter;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCentroidOfElement(unsigned index)
{
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2);
    #undef COVERAGE_IGNORE

    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);
    c_vector<double, SPACE_DIM> current_node;
    c_vector<double, SPACE_DIM> anticlockwise_node;

    double temp_centroid_x = 0;
    double temp_centroid_y = 0;

    unsigned num_nodes_in_element = p_element->GetNumNodes();

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        // Find locations of current node and anticlockwise node
        current_node = p_element->GetNodeLocation(local_index);
        anticlockwise_node = p_element->GetNodeLocation((local_index+1)%num_nodes_in_element);

        temp_centroid_x += (current_node[0]+anticlockwise_node[0])*(current_node[0]*anticlockwise_node[1]-current_node[1]*anticlockwise_node[0]);
        temp_centroid_y += (current_node[1]+anticlockwise_node[1])*(current_node[0]*anticlockwise_node[1]-current_node[1]*anticlockwise_node[0]);
    }

    double vertex_area = GetAreaOfElement(index);
    double centroid_coefficient = 1.0/(6.0*vertex_area);

    centroid(0) = centroid_coefficient*temp_centroid_x;
    centroid(1) = centroid_coefficient*temp_centroid_y;

    return centroid;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetAreaGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex)
{
    #define COVERAGE_IGNORE
    assert(SPACE_DIM==2);
    #undef COVERAGE_IGNORE

    unsigned num_nodes_in_element = pElement->GetNumNodes();
    unsigned next_local_index = (localIndex+1)%num_nodes_in_element;

    // We add an extra num_nodes_in_element in the line below as otherwise this term can be negative, which breaks the % operator
    unsigned previous_local_index = (num_nodes_in_element+localIndex-1)%num_nodes_in_element;

    c_vector<double, SPACE_DIM> previous_node_location = pElement->GetNodeLocation(previous_local_index);
    c_vector<double, SPACE_DIM> next_node_location = pElement->GetNodeLocation(next_local_index);
    c_vector<double, SPACE_DIM> difference_vector = this->GetVectorFromAtoB(previous_node_location, next_node_location);

    c_vector<double, SPACE_DIM> area_gradient;

    area_gradient[0] = 0.5*difference_vector[1];
    area_gradient[1] = -0.5*difference_vector[0];

    return area_gradient;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetPreviousEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex)
{
    #define COVERAGE_IGNORE
    assert(SPACE_DIM==2);
    #undef COVERAGE_IGNORE

    unsigned num_nodes_in_element = pElement->GetNumNodes();

    // We add an extra localIndex-1 in the line below as otherwise this term can be negative, which breaks the % operator
    unsigned previous_local_index = (num_nodes_in_element+localIndex-1)%num_nodes_in_element;

    unsigned current_global_index = pElement->GetNodeGlobalIndex(localIndex);
    unsigned previous_global_index = pElement->GetNodeGlobalIndex(previous_local_index);

    double previous_edge_length = this->GetDistanceBetweenNodes(current_global_index, previous_global_index);
    assert(previous_edge_length > DBL_EPSILON);

    c_vector<double, SPACE_DIM> previous_edge_gradient = this->GetVectorFromAtoB(pElement->GetNodeLocation(previous_local_index), pElement->GetNodeLocation(localIndex))/previous_edge_length;

    return previous_edge_gradient;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNextEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex)
{
    #define COVERAGE_IGNORE
    assert(SPACE_DIM==2);
    #undef COVERAGE_IGNORE

    unsigned next_local_index = (localIndex+1)%(pElement->GetNumNodes());

    unsigned current_global_index = pElement->GetNodeGlobalIndex(localIndex);
    unsigned next_global_index = pElement->GetNodeGlobalIndex(next_local_index);

    double next_edge_length = this->GetDistanceBetweenNodes(current_global_index, next_global_index);
    assert(next_edge_length > DBL_EPSILON);

    c_vector<double, SPACE_DIM> next_edge_gradient = this->GetVectorFromAtoB(pElement->GetNodeLocation(next_local_index), pElement->GetNodeLocation(localIndex))/next_edge_length;

    return next_edge_gradient;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetPerimeterGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex)
{
    #define COVERAGE_IGNORE
    assert(SPACE_DIM==2);
    #undef COVERAGE_IGNORE

    c_vector<double, SPACE_DIM> previous_edge_gradient = GetPreviousEdgeGradientOfElementAtNode(pElement, localIndex);
    c_vector<double, SPACE_DIM> next_edge_gradient = GetNextEdgeGradientOfElementAtNode(pElement, localIndex);

    return previous_edge_gradient + next_edge_gradient;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> VertexMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMomentsOfElement(unsigned index)
{
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2);
    #undef COVERAGE_IGNORE

    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);
    unsigned num_nodes_in_element = p_element->GetNumNodes();
    c_vector<double, 2> centroid = GetCentroidOfElement(index);

    c_vector<double, 3> moments = zero_vector<double>(3);

    unsigned node_1;
    unsigned node_2;

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        node_1 = local_index;
        node_2 = (local_index+1)%num_nodes_in_element;

        // Original position of nodes
        c_vector<double, 2> original_pos_1 = p_element->GetNodeLocation(node_1);
        c_vector<double, 2> original_pos_2 = p_element->GetNodeLocation(node_2);

        // Node position so centerd on origin
        c_vector<double, 2> pos_1 = this->GetVectorFromAtoB(centroid, original_pos_1);
        c_vector<double, 2> pos_2 = this->GetVectorFromAtoB(centroid, original_pos_2);

        /*
         * Note these formulae require the polygon to be centered on the origin
         */
        double a = pos_1(0)*pos_2(1)-pos_2(0)*pos_1(1);

        // Ixx
        moments(0) += (  pos_1(1)*pos_1(1)
                       + pos_1(1)*pos_2(1)
                       + pos_2(1)*pos_2(1) ) * a;

        // Iyy
        moments(1) += (  pos_1(0)*pos_1(0)
                       + pos_1(0)*pos_2(0)
                       + pos_2(0)*pos_2(0) ) * a;

        // Ixy
        moments(2) += (  pos_1(0)*pos_2(1)
                       + 2*pos_1(0)*pos_1(1)
                       + 2*pos_2(0)*pos_2(1)
                       + pos_2(0)*pos_1(1) ) * a;
    }

    moments(0) /= 12;
    moments(1) /= 12;
    moments(2) /= 24;

    return moments;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetShortAxisOfElement(unsigned index)
{
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2);
    #undef COVERAGE_IGNORE

    c_vector<double, SPACE_DIM> short_axis = zero_vector<double>(SPACE_DIM);
    c_vector<double, 3> moments = CalculateMomentsOfElement(index);

    double largest_eigenvalue, discriminant;

    discriminant = sqrt((moments(0) - moments(1))*(moments(0) - moments(1)) + 4.0*moments(2)*moments(2));
    // This is always the largest eigenvalue as both eigenvalues are real as it is a
    // symmetric matrix
    largest_eigenvalue = (moments(0) + moments(1) + discriminant)*0.5;
    if (fabs(discriminant) < 1e-10)
    {
        // Return a random unit vector
        short_axis(0) = RandomNumberGenerator::Instance()->ranf();
        short_axis(1) = sqrt(1.0 - short_axis(0)*short_axis(0));
    }
    else
    {
        if (moments(2) == 0.0)
        {
            short_axis(0) = (moments(0) < moments(1)) ? 0.0 : 1.0;
            short_axis(1) = (moments(0) < moments(1)) ? 1.0 : 0.0;
        }
        else
        {
            short_axis(0) = 1.0;
            short_axis(1) = (moments(0) - largest_eigenvalue)/moments(2);

            double length_short_axis = norm_2(short_axis);

            short_axis /= length_short_axis;
        }
    }
    return short_axis;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNewNode)
{
    if (mDeletedNodeIndices.empty())
    {
        pNewNode->SetIndex(this->mNodes.size());
        this->mNodes.push_back(pNewNode);
    }
    else
    {
        unsigned index = mDeletedNodeIndices.back();
        pNewNode->SetIndex(index);
        mDeletedNodeIndices.pop_back();
        delete this->mNodes[index];
        this->mNodes[index] = pNewNode;
    }
    return pNewNode->GetIndex();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::AddElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pNewElement)
{
    unsigned new_element_index = pNewElement->GetIndex();

    if (new_element_index == mElements.size())
    {
        mElements.push_back(pNewElement);
    }
    else
    {
        mElements[new_element_index] = pNewElement;
    }
    pNewElement->RegisterWithNodes();
    return pNewElement->GetIndex();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point)
{
    this->mNodes[nodeIndex]->SetPoint(point);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::DeleteElementPriorToReMesh(unsigned index)
{
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2);
    #undef COVERAGE_IGNORE

    // Mark any nodes that are contained only in this element as deleted
    for (unsigned i=0; i<this->mElements[index]->GetNumNodes(); i++)
    {
        Node<SPACE_DIM>* p_node = this->mElements[index]->GetNode(i);

        if (p_node->rGetContainingElementIndices().size()==1)
        {
            p_node->MarkAsDeleted();
            mDeletedNodeIndices.push_back(p_node->GetIndex());
        }

        // Mark all the nodes contained in the removed element as boundary nodes
        p_node->SetAsBoundaryNode(true);

    }

    // Mark this element as deleted
    this->mElements[index]->MarkAsDeleted();
    mDeletedElementIndices.push_back(index);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideEdge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(shared_elements.size() > 0);

    // Create a new node (position is not important as it will be changed)
    Node<SPACE_DIM>* p_new_node = new Node<SPACE_DIM>(GetNumNodes(), false, 0.0, 0.0);

    // Update the node location
    c_vector<double, SPACE_DIM> new_node_position = pNodeA->rGetLocation() + 0.5*GetVectorFromAtoB(pNodeA->rGetLocation(), pNodeB->rGetLocation());
    ChastePoint<SPACE_DIM> point(new_node_position);
    p_new_node->SetPoint(new_node_position);

    // Add node to mesh
    this->mNodes.push_back(p_new_node);

    // Iterate over common elements
    for (std::set<unsigned>::iterator iter = shared_elements.begin();
         iter != shared_elements.end();
         ++iter)
    {
        // Find which node has the lower local index in this element
        unsigned local_indexA = GetElement(*iter)->GetNodeLocalIndex(pNodeA->GetIndex());
        unsigned local_indexB = GetElement(*iter)->GetNodeLocalIndex(pNodeB->GetIndex());

        unsigned index = local_indexB;
        if ( local_indexB > local_indexA )
        {
            index = local_indexA;
        }
        if ( (local_indexA == 0) && (local_indexB == GetElement(*iter)->GetNumNodes()-1))
        {
            index = local_indexB;
        }
        if ( (local_indexB == 0) && (local_indexA == GetElement(*iter)->GetNumNodes()-1))
        {
            index = local_indexA;
        }
        // Add new node to this element
        GetElement(*iter)->AddNode(index, p_new_node);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::RemoveDeletedNodesAndElements(VertexElementMap& rElementMap)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM==2 || SPACE_DIM==3);
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Make sure the map is big enough
    rElementMap.Resize(GetNumAllElements());

    // Remove deleted elements
    std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> live_elements;
    for (unsigned i=0; i<mElements.size(); i++)
    {
        if (mElements[i]->IsDeleted())
        {
            delete this->mElements[i];
        }
        else
        {
            live_elements.push_back(mElements[i]);
            rElementMap.SetNewIndex(i, (unsigned)(live_elements.size()-1));
        }
    }

    assert(mDeletedElementIndices.size() == mElements.size() - live_elements.size());
    mDeletedElementIndices.clear();
    mElements = live_elements;

    for (unsigned i=0; i<mElements.size(); i++)
    {
        mElements[i]->ResetIndex(i);
    }

    // Remove deleted nodes
    RemoveDeletedNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::RemoveDeletedNodes()
{
    std::vector<Node<SPACE_DIM>*> live_nodes;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        if (this->mNodes[i]->IsDeleted())
        {
            delete this->mNodes[i];
        }
        else
        {
            live_nodes.push_back(this->mNodes[i]);
        }
    }

    assert(mDeletedNodeIndices.size() == this->mNodes.size() - live_nodes.size());
    this->mNodes = live_nodes;
    mDeletedNodeIndices.clear();

    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->SetIndex(i);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh(VertexElementMap& rElementMap)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM==2 || SPACE_DIM==3);
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    if (SPACE_DIM==2)
    {
        // Remove deleted nodes and elements
        RemoveDeletedNodesAndElements(rElementMap);

        /*
         * We do not need to call Clear() and remove all current data, since
         * cell birth, rearrangement and death result only in local remeshing
         * of a vertex-based mesh.
         */

        // Restart check after each T1/T2Swap as elements are changed
        bool recheck_mesh = true;
        while (recheck_mesh == true)
        {
            recheck_mesh = false;

            // Loop over elements to check for T2Swaps
            // Separate loops as need to check for T2Swaps first
            for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_iter = GetElementIteratorBegin();
                 elem_iter != GetElementIteratorEnd();
                 ++elem_iter)
            {
                if (!recheck_mesh)
                {
                    if (elem_iter->GetNumNodes() == 3u)
                    {
                        /*
                         * Perform T2 swaps where necesary
                         * Check there are only 3 nodes and the element is small enough
                         */
                        if (GetAreaOfElement(elem_iter->GetIndex()) < GetT2Threshold())
                        {
                            PerformT2Swap(*elem_iter);
                            // Now remove the deleted nodes (if we don't do this then the search for T1Swap causes errors)
                            RemoveDeletedNodesAndElements(rElementMap);
                            recheck_mesh = true;
                            break;
                        }
                    }
                }
            }
            // Loop over elements to check for T1Swaps
            for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_iter = GetElementIteratorBegin();
                 elem_iter != GetElementIteratorEnd();
                 ++elem_iter)
            {
                if (!recheck_mesh)
                {
                    unsigned num_nodes = elem_iter->GetNumNodes();
                    assert(num_nodes > 0); // if not element should be deleted

                    unsigned new_num_nodes = num_nodes;

                    /*
                     *  Perform T1 swaps and divide edges where necesary
                     *  Check there are > 3 nodes in both elements that contain the pair of nodes
                     *  and the edges are small enough
                     */

                    // Loop over element vertices
                    for (unsigned local_index=0; local_index<num_nodes; local_index++)
                    {
                        // Find locations of current node and anticlockwise node
                        Node<SPACE_DIM>* p_current_node = elem_iter->GetNode(local_index);
                        unsigned local_index_plus_one = (local_index+1)%new_num_nodes; /// \todo Should use iterators to tidy this up
                        Node<SPACE_DIM>* p_anticlockwise_node = elem_iter->GetNode(local_index_plus_one);

                        // Find distance between nodes
                        double distance_between_nodes = this->GetDistanceBetweenNodes(p_current_node->GetIndex(), p_anticlockwise_node->GetIndex());

                        std::set<unsigned> elements_of_node_a = p_current_node->rGetContainingElementIndices();
                        std::set<unsigned> elements_of_node_b = p_anticlockwise_node->rGetContainingElementIndices();

                        std::set<unsigned> all_elements;
                        std::set_union(elements_of_node_a.begin(), elements_of_node_a.end(),
                                       elements_of_node_b.begin(), elements_of_node_b.end(),
                                       std::inserter(all_elements, all_elements.begin()));

                        // Track if either node is in a trinagular element.
                        bool triangular_element = false;

                        for (std::set<unsigned>::const_iterator it = all_elements.begin();
                             it != all_elements.end();
                             ++it)
                        {
                            ///\todo this magic number needs to be investigated.
                            if (this->GetElement(*it)->GetNumNodes() <= 3)//&&(this->GetAreaOfElement(*it) < 2.0*GetT2Threshold()))
                            {
                                triangular_element = true;
                            }
                        }

                        // If the nodes are too close together and we don't have any triangular elements connected to the nodes, perform a swap
                        if ((!triangular_element) && (distance_between_nodes < mCellRearrangementThreshold))
                        {
                            // Identify the type of node swap/merge needed then call method to perform swap/merge
                            IdentifySwapType(p_current_node, p_anticlockwise_node, rElementMap);
                            recheck_mesh = true;
                            break;
                        }

                        if (distance_between_nodes > mEdgeDivisionThreshold)
                        {
                            // If the nodes are too far apart, divide the edge
                            DivideEdge(p_current_node, p_anticlockwise_node);
                            new_num_nodes++;
                        }
                    }
                }
                else
                {
                    break;
                }
            }
        }

        // ... end of element rearrangement code

        // areas and perimeters of elements are sorted in PerformT1Swap() method

        // Check that no nodes have overlapped elements
        /// \todo Only need to check this next bit if the element/node is on the boundary (see #933 and #943)
        for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_iter = GetElementIteratorBegin();
             elem_iter != GetElementIteratorEnd();
             ++elem_iter)
        {
            unsigned num_nodes = elem_iter->GetNumNodes();

            // Loop over element vertices
            for (unsigned local_index=0; local_index<num_nodes; local_index++)
            {
                // Find locations of current node and anticlockwise node
                Node<SPACE_DIM>* p_current_node = elem_iter->GetNode(local_index);

                if (p_current_node->IsBoundaryNode())
                {
                    for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator other_iter = GetElementIteratorBegin();
                         other_iter != GetElementIteratorEnd();
                         ++other_iter)
                    {
                        if (other_iter != elem_iter)
                        {
                            if (ElementIncludesPoint(p_current_node->rGetLocation(), other_iter->GetIndex()))
                            {
//                                TRACE("*****************************Overlap******************************");
                                MoveOverlappingNodeOntoEdgeOfElement(p_current_node, other_iter->GetIndex());
                            }
                        }
                    }
                }
            }
        }
    }
    else // 3D
    {
        #define COVERAGE_IGNORE
        EXCEPTION("Remeshing has not been implemented in 3D (see #827 and #860)\n");
        #undef COVERAGE_IGNORE
        /// \todo put code for remeshing in 3D here - see #866 and the paper doi:10.1016/j.jtbi.2003.10.001
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh()
{
    VertexElementMap map(GetNumElements());
    ReMesh(map);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
    // Create a set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices;

    // Find the indices of the elements owned by this node
    std::set<unsigned> containing_elem_indices = this->GetNode(nodeIndex)->rGetContainingElementIndices();

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter = containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
         ++elem_iter)
    {
        // Find the local index of this node in this element
        unsigned local_index = GetElement(*elem_iter)->GetNodeLocalIndex(nodeIndex);

        // Find the global indices of the preceding and successive nodes in this element
        unsigned num_nodes = GetElement(*elem_iter)->GetNumNodes();
        unsigned previous_local_index = (local_index + num_nodes - 1)%num_nodes;
        unsigned next_local_index = (local_index + 1)%num_nodes;

        // Add the global indices of these two nodes to the set of neighbouring node indices
        neighbouring_node_indices.insert(GetElement(*elem_iter)->GetNodeGlobalIndex(previous_local_index));
        neighbouring_node_indices.insert(GetElement(*elem_iter)->GetNodeGlobalIndex(next_local_index));
    }

    return neighbouring_node_indices;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodeNotAlsoInElement(unsigned nodeIndex, unsigned elemIndex)
{
    /// \todo We should probably assert here that the node is in fact contained in the element

    // Create a set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices_not_in_this_element;

    // Get the indices of this node's neighbours
    std::set<unsigned> node_neighbours = GetNeighbouringNodeIndices(nodeIndex);

    // Get the indices of the nodes contained in this element
    std::set<unsigned> node_indices_in_this_element;
    for (unsigned local_index=0; local_index<GetElement(elemIndex)->GetNumNodes(); local_index++)
    {
        unsigned global_index = GetElement(elemIndex)->GetNodeGlobalIndex(local_index);
        node_indices_in_this_element.insert(global_index);
    }

    // Check if each neighbour is also in this element; if not, add it to the set
    for (std::set<unsigned>::iterator iter = node_neighbours.begin();
         iter != node_neighbours.end();
         ++iter)
    {
        if (node_indices_in_this_element.find(*iter) == node_indices_in_this_element.end())
        {
            neighbouring_node_indices_not_in_this_element.insert(*iter);
        }
    }

    return neighbouring_node_indices_not_in_this_element;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::IdentifySwapType(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, VertexElementMap& rElementMap)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // this method only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Find the sets of elements containing nodes A and B
    std::set<unsigned> nodeA_elem_indices = pNodeA->rGetContainingElementIndices();

    std::set<unsigned> nodeB_elem_indices = pNodeB->rGetContainingElementIndices();

    // Form the set union
    std::set<unsigned> all_indices, temp_set;
    std::set_union(nodeA_elem_indices.begin(), nodeA_elem_indices.end(),
                   nodeB_elem_indices.begin(), nodeB_elem_indices.end(),
                   std::inserter(temp_set, temp_set.begin()));
    all_indices.swap(temp_set); // temp_set will be deleted

    if ((nodeA_elem_indices.size()>3) || (nodeB_elem_indices.size()>3))
    {
        /*
         * Looks like
         *
         *  \
         *   \ A   B
         * ---o---o---
         *   /
         *  /
         *
         */
        EXCEPTION("A node is contained in more than three elements"); // the code can't handle this case
    }
    else // each node is contained in at most three elements
    {
        switch (all_indices.size())
        {
            case 1:
            {
                /*
                 * In this case, each node is contained in a single element, so the nodes
                 * lie on the boundary of the mesh:
                 *
                 *    A   B
                 * ---o---o---
                 *
                 * We merge the nodes.
                 */
                PerformNodeMerge(pNodeA, pNodeB);

                // Remove the deleted node and re-index
                RemoveDeletedNodes();
                break;
            }
            case 2:
            {
                if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)
                {
                    /*
                     * In this case, each node is contained in two elements, so the nodes
                     * lie on an internal edge:
                     *
                     *    A   B
                     * ---o---o---
                     *
                     * We merge the nodes in this case.
                     */
                     PerformNodeMerge(pNodeA, pNodeB);

                     // Remove the deleted node and re-index
                     RemoveDeletedNodes();
                }
                else
                {
                    /*
                     * In this case, the node configuration looks like:
                     *
                     * Outside
                     *         /
                     *   --o--o (2)
                     *     (1) \
                     *
                     * We merge the nodes in this case.
                     */
                     PerformNodeMerge(pNodeA, pNodeB);

                     // Remove the deleted node and re-index
                     RemoveDeletedNodes();
                }
                break;
            }
            case 3:
            {
                if (nodeA_elem_indices.size()==1 || nodeB_elem_indices.size()==1)
                {
                    /*
                     * In this case, one node is contained in one element and
                     * the other node is contained in three elements:
                     *
                     *    A   B
                     *
                     *  empty   /
                     *         / (3)
                     * ---o---o-----   (element number in brackets)
                     *  (1)    \ (2)
                     *          \
                     *
                     * We perform a node merge in this case.
                     */
                    PerformNodeMerge(pNodeA, pNodeB);

                    // Remove the deleted node and re-index
                    RemoveDeletedNodes();
                }
                else if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)
                {
                    /*
                     * In this case, the node configuration looks like:
                     *
                     *     A  B                  A  B
                     *   \ empty/              \      /
                     *    \    /                \(1) /
                     * (3) o--o (1)  or      (2) o--o (3)    (element number in brackets)
                     *    / (2)\                /    \
                     *   /      \              /empty \
                     *
                     * We perform a Type 1 swap in this case.
                     */
                    PerformT1Swap(pNodeA, pNodeB, all_indices);
                }
                else
                {
                    /*
                     * In this case, one of the nodes is contained in two elements
                     * and the other node is contained in three elements.
                     */
                    assert (   (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==3)
                            || (nodeA_elem_indices.size()==3 && nodeB_elem_indices.size()==2) );

                    /*
                     * Let node alpha be the node contained in two elements and node beta
                     * the node contained in three elements.
                     */
                    Node<SPACE_DIM>* p_node_alpha = (nodeA_elem_indices.size()==2) ? pNodeA : pNodeB;
                    Node<SPACE_DIM>* p_node_beta = (nodeA_elem_indices.size()==2) ? pNodeB : pNodeA;

                    // Get the set of elements containing node alpha and assert there are two such elements
                    std::set<unsigned> node_alpha_elem_indices = p_node_alpha->rGetContainingElementIndices();
                    assert(node_alpha_elem_indices.size() == 2u);

                    // Get a pointer to the first of these elements and assert it contains node alpha
                    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = mElements[*node_alpha_elem_indices.begin()];

                    unsigned node_alpha_local_index = p_element->GetNodeLocalIndex(p_node_alpha->GetIndex());
                    assert(node_alpha_local_index < UINT_MAX);

                    /*
                     * Get the local indices of the nodes neighbouring node alpha in this element.
                     * We add an extra temp in the second call to the % operator term to ensure we
                     * are operating on a positive number (otherwise the % operator may break).
                     */
                    unsigned temp = p_element->GetNumNodes();
                    unsigned node_alpha_local_index_before = (node_alpha_local_index+1)%temp;
                    unsigned node_alpha_local_index_after = (node_alpha_local_index+temp-1)%temp;

                    // Get pointers to these nodes and assert one of them is p_node_beta
                    Node<SPACE_DIM>* p_node1 = p_element->GetNode(node_alpha_local_index_before);
                    Node<SPACE_DIM>* p_node2 = p_element->GetNode(node_alpha_local_index_after);
                    assert(p_node1 == p_node_beta || p_node2 == p_node_beta);

                    // Get whichever of the nodes is NOT p_node_beta, call it node gamma
                    Node<SPACE_DIM>* p_node_gamma = (p_node1 == p_node_beta) ? p_node2 : p_node1;

                    // Get the set of elements containing nodes beta and gamma
                    std::set<unsigned> node_beta_elem_indices = p_node_beta->rGetContainingElementIndices();
                    std::set<unsigned> node_gamma_elem_indices = p_node_gamma->rGetContainingElementIndices();

                    // Form the set intersection
                    std::set<unsigned> intersection_indices, temp_set2;
                    std::set_intersection(node_beta_elem_indices.begin(), node_beta_elem_indices.end(),
                                   node_gamma_elem_indices.begin(), node_gamma_elem_indices.end(),
                                   std::inserter(temp_set2, temp_set2.begin()));
                    intersection_indices.swap(temp_set2); // temp_set2 will be deleted

                    // The number of intersecting indices must be 1 or 2
                    switch (intersection_indices.size())
                    {
                        case 1:
                        {
                            /*
                             * In this case, the node configuration looks like:
                             *
                             *     A  B                      A  B
                             *   \      /                  \      /
                             *    \ (1)/                    \(1) /
                             * (3) o--o (empty)  or  (empty) o--o (3)    (element number in brackets)
                             *    / (2)\                    /(2) \
                             *   /      \                  /      \
                             *
                             * We perform a Type 1 swap in this case.
                             */
                            PerformT1Swap(pNodeA, pNodeB, all_indices);
                            break;
                        }
                        case 2:
                        {
                            /*
                             * In this case, the node configuration looks like:
                             *
                             *     A  B             A  B
                             *   \                       /
                             *    \  (1)           (1)  /
                             * (3) o--o---   or  ---o--o (3)    (element number in brackets)
                             *    /  (2)           (2)  \
                             *   /                       \
                             *
                             * We perform a node merge in this case.
                             */
                            PerformNodeMerge(pNodeA, pNodeB);

                            // Remove the deleted node and re-index
                            RemoveDeletedNodes();
                            break;
                        }
                        default:
                            NEVER_REACHED;
                    }
                }
                break;
            }
            case 4:
            {
                /*
                 * In this case, the node configuration looks like:
                 *
                 *   \(1)/
                 *    \ / Node A
                 * (2) |   (4)      (element number in brackets)
                 *    / \ Node B
                 *   /(3)\
                 *
                 * We perform a Type 1 swap in this case.
                 */
                PerformT1Swap(pNodeA, pNodeB, all_indices);
                break;
            }
            default:
                // This can't happen
                NEVER_REACHED;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformNodeMerge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    // Sort the nodes by index
    unsigned nodeA_index = pNodeA->GetIndex();
    unsigned nodeB_index = pNodeB->GetIndex();

    unsigned lo_node_index = (nodeA_index < nodeB_index) ? nodeA_index : nodeB_index; // low index
    unsigned hi_node_index = (nodeA_index < nodeB_index) ? nodeB_index : nodeA_index; // high index

    // Get pointers to the nodes, sorted by index
    Node<SPACE_DIM>* p_lo_node = this->GetNode(lo_node_index);
    Node<SPACE_DIM>* p_hi_node = this->GetNode(hi_node_index);

    // Find the sets of elements containing each of the nodes, sorted by index
    std::set<unsigned> lo_node_elem_indices = p_lo_node->rGetContainingElementIndices();
    std::set<unsigned> hi_node_elem_indices = p_hi_node->rGetContainingElementIndices();

    // Move the low-index node to the mid-point
    c_vector<double, SPACE_DIM> node_midpoint = p_lo_node->rGetLocation() + 0.5*this->GetVectorFromAtoB(p_lo_node->rGetLocation(), p_hi_node->rGetLocation());
    c_vector<double, SPACE_DIM>& r_lo_node_location = p_lo_node->rGetModifiableLocation();
    r_lo_node_location = node_midpoint;

    // Update the elements previously containing the high-index node to contain the low-index node
    for (std::set<unsigned>::const_iterator it = hi_node_elem_indices.begin();
         it != hi_node_elem_indices.end();
         ++it)
    {
        // Find the local index of the high-index node in this element
        unsigned hi_node_local_index =  mElements[*it]->GetNodeLocalIndex(hi_node_index);
        assert(hi_node_local_index < UINT_MAX); // this element should contain the high-index node

        /*
         * If this element already contains the low-index node, then just remove the high-index node.
         * Otherwise replace it with the low-index node in the element and remove it from mNodes.
         */
        if (lo_node_elem_indices.count(*it) > 0)
        {
            mElements[*it]->DeleteNode(hi_node_local_index); // think this method removes the high-index node from mNodes
        }
        else
        {
            // Replace the high-index node with the low-index node in this element
            mElements[*it]->UpdateNode(hi_node_local_index, p_lo_node);
        }
    }
    assert(!(this->mNodes[hi_node_index]->IsDeleted()));
    this->mNodes[hi_node_index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(hi_node_index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformT1Swap(Node<SPACE_DIM>* pNodeA,
                                                       Node<SPACE_DIM>* pNodeB,
                                                       std::set<unsigned>& rElementsContainingNodes)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    /*
     * Restructure elements - remember to update nodes and elements.
     *
     * We need to implement the following changes:
     *
     * The element whose index was in nodeA_elem_indices but not nodeB_elem_indices,
     * and the element whose index was in nodeB_elem_indices but not nodeA_elem_indices,
     * should now both contain nodes A and B.
     *
     * The element whose index was in nodeA_elem_indices and nodeB_elem_indices, and which
     * node C lies inside, should now only contain node A.
     *
     * The element whose index was in nodeA_elem_indices and nodeB_elem_indices, and which
     * node D lies inside, should now only contain node B.
     *
     * Iterate over all elements involved and identify which element they are
     * in the diagram then update the nodes as necessary.
     *
     *   \(1)/
     *    \ / Node A
     * (2) |   (4)     elements in brackets
     *    / \ Node B
     *   /(3)\
     *
     */

    /*
     * Compute the locations of two new nodes C, D, placed on either side of the
     * edge E_old formed by nodes current_node and anticlockwise_node, such
     * that the edge E_new formed by the new nodes is the perpendicular bisector
     * of E_old, with |E_new| 'just larger' than mThresholdDistance.
     */

    // Find the sets of elements containing nodes A and B
    std::set<unsigned> nodeA_elem_indices = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> nodeB_elem_indices = pNodeB->rGetContainingElementIndices();

    double distance_between_nodes_CD = 1.5*mCellRearrangementThreshold; //this was 2*mCellRearrangementThreshold;

    c_vector<double, SPACE_DIM> nodeA_location = pNodeA->rGetLocation();
    c_vector<double, SPACE_DIM> nodeB_location = pNodeB->rGetLocation();

    c_vector<double, SPACE_DIM> a_to_b = this->GetVectorFromAtoB(nodeA_location, nodeB_location);
    c_vector<double, SPACE_DIM> perpendicular_vector;
    perpendicular_vector(0) = -a_to_b(1);
    perpendicular_vector(1) = a_to_b(0);

    c_vector<double, SPACE_DIM> c_to_d = distance_between_nodes_CD / norm_2(a_to_b) * perpendicular_vector;
    c_vector<double, SPACE_DIM> nodeC_location = nodeA_location + 0.5*a_to_b - 0.5*c_to_d;
    c_vector<double, SPACE_DIM> nodeD_location = nodeC_location + c_to_d;

    /*
     * Move node A to C and node B to D
     */

    c_vector<double, SPACE_DIM>& r_nodeA_location = pNodeA->rGetModifiableLocation();
    r_nodeA_location = nodeC_location;

    c_vector<double, SPACE_DIM>& r_nodeB_location = pNodeB->rGetModifiableLocation();
    r_nodeB_location = nodeD_location;

    for (std::set<unsigned>::const_iterator it = rElementsContainingNodes.begin();
         it != rElementsContainingNodes.end();
         ++it)
    {
        if (nodeA_elem_indices.find(*it) == nodeA_elem_indices.end()) // not in nodeA_elem_indices so element 3
        {
            /*
             * In this case the element index was not in
             * nodeA_elem_indices, so this element
             * does not contain node A. Therefore we must add node A
             * (which has been moved to node C) to this element.
             *
             * Locate local index of node B in element then add node A after
             * in anticlockwise direction.
             */

            unsigned nodeB_local_index =  mElements[*it]->GetNodeLocalIndex(pNodeB->GetIndex());
            assert(nodeB_local_index < UINT_MAX); // this element should contain node B

            mElements[*it]->AddNode(nodeB_local_index, pNodeA);
        }
        else if (nodeB_elem_indices.find(*it) == nodeB_elem_indices.end()) // not in nodeB_elem_indices so element 1
        {
            /*
             * In this case the element index was not in
             * nodeB_elem_indices, so this element
             * does not contain node B. Therefore we must add node B
             * (which has been moved to node D) to this element.
             *
             * Locate local index of node A in element then add node B after
             * in anticlockwise direction.
             */
            unsigned nodeA_local_index =  mElements[*it]->GetNodeLocalIndex(pNodeA->GetIndex());
            assert(nodeA_local_index < UINT_MAX); // this element should contain node A
            mElements[*it]->AddNode(nodeA_local_index, pNodeB);
        }
        else
        {
            /*
             * In this case the element index was in both nodeB_elem_indices and nodeB_elem_indices
             * so is element 2 or 4
             */

            /*
             * Locate local index of nodeA and nodeB and use the oredering to
             * identify the element, if nodeB_index > nodeA_index then element 4
             * and if nodeA_index > nodeB_index then element 2
             */
            unsigned nodeA_local_index =  mElements[*it]->GetNodeLocalIndex(pNodeA->GetIndex());
            assert(nodeA_local_index < UINT_MAX); // this element should contain node A

            unsigned nodeB_local_index =  mElements[*it]->GetNodeLocalIndex(pNodeB->GetIndex());
            assert(nodeB_local_index < UINT_MAX); // this element should contain node B

            unsigned nodeB_local_index_plus_one = (nodeB_local_index + 1)%(mElements[*it]->GetNumNodes());

            if (nodeA_local_index == nodeB_local_index_plus_one)
            {
                /*
                 * In this case the local index of nodeA is the local index of
                 * nodeB plus one so we are in element 2 so we remove nodeB
                 */
                 mElements[*it]->DeleteNode(nodeB_local_index);
            }
            else
            {
                assert(nodeB_local_index == (nodeA_local_index + 1)%(mElements[*it]->GetNumNodes())); // as A and B are next to each other
                /*
                 * In this case the local index of nodeA is the local index of
                 * nodeB minus one so we are in element 4 so we remove nodeA
                 */
                 mElements[*it]->DeleteNode(nodeA_local_index);
            }
        }
    }

    // Sort out boundary nodes
    if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
    {
        if (pNodeA->GetNumContainingElements()==3)
        {
            pNodeA->SetAsBoundaryNode(false);
        }
        else
        {
            pNodeA->SetAsBoundaryNode(true);
        }
        if (pNodeB->GetNumContainingElements()==3)
        {
            pNodeB->SetAsBoundaryNode(false);
        }
        else
        {
            pNodeB->SetAsBoundaryNode(true);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformT2Swap(VertexElement<ELEMENT_DIM,SPACE_DIM>& rElement)
{
   // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    /* For removing a small triangle
     *
     *  \
     *   \
     *   |\
     *   | \
     *   |  \_ _ _ _ _
     *   |  /
     *   | /
     *   |/
     *   /
     *  /
     *
     * To

     *  \
     *   \
     *    \
     *     \
     *      \_ _ _ _ _
     *      /
     *     /
     *    /
     *   /
     *  /
     *
     */

     // Assert that the triangle element has only three nodes (!)

     assert(rElement.GetNumNodes() == 3u);

     c_vector<double, SPACE_DIM>& new_node_location = rElement.GetNode(0)->rGetModifiableLocation();
     new_node_location = GetCentroidOfElement(rElement.GetIndex());

     c_vector<unsigned, 3> neighbouring_elem_nums;

     for (unsigned i=0; i<3; i++)
     {
         std::set<unsigned> elements_of_node_a = rElement.GetNode((i+1)%3)->rGetContainingElementIndices();
         std::set<unsigned> elements_of_node_b = rElement.GetNode((i+2)%3)->rGetContainingElementIndices();

         std::set<unsigned> common_elements;
         std::set_intersection(elements_of_node_a.begin(), elements_of_node_a.end(),
                                elements_of_node_b.begin(), elements_of_node_b.end(),
                                std::inserter(common_elements, common_elements.begin()));

         assert(common_elements.size() == 2u);
         common_elements.erase(rElement.GetIndex());
         assert(common_elements.size() == 1u);

         neighbouring_elem_nums(i) = *(common_elements.begin());
     }

     // Extract the neighbouring elements
     VertexElement<ELEMENT_DIM,SPACE_DIM>* p_neighbouring_element_0 = this->GetElement(neighbouring_elem_nums(0));
     VertexElement<ELEMENT_DIM,SPACE_DIM>* p_neighbouring_element_1 = this->GetElement(neighbouring_elem_nums(1));
     VertexElement<ELEMENT_DIM,SPACE_DIM>* p_neighbouring_element_2 = this->GetElement(neighbouring_elem_nums(2));

     // Need to check that none of the neighbouring elements are triangles
     if (    (p_neighbouring_element_0->GetNumNodes() > 3u)
         && (p_neighbouring_element_1->GetNumNodes() > 3u)
         && (p_neighbouring_element_2->GetNumNodes() > 3u) )
     {
         // Neighbour 0 - replace node 1 with node 0, delete node 2
         p_neighbouring_element_0->ReplaceNode(rElement.GetNode(1), rElement.GetNode(0));
         p_neighbouring_element_0->DeleteNode(p_neighbouring_element_0->GetNodeLocalIndex(rElement.GetNodeGlobalIndex(2)));

         // Neighbour 1 - delete node 2
         p_neighbouring_element_1->DeleteNode(p_neighbouring_element_1->GetNodeLocalIndex(rElement.GetNodeGlobalIndex(2)));

         // Neighbour 2 - delete node 1
         p_neighbouring_element_2->DeleteNode(p_neighbouring_element_2->GetNodeLocalIndex(rElement.GetNodeGlobalIndex(1)));

         // Also have to mark pElement, pElement->GetNode(1), pElement->GetNode(2) as deleted.
         mDeletedNodeIndices.push_back(rElement.GetNodeGlobalIndex(1));
         mDeletedNodeIndices.push_back(rElement.GetNodeGlobalIndex(2));
         rElement.GetNode(1)->MarkAsDeleted();
         rElement.GetNode(2)->MarkAsDeleted();

         mDeletedElementIndices.push_back(rElement.GetIndex());
         rElement.MarkAsDeleted();
     }
     else
     {
        EXCEPTION("One of the neighbours of a small triangular element is also a triangle - dealing with this has not been implemented yet");
     }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElementAlongGivenAxis(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, c_vector<double, SPACE_DIM> axisOfDivision)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Get the centroid of the element
    c_vector<double, SPACE_DIM> centroid = GetCentroidOfElement(pElement->GetIndex());

    // Create a vector perpendicular to the axis of division
    c_vector<double, SPACE_DIM> perp_axis;
    perp_axis(0) = -axisOfDivision(1);
    perp_axis(1) = axisOfDivision(0);

    /*
     * Find which edges the axis of division crosses by finding any node
     * that lies on the opposite side of the axis of division to its next
     * neighbour.
     */
    unsigned num_nodes = pElement->GetNumNodes();
    std::vector<unsigned> intersecting_nodes;
    for (unsigned i=0; i<num_nodes; i++)
    {
        bool is_current_node_on_left = (inner_prod(GetVectorFromAtoB(pElement->GetNodeLocation(i), centroid), perp_axis) >= 0);
        bool is_next_node_on_left = (inner_prod(GetVectorFromAtoB(pElement->GetNodeLocation((i+1)%num_nodes), centroid), perp_axis) >= 0);

        if (is_current_node_on_left != is_next_node_on_left)
        {
            intersecting_nodes.push_back(i);
        }
    }

    // If the axis of division does not cross two edges then we cannot proceed
    if (intersecting_nodes.size() != 2)
    {
        EXCEPTION("Cannot proceed with element division: the given axis of division does not cross two edges of the element");
    }

    std::vector<unsigned> division_node_global_indices;
    unsigned nodes_added = 0;

    // Find the intersections between the axis of division and the element edges
    for (unsigned i=0; i<intersecting_nodes.size(); i++)
    {
        /*
         * Get pointers to the nodes forming the edge into which one new node will be inserted.
         *
         * Note that when we use the first entry of intersecting_nodes to add a node,
         * we change the local index of the second entry of intersecting_nodes in
         * pElement, so must account for this by moving one entry further on.
         */
        Node<SPACE_DIM>* p_node_A = pElement->GetNode((intersecting_nodes[i]+nodes_added)%pElement->GetNumNodes());
        Node<SPACE_DIM>* p_node_B = pElement->GetNode((intersecting_nodes[i]+nodes_added+1)%pElement->GetNumNodes());

        c_vector<double, SPACE_DIM> position_a = p_node_A->rGetLocation();
        c_vector<double, SPACE_DIM> position_b = p_node_B->rGetLocation();

        c_vector<double, SPACE_DIM> a_to_b = GetVectorFromAtoB(position_a, position_b);

        // Find the location of the intersection
        double determinant = a_to_b[0]*axisOfDivision[1] - a_to_b[1]*axisOfDivision[0];

        c_vector<double, SPACE_DIM> moved_centroid = position_a + GetVectorFromAtoB(position_a, centroid); // allow for periodicity

        double alpha = (moved_centroid[0]*a_to_b[1] - position_a[0]*a_to_b[1]
                        -moved_centroid[1]*a_to_b[0] + position_a[1]*a_to_b[0])/determinant;

        c_vector<double, SPACE_DIM> intersection = moved_centroid + alpha*axisOfDivision;

        /*
         * If then new node is too close to one of the edge nodes, then reposition it
         * a distance 2*mCellRearrangementThreshold further along the edge.
         */
        c_vector<double, SPACE_DIM> a_to_intersection = this->GetVectorFromAtoB(position_a, intersection);
        if (norm_2(a_to_intersection) < mCellRearrangementThreshold)
        {
            intersection = position_a + 2*mCellRearrangementThreshold*a_to_b/norm_2(a_to_b);
        }

        c_vector<double, SPACE_DIM> b_to_intersection = this->GetVectorFromAtoB(position_b, intersection);
        if (norm_2(b_to_intersection) < mCellRearrangementThreshold)
        {
            intersection = position_b - 2*mCellRearrangementThreshold*a_to_b/norm_2(a_to_b);
        }

        bool is_boundary = false;
        if (p_node_A->IsBoundaryNode() && p_node_B->IsBoundaryNode())
        {
            is_boundary = true;
        }

        // Add a new node to the mesh at the location of the intersection
        unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(0, is_boundary, intersection[0], intersection[1]));
        nodes_added++;

        // Now make sure node is added to neighbouring elements

        // Find the indices of the elements owned by each node on the edge into which one new node will be inserted
        std::set<unsigned> elems_containing_node1 = p_node_A->rGetContainingElementIndices();
        std::set<unsigned> elems_containing_node2 = p_node_B->rGetContainingElementIndices();

        // Find common elements
        std::set<unsigned> shared_elements;
        std::set_intersection(elems_containing_node1.begin(),
                              elems_containing_node1.end(),
                              elems_containing_node2.begin(),
                              elems_containing_node2.end(),
                              std::inserter(shared_elements, shared_elements.begin()));

        // Iterate over common elements
        for (std::set<unsigned>::iterator iter = shared_elements.begin();
             iter != shared_elements.end();
             ++iter)
        {
            // Find which node has the lower local index in this element
            unsigned local_indexA = GetElement(*iter)->GetNodeLocalIndex(p_node_A->GetIndex());
            unsigned local_indexB = GetElement(*iter)->GetNodeLocalIndex(p_node_B->GetIndex());

            unsigned index = local_indexB;
            if (local_indexB > local_indexA)
            {
                index = local_indexA;
            }
            if ( (local_indexA == 0) && (local_indexB == GetElement(*iter)->GetNumNodes()-1))
            {
                index = local_indexB;
            }
            if ( (local_indexB == 0) && (local_indexA == GetElement(*iter)->GetNumNodes()-1))
            {
                index = local_indexA;
            }
            // Add new node to this element
            GetElement(*iter)->AddNode(index, this->GetNode(new_node_global_index));
        }
        // Store index of new node
        division_node_global_indices.push_back(new_node_global_index);

    }

    // Now call DivideElement() to divide the element using the new nodes
    unsigned new_element_index = DivideElement(pElement, pElement->GetNodeLocalIndex(division_node_global_indices[0]), pElement->GetNodeLocalIndex(division_node_global_indices[1]));
    return new_element_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElementAlongShortAxis(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Find the short axis of the element
    c_vector<double, SPACE_DIM> short_axis = GetShortAxisOfElement(pElement->GetIndex());

    unsigned new_element_index = DivideElementAlongGivenAxis(pElement, short_axis);
    return new_element_index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned nodeAIndex, unsigned nodeBIndex)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Sort nodeA and nodeB such that nodeBIndex > nodeAindex
    assert(nodeBIndex != nodeAIndex);

    unsigned node1_index = (nodeAIndex < nodeBIndex) ? nodeAIndex : nodeBIndex; // low index
    unsigned node2_index = (nodeAIndex < nodeBIndex) ? nodeBIndex : nodeAIndex; // high index

    // Store the number of nodes in the element (this changes when nodes are deleted from the element)
    unsigned num_nodes = pElement->GetNumNodes();

    // Copy the nodes in this element
    std::vector<Node<SPACE_DIM>*> nodes_elem;
    for (unsigned i=0; i<num_nodes; i++)
    {
        nodes_elem.push_back(pElement->GetNode(i));
    }

    // Get the index of the new element
    unsigned new_element_index;
    if (mDeletedElementIndices.empty())
    {
        new_element_index = mElements.size();
    }
    else
    {
        new_element_index = mDeletedElementIndices.back();
        mDeletedElementIndices.pop_back();
        delete mElements[new_element_index];
    }

    // Add the new element to the mesh
    AddElement(new VertexElement<ELEMENT_DIM,SPACE_DIM>(new_element_index, nodes_elem));

    /**
     * Remove the correct nodes from each element. Choose the
     * original element to be below (in the y direction) the
     * new element, so as to keep stem cells at the bottom.
     *
     * \todo Woah! Why does this comment refer to stem cells?
     * Is this code crypt-dependent? It shouldn't be... see also #1099
     */

    /// Find lowest element \todo this could be more efficient
    double height_midpoint_1 = 0.0;
    double height_midpoint_2 = 0.0;
    unsigned counter_1 = 0;
    unsigned counter_2 = 0;

    for (unsigned i=0; i<num_nodes; i++)
    {
        if (i>=node1_index && i<=node2_index)
        {
            height_midpoint_1 += pElement->GetNode(i)->rGetLocation()[1];
            counter_1++;
        }
        if (i<=node1_index || i>=node2_index)
        {
            height_midpoint_2 += pElement->GetNode(i)->rGetLocation()[1];
            counter_2++;
        }
    }
    height_midpoint_1 /= (double)counter_1;
    height_midpoint_2 /= (double)counter_2;

    for (unsigned i=num_nodes; i>0; i--)
    {
        if (i-1 < node1_index || i-1 > node2_index)
        {
            if (height_midpoint_1 < height_midpoint_2)
            {
                pElement->DeleteNode(i-1);
            }
            else
            {
                mElements[new_element_index]->DeleteNode(i-1);
            }
        }
        else if (i-1 > node1_index && i-1 < node2_index)
        {
            if (height_midpoint_1 < height_midpoint_2)
            {
                mElements[new_element_index]->DeleteNode(i-1);
            }
            else
            {
                pElement->DeleteNode(i-1);
            }
        }
    }

    return new_element_index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader)
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
        unsigned is_boundary_node = (unsigned) node_data[2];
        node_data.pop_back();
        this->mNodes.push_back(new Node<SPACE_DIM>(i, node_data, is_boundary_node));
    }

    rMeshReader.Reset();

    // Reserve memory for nodes
    mElements.reserve(rMeshReader.GetNumElements());

    // Add elements
    for (unsigned elem_index=0; elem_index<num_elements; elem_index++)
    {
        ElementData element_data = rMeshReader.GetNextElementData();
        std::vector<Node<SPACE_DIM>*> nodes;

        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j=0; j<num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        VertexElement<ELEMENT_DIM,SPACE_DIM>* p_element = new VertexElement<ELEMENT_DIM,SPACE_DIM>(elem_index, nodes);
        mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = element_data.AttributeValue;
            p_element->SetRegion(attribute_value);
        }
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::ElementIncludesPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Initialise boolean
    bool element_includes_point = false;

    // Get the element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(elementIndex);
    unsigned num_nodes = p_element->GetNumNodes();

    // Remap the origin to the first vertex to allow alternative distance metrics to be used in subclasses
    c_vector<double, SPACE_DIM> first_vertex = p_element->GetNodeLocation(0);

    c_vector<double, SPACE_DIM> test_point =  GetVectorFromAtoB(first_vertex,rTestPoint);

    // Loop over edges of the element
    for (unsigned local_index=0; local_index<num_nodes; local_index++)
    {

        // Get the end points of this edge
        // Remap to the origin to allow alternative distance metrics to be used in subclasses
        c_vector<double, SPACE_DIM> vertexA = GetVectorFromAtoB(first_vertex, p_element->GetNodeLocation(local_index));
        c_vector<double, SPACE_DIM> vertexB = GetVectorFromAtoB(first_vertex, p_element->GetNodeLocation((local_index+1)%num_nodes));


        // Check if this edge crosses the ray running out horizontally (increasing x, fixed y) from the test point
        c_vector<double, SPACE_DIM> vector_a_to_point = GetVectorFromAtoB(vertexA, test_point);
        c_vector<double, SPACE_DIM> vector_b_to_point = GetVectorFromAtoB(vertexB, test_point);
        c_vector<double, SPACE_DIM> vector_a_to_b = GetVectorFromAtoB(vertexA, vertexB);

        // Pathological case - test point coincides with vertexA or vertexB
        if (    (norm_2(vector_a_to_point) < DBL_EPSILON)
             || (norm_2(vector_b_to_point) < DBL_EPSILON) )
        {
            return false;
        }

        // Pathological case - ray coincides with horizontal edge
        if ( (fabs(vector_a_to_b[1]) < DBL_EPSILON) &&
             (fabs(vector_a_to_point[1]) < DBL_EPSILON) &&
             (fabs(vector_b_to_point[1]) < DBL_EPSILON) )
        {
            if ( (vector_a_to_point[0]>0) != (vector_b_to_point[0]>0) )
            {
                return false;
            }
        }

        /// \todo Need to carefully check all pathological cases (see #933)

        // Non-pathological case
        // A and B on different sides of the line y = test_point[1]
        if ( (vertexA[1] > test_point[1]) != (vertexB[1] > test_point[1]) )
        {
            // intersection of y=test_point[1] and vector_a_to_b is on the Right of test_point
            if (test_point[0] < vertexA[0] + vector_a_to_b[0]*vector_a_to_point[1]/vector_a_to_b[1])
            {
                element_includes_point = !element_includes_point;
            }
        }
    }
    return element_includes_point;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetLocalIndexForElementEdgeClosestToPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Get the element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(elementIndex);
    unsigned num_nodes = p_element->GetNumNodes();

    double min_squared_distance = DBL_MAX;
    unsigned min_distance_edge_index = UINT_MAX;

    // Loop over edges of the element
    for (unsigned local_index=0; local_index<num_nodes; local_index++)
    {
        // Get the end points of this edge
        c_vector<double, SPACE_DIM> vertexA = p_element->GetNodeLocation(local_index);
        c_vector<double, SPACE_DIM> vertexB = p_element->GetNodeLocation((local_index+1)%num_nodes);

        c_vector<double, SPACE_DIM> vector_a_to_point = this->GetVectorFromAtoB(vertexA, rTestPoint);
        c_vector<double, SPACE_DIM> vector_a_to_b = this->GetVectorFromAtoB(vertexA, vertexB);

        c_vector<double, SPACE_DIM> edge_ab_unit_vector = vector_a_to_b/norm_2(vector_a_to_b);
        double distance_parallel_to_edge = inner_prod(vector_a_to_point, edge_ab_unit_vector);

        double squared_distance_normal_to_edge = SmallPow(norm_2(vector_a_to_point), 2) - SmallPow(distance_parallel_to_edge, 2);

        if (squared_distance_normal_to_edge < min_squared_distance)
        {
            min_squared_distance = squared_distance_normal_to_edge;
            min_distance_edge_index = local_index;
        }
    }

    return min_distance_edge_index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::MoveOverlappingNodeOntoEdgeOfElement(Node<SPACE_DIM>* pNode, unsigned elementIndex)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Get the local index of the node in the element after which the new node is to be added
    unsigned local_index = GetLocalIndexForElementEdgeClosestToPoint(pNode->rGetLocation(), elementIndex);

    // Move the new node back onto the edge
    c_vector<double, SPACE_DIM> node_location = pNode->rGetModifiableLocation();

    // Get element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(elementIndex);
    unsigned num_nodes = p_element->GetNumNodes();

    c_vector<double, SPACE_DIM> vertexA = p_element->GetNodeLocation(local_index);
    c_vector<double, SPACE_DIM> vertexB = p_element->GetNodeLocation((local_index+1)%num_nodes);

    c_vector<double, SPACE_DIM> vector_a_to_point = this->GetVectorFromAtoB(vertexA, node_location);
    c_vector<double, SPACE_DIM> vector_a_to_b = this->GetVectorFromAtoB(vertexA, vertexB);

    c_vector<double, SPACE_DIM> edge_ab_unit_vector = vector_a_to_b/norm_2(vector_a_to_b);

    pNode->rGetModifiableLocation() = vertexA + edge_ab_unit_vector*inner_prod(vector_a_to_point, edge_ab_unit_vector);

    // Add the new node to the element (this also updates the node)
    GetElement(elementIndex)->AddNode(local_index, pNode);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::Translate(
    const double xMovement,
    const double yMovement,
    const double zMovement)
{
    c_vector<double, SPACE_DIM> displacement;

    switch (SPACE_DIM)
    {
        case 3:
            displacement[2] = zMovement;
        case 2:
            displacement[1] = yMovement;
        case 1:
            displacement[0] = xMovement;
    }

    Translate(displacement);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::Translate(c_vector<double, SPACE_DIM>& rDisplacement)
{
    unsigned num_nodes = this->GetNumAllNodes();

    for (unsigned i=0; i<num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = this->mNodes[i]->rGetModifiableLocation();
        r_location += rDisplacement;
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class VertexMesh<1,1>;
template class VertexMesh<1,2>;
template class VertexMesh<1,3>;
template class VertexMesh<2,2>;
template class VertexMesh<2,3>;
template class VertexMesh<3,3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexMesh);
