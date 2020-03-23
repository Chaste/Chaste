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

#include "MutableVertexMesh.hpp"

#include "LogFile.hpp"
#include "UblasCustomFunctions.hpp"
#include "Warnings.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::MutableVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                                             std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertexElements,
                                                             double cellRearrangementThreshold,
                                                             double t2Threshold,
                                                             double cellRearrangementRatio,
                                                             double protorosetteFormationProbability,
                                                             double protorosetteResolutionProbabilityPerTimestep,
                                                             double rosetteResolutionProbabilityPerTimestep)
        : mCellRearrangementThreshold(cellRearrangementThreshold),
          mCellRearrangementRatio(cellRearrangementRatio),
          mT2Threshold(t2Threshold),
          mProtorosetteFormationProbability(protorosetteFormationProbability),
          mProtorosetteResolutionProbabilityPerTimestep(protorosetteResolutionProbabilityPerTimestep),
          mRosetteResolutionProbabilityPerTimestep(rosetteResolutionProbabilityPerTimestep),
          mCheckForInternalIntersections(false),
          mDistanceForT3SwapChecking(5.0)
{
    // Threshold parameters must be strictly positive
    assert(cellRearrangementThreshold > 0.0);
    assert(t2Threshold > 0.0);
    assert(protorosetteFormationProbability >= 0.0);
    assert(protorosetteFormationProbability <= 1.0);
    assert(protorosetteResolutionProbabilityPerTimestep >= 0.0);
    assert(protorosetteResolutionProbabilityPerTimestep <= 1.0);
    assert(rosetteResolutionProbabilityPerTimestep >= 0.0);
    assert(rosetteResolutionProbabilityPerTimestep <= 1.0);

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
        this->mElements.push_back(p_temp_vertex_element);
    }

    // If in 3D, then also populate mFaces
    if (SPACE_DIM == 3)
    {
        // Use a std::set to keep track of which faces have been added to mFaces
        std::set<unsigned> faces_counted;

        // Loop over mElements
        for (unsigned elem_index=0; elem_index<this->mElements.size(); elem_index++)
        {
            // Loop over faces of this element
            for (unsigned face_index=0; face_index<this->mElements[elem_index]->GetNumFaces(); face_index++)
            {
                VertexElement<ELEMENT_DIM-1, SPACE_DIM>* p_face = this->mElements[elem_index]->GetFace(face_index);

                // If this face is not already contained in mFaces, then add it and update faces_counted
                if (faces_counted.find(p_face->GetIndex()) == faces_counted.end())
                {
                    this->mFaces.push_back(p_face);
                    faces_counted.insert(p_face->GetIndex());
                }
            }
        }
    }

    // Register elements with nodes
    for (unsigned index=0; index<this->mElements.size(); index++)
    {
        VertexElement<ELEMENT_DIM,SPACE_DIM>* p_temp_vertex_element = this->mElements[index];
        for (unsigned node_index=0; node_index<p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<SPACE_DIM>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }

    this->GenerateEdgesFromElements(vertexElements);

    this->mMeshChangesDuringSimulation = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::MutableVertexMesh()
    : mCellRearrangementThreshold(0.01),
      mCellRearrangementRatio(1.5),
      mT2Threshold(0.001),
      mProtorosetteFormationProbability(0.0),
      mProtorosetteResolutionProbabilityPerTimestep(0.0),
      mRosetteResolutionProbabilityPerTimestep(0.0),
      mCheckForInternalIntersections(false),
      mDistanceForT3SwapChecking(5.0)
{
    // Note that the member variables initialised above will be overwritten as soon as archiving is complete
    this->mMeshChangesDuringSimulation = true;
    Clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::~MutableVertexMesh()
{
    Clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCellRearrangementThreshold() const
{
    return mCellRearrangementThreshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetT2Threshold() const
{
    return mT2Threshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCellRearrangementRatio() const
{
    return mCellRearrangementRatio;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetProtorosetteFormationProbability() const
{
    return this->mProtorosetteFormationProbability;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetProtorosetteResolutionProbabilityPerTimestep() const
{
    return this->mProtorosetteResolutionProbabilityPerTimestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetRosetteResolutionProbabilityPerTimestep() const
{
    return this->mRosetteResolutionProbabilityPerTimestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetDistanceForT3SwapChecking(double distanceForT3SwapChecking)
{
    mDistanceForT3SwapChecking = distanceForT3SwapChecking;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetDistanceForT3SwapChecking() const
{
    return mDistanceForT3SwapChecking;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCheckForInternalIntersections() const
{
    return mCheckForInternalIntersections;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetCellRearrangementThreshold(double cellRearrangementThreshold)
{
    mCellRearrangementThreshold = cellRearrangementThreshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetT2Threshold(double t2Threshold)
{
    mT2Threshold = t2Threshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetCellRearrangementRatio(double cellRearrangementRatio)
{
    mCellRearrangementRatio = cellRearrangementRatio;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetProtorosetteFormationProbability(double protorosetteFormationProbability)
{
    // Check that the new value is in [0,1]
    if (protorosetteFormationProbability < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }
    if (protorosetteFormationProbability > 1.0)
    {
        EXCEPTION("Attempting to assign a probability greater than one.");
    }

    // Assign the new value
    mProtorosetteFormationProbability = protorosetteFormationProbability;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetProtorosetteResolutionProbabilityPerTimestep(double protorosetteResolutionProbabilityPerTimestep)
{
    // Check that the new value is in [0,1]
    if (protorosetteResolutionProbabilityPerTimestep < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }
    if (protorosetteResolutionProbabilityPerTimestep > 1.0)
    {
        EXCEPTION("Attempting to assign a probability greater than one.");
    }

    // Assign the new value
    mProtorosetteResolutionProbabilityPerTimestep = protorosetteResolutionProbabilityPerTimestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetRosetteResolutionProbabilityPerTimestep(double rosetteResolutionProbabilityPerTimestep)
{
    // Check that the new value is in [0,1]
    if (rosetteResolutionProbabilityPerTimestep < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }
    if (rosetteResolutionProbabilityPerTimestep > 1.0)
    {
        EXCEPTION("Attempting to assign a probability greater than one.");
    }

    // Assign the new value
    mRosetteResolutionProbabilityPerTimestep = rosetteResolutionProbabilityPerTimestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetCheckForInternalIntersections(bool checkForInternalIntersections)
{
    mCheckForInternalIntersections = checkForInternalIntersections;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    mDeletedNodeIndices.clear();
    mDeletedElementIndices.clear();

    VertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size() - mDeletedNodeIndices.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return this->mElements.size() - mDeletedElementIndices.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetLastT2SwapLocation()
{
    return mLastT2SwapLocation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNewNode)
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
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::AddElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pNewElement)
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
    pNewElement->SetEdgeHelper(&this->mEdges);
    pNewElement->BuildEdges();
    return pNewElement->GetIndex();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point)
{
    this->mNodes[nodeIndex]->SetPoint(point);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElementAlongGivenAxis(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement,
                                                                                c_vector<double, SPACE_DIM> axisOfDivision,
                                                                                bool placeOriginalElementBelow)
{
    assert(SPACE_DIM == 2);                // LCOV_EXCL_LINE
    assert(ELEMENT_DIM == SPACE_DIM);    // LCOV_EXCL_LINE

    //Storing element edge map prior to division. This is needed for division recording.
    std::vector<unsigned int> edgeIds;
    for (unsigned i = 0; i < pElement->GetNumEdges(); i++)
    {
        edgeIds.push_back(pElement->GetEdge(i)->GetIndex());
    }

    // Get the centroid of the element
    c_vector<double, SPACE_DIM> centroid = this->GetCentroidOfElement(pElement->GetIndex());
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
    bool is_current_node_on_left = (inner_prod(this->GetVectorFromAtoB(pElement->GetNodeLocation(0), centroid), perp_axis) >= 0);
    for (unsigned i=0; i<num_nodes; i++)
    {
        bool is_next_node_on_left = (inner_prod(this->GetVectorFromAtoB(pElement->GetNodeLocation((i+1)%num_nodes), centroid), perp_axis) >= 0);
        if (is_current_node_on_left != is_next_node_on_left)
        {
            intersecting_nodes.push_back(i);
        }
        is_current_node_on_left = is_next_node_on_left;
    }

    // If the axis of division does not cross two edges then we cannot proceed
    if (intersecting_nodes.size() != 2)
    {
        EXCEPTION("Cannot proceed with element division: the given axis of division does not cross two edges of the element");
    }

    std::vector<unsigned> division_node_global_indices;
    unsigned nodes_added = 0;

    //Keeps track of elements (first) whose edge (second) had been split
    std::vector<std::pair<VertexElement<ELEMENT_DIM, SPACE_DIM>*, unsigned int> > edge_split_pairs;
    std::vector<double> relative_new_node;
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

        // Find the indices of the elements owned by each node on the edge into which one new node will be inserted
        std::set<unsigned> elems_containing_node_A = p_node_A->rGetContainingElementIndices();
        std::set<unsigned> elems_containing_node_B = p_node_B->rGetContainingElementIndices();

        c_vector<double, SPACE_DIM> position_a = p_node_A->rGetLocation();
        c_vector<double, SPACE_DIM> position_b = p_node_B->rGetLocation();
        c_vector<double, SPACE_DIM> a_to_b = this->GetVectorFromAtoB(position_a, position_b);

        c_vector<double, SPACE_DIM> intersection;

        if (norm_2(a_to_b) < 2.0*mCellRearrangementRatio*mCellRearrangementThreshold)
        {
            WARNING("Edge is too small for normal division; putting node in the middle of a and b. There may be T1 swaps straight away.");
            ///\todo or should we move a and b apart, it may interfere with neighbouring edges? (see #1399 and #2401)
            intersection = position_a + 0.5*a_to_b;
        }
        else
        {
            // Find the location of the intersection
            double determinant = a_to_b[0]*axisOfDivision[1] - a_to_b[1]*axisOfDivision[0];

            // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
            c_vector<double, SPACE_DIM> moved_centroid;
            moved_centroid = position_a + this->GetVectorFromAtoB(position_a, centroid);

            double alpha = (moved_centroid[0]*a_to_b[1] - position_a[0]*a_to_b[1]
                            -moved_centroid[1]*a_to_b[0] + position_a[1]*a_to_b[0])/determinant;

            intersection = moved_centroid + alpha*axisOfDivision;

            /*
             * If then new node is too close to one of the edge nodes, then reposition it
             * a distance mCellRearrangementRatio*mCellRearrangementThreshold further along the edge.
             */
            c_vector<double, SPACE_DIM> a_to_intersection = this->GetVectorFromAtoB(position_a, intersection);
            if (norm_2(a_to_intersection) < mCellRearrangementThreshold)
            {
                intersection = position_a + mCellRearrangementRatio*mCellRearrangementThreshold*a_to_b/norm_2(a_to_b);
            }

            c_vector<double, SPACE_DIM> b_to_intersection = this->GetVectorFromAtoB(position_b, intersection);
            if (norm_2(b_to_intersection) < mCellRearrangementThreshold)
            {
                assert(norm_2(a_to_intersection) > mCellRearrangementThreshold); // to prevent moving intersection back to original position

                intersection = position_b - mCellRearrangementRatio*mCellRearrangementThreshold*a_to_b/norm_2(a_to_b);
            }
        }

        /*
         * The new node is boundary node if the 2 nodes are boundary nodes and the elements don't look like
         *   ___A___
         *  |   |   |
         *  |___|___|
         *      B
         */
        bool is_boundary = false;
        if (p_node_A->IsBoundaryNode() && p_node_B->IsBoundaryNode())
        {
            if (elems_containing_node_A.size() != 2 ||
                elems_containing_node_B.size() != 2 ||
                elems_containing_node_A != elems_containing_node_B)
            {
                is_boundary = true;
            }
        }

        // Add a new node to the mesh at the location of the intersection
        unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(0, is_boundary, intersection[0], intersection[1]));
        nodes_added++;

        // Now make sure the new node is added to all neighbouring elements

        // Find common elements
        std::set<unsigned> shared_elements;
        std::set_intersection(elems_containing_node_A.begin(),
                              elems_containing_node_A.end(),
                              elems_containing_node_B.begin(),
                              elems_containing_node_B.end(),
                              std::inserter(shared_elements, shared_elements.begin()));

        // Iterate over common elements, including the element to be divided
        unsigned node_A_index = p_node_A->GetIndex();
        unsigned node_B_index = p_node_B->GetIndex();
        bool original_element = false;
        for (std::set<unsigned>::iterator iter = shared_elements.begin();
             iter != shared_elements.end();
             ++iter)
        {
            VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(*iter);
            original_element = p_element->GetIndex()==pElement->GetIndex();

            // Find which node has the lower local index in this element
            unsigned local_indexA = p_element->GetNodeLocalIndex(node_A_index);
            unsigned local_indexB = p_element->GetNodeLocalIndex(node_B_index);

            unsigned index = local_indexB;

            // If node B has a higher index then use node A's index...
            if (local_indexB > local_indexA)
            {
                index = local_indexA;

                // ...unless nodes A and B share the element's last edge
                if ((local_indexA == 0) && (local_indexB == p_element->GetNumNodes()-1))
                {
                    index = local_indexB;
                }
            }
            else if ((local_indexB == 0) && (local_indexA == p_element->GetNumNodes()-1))
            {
                // ...otherwise use node B's index, unless nodes A and B share the element's last edge
                index = local_indexA;
            }

            // Record the edge split in neighbouring elements
            if (!original_element&&mTrackMeshOperations)
            {
                const unsigned int n_edges = p_element->GetNumEdges();
                const unsigned int nextIndex = (index+1)%n_edges;

                auto prev_node = p_element->GetNode(index)->rGetLocation();
                auto next_node = p_element->GetNode(nextIndex)->rGetLocation();
                auto curr_node = this->GetNode(new_node_global_index)->rGetLocation();

                c_vector<double, SPACE_DIM> last_to_next = this->GetVectorFromAtoB(prev_node, next_node);
                c_vector<double, SPACE_DIM> last_to_curr = this->GetVectorFromAtoB(prev_node, curr_node);

                double old_distance = norm_2(last_to_next);
                double prev_curr_distance = norm_2(last_to_curr);
                if (old_distance<prev_curr_distance)
                {
                    EXCEPTION("New node in split edge is not between a pair old nodes");
                }
                double theta = prev_curr_distance/old_distance;
                relative_new_node.push_back(theta);
                edge_split_pairs.push_back(std::pair<VertexElement<ELEMENT_DIM,SPACE_DIM>*, unsigned int>(p_element,index));
            }
            // Add new node to this element
            p_element->AddNode(this->GetNode(new_node_global_index), index);
        }

        // Store index of new node
        division_node_global_indices.push_back(new_node_global_index);
    }

    // Now call DivideElement() to divide the element using the new nodes
    unsigned new_element_index = DivideElement(pElement,
                                               pElement->GetNodeLocalIndex(division_node_global_indices[0]),
                                               pElement->GetNodeLocalIndex(division_node_global_indices[1]),
                                               placeOriginalElementBelow);

    //Record cell division info
    CellDivisionInfo<SPACE_DIM> division_info;
    division_info.mLocation = centroid;
    division_info.mDaughterLocation1 = this->GetCentroidOfElement(pElement->GetIndex());
    c_vector<double, SPACE_DIM> long_axis = this->GetShortAxisOfElement(pElement->GetIndex());
    division_info.mDaughterLongAxis1(0) = -long_axis(1);
    division_info.mDaughterLongAxis1(1) = long_axis(0);

    division_info.mDaughterLocation2 = this->GetCentroidOfElement(new_element_index);
    long_axis = this->GetShortAxisOfElement(new_element_index);
    division_info.mDaughterLongAxis2(0) = -long_axis(1);
    division_info.mDaughterLongAxis2(1) = long_axis(0);

    division_info.mDivisionAxis = axisOfDivision;
    mOperationRecorder.RecordCellDivisionInfo(division_info);

    if (mTrackMeshOperations)
    {
        // Record edge rearrangements in the daughter cells ...
        mOperationRecorder.RecordCellDivideOperation(edgeIds, pElement, this->mElements[new_element_index]);
        // ... and in each neighbouring cell
        for (unsigned int i=0; i<edge_split_pairs.size(); ++i)
        {
            mOperationRecorder.RecordEdgeSplitOperation(edge_split_pairs[i].first,
                                                        edge_split_pairs[i].second,
                                                        relative_new_node[i], true);
        }
    }
    return new_element_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElementAlongShortAxis(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement,
                                                                                bool placeOriginalElementBelow)
{
    assert(SPACE_DIM == 2);                // LCOV_EXCL_LINE
    assert(ELEMENT_DIM == SPACE_DIM);    // LCOV_EXCL_LINE

    c_vector<double, SPACE_DIM> short_axis = this->GetShortAxisOfElement(pElement->GetIndex());


    unsigned new_element_index = DivideElementAlongGivenAxis(pElement, short_axis, placeOriginalElementBelow);
    return new_element_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement,
                                                                  unsigned nodeAIndex,
                                                                  unsigned nodeBIndex,
                                                                  bool placeOriginalElementBelow)
{
    assert(SPACE_DIM == 2);                // LCOV_EXCL_LINE
    assert(ELEMENT_DIM == SPACE_DIM);    // LCOV_EXCL_LINE

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
        new_element_index = this->mElements.size();
    }
    else
    {
        new_element_index = mDeletedElementIndices.back();
        mDeletedElementIndices.pop_back();
        delete this->mElements[new_element_index];
    }

    // Add the new element to the mesh
    AddElement(new VertexElement<ELEMENT_DIM,SPACE_DIM>(new_element_index, nodes_elem));

    /**
     * Remove the correct nodes from each element. If placeOriginalElementBelow is true,
     * place the original element below (in the y direction) the new element; otherwise,
     * place it above.
     */

    // Find lowest element
    ///\todo this could be more efficient (see #2401)
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
        else if (i-1 > node1_index && i-1 < node2_index)
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
    // Re-build edges when division is performed
    this->mElements[new_element_index]->RebuildEdges();
    pElement->RebuildEdges();
    return new_element_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::DeleteElementPriorToReMesh(unsigned index)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE

    // Mark any nodes that are contained only in this element as deleted
    for (unsigned i=0; i<this->mElements[index]->GetNumNodes(); i++)
    {
        Node<SPACE_DIM>* p_node = this->mElements[index]->GetNode(i);

        if (p_node->rGetContainingElementIndices().size() == 1)
        {
            DeleteNodePriorToReMesh(p_node->GetIndex());
        }

        // Mark all the nodes contained in the removed element as boundary nodes
        p_node->SetAsBoundaryNode(true);
    }

    // Mark this element as deleted
    this->mElements[index]->MarkAsDeleted();
    mDeletedElementIndices.push_back(index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::DeleteNodePriorToReMesh(unsigned index)
{
    this->mNodes[index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideEdge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
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

    // Check that the nodes have a common edge and not more than 2
    assert(!shared_elements.empty());
    assert(shared_elements.size()<=2u);

    // Specify if it's a boundary node
    bool is_boundary_node = false;
    if (shared_elements.size()==1u)
    {
        // If only one shared element then must be on the boundary.
        assert((pNodeA->IsBoundaryNode()) && (pNodeB->IsBoundaryNode()));
        is_boundary_node = true;
    }

    // Create a new node (position is not important as it will be changed)
    Node<SPACE_DIM>* p_new_node = new Node<SPACE_DIM>(GetNumNodes(), is_boundary_node, 0.0, 0.0);

    // Update the node location
    c_vector<double, SPACE_DIM> new_node_position = pNodeA->rGetLocation() + 0.5*this->GetVectorFromAtoB(pNodeA->rGetLocation(), pNodeB->rGetLocation());
    ChastePoint<SPACE_DIM> point(new_node_position);
    p_new_node->SetPoint(new_node_position);

    // Add node to mesh
    this->mNodes.push_back(p_new_node);

    // Iterate over common elements
    unsigned node_A_index = pNodeA->GetIndex();
    unsigned node_B_index = pNodeB->GetIndex();
    for (std::set<unsigned>::iterator iter = shared_elements.begin();
         iter != shared_elements.end();
         ++iter)
    {
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(*iter);

        // Find which node has the lower local index in this element
        unsigned local_indexA = p_element->GetNodeLocalIndex(node_A_index);
        unsigned local_indexB = p_element->GetNodeLocalIndex(node_B_index);

        unsigned index = local_indexB;

        // If node B has a higher index then use node A's index...
        if (local_indexB > local_indexA)
        {
            index = local_indexA;

            // ...unless nodes A and B share the element's last edge
            if ((local_indexA == 0) && (local_indexB == p_element->GetNumNodes()-1))
            {
                index = local_indexB;
            }
        }
        else if ((local_indexB == 0) && (local_indexA == p_element->GetNumNodes()-1))
        {
            // ...otherwise use node B's index, unless nodes A and B share the element's last edge
            index = local_indexA;
        }

        // Add new node to this element
        this->GetElement(*iter)->AddNode(p_new_node, index);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::RemoveDeletedNodesAndElements(VertexElementMap& rElementMap)
{
    // Make sure the map is big enough.  Each entry will be set in the loop below.
    rElementMap.Resize(this->GetNumAllElements());
    // Remove any elements that have been marked for deletion and store all other elements in a temporary structure
    std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> live_elements;
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        if (this->mElements[i]->IsDeleted())
        {
            delete this->mElements[i];
            rElementMap.SetDeleted(i);
        }
        else
        {
            live_elements.push_back(this->mElements[i]);
            rElementMap.SetNewIndex(i, (unsigned)(live_elements.size()-1));
        }
    }

    // Sanity check
    assert(mDeletedElementIndices.size() == this->mElements.size() - live_elements.size());

    // Repopulate the elements vector and reset the list of deleted element indices
    mDeletedElementIndices.clear();
    this->mElements = live_elements;

    // Finally, reset the element indices to run from zero
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        this->mElements[i]->ResetIndex(i);
    }

    // Remove deleted edges and update the node-edge mapping
    this->mEdges.RemoveDeletedEdges();

    // Remove deleted nodes
    RemoveDeletedNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::RemoveDeletedNodes()
{
    // Remove any nodes that have been marked for deletion and store all other nodes in a temporary structure
    // Also mark edges associated with the deleted nodes
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

    // Sanity check
    assert(mDeletedNodeIndices.size() == this->mNodes.size() - live_nodes.size());
    // Repopulate the nodes vector and reset the list of deleted node indices
    this->mNodes = live_nodes;
    mDeletedNodeIndices.clear();

    // Finally, reset the node indices to run from zero
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->SetIndex(i);
    }

    // Remove deleted edges and update the node-edge mapping
    this->mEdges.RemoveDeletedEdges();
    // Update the node-edge mapping
    this->mEdges.UpdateEdgesMapKey();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh(VertexElementMap& rElementMap)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    assert(SPACE_DIM==2 || SPACE_DIM==3);     // LCOV_EXCL_LINE
    assert(ELEMENT_DIM == SPACE_DIM);         // LCOV_EXCL_LINE


    if (SPACE_DIM == 2)
    {
        // Make sure the map is big enough
        rElementMap.Resize(this->GetNumAllElements());

        /*
         * To begin the remeshing process, we do not need to call Clear() and remove all current data,
         * since cell birth, rearrangement and death result only in local remeshing of a vertex-based
         * mesh. Instead, we just remove any deleted elements and nodes.
         */
        RemoveDeletedNodesAndElements(rElementMap);
        bool recheck_mesh = true;
        while (recheck_mesh == true)
        {
            // We check for any short edges and perform swaps if necessary and possible.
            recheck_mesh = CheckForSwapsFromShortEdges();
        }

        // Check for element intersections
        recheck_mesh = true;
        while (recheck_mesh == true)
        {
            // Check mesh for intersections, and perform T3 swaps where required
            recheck_mesh = CheckForIntersections();
        }

        RemoveDeletedNodes();

        /*
         * This is handled in a separate method to allow child classes to implement additional ReMeshing functionality
         * (see #2664).
         */
        this->CheckForRosettes();


    }
    else // 3D
    {
// LCOV_EXCL_START
        EXCEPTION("Remeshing has not been implemented in 3D (see #827 and #860)\n");
// LCOV_EXCL_STOP
        ///\todo Implement ReMesh() in 3D (see #1422)
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh()
{
    VertexElementMap map(GetNumElements());
    ReMesh(map);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::CheckForSwapsFromShortEdges()
{
    // Loop over elements to check for T1 swaps
    for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_iter = this->GetElementIteratorBegin();
         elem_iter != this->GetElementIteratorEnd();
         ++elem_iter)
    {
        ///\todo Could we search more efficiently by just iterating over edges? (see #2401)

        unsigned num_nodes = elem_iter->GetNumNodes();
        assert(num_nodes > 0);

        // Loop over the nodes contained in this element
        for (unsigned local_index=0; local_index<num_nodes; local_index++)
        {
            // Find locations of the current node and anticlockwise node
            Node<SPACE_DIM>* p_current_node = elem_iter->GetNode(local_index);
            unsigned local_index_plus_one = (local_index+1)%num_nodes;    ///\todo Use iterators to tidy this up (see #2401)
            Node<SPACE_DIM>* p_anticlockwise_node = elem_iter->GetNode(local_index_plus_one);

            // Find distance between nodes
            double distance_between_nodes = this->GetDistanceBetweenNodes(p_current_node->GetIndex(), p_anticlockwise_node->GetIndex());

            // If the nodes are too close together...
            if (distance_between_nodes < mCellRearrangementThreshold)
            {
                // ...then check if any triangular elements are shared by these nodes...
                std::set<unsigned> elements_of_node_a = p_current_node->rGetContainingElementIndices();
                std::set<unsigned> elements_of_node_b = p_anticlockwise_node->rGetContainingElementIndices();

                std::set<unsigned> shared_elements;
                std::set_intersection(elements_of_node_a.begin(), elements_of_node_a.end(),
                               elements_of_node_b.begin(), elements_of_node_b.end(),
                               std::inserter(shared_elements, shared_elements.begin()));

                bool both_nodes_share_triangular_element = false;
                for (std::set<unsigned>::const_iterator it = shared_elements.begin();
                     it != shared_elements.end();
                     ++it)
                {
                    if (this->GetElement(*it)->GetNumNodes() <= 3)
                    {
                        both_nodes_share_triangular_element = true;
                        break;
                    }
                }

                // ...and if none are, then perform the required type of swap and halt the search, returning true
                if (!both_nodes_share_triangular_element)
                {
                    IdentifySwapType(p_current_node, p_anticlockwise_node);
                    return true;
                }
            }
        }
    }

    return false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::CheckForT2Swaps(VertexElementMap& rElementMap)
{
    // Loop over elements to check for T2 swaps
    for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_iter = this->GetElementIteratorBegin();
         elem_iter != this->GetElementIteratorEnd();
         ++elem_iter)
    {
        // If this element is triangular...
        if (elem_iter->GetNumNodes() == 3)
        {
            // ...and smaller than the threshold area...
            if (this->GetVolumeOfElement(elem_iter->GetIndex()) < GetT2Threshold())
            {
                // ...then perform a T2 swap and break out of the loop
                PerformT2Swap(*elem_iter);
                ///\todo: cover this line in a test
                rElementMap.SetDeleted(elem_iter->GetIndex());
                return true;
            }
        }
    }
    return false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::CheckForIntersections()
{
    // If checking for internal intersections as well as on the boundary, then check that no nodes have overlapped any elements...
    if (mCheckForInternalIntersections)
    {
        ///\todo Change to only loop over neighbouring elements (see #2401)
        for (typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator node_iter = this->GetNodeIteratorBegin();
             node_iter != this->GetNodeIteratorEnd();
             ++node_iter)
        {
            assert(!(node_iter->IsDeleted()));

            for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_iter = this->GetElementIteratorBegin();
                 elem_iter != this->GetElementIteratorEnd();
                 ++elem_iter)
            {
                unsigned elem_index = elem_iter->GetIndex();

                // Check that the node is not part of this element
                if (node_iter->rGetContainingElementIndices().count(elem_index) == 0)
                {
                    if (this->ElementIncludesPoint(node_iter->rGetLocation(), elem_index))
                    {
                        PerformIntersectionSwap(&(*node_iter), elem_index);
                        return true;
                    }
                }
            }
        }
    }
    else
    {
        // ...otherwise, just check that no boundary nodes have overlapped any boundary elements
        // First: find all boundary element and calculate their centroid only once
        std::vector<unsigned> boundary_element_indices;
        std::vector< c_vector<double, SPACE_DIM> > boundary_element_centroids;
        for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_iter = this->GetElementIteratorBegin();
                elem_iter != this->GetElementIteratorEnd();
                ++elem_iter)
        {
            if (elem_iter->IsElementOnBoundary())
            {
                unsigned element_index = elem_iter->GetIndex();
                boundary_element_indices.push_back(element_index);
                // should be a map but I am too lazy to look up the syntax
                boundary_element_centroids.push_back(this->GetCentroidOfElement(element_index));
            }
        }

        // Second: Check intersections only for those nodes and elements within
        // mDistanceForT3SwapChecking within each other (node<-->element centroid)
        for (typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator node_iter = this->GetNodeIteratorBegin();
                node_iter != this->GetNodeIteratorEnd();
                ++node_iter)
        {
            if (node_iter->IsBoundaryNode())
            {
                assert(!(node_iter->IsDeleted()));

                // index in boundary_element_centroids and boundary_element_indices
                unsigned boundary_element_index = 0;
                for (std::vector<unsigned>::iterator elem_iter = boundary_element_indices.begin();
                        elem_iter != boundary_element_indices.end();
                        ++elem_iter)
                {
                    // Check that the node is not part of this element
                    if (node_iter->rGetContainingElementIndices().count(*elem_iter) == 0)
                    {
                        c_vector<double, SPACE_DIM> node_location = node_iter->rGetLocation();
                        c_vector<double, SPACE_DIM> element_centroid = boundary_element_centroids[boundary_element_index];
                        double node_element_distance = norm_2(this->GetVectorFromAtoB(node_location, element_centroid));

                        if ( node_element_distance < mDistanceForT3SwapChecking )
                        {
                            if (this->ElementIncludesPoint(node_iter->rGetLocation(), *elem_iter))
                            {
                                this->PerformT3Swap(&(*node_iter), *elem_iter);
                                return true;
                            }
                        }
                    }
                    // increment the boundary element index
                    boundary_element_index +=1u;
                }
            }
        }

    }
    return false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::IdentifySwapType(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    // Find the sets of elements containing nodes A and B
    std::set<unsigned> nodeA_elem_indices = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> nodeB_elem_indices = pNodeB->rGetContainingElementIndices();

    // Form the set union
    std::set<unsigned> all_indices, temp_union_set;
    std::set_union(nodeA_elem_indices.begin(), nodeA_elem_indices.end(),
                   nodeB_elem_indices.begin(), nodeB_elem_indices.end(),
                   std::inserter(temp_union_set, temp_union_set.begin()));
    all_indices.swap(temp_union_set); // temp_set will be deleted, all_indices now contains all the indices of elements
                                      // that touch the potentially swapping nodes

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

        /*
         * This case is handled in a separate method to allow child classes to implement different
         * functionality for high-order-junction remodelling events (see #2664).
         */
        this->HandleHighOrderJunctions(pNodeA, pNodeB);
    }
    else // each node is contained in at most three elements
    {
        switch (all_indices.size())
        {
            case 1:
            {
                /*
                 * Each node is contained in a single element, so the nodes must lie on the boundary
                 * of the mesh, as shown below. In this case, we merge the nodes and tidy up node
                 * indices through calls to PerformNodeMerge() and RemoveDeletedNodes().
                 *
                 *    A   B
                 * ---o---o---
                 */
                assert(pNodeA->IsBoundaryNode());
                assert(pNodeB->IsBoundaryNode());
                PerformNodeMerge(pNodeA, pNodeB);
                RemoveDeletedNodes();
                break;
            }
            case 2:
            {
                if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)
                {
                    if (pNodeA->IsBoundaryNode() && pNodeB->IsBoundaryNode())
                    {
                        /*
                         * The node configuration is as shown below, with voids on either side. In this case
                         * we perform a T1 swap, which separates the elements.
                         *
                         *   \   /
                         *    \ / Node A
                         * (1) |   (2)      (element number in brackets)
                         *    / \ Node B
                         *   /   \
                         */
                         PerformT1Swap(pNodeA, pNodeB,all_indices);
                    }
                    else if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
                    {
                        /*
                         * The node configuration is as shown below, with a void on one side. We should not
                         * be able to reach this case at present, since we allow only for three-way junctions
                         * or boundaries, so we throw an exception.
                         *
                         *   \   /
                         *    \ / Node A
                         * (1) |   (2)      (element number in brackets)
                         *     x Node B
                         *     |
                         */
                        EXCEPTION("There is a non-boundary node contained only in two elements; something has gone wrong.");
                    }
                    else
                    {
                        /*
                         * Each node is contained in two elements, so the nodes lie on an internal edge, as shown below.
                         * We should not be able to reach this case at present, since we allow only for three-way junctions
                         * or boundaries, so we throw an exception.
                         *
                         *    A   B
                         * ---o---o---
                         */
                        EXCEPTION("There are non-boundary nodes contained only in two elements; something has gone wrong.");
                    }
                }// from [if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)]
                else
                {
                    /*
                     * The node configuration looks like that shown below. In this case, we merge the nodes
                     * and tidy up node indices through calls to PerformNodeMerge() and  RemoveDeletedNodes().
                     *
                     * Outside
                     *         /
                     *   --o--o (2)
                     *     (1) \
                     *
                     * ///\todo this should be a T1 swap (see #1263 and #2401)
                     * Referring to the todo: this should probably stay a node-merge. If this is a T1 swap then
                     * the single boundary node will travel from element 1 to element 2, but still remain a single node.
                     * I.e. we would not reduce the total number of nodes in this situation.
                     */
                    PerformNodeMerge(pNodeA, pNodeB);
                    RemoveDeletedNodes();
                }
                break;
            }
            case 3:
            {
                if (nodeA_elem_indices.size()==1 || nodeB_elem_indices.size()==1)
                {
                    /*
                     * One node is contained in one element and the other node is contained in three elements.
                     * We should not be able to reach this case at present, since we allow each boundary node
                     * to be contained in at most two elements, so we throw an exception.
                     *
                     *    A   B
                     *
                     *  empty   /
                     *         / (3)
                     * ---o---o-----   (element number in brackets)
                     *  (1)    \ (2)
                     *          \
                     */
                    assert(pNodeA->IsBoundaryNode());
                    assert(pNodeB->IsBoundaryNode());

                    EXCEPTION("There is a boundary node contained in three elements something has gone wrong.");
                }
                else if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)
                {
                    // The short edge must be at the boundary. We need to check whether this edge is
                    // adjacent to a triangular void before we swap. If it is a triangular void, we perform a T2-type swap.
                    // If not, then we perform a normal T1 swap. I.e. in detail we need to check whether the
                    // element in nodeA_elem_indices which is not in nodeB_elem_indices contains a shared node
                    // with the element in nodeB_elem_indices which is not in nodeA_elem_indices.

                    std::set<unsigned> element_A_not_B, temp_set;
                    std::set_difference(all_indices.begin(), all_indices.end(), nodeB_elem_indices.begin(),
                            nodeB_elem_indices.end(), std::inserter(temp_set, temp_set.begin()));
                    element_A_not_B.swap(temp_set);

                    // There must be only one such element
                    assert(element_A_not_B.size() == 1);

                    std::set<unsigned> element_B_not_A;
                    std::set_difference(all_indices.begin(), all_indices.end(), nodeA_elem_indices.begin(),
                            nodeA_elem_indices.end(), std::inserter(temp_set, temp_set.begin()));
                    element_B_not_A.swap(temp_set);

                    // There must be only one such element
                    assert(element_B_not_A.size() == 1);

                    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_A_not_B = this->mElements[*element_A_not_B.begin()];
                    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_B_not_A = this->mElements[*element_B_not_A.begin()];

                    unsigned local_index_1 = p_element_A_not_B->GetNodeLocalIndex(pNodeA->GetIndex());
                    unsigned next_node_1 = p_element_A_not_B->GetNodeGlobalIndex((local_index_1 + 1)%(p_element_A_not_B->GetNumNodes()));
                    unsigned previous_node_1 = p_element_A_not_B->GetNodeGlobalIndex(
                            (local_index_1 + p_element_A_not_B->GetNumNodes() - 1)%(p_element_A_not_B->GetNumNodes()));
                    unsigned local_index_2 = p_element_B_not_A->GetNodeLocalIndex(pNodeB->GetIndex());
                    unsigned next_node_2 = p_element_B_not_A->GetNodeGlobalIndex(
                            (local_index_2 + 1)%(p_element_B_not_A->GetNumNodes()));
                    unsigned previous_node_2 = p_element_B_not_A->GetNodeGlobalIndex(
                            (local_index_2 + p_element_B_not_A->GetNumNodes() - 1)%(p_element_B_not_A->GetNumNodes()));

                    if (next_node_1 == previous_node_2 || next_node_2 == previous_node_1)
                     {
                        /*
                         * The node configuration looks like that shown below, and both nodes must be on the boundary.
                         * In this case we remove the void through a call to PerformVoidRemoval().
                         *
                         *    A  C  B                A      B
                         *      /\                 \        /
                         *     /v \                 \  (1) /
                         * (3)o----o (1)  or     (2) o----o (3)    (element number in brackets, v is a void)
                         *   /  (2) \                 \v /
                         *  /        \                 \/
                         *                             C
                         */
                        assert(pNodeA->IsBoundaryNode());
                        assert(pNodeB->IsBoundaryNode());

                        // Get the third node in the triangular void

                        unsigned nodeC_index;
                        if (next_node_1 == previous_node_2 && next_node_2 != previous_node_1)
                        {
                            nodeC_index = next_node_1;
                        }
                        else if (next_node_2 == previous_node_1 && next_node_1 != previous_node_2)
                        {
                            nodeC_index = next_node_2;
                        }
                        else
                        {
                             assert(next_node_1 == previous_node_2 && next_node_2 == previous_node_1);
                             /**
                              * Here, the triangular element would be along the short edge. Since we
                              * are already checking in CheckForSwapsFromShortEdges() whether the element
                              * is triangular, this exception is redundant for simulations. We leave it in for
                              * clarity.
                              * ///\todo: consider removing the checking for this exception (see #2401)
                              */
                             EXCEPTION("Triangular element next to triangular void, not implemented yet.");
                        }

                        if (p_element_A_not_B->GetNumNodes() == 3u || p_element_B_not_A->GetNumNodes() == 3u)
                        {
                            /**
                             * If this is true then one of the elements adjacent to the triangular void
                             * is triangular. This element will then not share the short edge that is considered
                             * for a swap. Nevertheless, it would loose an edge during the swap. We are currently
                             * not able to deal with this situation.
                             * Related to #2533 and #2401.
                             */
                             EXCEPTION("Triangular element next to triangular void, not implemented yet.");
                        }
                        PerformVoidRemoval(pNodeA, pNodeB, this->mNodes[nodeC_index]);
                    }
                    else
                    {
                        /*
                         * The node configuration looks like that below, and both nodes must lie on the boundary.
                         * In this case we perform a T1 swap.
                         *
                         *     A  B                  A  B
                         *   \ empty/              \      /
                         *    \    /                \(1) /
                         * (3) o--o (1)  or      (2) o--o (3)    (element number in brackets)
                         *    / (2)\                /    \
                         *   /      \              /empty \
                         */
                        assert(pNodeA->IsBoundaryNode());
                        assert(pNodeB->IsBoundaryNode());

                        PerformT1Swap(pNodeA, pNodeB, all_indices);
                    }
                } // from else if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)
                else
                {
                    // In this case, one node must be contained in two elements and the other in three elements.
                    assert (   (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==3)
                            || (nodeA_elem_indices.size()==3 && nodeB_elem_indices.size()==2) );

                    // They can't both be boundary nodes
                    assert(!(pNodeA->IsBoundaryNode() && pNodeB->IsBoundaryNode()));

                    if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
                    {
                        /*
                         * The node configuration looks like that shown below. We perform a T1 swap in this case.
                         *
                         *     A  B                      A  B
                         *   \      /                  \      /
                         *    \ (1)/                    \(1) /
                         * (3) o--o (empty)  or  (empty) o--o (3)    (element number in brackets)
                         *    / (2)\                    /(2) \
                         *   /      \                  /      \
                         */
                        PerformT1Swap(pNodeA, pNodeB, all_indices);
                    }
                    else
                    {
                        /*
                         * The node configuration looks like that shown below. We should not be able to reach this case
                         * at present, since we allow only for three-way junctions or boundaries, so we throw an exception.
                         *
                         *     A  B             A  B
                         *   \                       /
                         *    \  (1)           (1)  /
                         * (3) o--o---   or  ---o--o (3)    (element number in brackets)
                         *    /  (2)           (2)  \
                         *   /                       \
                         */
                        EXCEPTION("There are non-boundary nodes contained only in two elements; something has gone wrong.");
                    }
                }
                break;
            }
            case 4:
            {
                /*
                 * The node configuration looks like that shown below. We perform a T1 swap in this case.
                 *
                 *   \(1)/
                 *    \ / Node A
                 * (2) |   (4)      (element number in brackets)
                 *    / \ Node B
                 *   /(3)\
                 */

                /*
                 * This case is handled in a separate method to allow child classes to implement different
                 * functionality for junction remodelling events (see #2664).
                 */
                if (mProtorosetteFormationProbability > RandomNumberGenerator::Instance()->ranf())
                {
                    this->PerformNodeMerge(pNodeA, pNodeB);
                    this->RemoveDeletedNodes();
                }
                else
                {
                    this->PerformT1Swap(pNodeA, pNodeB, all_indices);
                }
                break;
            }
            default:
                // This can't happen
                NEVER_REACHED;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformNodeMerge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    // Find the sets of elements containing each of the nodes, sorted by index
    std::set<unsigned> nodeA_elem_indices = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> nodeB_elem_indices = pNodeB->rGetContainingElementIndices();

    // Move node A to the mid-point
    pNodeA->rGetModifiableLocation() += 0.5 * this->GetVectorFromAtoB(pNodeA->rGetLocation(), pNodeB->rGetLocation());

    // Update the elements previously containing node B to contain node A
    unsigned node_B_index = pNodeB->GetIndex();
    // For rebuilding edges affected elements
    std::set<unsigned> rebuilt_elements;
    for (std::set<unsigned>::const_iterator it = nodeB_elem_indices.begin(); it != nodeB_elem_indices.end(); ++it)
    {
        // Find the local index of node B in this element
        unsigned node_B_local_index = this->mElements[*it]->GetNodeLocalIndex(node_B_index);
        assert(node_B_local_index < UINT_MAX); // this element contains node B

        /*
         * If this element already contains node A, then just remove node B.
         * Otherwise replace it with node A in the element and remove it from mNodes.
         */
        if (nodeA_elem_indices.count(*it) != 0)
        {
            std::vector<unsigned int> edgeIds;
            for (unsigned i = 0; i < this->mElements[*it]->GetNumEdges(); i++)
            {
                edgeIds.push_back(this->mElements[*it]->GetEdge(i)->GetIndex());
            }
            const unsigned int node_A_local_index = this->mElements[*it]->GetNodeLocalIndex(pNodeA->GetIndex());
            this->mElements[*it]->DeleteNode(node_B_local_index);

            if (mTrackMeshOperations)
            {
                mOperationRecorder.RecordNodeMergeOperation(edgeIds, this->mElements[*it],
                                                            std::pair<unsigned int, unsigned int>(node_A_local_index, node_B_local_index));
            }
        }
        else
        {
            // Replace node B with node A in this element
            this->mElements[*it]->UpdateNode(node_B_local_index, pNodeA);
        }
        rebuilt_elements.insert(this->mElements[*it]->GetIndex());
    }

    assert(!(this->mNodes[node_B_index]->IsDeleted()));
    this->mNodes[node_B_index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(node_B_index);
    for (unsigned i:rebuilt_elements)
        this->GetElement(i)->RebuildEdges();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformT1Swap(Node<SPACE_DIM>* pNodeA,
                                                              Node<SPACE_DIM>* pNodeB,
                                                              std::set<unsigned>& rElementsContainingNodes)
{
    // First compute and store the location of the T1 swap, which is at the midpoint of nodes A and B
    double distance_between_nodes_CD = mCellRearrangementRatio*mCellRearrangementThreshold;

    c_vector<double, SPACE_DIM> nodeA_location = pNodeA->rGetLocation();
    c_vector<double, SPACE_DIM> nodeB_location = pNodeB->rGetLocation();
    c_vector<double, SPACE_DIM> vector_AB = this->GetVectorFromAtoB(nodeA_location, nodeB_location);


    double distance_AB = norm_2(vector_AB);
    if (distance_AB < 1e-10) ///\todo remove magic number? (see #1884 and #2401)
    {
        EXCEPTION("Nodes are too close together, this shouldn't happen");
    }

    /*
     * Compute the locations of two new nodes C, D, placed on either side of the
     * edge E_old formed by nodes A and B, such that the edge E_new formed by the
     * new nodes is the perpendicular bisector of E_old, with |E_new| 'just larger'
     * (mCellRearrangementRatio) than mThresholdDistance.
     *
     * We implement the following changes to the mesh:
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
     */

    // Move nodes A and B to C and D respectively
    c_vector<double, SPACE_DIM> vector_CD;
    vector_CD(0) = -vector_AB(1) * distance_between_nodes_CD / distance_AB;
    vector_CD(1) =  vector_AB(0) * distance_between_nodes_CD / distance_AB;

    //Record T1Swap
    T1SwapInfo<SPACE_DIM> swap_info;
    swap_info.mLocation = nodeA_location+0.5*vector_AB;
    swap_info.mPreSwapEdge = vector_AB;
    swap_info.mPostSwapEdge = vector_CD;
    mOperationRecorder.RecordT1Swap(swap_info);

    c_vector<double, SPACE_DIM> nodeC_location = nodeA_location + 0.5*vector_AB - 0.5*vector_CD;
    c_vector<double, SPACE_DIM> nodeD_location = nodeC_location + vector_CD;

    pNodeA->rGetModifiableLocation() = nodeC_location;
    pNodeB->rGetModifiableLocation() = nodeD_location;

    // Find the sets of elements containing nodes A and B
    std::set<unsigned> nodeA_elem_indices = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> nodeB_elem_indices = pNodeB->rGetContainingElementIndices();

    //For rebuilding affected elements
    std::set<unsigned> rebuilt_elements;
    for (std::set<unsigned>::const_iterator it = rElementsContainingNodes.begin();
         it != rElementsContainingNodes.end();
         ++it)
    {
        std::vector<unsigned int> oldIds;
        for (unsigned int i=0; i<this->mElements[*it]->GetNumEdges(); ++i)
        {
            oldIds.push_back(this->mElements[*it]->GetEdge(i)->GetIndex());
        }
        // If, as in element 3 above, this element does not contain node A (now C)...
        if (nodeA_elem_indices.find(*it) == nodeA_elem_indices.end())
        {
            // ...then add it to the element just after node B (now D), going anticlockwise
            unsigned nodeB_local_index = this->mElements[*it]->GetNodeLocalIndex(pNodeB->GetIndex());
            assert(nodeB_local_index < UINT_MAX);

            std::vector<unsigned int> oldIds;
            for (unsigned int i=0; i<this->mElements[*it]->GetNumEdges(); ++i)
            {
                oldIds.push_back(this->mElements[*it]->GetEdge(i)->GetIndex());
            }

            this->mElements[*it]->AddNode(pNodeA, nodeB_local_index);
            if (mTrackMeshOperations)
                mOperationRecorder.RecordNewEdgeOperation(this->mElements[*it], nodeB_local_index);

        }
        else if (nodeB_elem_indices.find(*it) == nodeB_elem_indices.end())
        {
            // Do similarly if the element does not contain node B (now D), as in element 1 above
            unsigned nodeA_local_index = this->mElements[*it]->GetNodeLocalIndex(pNodeA->GetIndex());
            assert(nodeA_local_index < UINT_MAX);

            std::vector<unsigned int> oldIds;
            for (unsigned int i=0; i<this->mElements[*it]->GetNumEdges(); ++i)
            {
                oldIds.push_back(this->mElements[*it]->GetEdge(i)->GetIndex());
            }

            this->mElements[*it]->AddNode(pNodeB, nodeA_local_index);
            if (mTrackMeshOperations)
                mOperationRecorder.RecordNewEdgeOperation(this->mElements[*it], nodeA_local_index);
        }
        else
        {
            // If the element contains both nodes A and B (now C and D respectively)...
            unsigned nodeA_local_index = this->mElements[*it]->GetNodeLocalIndex(pNodeA->GetIndex());
            unsigned nodeB_local_index = this->mElements[*it]->GetNodeLocalIndex(pNodeB->GetIndex());

            assert(nodeA_local_index < UINT_MAX);
            assert(nodeB_local_index < UINT_MAX);

            /*
             * Locate local index of nodeA and nodeB and use the ordering to
             * identify the element, if nodeB_index > nodeA_index then element 4
             * and if nodeA_index > nodeB_index then element 2
             */
            unsigned nodeB_local_index_plus_one = (nodeB_local_index + 1)%(this->mElements[*it]->GetNumNodes());
            /*
             * T1 swap and subsequent A-B edge shrinkage is recorded as node merging
             */
            std::pair<unsigned int, unsigned int> deleted_node_indices;
            if (nodeA_local_index == nodeB_local_index_plus_one)
            {
                /*
                 * In this case the local index of nodeA is the local index of
                 * nodeB plus one so we are in element 2 so we remove nodeB
                 */
                this->mElements[*it]->DeleteNode(nodeB_local_index);
                deleted_node_indices.first = nodeA_local_index;
                deleted_node_indices.second = nodeB_local_index;
            }
            else
            {
                assert(nodeB_local_index == (nodeA_local_index + 1)%(this->mElements[*it]->GetNumNodes())); // as A and B are next to each other
                /*
                 * In this case the local index of nodeA is the local index of
                 * nodeB minus one so we are in element 4 so we remove nodeA
                 */
                this->mElements[*it]->DeleteNode(nodeA_local_index);
                deleted_node_indices.first = nodeB_local_index;
                deleted_node_indices.second = nodeA_local_index;
            }
            if (mTrackMeshOperations)
                mOperationRecorder.RecordNodeMergeOperation(oldIds, this->mElements[*it],
                                                            deleted_node_indices);
        }
        rebuilt_elements.insert(this->mElements[*it]->GetIndex());
    }

    // Sort out boundary nodes
    if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
    {
        if (pNodeA->GetNumContainingElements() == 3)
        {
            pNodeA->SetAsBoundaryNode(false);
        }
        else
        {
            pNodeA->SetAsBoundaryNode(true);
        }
        if (pNodeB->GetNumContainingElements() == 3)
        {
            pNodeB->SetAsBoundaryNode(false);
        }
        else
        {
            pNodeB->SetAsBoundaryNode(true);
        }
    }

    for (unsigned i:rebuilt_elements)
        this->GetElement(i)->RebuildEdges();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformIntersectionSwap(Node<SPACE_DIM>* pNode, unsigned elementIndex)
{
    assert(SPACE_DIM == 2);                    // LCOV_EXCL_LINE
    assert(ELEMENT_DIM == SPACE_DIM);        // LCOV_EXCL_LINE

    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(elementIndex);
    unsigned num_nodes = p_element->GetNumNodes();

    std::set<unsigned> elements_containing_intersecting_node;

    for (unsigned node_local_index=0; node_local_index<num_nodes; node_local_index++)
    {
        unsigned node_global_index = p_element->GetNodeGlobalIndex(node_local_index);

        std::set<unsigned> node_elem_indices = this->GetNode(node_global_index)->rGetContainingElementIndices();

        for (std::set<unsigned>::const_iterator elem_iter = node_elem_indices.begin();
             elem_iter != node_elem_indices.end();
             ++elem_iter)
        {
            VertexElement<ELEMENT_DIM, SPACE_DIM>* p_neighbouring_element = this->GetElement(*elem_iter);
            unsigned num_nodes_in_neighbouring_element = p_neighbouring_element->GetNumNodes();

            // Check if element contains the intersecting node
            for (unsigned node_index_2 = 0; node_index_2 < num_nodes_in_neighbouring_element; node_index_2++)
            {
                if (p_neighbouring_element->GetNodeGlobalIndex(node_index_2) == pNode->GetIndex())
                {
                    elements_containing_intersecting_node.insert(p_neighbouring_element->GetIndex());
                }
            }
        }
    }
    /*
     * If there are not two elements containing the intersecting node then the node is coming from the other side of the element
     * and there is no way to fix it unless you want to make two new elements.
     */
    assert(elements_containing_intersecting_node.size() == 2);

    std::set<unsigned> all_elements_containing_intersecting_node = pNode->rGetContainingElementIndices();

    std::set<unsigned> intersecting_element;

    std::set_difference(all_elements_containing_intersecting_node.begin(), all_elements_containing_intersecting_node.end(),
                        elements_containing_intersecting_node.begin(), elements_containing_intersecting_node.end(),
                        std::inserter(intersecting_element, intersecting_element.begin()));

    /*
     * Identify nodes and elements to perform switch on
     * Intersecting node is node A
     * Other node is node B
     *
     * Element 1 only contains node A
     * Element 2 has nodes B and A (in that order)
     * Element 3 only contains node B
     * Element 4 has nodes A and B (in that order)
     */
    unsigned node_A_index = pNode->GetIndex();
    unsigned node_B_index;
    unsigned element_1_index = *(intersecting_element.begin());
    unsigned element_2_index;
    unsigned element_3_index = elementIndex;
    unsigned element_4_index;

    std::set<unsigned>::iterator iter = elements_containing_intersecting_node.begin();
    unsigned element_a_index = *(iter);
    iter++;
    unsigned element_b_index = *(iter);

    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_a = this->GetElement(element_a_index);
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_b = this->GetElement(element_b_index);

    std::set<unsigned> element_a_nodes;
    for (unsigned node_index = 0;  node_index < p_element_a->GetNumNodes(); node_index++)
    {
        element_a_nodes.insert(p_element_a->GetNodeGlobalIndex(node_index));
    }

    std::set<unsigned> element_b_nodes;
    for (unsigned node_index = 0;  node_index < p_element_b->GetNumNodes(); node_index++)
    {
        element_b_nodes.insert(p_element_b->GetNodeGlobalIndex(node_index));
    }

    std::set<unsigned> switching_nodes;
    std::set_intersection(element_a_nodes.begin(), element_a_nodes.end(),
                          element_b_nodes.begin(), element_b_nodes.end(),
                          std::inserter(switching_nodes, switching_nodes.begin()));

    assert(switching_nodes.size() == 2);

    // Check intersecting node is this set
    assert(switching_nodes.find(node_A_index) != switching_nodes.end());
    switching_nodes.erase(node_A_index);

    assert(switching_nodes.size() == 1);

    node_B_index = *(switching_nodes.begin());

    // Now identify elements 2 and 4
    unsigned node_A_local_index_in_a = p_element_a->GetNodeLocalIndex(node_A_index);
    unsigned node_B_local_index_in_a = p_element_a->GetNodeLocalIndex(node_B_index);

    if ((node_B_local_index_in_a+1)%p_element_a->GetNumNodes() == node_A_local_index_in_a)
    {
        assert((p_element_b->GetNodeLocalIndex(node_A_index)+1)%p_element_b->GetNumNodes()
               == p_element_b->GetNodeLocalIndex(node_B_index));

        // Element 2 is element a, element 4 is element b
        element_2_index = element_a_index;
        element_4_index = element_b_index;
    }
    else
    {
        assert((p_element_b->GetNodeLocalIndex(node_B_index)+1)%p_element_b->GetNumNodes()
               == p_element_b->GetNodeLocalIndex(node_A_index));

        // Element 2 is element b, element 4 is element a
        element_2_index = element_b_index;
        element_4_index = element_a_index;
    }

    unsigned intersected_edge = this->GetLocalIndexForElementEdgeClosestToPoint(pNode->rGetLocation(), elementIndex);

    unsigned node_A_local_index_in_1 = this->GetElement(element_1_index)->GetNodeLocalIndex(node_A_index);

    unsigned node_A_local_index_in_2 = this->GetElement(element_2_index)->GetNodeLocalIndex(node_A_index);
    unsigned node_B_local_index_in_2 = this->GetElement(element_2_index)->GetNodeLocalIndex(node_B_index);

    unsigned node_B_local_index_in_3 = this->GetElement(elementIndex)->GetNodeLocalIndex(node_B_index);

    unsigned node_A_local_index_in_4 = this->GetElement(element_4_index)->GetNodeLocalIndex(node_A_index);
    unsigned node_B_local_index_in_4 = this->GetElement(element_4_index)->GetNodeLocalIndex(node_B_index);

    if (intersected_edge==node_B_local_index_in_3)
    {
        /*
         * Add node B to element 1 after node A
         * Add node A to element 3 after node B
         *
         * Remove node B from element 2
         * Remove node A from element 4
         */
        this->mElements[element_1_index]->AddNode(this->mNodes[node_B_index], node_A_local_index_in_1);
        this->mElements[element_3_index]->AddNode(this->mNodes[node_A_index], node_B_local_index_in_3);

        this->mElements[element_2_index]->DeleteNode(node_B_local_index_in_2);
        this->mElements[element_4_index]->DeleteNode(node_A_local_index_in_4);
    }
    else
    {
        assert((intersected_edge+1)%num_nodes==node_B_local_index_in_3);

        // Add node B to element 1 before node A and add node A to element 3 before node B
        unsigned node_before_A_in_1 = (node_A_local_index_in_1 - 1)%this->GetElement(element_1_index)->GetNumNodes();
        unsigned node_before_B_in_3 = (node_B_local_index_in_3 - 1)%this->GetElement(element_3_index)->GetNumNodes();
        this->mElements[element_1_index]->AddNode(this->mNodes[node_B_index], node_before_A_in_1);
        this->mElements[element_3_index]->AddNode(this->mNodes[node_A_index], node_before_B_in_3);

        // Remove node A from element 2 and remove node B from element 4
        this->mElements[element_2_index]->DeleteNode(node_A_local_index_in_2);
        this->mElements[element_4_index]->DeleteNode(node_B_local_index_in_4);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformT2Swap(VertexElement<ELEMENT_DIM,SPACE_DIM>& rElement)
{
    // The given element must be triangular for us to be able to perform a T2 swap on it
    assert(rElement.GetNumNodes() == 3);
    // Note that we define this vector before setting it, as otherwise the profiling build will break (see #2367)
    c_vector<double, SPACE_DIM> new_node_location;
    new_node_location = this->GetCentroidOfElement(rElement.GetIndex());
    mLastT2SwapLocation = new_node_location;

    // Create a new node at the element's centroid; this will be a boundary node if any existing nodes were on the boundary
    bool is_node_on_boundary = false;
    for (unsigned i=0; i<3; i++)
    {
        if (rElement.GetNode(i)->IsBoundaryNode())
        {
            is_node_on_boundary = true;
            break;
        }
    }
    unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(GetNumNodes(), new_node_location, is_node_on_boundary));
    Node<SPACE_DIM>* p_new_node = this->GetNode(new_node_global_index);

    std::set<unsigned int> neigh_indices;

    // Loop over each of the three nodes contained in rElement
    for (unsigned i=0; i<3; i++)
    {
        // For each node, find the set of other elements containing it
        Node<SPACE_DIM>* p_node = rElement.GetNode(i);

        std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
        containing_elements.erase(rElement.GetIndex());
        // For each of these elements...
        for (std::set<unsigned>::iterator elem_iter = containing_elements.begin(); elem_iter != containing_elements.end(); ++elem_iter)
        {
            VertexElement<ELEMENT_DIM,SPACE_DIM>* p_this_elem = this->GetElement(*elem_iter);

            neigh_indices.insert(p_this_elem->GetIndex());

            // ...throw an exception if the element is triangular...
            if (p_this_elem->GetNumNodes() < 4)
            {
                EXCEPTION("One of the neighbours of a small triangular element is also a triangle - dealing with this has not been implemented yet");
            }

            // ...otherwise, replace p_node with p_new_node unless this has already happened (in which case, delete p_node from the element)
            if (p_this_elem->GetNodeLocalIndex(new_node_global_index) == UINT_MAX)
            {
                p_this_elem->ReplaceNode(p_node, p_new_node);
            }
            else
            {
                std::vector<unsigned int> oldIds;
                std::pair<unsigned int, unsigned int> node_pair;
                if (mTrackMeshOperations)
                {
                    for (unsigned int k=0; k<p_this_elem->GetNumEdges(); ++k)
                    {
                        oldIds.push_back(p_this_elem->GetEdge(k)->GetIndex());
                    }
                    node_pair.first= p_this_elem->GetNodeLocalIndex(new_node_global_index);
                    node_pair.second = p_this_elem->GetNodeLocalIndex(p_node->GetIndex());
                }
                p_this_elem->DeleteNode(p_this_elem->GetNodeLocalIndex(p_node->GetIndex()));
                if (mTrackMeshOperations)
                    mOperationRecorder.RecordNodeMergeOperation(oldIds, p_this_elem, node_pair, true);
            }
        }
    }

    // We also have to mark pElement, pElement->GetNode(0), pElement->GetNode(1), and pElement->GetNode(2) as deleted
    mDeletedNodeIndices.push_back(rElement.GetNodeGlobalIndex(0));
    mDeletedNodeIndices.push_back(rElement.GetNodeGlobalIndex(1));
    mDeletedNodeIndices.push_back(rElement.GetNodeGlobalIndex(2));

    rElement.GetNode(0)->MarkAsDeleted();
    rElement.GetNode(1)->MarkAsDeleted();
    rElement.GetNode(2)->MarkAsDeleted();

    mDeletedElementIndices.push_back(rElement.GetIndex());
    rElement.MarkAsDeleted();

    for (unsigned i:neigh_indices)
    {
        this->GetElement(i)->RebuildEdges();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformT3Swap(Node<SPACE_DIM>* pNode, unsigned elementIndex)
{
    assert(SPACE_DIM == 2);                 // LCOV_EXCL_LINE - code will be removed at compile time
    assert(ELEMENT_DIM == SPACE_DIM);    // LCOV_EXCL_LINE - code will be removed at compile time

    assert(pNode->IsBoundaryNode());

    // Store the index of the elements containing the intersecting node
    std::set<unsigned> elements_containing_intersecting_node = pNode->rGetContainingElementIndices();

    // Get the local index of the node in the intersected element after which the new node is to be added
    unsigned node_A_local_index = this->GetLocalIndexForElementEdgeClosestToPoint(pNode->rGetLocation(), elementIndex);

    // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
    c_vector<double, SPACE_DIM> node_location;
    node_location = pNode->rGetModifiableLocation();

    // Get element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(elementIndex);
    unsigned num_nodes = p_element->GetNumNodes();

    // Get the nodes at either end of the edge to be divided
    unsigned vertexA_index = p_element->GetNodeGlobalIndex(node_A_local_index);
    unsigned vertexB_index = p_element->GetNodeGlobalIndex((node_A_local_index+1)%num_nodes);

    // Check these nodes are also boundary nodes if this fails then the elements have become concave and you need a smaller timestep
    if (!this->mNodes[vertexA_index]->IsBoundaryNode() || !this->mNodes[vertexB_index]->IsBoundaryNode())
    {
        EXCEPTION("A boundary node has intersected a non-boundary edge; this is because the boundary element has become concave. You need to rerun the simulation with a smaller time step to prevent this.");
    }

    // Get the nodes at either end of the edge to be divided and calculate intersection
    c_vector<double, SPACE_DIM> vertexA = p_element->GetNodeLocation(node_A_local_index);
    c_vector<double, SPACE_DIM> vertexB = p_element->GetNodeLocation((node_A_local_index+1)%num_nodes);
    c_vector<double, SPACE_DIM> vector_a_to_point = this->GetVectorFromAtoB(vertexA, node_location);

    c_vector<double, SPACE_DIM> vector_a_to_b = this->GetVectorFromAtoB(vertexA, vertexB);

    c_vector<double, SPACE_DIM> edge_ab_unit_vector = vector_a_to_b/norm_2(vector_a_to_b);
    c_vector<double, SPACE_DIM> intersection = vertexA + edge_ab_unit_vector*inner_prod(vector_a_to_point, edge_ab_unit_vector);

    // Store the location of the T3 swap, the location of the intersection with the edge
    ///\todo the intersection location is sometimes overwritten when WidenEdgeOrCorrectIntersectionLocationIfNecessary
    // is called (see #2401) - we should correct this in these cases!
    T3SwapInfo<SPACE_DIM> swap_info;
    swap_info.mLocation = intersection;
    mOperationRecorder.RecordT3Swap(swap_info);

    // For rebuilding edges
    std::set<unsigned> rebuilt_elements;
    if (pNode->GetNumContainingElements() == 1)
    {
        // Get the index of the element containing the intersecting node
        unsigned intersecting_element_index = *elements_containing_intersecting_node.begin();

        // Get element
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_intersecting_element = this->GetElement(intersecting_element_index);

        //Edge map before swap
        const unsigned int n_edges = p_intersecting_element->GetNumEdges();
        std::vector<unsigned int> oldIds(n_edges);
        for (unsigned int i=0; i<n_edges; ++i)
        {
            oldIds[i] = p_intersecting_element->GetEdge(i)->GetIndex();
        }

        unsigned local_index = p_intersecting_element->GetNodeLocalIndex(pNode->GetIndex());
        unsigned next_node = p_intersecting_element->GetNodeGlobalIndex((local_index + 1)%(p_intersecting_element->GetNumNodes()));
        unsigned previous_node = p_intersecting_element->GetNodeGlobalIndex((local_index + p_intersecting_element->GetNumNodes() - 1)%(p_intersecting_element->GetNumNodes()));

        // Check to see if the nodes adjacent to the intersecting node are contained in the intersected element between vertices A and B
        if (next_node == vertexA_index || previous_node == vertexA_index || next_node == vertexB_index || previous_node == vertexB_index)
        {
            unsigned common_vertex_index;

            if (next_node == vertexA_index || previous_node == vertexA_index)
            {
                common_vertex_index = vertexA_index;
            }
            else
            {
                common_vertex_index = vertexB_index;
            }

            assert(this->mNodes[common_vertex_index]->GetNumContainingElements()>1);

            std::set<unsigned> elements_containing_common_vertex = this->mNodes[common_vertex_index]->rGetContainingElementIndices();
            std::set<unsigned>::const_iterator it = elements_containing_common_vertex.begin();
            VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_common_1 = this->GetElement(*it);
            it++;
            VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_common_2 = this->GetElement(*it);

            // Find the number and indices of common vertices between element_1 and element_2
            unsigned num_common_vertices = 0;
            std::vector<unsigned> common_vertex_indices;
            for (unsigned i=0; i<p_element_common_1->GetNumNodes(); i++)
            {
                for (unsigned j=0; j<p_element_common_2->GetNumNodes(); j++)
                {
                    if (p_element_common_1->GetNodeGlobalIndex(i)==p_element_common_2->GetNodeGlobalIndex(j))
                    {
                        num_common_vertices++;
                        common_vertex_indices.push_back(p_element_common_1->GetNodeGlobalIndex(i));
                    }
                }
            }

            if (num_common_vertices == 1 || this->mNodes[common_vertex_index]->GetNumContainingElements() > 2)
            {
                /*
                 * This is the situation here.
                 *
                 *  From          To
                 *   _             _
                 *    |  <---       |
                 *    |  /\         |\
                 *    | /  \        | \
                 *   _|/____\      _|__\
                 *
                 * The edge goes from vertexA--vertexB to vertexA--pNode--vertexB
                 */

                // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
                intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(vertexA_index, vertexB_index, intersection);

                // Move original node
                pNode->rGetModifiableLocation() = intersection;

                //Record the edge split data
                const double a_to_b_length = norm_2(vector_a_to_b);
                c_vector<double, SPACE_DIM> vector_a_to_node = this->GetVectorFromAtoB(vertexA, intersection);
                const double a_to_node_length = norm_2(vector_a_to_node);

                // Add the moved nodes to the element (this also updates the node)
                this->GetElement(elementIndex)->AddNode(pNode, node_A_local_index);

                if (mTrackMeshOperations)
                    mOperationRecorder.RecordEdgeSplitOperation(p_element, node_A_local_index, a_to_node_length/a_to_b_length);
                rebuilt_elements.insert(p_element->GetIndex());
                // Check the nodes are updated correctly
                assert(pNode->GetNumContainingElements() == 2);
            }
            else if (num_common_vertices == 2)
            {
                // The two elements must have an edge in common.  Find whether the common edge is the same as the
                // edge that is merged onto.

                if ((common_vertex_indices[0]==vertexA_index && common_vertex_indices[1]==vertexB_index) ||
                    (common_vertex_indices[1]==vertexA_index && common_vertex_indices[0]==vertexB_index))
                {
                    /*
                     * Due to a previous T3 swap the situation looks like this.
                     *
                     *              pNode
                     *     \         |\    /
                     *      \        | \  /
                     *       \_______|__\/
                     *       /A      |     B
                     *      /         \
                     *
                     * A T3 Swap would merge pNode onto an edge of its own element.
                     * We prevent this by just removing pNode. By doing this we also avoid the
                     * intersecting element to be concave.
                     */

                    //Obtain necessary data for event recording
                    const unsigned int downstream_index = (local_index + p_intersecting_element->GetNumNodes() - 1)%(p_intersecting_element->GetNumNodes());

                    // Delete pNode in the intersecting element
                    p_intersecting_element->DeleteNode(local_index);

                    //We record node deletion as node merging.
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordNodeMergeOperation(oldIds, p_intersecting_element, std::pair<unsigned int, unsigned int>(downstream_index, local_index));
                    rebuilt_elements.insert(p_intersecting_element->GetIndex());

                    // Mark all three nodes as deleted
                    pNode->MarkAsDeleted();
                    mDeletedNodeIndices.push_back(pNode->GetIndex());
                }
                else
                {
                    /*
                     * This is the situation here.
                     *
                     * C is common_vertex D is the other one.
                     *
                     *  From          To
                     *   _ D          _
                     *    | <---       |
                     *    | /\         |\
                     *   C|/  \        | \
                     *   _|____\      _|__\
                     *
                     *  The edge goes from vertexC--vertexB to vertexC--pNode--vertexD
                     *  then vertex B is removed as it is no longer needed.
                     */

                    // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
                    intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(vertexA_index, vertexB_index, intersection);

                    // Move original node
                    pNode->rGetModifiableLocation() = intersection;

                    // Replace common_vertex with the the moved node (this also updates the nodes)
                    this->GetElement(elementIndex)->ReplaceNode(this->mNodes[common_vertex_index], pNode);

                    // Remove common_vertex
                    unsigned common_vertex_local_index = this->GetElement(intersecting_element_index)->GetNodeLocalIndex(common_vertex_index);
                    this->GetElement(intersecting_element_index)->DeleteNode(common_vertex_local_index);
                    assert(this->mNodes[common_vertex_index]->GetNumContainingElements() == 0);

                    //Record edge merging in the intersecting element
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeMergeOperation(p_intersecting_element, common_vertex_local_index);
                    rebuilt_elements.insert(p_intersecting_element->GetIndex());

                    this->mNodes[common_vertex_index]->MarkAsDeleted();
                    mDeletedNodeIndices.push_back(common_vertex_index);

                    // Check the nodes are updated correctly
                    assert(pNode->GetNumContainingElements() == 2);
                }
            }
            else if (num_common_vertices == 4)
            {
                /*
                 * The two elements share edges CA and BD due to previous swaps but not the edge AB
                 *
                 *  From          To
                 *  D___         D___
                 *    |            |
                 *   B|\           |
                 *    | \          |
                 *    | /          |
                 *   A|/           |
                 *  C_|__        C_|__
                 *
                 *  We just remove the intersecting node as well as vertices A and B.
                 */

                // Delete node A and B in the intersected element
                this->GetElement(elementIndex)->DeleteNode(node_A_local_index);
                if (mTrackMeshOperations)
                    mOperationRecorder.RecordEdgeMergeOperation(this->GetElement(elementIndex), node_A_local_index);
                unsigned node_B_local_index = this->
                        GetElement(elementIndex)->GetNodeLocalIndex(vertexB_index);
                this->GetElement(elementIndex)->DeleteNode(node_B_local_index);
                if (mTrackMeshOperations)
                    mOperationRecorder.RecordEdgeMergeOperation(this->GetElement(elementIndex), node_B_local_index);

                // Delete nodes A and B in the intersecting element
                unsigned node_A_local_index_intersecting_element = this->
                        GetElement(intersecting_element_index)->GetNodeLocalIndex(vertexA_index);
                this->GetElement(intersecting_element_index)->DeleteNode(node_A_local_index_intersecting_element);
                if (mTrackMeshOperations)
                    mOperationRecorder.RecordEdgeMergeOperation(this->GetElement(intersecting_element_index), node_A_local_index_intersecting_element);
                unsigned node_B_local_index_intersecting_element = this->
                        GetElement(intersecting_element_index)->GetNodeLocalIndex(vertexB_index);
                this->GetElement(intersecting_element_index)->DeleteNode(node_B_local_index_intersecting_element);
                if (mTrackMeshOperations)
                    mOperationRecorder.RecordEdgeMergeOperation(this->GetElement(intersecting_element_index), node_B_local_index_intersecting_element);

                // Delete pNode in the intersecting element
                unsigned p_node_local_index = this->
                        GetElement(intersecting_element_index)->GetNodeLocalIndex(pNode->GetIndex());
                this->GetElement(intersecting_element_index)->DeleteNode(p_node_local_index);
                if (mTrackMeshOperations)
                    mOperationRecorder.RecordEdgeMergeOperation(this->GetElement(intersecting_element_index), p_node_local_index);

                rebuilt_elements.insert(elementIndex);
                rebuilt_elements.insert(intersecting_element_index);
                // Mark all three nodes as deleted
                pNode->MarkAsDeleted();
                mDeletedNodeIndices.push_back(pNode->GetIndex());
                this->mNodes[vertexA_index]->MarkAsDeleted();
                mDeletedNodeIndices.push_back(vertexA_index);
                this->mNodes[vertexB_index]->MarkAsDeleted();
                mDeletedNodeIndices.push_back(vertexB_index);
            }
            else
            {
                // This can't happen as nodes can't be on the internal edge of 2 elements.
                NEVER_REACHED;
            }
        }
        else
        {
            /*
             *  From          To
             *   ____        _______
             *                 / \
             *    /\   ^      /   \
             *   /  \  |
             *
             *  The edge goes from vertexA--vertexB to vertexA--new_node--pNode--vertexB
             */

            // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
            intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(vertexA_index, vertexB_index, intersection);
            edge_ab_unit_vector = this->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

            // Move original node
            pNode->rGetModifiableLocation() = intersection + 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

            // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
            c_vector<double, SPACE_DIM> new_node_location;
            new_node_location = intersection - 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

            // Add new node which will always be a boundary node
            unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(0, true, new_node_location[0], new_node_location[1]));

            // Add the moved node (this also updates the node)
            //and record the split of the edge vertexA--vertexB due to insertion of pNode
            this->GetElement(elementIndex)->AddNode(pNode, node_A_local_index);
            //Recording edge split
            c_vector<double, SPACE_DIM> vector_a_to_node = this->GetVectorFromAtoB(pNode->rGetLocation(), vertexA);
            if (mTrackMeshOperations)
                mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex), node_A_local_index,
                                                            norm_2(vector_a_to_node)/norm_2(vector_a_to_b));

            //Add the new node to the element (this also updates the node)
            //and record the split of the edge vertexA--pNode by insertion of the new node
            this->GetElement(elementIndex)->AddNode(this->mNodes[new_node_global_index], node_A_local_index);

            //Recording edge split
            c_vector<double, SPACE_DIM> vector_a_to_new_node = this->GetVectorFromAtoB(new_node_location, vertexA);
            if (mTrackMeshOperations)
                mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex), node_A_local_index,
                                                            norm_2(vector_a_to_new_node)/norm_2(vector_a_to_node));
            //New node must be between vertexA and pNode
            assert(norm_2(vector_a_to_new_node)/norm_2(vector_a_to_node)<1);

            // Add the new node to the original element containing pNode (this also updates the node)
            const unsigned int new_pNode_local_index = this->GetElement(intersecting_element_index)->GetNodeLocalIndex(pNode->GetIndex());
            this->GetElement(intersecting_element_index)->AddNode(this->mNodes[new_node_global_index],
                                                                  new_pNode_local_index);
            if (mTrackMeshOperations)
                mOperationRecorder.RecordNewEdgeOperation(this->GetElement(intersecting_element_index), new_pNode_local_index);
            // The nodes must have been updated correctly
            assert(pNode->GetNumContainingElements() == 2);
            assert(this->mNodes[new_node_global_index]->GetNumContainingElements() == 2);

            rebuilt_elements.insert(elementIndex);
            rebuilt_elements.insert(intersecting_element_index);
        }
    }
    else if (pNode->GetNumContainingElements() == 2)
    {
        // Find the nodes contained in elements containing the intersecting node
        std::set<unsigned>::const_iterator it = elements_containing_intersecting_node.begin();

        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_1 = this->GetElement(*it);
        unsigned num_nodes_elem_1 = p_element_1->GetNumNodes();
        it++;

        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_2 = this->GetElement(*it);
        unsigned num_nodes_elem_2 = p_element_2->GetNumNodes();

        unsigned node_global_index = pNode->GetIndex();

        unsigned local_index_1 = p_element_1->GetNodeLocalIndex(node_global_index);
        unsigned next_node_1 = p_element_1->GetNodeGlobalIndex((local_index_1 + 1)%num_nodes_elem_1);
        unsigned previous_node_1 = p_element_1->GetNodeGlobalIndex((local_index_1 + num_nodes_elem_1 - 1)%num_nodes_elem_1);

        unsigned local_index_2 = p_element_2->GetNodeLocalIndex(node_global_index);
        unsigned next_node_2 = p_element_2->GetNodeGlobalIndex((local_index_2 + 1)%num_nodes_elem_2);
        unsigned previous_node_2 = p_element_2->GetNodeGlobalIndex((local_index_2 + num_nodes_elem_2 - 1)%num_nodes_elem_2);

        // Check to see if the nodes adjacent to the intersecting node are contained in the intersected element between vertices A and B
        if ((next_node_1 == vertexA_index || previous_node_1 == vertexA_index || next_node_2 == vertexA_index || previous_node_2 == vertexA_index) &&
                (next_node_1 == vertexB_index || previous_node_1 == vertexB_index || next_node_2 == vertexB_index || previous_node_2 == vertexB_index))
        {
            /*
             * Here we have
             *        __
             *      /|             /
             *  __ / |     --> ___/
             *     \ |            \
             *      \|__           \
             *
             * Where the node on the left has overlapped the edge A B
             *
             * Move p_node to the intersection on A B and merge AB and p_node
             */

            // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
            intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(vertexA_index, vertexB_index, intersection);
            edge_ab_unit_vector = this->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

            // Check they are all boundary nodes
            assert(pNode->IsBoundaryNode());
            assert(this->mNodes[vertexA_index]->IsBoundaryNode());
            assert(this->mNodes[vertexB_index]->IsBoundaryNode());

            // Move p_node to the intersection with the edge AB
            pNode->rGetModifiableLocation() = intersection;
            pNode->SetAsBoundaryNode(false);

            // Add pNode to the intersected element
            this->GetElement(elementIndex)->AddNode(pNode, node_A_local_index);

            //Record this as edge split
            const double a_to_b_length = norm_2(vector_a_to_b);
            c_vector<double, SPACE_DIM> vector_a_to_p = this->GetVectorFromAtoB(intersection, vertexA);
            const double split_ratio = norm_2(vector_a_to_p)/a_to_b_length;
            if (mTrackMeshOperations)
                mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex), node_A_local_index, split_ratio);
            rebuilt_elements.insert(elementIndex);

            // Remove vertex A from elements and record this as edge merge operation
            std::set<unsigned> elements_containing_vertex_A = this->mNodes[vertexA_index]->rGetContainingElementIndices();
            for (std::set<unsigned>::const_iterator iter = elements_containing_vertex_A.begin();
                    iter != elements_containing_vertex_A.end();
                    iter++)
            {
                const unsigned int this_vertexA_local_index = this->GetElement(*iter)->GetNodeLocalIndex(vertexA_index);
                this->GetElement(*iter)->DeleteNode(this_vertexA_local_index);
                if (mTrackMeshOperations)
                    mOperationRecorder.RecordEdgeMergeOperation(this->GetElement(*iter), this_vertexA_local_index);
                rebuilt_elements.insert(this->GetElement(*iter)->GetIndex());
            }

            // Remove vertex A from the mesh
            assert(this->mNodes[vertexA_index]->GetNumContainingElements() == 0);
            this->mNodes[vertexA_index]->MarkAsDeleted();
            mDeletedNodeIndices.push_back(vertexA_index);

            // Remove vertex B from elements and record this as edge merge operation
            std::set<unsigned> elements_containing_vertex_B = this->mNodes[vertexB_index]->rGetContainingElementIndices();
            for (std::set<unsigned>::const_iterator iter = elements_containing_vertex_B.begin();
                 iter != elements_containing_vertex_B.end();
                 iter++)
            {
                const unsigned int this_vertexB_local_index = this->GetElement(*iter)->GetNodeLocalIndex(vertexB_index);
                this->GetElement(*iter)->DeleteNode(this_vertexB_local_index);
                if (mTrackMeshOperations)
                    mOperationRecorder.RecordEdgeMergeOperation(this->GetElement(*iter), this_vertexB_local_index);
                rebuilt_elements.insert(this->GetElement(*iter)->GetIndex());
            }

            // Remove vertex B from the mesh
            assert(this->mNodes[vertexB_index]->GetNumContainingElements()==0);
            this->mNodes[vertexB_index]->MarkAsDeleted();
            mDeletedNodeIndices.push_back(vertexB_index);
        }
        else
        {
            if (next_node_1 == vertexA_index || previous_node_1 == vertexA_index || next_node_2 == vertexA_index || previous_node_2 == vertexA_index)
            {
                // Get elements containing vertexA_index (the common vertex)

                assert(this->mNodes[vertexA_index]->GetNumContainingElements() > 1);

                std::set<unsigned> elements_containing_vertex_A = this->mNodes[vertexA_index]->rGetContainingElementIndices();
                std::set<unsigned>::const_iterator iter = elements_containing_vertex_A.begin();
                VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_common_1 = this->GetElement(*iter);
                iter++;
                VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_common_2 = this->GetElement(*iter);

                // Calculate the number of common vertices between element_1 and element_2
                unsigned num_common_vertices = 0;
                for (unsigned i=0; i<p_element_common_1->GetNumNodes(); i++)
                {
                    for (unsigned j=0; j<p_element_common_2->GetNumNodes(); j++)
                    {
                        if (p_element_common_1->GetNodeGlobalIndex(i) == p_element_common_2->GetNodeGlobalIndex(j))
                        {
                            num_common_vertices++;
                        }
                    }
                }

                if (num_common_vertices == 1 || this->mNodes[vertexA_index]->GetNumContainingElements() > 2)
                {
                    /*
                     *  From          To
                     *   _ B              _ B
                     *    |  <---          |
                     *    |   /|\          |\
                     *    |  / | \         | \
                     *    | /  |  \        |\ \
                     *   _|/___|___\      _|_\_\
                     *     A                A
                     *
                     * The edge goes from vertexA--vertexB to vertexA--pNode--new_node--vertexB
                     */

                    // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
                    intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(vertexA_index, vertexB_index, intersection);
                    edge_ab_unit_vector = this->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

                    // Move original node and change to non-boundary node
                    pNode->rGetModifiableLocation() = intersection - 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;
                    pNode->SetAsBoundaryNode(false);

                    // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
                    c_vector<double, SPACE_DIM> new_node_location;
                    new_node_location = intersection + 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

                    // Add new node, which will always be a boundary node
                    unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(0, true, new_node_location[0], new_node_location[1]));

                    // Add the moved nodes to the element (this also updates the node)
                    this->GetElement(elementIndex)->AddNode(this->mNodes[new_node_global_index], node_A_local_index);

                    const double a_to_b_length = norm_2(vector_a_to_b);
                    c_vector<double, SPACE_DIM> a_to_new_node_vector = this->GetVectorFromAtoB(vertexA, new_node_location);
                    const double a_to_new_length = norm_2(a_to_new_node_vector);
                    const double ratio_new_node = a_to_new_length/a_to_b_length;
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex), node_A_local_index, ratio_new_node);
                    rebuilt_elements.insert(elementIndex);

                    this->GetElement(elementIndex)->AddNode(pNode, node_A_local_index);

                    c_vector<double, SPACE_DIM> a_to_p_vector = this->GetVectorFromAtoB(vertexA, pNode->rGetLocation());
                    const double ratio_p_node = norm_2(a_to_p_vector)/a_to_new_length;
                    assert(ratio_p_node<=1);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex), node_A_local_index, ratio_p_node);

                    // Add the new nodes to the original elements containing pNode (this also updates the node)
                    if (next_node_1 == previous_node_2)
                    {
                        const unsigned int insertion_local_index = (local_index_1 + p_element_1->GetNumNodes() - 1)%(p_element_1->GetNumNodes());
                        p_element_1->AddNode(this->mNodes[new_node_global_index], insertion_local_index);
                        if (mTrackMeshOperations)
                            mOperationRecorder.RecordNewEdgeOperation(p_element_1, insertion_local_index);
                        rebuilt_elements.insert(p_element_1->GetIndex());
                    }
                    else
                    {
                        assert(next_node_2 == previous_node_1);
                        const unsigned int insertion_local_index = (local_index_2 + p_element_2->GetNumNodes() - 1)%(p_element_2->GetNumNodes());
                        p_element_2->AddNode(this->mNodes[new_node_global_index], insertion_local_index);
                        if (mTrackMeshOperations)
                            mOperationRecorder.RecordNewEdgeOperation(p_element_2, insertion_local_index);
                        rebuilt_elements.insert(p_element_2->GetIndex());
                    }

                    // Check the nodes are updated correctly
                    assert(pNode->GetNumContainingElements() == 3);
                    assert(this->mNodes[new_node_global_index]->GetNumContainingElements() == 2);
                }
                else if (num_common_vertices == 2)
                {
                    /*
                     *  From          To
                     *   _ B              _ B
                     *    |<---          |
                     *    | /|\          |\
                     *    |/ | \         | \
                     *    |  |  \        |\ \
                     *   _|__|___\      _|_\_\
                     *     A              A
                     *
                     * The edge goes from vertexA--vertexB to vertexA--pNode--new_node--vertexB
                     * then vertexA is removed
                     */

                    // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
                    intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(vertexA_index, vertexB_index, intersection);
                    edge_ab_unit_vector = this->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

                    // Move original node and change to non-boundary node
                    pNode->rGetModifiableLocation() = intersection - 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;
                    pNode->SetAsBoundaryNode(false);

                    // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
                    c_vector<double, SPACE_DIM> new_node_location;
                    new_node_location = intersection + 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

                    // Add new node, which will always be a boundary node
                    unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(0, true, new_node_location[0], new_node_location[1]));

                    // Add the moved nodes to the element (this also updates the node)
                    this->GetElement(elementIndex)->AddNode(this->mNodes[new_node_global_index], node_A_local_index);

                    const double a_to_b_length = norm_2(vector_a_to_b);
                    c_vector<double, SPACE_DIM> a_to_new_node_vector = this->GetVectorFromAtoB(vertexA, new_node_location);
                    const double a_to_new_length = norm_2(a_to_new_node_vector);
                    const double ratio_new_node = a_to_new_length/a_to_b_length;
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex), node_A_local_index, ratio_new_node);

                    this->GetElement(elementIndex)->AddNode(pNode, node_A_local_index);

                    c_vector<double, SPACE_DIM> a_to_p_vector = this->GetVectorFromAtoB(vertexA, pNode->rGetLocation());
                    const double ratio_p_node = norm_2(a_to_p_vector)/a_to_new_length;
                    assert(ratio_p_node<=1);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex), node_A_local_index, ratio_p_node);

                    rebuilt_elements.insert(elementIndex);
                    // Add the new nodes to the original elements containing pNode (this also updates the node)
                    if (next_node_1 == previous_node_2)
                    {
                        const unsigned int insertion_local_index = (local_index_1 + p_element_1->GetNumNodes() - 1)%(p_element_1->GetNumNodes());
                        p_element_1->AddNode(this->mNodes[new_node_global_index], insertion_local_index);
                        if (mTrackMeshOperations)
                            mOperationRecorder.RecordNewEdgeOperation(p_element_1, insertion_local_index);
                        rebuilt_elements.insert(p_element_1->GetIndex());
                    }
                    else
                    {
                        assert(next_node_2 == previous_node_1);
                        const unsigned int insertion_local_index = (local_index_2 + p_element_2->GetNumNodes() - 1)%(p_element_2->GetNumNodes());
                        p_element_2->AddNode(this->mNodes[new_node_global_index], insertion_local_index);
                        if (mTrackMeshOperations)
                            mOperationRecorder.RecordNewEdgeOperation(p_element_2, insertion_local_index);
                        rebuilt_elements.insert(p_element_2->GetIndex());
                    }

                    // Remove vertex A from the mesh
                    const unsigned int local_index_1 = p_element_common_1->GetNodeLocalIndex(vertexA_index);
                    const unsigned int local_index_2 = p_element_common_2->GetNodeLocalIndex(vertexA_index);
                    p_element_common_1->DeleteNode(local_index_1);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeMergeOperation(p_element_common_1, local_index_1);
                    p_element_common_2->DeleteNode(local_index_2);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeMergeOperation(p_element_common_2, local_index_2);
                    assert(this->mNodes[vertexA_index]->GetNumContainingElements()==0);
                    rebuilt_elements.insert(p_element_common_1->GetIndex());
                    rebuilt_elements.insert(p_element_common_2->GetIndex());

                    this->mNodes[vertexA_index]->MarkAsDeleted();
                    mDeletedNodeIndices.push_back(vertexA_index);

                    // Check the nodes are updated correctly
                    assert(pNode->GetNumContainingElements() == 3);
                    assert(this->mNodes[new_node_global_index]->GetNumContainingElements() == 2);
                }
                else
                {
                    // This can't happen as nodes can't be on the internal edge of two elements
                    NEVER_REACHED;
                }
            }
            else if (next_node_1 == vertexB_index || previous_node_1 == vertexB_index || next_node_2 == vertexB_index || previous_node_2 == vertexB_index)
            {
                // Get elements containing vertexB_index (the common vertex)

                assert(this->mNodes[vertexB_index]->GetNumContainingElements()>1);

                std::set<unsigned> elements_containing_vertex_B = this->mNodes[vertexB_index]->rGetContainingElementIndices();
                std::set<unsigned>::const_iterator iter = elements_containing_vertex_B.begin();
                VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_common_1 = this->GetElement(*iter);
                iter++;
                VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_common_2 = this->GetElement(*iter);

                // Calculate the number of common vertices between element_1 and element_2
                unsigned num_common_vertices = 0;
                for (unsigned i=0; i<p_element_common_1->GetNumNodes(); i++)
                {
                    for (unsigned j=0; j<p_element_common_2->GetNumNodes(); j++)
                    {
                        if (p_element_common_1->GetNodeGlobalIndex(i) == p_element_common_2->GetNodeGlobalIndex(j))
                        {
                            num_common_vertices++;
                        }
                    }
                }

                if (num_common_vertices == 1 || this->mNodes[vertexB_index]->GetNumContainingElements() > 2)
                {
                    /*
                     *  From          To
                     *   _B_________      _B____
                     *    |\   |   /       | / /
                     *    | \  |  /        |/ /
                     *    |  \ | /         | /
                     *    |   \|/          |/
                     *   _|   <---        _|
                     *    A
                     *
                     * The edge goes from vertexA--vertexB to vertexA--new_node--pNode--vertexB
                     */

                    // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
                    intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(vertexA_index, vertexB_index, intersection);
                    edge_ab_unit_vector = this->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

                    // Move original node and change to non-boundary node
                    pNode->rGetModifiableLocation() = intersection + 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;
                    pNode->SetAsBoundaryNode(false);

                    // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
                    c_vector<double, SPACE_DIM> new_node_location;
                    new_node_location = intersection - 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

                    // Add new node which will always be a boundary node
                    unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(0, true, new_node_location[0], new_node_location[1]));

                    // Add the moved nodes to the element (this also updates the node)
                    this->GetElement(elementIndex)->AddNode(pNode, node_A_local_index);
                    const double a_to_b_length = norm_2(vector_a_to_b);
                    c_vector<double, SPACE_DIM> a_to_p_vector = this->GetVectorFromAtoB(vertexA, pNode->rGetLocation());
                    const double ratio_p_node = norm_2(a_to_p_vector)/a_to_b_length;
                    assert(ratio_p_node<=1);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex), node_A_local_index, ratio_p_node);

                    this->GetElement(elementIndex)->AddNode(this->mNodes[new_node_global_index], node_A_local_index);
                    c_vector<double, SPACE_DIM> a_to_new_node_vector = this->GetVectorFromAtoB(vertexA, new_node_location);
                    const double a_to_new_length = norm_2(a_to_new_node_vector);
                    const double ratio_new_node = a_to_new_length/norm_2(a_to_p_vector);
                    assert(ratio_new_node<=1);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex), node_A_local_index, ratio_new_node);
                    rebuilt_elements.insert(elementIndex);
                    // Add the new nodes to the original elements containing pNode (this also updates the node)
                    if (next_node_1 == previous_node_2)
                    {
                        p_element_2->AddNode(this->mNodes[new_node_global_index], local_index_2);
                        if (mTrackMeshOperations)
                            mOperationRecorder.RecordNewEdgeOperation(p_element_2, local_index_2);
                        rebuilt_elements.insert(p_element_2->GetIndex());
                    }
                    else
                    {
                        assert(next_node_2 == previous_node_1);
                        p_element_1->AddNode(this->mNodes[new_node_global_index], local_index_1);
                        if (mTrackMeshOperations)
                            mOperationRecorder.RecordNewEdgeOperation(p_element_1, local_index_1);
                        rebuilt_elements.insert(p_element_1->GetIndex());
                    }

                    // Check the nodes are updated correctly
                    assert(pNode->GetNumContainingElements() == 3);
                    assert(this->mNodes[new_node_global_index]->GetNumContainingElements() == 2);
                }
                else if (num_common_vertices == 2)
                {
                    /*
                     *  From          To
                     *   _B_______      _B____
                     *    |  |   /       | / /
                     *    |  |  /        |/ /
                     *    |\ | /         | /
                     *    | \|/          |/
                     *   _| <---        _|
                     *    A
                     *
                     * The edge goes from vertexA--vertexB to vertexA--new_node--pNode--vertexB
                     * then vertexB is removed
                     */

                    // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
                    intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(vertexA_index, vertexB_index, intersection);
                    edge_ab_unit_vector = this->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

                    // Move original node and change to non-boundary node
                    pNode->rGetModifiableLocation() = intersection + 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;
                    pNode->SetAsBoundaryNode(false);

                    // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
                    c_vector<double, SPACE_DIM> new_node_location;
                    new_node_location = intersection - 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

                    // Add new node which will always be a boundary node
                    unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(0, true, new_node_location[0], new_node_location[1]));

                    // Add the moved nodes to the element (this also updates the node)
                    this->GetElement(elementIndex)->AddNode(pNode, node_A_local_index);
                    const double a_to_b_length = norm_2(vector_a_to_b);
                    c_vector<double, SPACE_DIM> a_to_p_vector = this->GetVectorFromAtoB(vertexA, pNode->rGetLocation());
                    const double ratio_p_node = norm_2(a_to_p_vector)/a_to_b_length;
                    assert(ratio_p_node<=1);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex), node_A_local_index, ratio_p_node);

                    this->GetElement(elementIndex)->AddNode(this->mNodes[new_node_global_index], node_A_local_index);
                    c_vector<double, SPACE_DIM> a_to_new_node_vector = this->GetVectorFromAtoB(vertexA, new_node_location);
                    const double a_to_new_length = norm_2(a_to_new_node_vector);
                    const double ratio_new_node = a_to_new_length/norm_2(a_to_p_vector);
                    assert(ratio_new_node<=1);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex), node_A_local_index, ratio_new_node);
                    rebuilt_elements.insert(elementIndex);
                    // Add the new nodes to the original elements containing pNode (this also updates the node)
                    if (next_node_1 == previous_node_2)
                    {
                        p_element_2->AddNode(this->mNodes[new_node_global_index], local_index_2);
                        if (mTrackMeshOperations)
                            mOperationRecorder.RecordNewEdgeOperation(p_element_2, local_index_2);
                        rebuilt_elements.insert(p_element_2->GetIndex());
                    }
                    else
                    {
                        assert(next_node_2 == previous_node_1);
                        p_element_1->AddNode(this->mNodes[new_node_global_index], local_index_1);
                        if (mTrackMeshOperations)
                            mOperationRecorder.RecordNewEdgeOperation(p_element_1, local_index_1);
                        rebuilt_elements.insert(p_element_1->GetIndex());
                    }

                    // Remove vertex B from the mesh
                    const unsigned int local_index_1 = p_element_common_1->GetNodeLocalIndex(vertexB_index);
                    const unsigned int local_index_2 = p_element_common_2->GetNodeLocalIndex(vertexB_index);
                    p_element_common_1->DeleteNode(local_index_1);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeMergeOperation(p_element_common_1, local_index_1);
                    p_element_common_2->DeleteNode(local_index_2);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordEdgeMergeOperation(p_element_common_2, local_index_2);
                    rebuilt_elements.insert(p_element_common_1->GetIndex());
                    rebuilt_elements.insert(p_element_common_2->GetIndex());

                    assert(this->mNodes[vertexB_index]->GetNumContainingElements()==0);

                    this->mNodes[vertexB_index]->MarkAsDeleted();
                    mDeletedNodeIndices.push_back(vertexB_index);

                    // Check the nodes are updated correctly
                    assert(pNode->GetNumContainingElements() == 3);
                    assert(this->mNodes[new_node_global_index]->GetNumContainingElements() == 2);
                }
                else
                {
                    // This can't happen as nodes can't be on the internal edge of two elements
                    NEVER_REACHED;
                }
            }
            else
            {
                /*
                 *  From          To
                 *   _____         _______
                 *                  / | \
                 *    /|\   ^      /  |  \
                 *   / | \  |
                 *
                 * The edge goes from vertexA--vertexB to vertexA--new_node_1--pNode--new_node_2--vertexB
                 */

                // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
                intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(vertexA_index, vertexB_index, intersection);
                edge_ab_unit_vector = this->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

                // Move original node and change to non-boundary node
                pNode->rGetModifiableLocation() = intersection;
                pNode->SetAsBoundaryNode(false);

                c_vector<double, SPACE_DIM> new_node_1_location;
                new_node_1_location = intersection - mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;
                c_vector<double, SPACE_DIM> new_node_2_location;
                new_node_2_location = intersection + mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

                // Add new nodes which will always be boundary nodes
                unsigned new_node_1_global_index = this->AddNode(new Node<SPACE_DIM>(0, true, new_node_1_location[0], new_node_1_location[1]));
                unsigned new_node_2_global_index = this->AddNode(new Node<SPACE_DIM>(0, true, new_node_2_location[0], new_node_2_location[1]));

                // Add the moved and new nodes to the element (this also updates the node)
                const double a_to_b_length = norm_2(vector_a_to_b);
                const c_vector<double, SPACE_DIM> a_to_node2_vector
                = this->GetVectorFromAtoB(new_node_2_location, vertexA);
                const double a_to_node2_length = norm_2(a_to_node2_vector);
                const c_vector<double, SPACE_DIM> a_to_nodeP_vector
                = this->GetVectorFromAtoB(intersection, vertexA);
                const double a_to_nodeP_length = norm_2(a_to_nodeP_vector);
                const c_vector<double, SPACE_DIM> a_to_node1_vector
                = this->GetVectorFromAtoB(new_node_1_location, vertexA);
                const double a_to_node1_length = norm_2(a_to_node1_vector);

                this->GetElement(elementIndex)->AddNode(this->mNodes[new_node_2_global_index], node_A_local_index);
                if (mTrackMeshOperations)
                    mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex),node_A_local_index,a_to_node2_length/a_to_b_length);
                this->GetElement(elementIndex)->AddNode(pNode, node_A_local_index);
                if (mTrackMeshOperations)
                    mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex),node_A_local_index,a_to_nodeP_length/a_to_node2_length);
                this->GetElement(elementIndex)->AddNode(this->mNodes[new_node_1_global_index], node_A_local_index);
                if (mTrackMeshOperations)
                    mOperationRecorder.RecordEdgeSplitOperation(this->GetElement(elementIndex),node_A_local_index,a_to_node1_length/a_to_nodeP_length);
                rebuilt_elements.insert(elementIndex);
                // Add the new nodes to the original elements containing pNode (this also updates the node)
                if (next_node_1 == previous_node_2)
                {
                    const unsigned int inserted_local_index_1 = (local_index_1 + p_element_1->GetNumNodes() - 1)%(p_element_1->GetNumNodes());
                    p_element_1->AddNode(this->mNodes[new_node_2_global_index], inserted_local_index_1);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordNewEdgeOperation(p_element_1, inserted_local_index_1);
                    p_element_2->AddNode(this->mNodes[new_node_1_global_index], local_index_2);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordNewEdgeOperation(p_element_2, local_index_2);

                }
                else
                {
                    assert(next_node_2 == previous_node_1);
                    p_element_1->AddNode(this->mNodes[new_node_1_global_index], local_index_1);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordNewEdgeOperation(p_element_1, local_index_1);
                    const unsigned int inserted_local_index_2 = (local_index_2 + p_element_2->GetNumNodes() - 1)%(p_element_2->GetNumNodes());
                    p_element_2->AddNode(this->mNodes[new_node_2_global_index], inserted_local_index_2);
                    if (mTrackMeshOperations)
                        mOperationRecorder.RecordNewEdgeOperation(p_element_2, local_index_2);
                }
                rebuilt_elements.insert(p_element_1->GetIndex());
                rebuilt_elements.insert(p_element_2->GetIndex());
                // Check the nodes are updated correctly
                assert(pNode->GetNumContainingElements() == 3);
                assert(this->mNodes[new_node_1_global_index]->GetNumContainingElements() == 2);
                assert(this->mNodes[new_node_2_global_index]->GetNumContainingElements() == 2);
            }
        }
    }
    else
    {
        EXCEPTION("Trying to merge a node, contained in more than 2 elements, into another element, this is not possible with the vertex mesh.");
    }
    for (unsigned i:rebuilt_elements)
        this->GetElement(i)->RebuildEdges();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformVoidRemoval(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, Node<SPACE_DIM>* pNodeC)
{
    // Calculate void centroid
    c_vector<double, SPACE_DIM> nodes_midpoint = pNodeA->rGetLocation()
            + this->GetVectorFromAtoB(pNodeA->rGetLocation(), pNodeB->rGetLocation()) / 3.0
            + this->GetVectorFromAtoB(pNodeA->rGetLocation(), pNodeC->rGetLocation()) / 3.0;
    /*
     * In two steps, merge nodes A, B and C into a single node.  This is implemented in such a way that
     * the ordering of their indices does not matter.
     */
    PerformNodeMerge(pNodeA, pNodeB);
    Node<SPACE_DIM>* p_merged_node = pNodeB;

    if (pNodeB->IsDeleted())
    {
        p_merged_node = pNodeA;
    }

    PerformNodeMerge(pNodeC, p_merged_node);
    if (p_merged_node->IsDeleted())
    {
        p_merged_node = pNodeC;
    }

    p_merged_node->rGetModifiableLocation() = nodes_midpoint;

    // Tag remaining node as non-boundary
    p_merged_node->SetAsBoundaryNode(false);

    // Remove the deleted nodes and re-index
    RemoveDeletedNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::HandleHighOrderJunctions(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    unsigned node_a_rank = pNodeA->rGetContainingElementIndices().size();
    unsigned node_b_rank = pNodeB->rGetContainingElementIndices().size();

    if ((node_a_rank > 3) && (node_b_rank > 3))
    {
        // The code can't handle this case
        EXCEPTION("Both nodes involved in a swap event are contained in more than three elements");
    }
    else // the rosette degree should increase in this case
    {
        assert(node_a_rank > 3 || node_b_rank > 3);
        this->PerformRosetteRankIncrease(pNodeA, pNodeB);
        this->RemoveDeletedNodes();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformRosetteRankIncrease(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    /*
     * One of the nodes will have 3 containing element indices, the other
     * will have at least four. We first identify which node is which.
     */

    unsigned node_a_index = pNodeA->GetIndex();
    unsigned node_b_index = pNodeB->GetIndex();

    unsigned node_a_rank = pNodeA->rGetContainingElementIndices().size();
    unsigned node_b_rank = pNodeB->rGetContainingElementIndices().size();

    unsigned lo_rank_index = (node_a_rank < node_b_rank) ? node_a_index : node_b_index;
    unsigned hi_rank_index = (node_a_rank < node_b_rank) ? node_b_index : node_a_index;

    // Get pointers to the nodes, sorted by index
    Node<SPACE_DIM>* p_lo_rank_node = this->GetNode(lo_rank_index);
    Node<SPACE_DIM>* p_hi_rank_node = this->GetNode(hi_rank_index);

    // Find the sets of elements containing each of the nodes, sorted by index
    std::set<unsigned> lo_rank_elem_indices = p_lo_rank_node->rGetContainingElementIndices();
    std::set<unsigned> hi_rank_elem_indices = p_hi_rank_node->rGetContainingElementIndices();

    /**
     * The picture below shows the situation we are in.  The central node (marked with an 'X')
     * is contained in (at least) four elements already (A, B, C and D), and this is the
     * rosette node which we have designated as hi_rank_node.
     *
     * The node shared by elements C, D and E has come within the cell rearrangement threshold
     * of the rosette node in order for this method to have been called.  We have designated
     * this node lo_rank_node.
     *
     * We now 'merge' hi_rank_node and lo_rank_node, but keep hi_rank_node where it is
     * (which is why we don't call PerformNodeMerge(), as that would move both nodes to
     * their average location).
     *
     * To accomplish this merge, we need do nothing to elements A and B.  We remove lo_rank_node
     * from elements C and D, and we replace lo_rank_node by hi_rank_node in element E.
     *
     *
     *      \  A  /
     *       \   /
     *        \ /
     *     B   X   C
     *        / \ ______
     *       /   \
     *      /  D  \  E
     */

    for (std::set<unsigned>::const_iterator it = lo_rank_elem_indices.begin();
         it != lo_rank_elem_indices.end();
         ++it)
    {
        // Find the local index of lo_rank_node in this element
        unsigned lo_rank_local_index = this->mElements[*it]->GetNodeLocalIndex(lo_rank_index);
        assert(lo_rank_local_index < UINT_MAX); // double check this element contains lo_rank_node

        /*
         * If this element already contains the hi_rank_node, we are in the situation of elements
         * C and D above, so we just remove lo_rank_node.
         *
         * Otherwise, we are in element E, so we must replace lo_rank_node with high-rank node,
         * and remove it from mNodes.
         *
         * We can check whether hi_rank_node is in this element using the set::count() method.
         */

        if (hi_rank_elem_indices.count(*it) > 0)
        {
            //For node merge recording
            std::vector<unsigned int> oldIds(this->mElements[*it]->GetNumEdges(),0);

            for (unsigned int i=0; i<oldIds.size(); ++i)
            {
                oldIds[i] = this->mElements[*it]->GetEdge(i)->GetIndex();
            }
            // Delete lo_rank_node from current element
            this->mElements[*it]->DeleteNode(lo_rank_local_index);
            //Record edge shrinkage
            const unsigned int hi_rank_local_index = this->mElements[*it]->GetNodeLocalIndex(hi_rank_index);
            if (mTrackMeshOperations)
                mOperationRecorder.RecordNodeMergeOperation(oldIds, this->mElements[*it],
                                                            std::pair<unsigned int, unsigned int>(hi_rank_local_index, lo_rank_local_index));
        }
        else
        {
            // Update lo_rank_node with all information (including index and location) of hi_rank_node
            this->mElements[*it]->UpdateNode(lo_rank_local_index, p_hi_rank_node);
        }
    }

    // Tidy up the mesh by ensuring the global instance of lo_rank_node is deleted
    assert(!(this->mNodes[lo_rank_index]->IsDeleted()));
    this->mNodes[lo_rank_index]->MarkAsDeleted();
    this->mDeletedNodeIndices.push_back(lo_rank_index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformProtorosetteResolution(Node<SPACE_DIM>* pProtorosetteNode)
{
    // Double check we are dealing with a protorosette
    assert(pProtorosetteNode->rGetContainingElementIndices().size() == 4);

    // Get random number (0, 1, 2 or 3), as the resolution axis is assumed to be random
    unsigned random_elem_increment = RandomNumberGenerator::Instance()->randMod(4);

    // Find global indices of elements around the protorosette node
    std::set<unsigned> protorosette_node_containing_elem_indices = pProtorosetteNode->rGetContainingElementIndices();

    // Select random element by advancing iterator a random number times
    std::set<unsigned>::const_iterator elem_index_iter(protorosette_node_containing_elem_indices.begin());
    advance(elem_index_iter, random_elem_increment);

    /**
     * Ordering elements as follows:
     *
     *      \  A  /
     *       \   /
     *        \ /
     *     B   X   D
     *        / \
     *       /   \
     *      /  C  \
     *
     * Element A is the randomly chosen element.  Element C, which is directly
     * opposite A, will end up separated from A, while the two elements B and
     * D which start adjacent to A will end up sharing a common edge:
     *
     *      \  A  /
     *       \   /
     *        \ /
     *         |
     *    B    |    D
     *         |
     *        / \
     *       /   \
     *      /  C  \
     *
     */

    /*
     * We need to find the global indices of elements B, C and D.  We do this with set intersections.
     */

    unsigned elem_a_idx = *elem_index_iter;
    unsigned elem_b_idx = UINT_MAX;
    unsigned elem_c_idx = UINT_MAX;
    unsigned elem_d_idx = UINT_MAX;

    // Get pointer to element we've chosen at random (element A)
    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_a = this->GetElement(elem_a_idx);

    // Get all necessary info about element A and the protorosette node
    unsigned num_nodes_elem_a = p_elem_a->GetNumNodes();
    unsigned protorosette_node_global_idx = pProtorosetteNode->GetIndex();
    unsigned protorosette_node_local_idx = p_elem_a->GetNodeLocalIndex(protorosette_node_global_idx);

    // Find global indices of previous (cw) and next (ccw) nodes, locally, from the protorosette node, in element A
    unsigned prev_node_global_idx = p_elem_a->GetNodeGlobalIndex((protorosette_node_local_idx + num_nodes_elem_a - 1) % num_nodes_elem_a);
    unsigned next_node_global_idx = p_elem_a->GetNodeGlobalIndex((protorosette_node_local_idx + 1) % num_nodes_elem_a);

    // Get the set of elements the previous and next nodes are contained in
    Node<SPACE_DIM>* p_prev_node = this->GetNode(prev_node_global_idx);
    Node<SPACE_DIM>* p_next_node = this->GetNode(next_node_global_idx);
    std::set<unsigned> prev_node_elem_indices = p_prev_node->rGetContainingElementIndices();
    std::set<unsigned> next_node_elem_indices = p_next_node->rGetContainingElementIndices();

    // Perform set intersections with the set of element indices which the protorosette node is contained in
    std::set<unsigned> intersection_with_prev;
    std::set<unsigned> intersection_with_next;

    // This intersection should contain just global indices for elements A and B
    std::set_intersection(protorosette_node_containing_elem_indices.begin(),
                          protorosette_node_containing_elem_indices.end(),
                          prev_node_elem_indices.begin(),
                          prev_node_elem_indices.end(),
                          std::inserter(intersection_with_prev, intersection_with_prev.begin()));

    // This intersection should contain just global indices for elements A and D
    std::set_intersection(protorosette_node_containing_elem_indices.begin(),
                          protorosette_node_containing_elem_indices.end(),
                          next_node_elem_indices.begin(),
                          next_node_elem_indices.end(),
                          std::inserter(intersection_with_next, intersection_with_next.begin()));

    assert(intersection_with_prev.size() == 2);
    assert(intersection_with_next.size() == 2);

    // Get global index of element B
    if (*intersection_with_prev.begin() != elem_a_idx)
    {
        elem_b_idx = *(intersection_with_prev.begin());
    }
    else
    {
        elem_b_idx = *(++(intersection_with_prev.begin()));
    }
    assert(elem_b_idx < UINT_MAX);

    // Get global index of element D
    if (*intersection_with_next.begin() != elem_a_idx)
    {
        elem_d_idx = *(intersection_with_next.begin());
    }
    else
    {
        elem_d_idx = *(++(intersection_with_next.begin()));
    }
    assert(elem_d_idx < UINT_MAX);

    // By elimination, the remaining unassigned index in the original set must be global index of element C
    for (elem_index_iter = protorosette_node_containing_elem_indices.begin();
         elem_index_iter != protorosette_node_containing_elem_indices.end();
         ++elem_index_iter)
    {
        if ((*elem_index_iter != elem_a_idx) && (*elem_index_iter != elem_b_idx) && (*elem_index_iter != elem_d_idx))
        {
            elem_c_idx = *elem_index_iter;
        }
    }
    assert(elem_c_idx < UINT_MAX);

    /**
     * Next, we compute where to place the two nodes which will replace the single protorosette node.
     *
     * We place each node along the line joining the protorosette node to the centroid of element which will contain it,
     * and the distance along this line is half of the swap distance.
     *
     * To do this, we will move the existing protorosette node in to element A, and create a new node in element C.  We
     * then need to tidy up the nodes by adding the new node to elements B, C and D, and removing the protorosette node
     * from element C.
     *
     * NOTE: as the protorosette node was not necessarily on the line joining the centroids of elements A and C, unlike
     * in a T1 swap, the new node locations will not necessarily be the full swap distance apart.
     */

    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_b = this->GetElement(elem_b_idx);
    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_c = this->GetElement(elem_c_idx);
    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_d = this->GetElement(elem_d_idx);

    double swap_distance = (this->mCellRearrangementRatio) * (this->mCellRearrangementThreshold);

    // Get normalized vectors to centre of elements A and B from protorosette node
    c_vector<double, SPACE_DIM> node_to_elem_a_centre = this->GetCentroidOfElement(elem_a_idx) - pProtorosetteNode->rGetLocation();
    node_to_elem_a_centre /= norm_2(node_to_elem_a_centre);

    c_vector<double, SPACE_DIM> node_to_elem_c_centre = this->GetCentroidOfElement(elem_c_idx) - pProtorosetteNode->rGetLocation();
    node_to_elem_c_centre /= norm_2(node_to_elem_c_centre);

    // Calculate new node locations
    c_vector<double, SPACE_DIM> new_location_of_protorosette_node = pProtorosetteNode->rGetLocation() + (0.5 * swap_distance) * node_to_elem_a_centre;
    c_vector<double, SPACE_DIM> location_of_new_node = pProtorosetteNode->rGetLocation() + (0.5 * swap_distance) * node_to_elem_c_centre;

    // Move protorosette node to new location
    pProtorosetteNode->rGetModifiableLocation() = new_location_of_protorosette_node;

    // Create new node in correct location
    unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(this->GetNumNodes(), location_of_new_node, false));
    Node<SPACE_DIM>* p_new_node = this->GetNode(new_node_global_index);

    /**
     * Here, we add the new node to elements B, C and D.
     *
     * The method AddNode() takes the local index of a node i, where the new node is to be inserted between nodes
     * i and i+1, so we need to find the local idx of the node directly before (clockwise from) where we wish to
     * insert it.
     *
     * For elements C and D, we just need the local index of the protorosette node, but for element B we need the one
     * before, modulo number of elements.
     */
    unsigned local_idx_elem_b = p_elem_b->GetNodeLocalIndex(protorosette_node_global_idx);
    local_idx_elem_b = (local_idx_elem_b + p_elem_b->GetNumNodes() - 1) % p_elem_b->GetNumNodes();
    unsigned local_idx_elem_c = p_elem_c->GetNodeLocalIndex(protorosette_node_global_idx);
    unsigned local_idx_elem_d = p_elem_d->GetNodeLocalIndex(protorosette_node_global_idx);

    p_elem_b->AddNode(p_new_node, local_idx_elem_b);
    if (mTrackMeshOperations)
        mOperationRecorder.RecordNewEdgeOperation(p_elem_b, local_idx_elem_b);
    p_elem_c->AddNode(p_new_node, local_idx_elem_c);
    p_elem_d->AddNode(p_new_node, local_idx_elem_d);
    if (mTrackMeshOperations)
        mOperationRecorder.RecordNewEdgeOperation(p_elem_d, local_idx_elem_d);

    // All that is left is to remove the original protorosette node from element C
    p_elem_c->DeleteNode(p_elem_c->GetNodeLocalIndex(protorosette_node_global_idx));
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformRosetteRankDecrease(Node<SPACE_DIM>* pRosetteNode)
{
    unsigned rosette_rank = pRosetteNode->rGetContainingElementIndices().size();

    // Double check we're dealing with a rosette
    assert(rosette_rank > 4);

    // Get random number in [0, 1, ..., n) where n is rank of rosette, as the resolution axis is assumed to be random
    unsigned random_elem_increment = RandomNumberGenerator::Instance()->randMod(rosette_rank);

    // Find global indices of elements around the protorosette node
    std::set<unsigned> rosette_node_containing_elem_indices = pRosetteNode->rGetContainingElementIndices();

    // Select random element by advancing iterator a random number times
    std::set<unsigned>::const_iterator elem_index_iter(rosette_node_containing_elem_indices.begin());
    advance(elem_index_iter, random_elem_increment);

    /**
     * We have now picked a vertex element at random from the rosette.
     * This element will now be disconnected from the rosette in a manner
     * analogous to performing a T1 swap.
     *
     * The node at the centre of the rosette will not move, and a new
     * node will be created a suitable distance away, along a line joining
     * the rosette node and the centroid of the element we have randomly
     * selected to move.
     *
     * Ordering of elements is as follows:
     *
     *      \  S  /
     *       \   /
     *     N  \ /  P
     *    -----X-----
     *        / \
     *       /   \
     *      /     \
     *
     * where element S is the selected element,
     * N is the next (counterclockwise) element from S, and
     * P is the previous (clockwise) element from S.
     *
     * Elements N and P will end up sharing an edge:
     *
     *      \  S  /
     *       \   /
     *        \ /
     *      N  |  P
     *         |
     *    -----*-----
     *        / \
     *       /   \
     *      /     \
     *
     */

    /*
     * We need to find the global indices of elements N and P.  We do this with set intersections.
     */

    // Get the vertex element S (which we randomly selected)
    unsigned elem_s_idx = *elem_index_iter;
    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_s = this->GetElement(elem_s_idx);

    unsigned elem_n_idx = UINT_MAX;
    unsigned elem_p_idx = UINT_MAX;

    // Get all necessary info about element S and the rosette node
    unsigned num_nodes_elem_s = p_elem_s->GetNumNodes();
    unsigned rosette_node_global_idx = pRosetteNode->GetIndex();
    unsigned rosette_node_local_idx = p_elem_s->GetNodeLocalIndex(rosette_node_global_idx);

    // Find global indices of previous (cw) and next (ccw) nodes, locally, from the rosette node, in element S
    unsigned prev_node_global_idx = p_elem_s->GetNodeGlobalIndex((rosette_node_local_idx + num_nodes_elem_s - 1) % num_nodes_elem_s);
    unsigned next_node_global_idx = p_elem_s->GetNodeGlobalIndex((rosette_node_local_idx + 1) % num_nodes_elem_s);

    // Get the set of elements that the previous and next nodes are contained in
    Node<SPACE_DIM>* p_prev_node = this->GetNode(prev_node_global_idx);
    Node<SPACE_DIM>* p_next_node = this->GetNode(next_node_global_idx);
    std::set<unsigned> prev_node_elem_indices = p_prev_node->rGetContainingElementIndices();
    std::set<unsigned> next_node_elem_indices = p_next_node->rGetContainingElementIndices();

    // Perform set intersections with the set of element indices that the rosette node is contained in
    std::set<unsigned> intersection_with_prev;
    std::set<unsigned> intersection_with_next;

    // This intersection should contain just global indices for elements S and N
    std::set_intersection(rosette_node_containing_elem_indices.begin(),
                          rosette_node_containing_elem_indices.end(),
                          prev_node_elem_indices.begin(),
                          prev_node_elem_indices.end(),
                          std::inserter(intersection_with_prev, intersection_with_prev.begin()));

    // This intersection should contain just global indices for elements S and P
    std::set_intersection(rosette_node_containing_elem_indices.begin(),
                          rosette_node_containing_elem_indices.end(),
                          next_node_elem_indices.begin(),
                          next_node_elem_indices.end(),
                          std::inserter(intersection_with_next, intersection_with_next.begin()));

    assert(intersection_with_prev.size() == 2);
    assert(intersection_with_next.size() == 2);

    // Get global index of element N
    if (*intersection_with_prev.begin() != elem_s_idx)
    {
        elem_n_idx = *intersection_with_prev.begin();
    }
    else
    {
        elem_n_idx = *(++(intersection_with_prev.begin()));
    }
    assert(elem_n_idx < UINT_MAX);

    // Get global index of element P
    if (*intersection_with_next.begin() != elem_s_idx)
    {
        elem_p_idx = *intersection_with_next.begin();
    }
    else
    {
        elem_p_idx = *(++(intersection_with_next.begin()));
    }
    assert(elem_p_idx < UINT_MAX);

    /**
     * Next, we compute where to place the new node which will separate the rosette node from element S.
     *
     * We place this node along the line joining the rosette node to the centroid of element S, at a distance of the
     * swap distance ( (rearrangement ratio) x (rearrangement threshold) )
     *
     * To do this, we create a new node in element S.  We then need to tidy up the nodes by adding the new node to
     * elements S, N and P, and by removing the rosette node from element S.
     */

    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_n = this->GetElement(elem_p_idx);
    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_p = this->GetElement(elem_n_idx);

    double swap_distance = (this->mCellRearrangementRatio) * (this->mCellRearrangementThreshold);

    // Calculate location of new node
    c_vector<double, 2> node_to_selected_elem = this->GetCentroidOfElement(elem_s_idx) - pRosetteNode->rGetLocation();
    node_to_selected_elem /= norm_2(node_to_selected_elem);
    c_vector<double, 2> new_node_location = pRosetteNode->rGetLocation() + (swap_distance * node_to_selected_elem);

    // Create new node in correct location
    unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(this->GetNumNodes(), new_node_location, false));
    Node<SPACE_DIM>* p_new_node = this->GetNode(new_node_global_index);

    /**
     * Here, we add the new node to elements S, N and P, and remove the rosette node from element S.
     *
     * The method AddNode() takes the local index of a node i, where the new node is to be inserted between nodes
     * i and i+1, so we need to find the local idx of the node directly before (clockwise from) where we wish to
     * insert it.
     *
     * For elements S and P, we just need the local index of the rosette node, but for element N we need the one
     * before, modulo number of elements.
     */

    // Add new node, and remove rosette node, from element S
    unsigned node_local_idx_in_elem_s = p_elem_s->GetNodeLocalIndex(rosette_node_global_idx);
    p_elem_s->AddNode(p_new_node, node_local_idx_in_elem_s);
    p_elem_s->DeleteNode(node_local_idx_in_elem_s);

    // Add new node to element N
    unsigned node_local_idx_in_elem_n = p_elem_n->GetNodeLocalIndex(rosette_node_global_idx);
    node_local_idx_in_elem_n = (node_local_idx_in_elem_n + p_elem_n->GetNumNodes() - 1) % p_elem_n->GetNumNodes();
    p_elem_n->AddNode(p_new_node, node_local_idx_in_elem_n);
    if (mTrackMeshOperations)
        mOperationRecorder.RecordNewEdgeOperation(p_elem_n, node_local_idx_in_elem_n);
    // Add new node to element P
    unsigned node_local_idx_in_elem_p = p_elem_p->GetNodeLocalIndex(rosette_node_global_idx);
    p_elem_p->AddNode(p_new_node, node_local_idx_in_elem_p);
    if (mTrackMeshOperations)
        mOperationRecorder.RecordNewEdgeOperation(p_elem_p, node_local_idx_in_elem_p);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::CheckForRosettes()
{
    /**
     * First, we loop over each node and populate vectors of protorosette and rosette nodes which need to undergo
     * resolution.
     *
     * We do not perform the resolution events in this initial loop because the resolution events involve changing
     * nodes in the mesh.
     */

    // Vectors to store the nodes that need resolution events
    std::vector<Node<SPACE_DIM>* > protorosette_nodes;
    std::vector<Node<SPACE_DIM>* > rosette_nodes;

    // First loop in which we populate these vectors
    unsigned num_nodes = this->GetNumAllNodes();
    for (unsigned node_idx = 0 ; node_idx < num_nodes ; node_idx++)
    {
        Node<SPACE_DIM>* current_node = this->GetNode(node_idx);
        unsigned node_rank = current_node->rGetContainingElementIndices().size();

        if (node_rank < 4)
        {
            // Nothing to do if the node is not high-rank
            continue;
        }
        else if (node_rank == 4)
        {
            // For protorosette nodes, we check against a random number to decide if resolution is necessary
            if (mProtorosetteResolutionProbabilityPerTimestep >= RandomNumberGenerator::Instance()->ranf())
            {
                protorosette_nodes.push_back(current_node);
            }
        }
        else // if (node_rank > 4)
        {
            // For rosette nodes, we check against a random number to decide if resolution is necessary
            if (mRosetteResolutionProbabilityPerTimestep >= RandomNumberGenerator::Instance()->ranf())
            {
                rosette_nodes.push_back(current_node);
            }
        }
    }

    /**
     * Finally, we loop over the contents of each node vector and perform the necessary resolution events.
     *
     * Because each resolution event changes nodes, we include several assertions to catch possible unconsidered
     * behaviour.
     */

    // First, resolve any protorosettes
    for (unsigned node_idx = 0 ; node_idx < protorosette_nodes.size() ; node_idx++)
    {
        Node<SPACE_DIM>* current_node = protorosette_nodes[node_idx];

        // Verify that node has not been marked for deletion, and that it is still contained in four elements
        assert( !(current_node->IsDeleted()) );
        assert( current_node->rGetContainingElementIndices().size() == 4 );

        // Perform protorosette resolution
        this->PerformProtorosetteResolution(current_node);
    }

    // Finally, resolve any rosettes
    for (unsigned node_idx = 0 ; node_idx < rosette_nodes.size() ; node_idx++)
    {
        Node<SPACE_DIM>* current_node = rosette_nodes[node_idx];

        // Verify that node has not been marked for deletion, and that it is still contained in at least four elements
        assert( !(current_node->IsDeleted()) );
        assert( current_node->rGetContainingElementIndices().size() > 4 );

        // Perform protorosette resolution
        this->PerformRosetteRankDecrease(current_node);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 2> MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::WidenEdgeOrCorrectIntersectionLocationIfNecessary(
        unsigned indexA, unsigned indexB, c_vector<double,2> intersection)
{
    /**
     * If the edge is shorter than 4.0*mCellRearrangementRatio*mCellRearrangementThreshold move vertexA and vertexB
     * 4.0*mCellRearrangementRatio*mCellRearrangementThreshold apart.
     * \todo investigate if moving A and B causes other issues with nearby nodes (see #2401)
     *
     * Note: this distance is so that there is always enough room for new nodes (if necessary)
     * \todo currently this assumes a worst case scenario of 3 nodes between A and B could be less movement for other cases
     *       (see #1399 and #2401)
     */
    c_vector<double, SPACE_DIM> vertexA = this->GetNode(indexA)->rGetLocation();
    c_vector<double, SPACE_DIM> vertexB = this->GetNode(indexB)->rGetLocation();
    c_vector<double, SPACE_DIM> vector_a_to_b = this->GetVectorFromAtoB(vertexA, vertexB);

    if (norm_2(vector_a_to_b) < 4.0*mCellRearrangementRatio*mCellRearrangementThreshold)
    {
        WARNING("Trying to merge a node onto an edge which is too small.");

        c_vector<double, SPACE_DIM> centre_a_and_b = vertexA + 0.5*vector_a_to_b;

        vertexA = centre_a_and_b  - 2.0*mCellRearrangementRatio*mCellRearrangementThreshold*vector_a_to_b/norm_2(vector_a_to_b);
        ChastePoint<SPACE_DIM> vertex_A_point(vertexA);
        SetNode(indexA, vertex_A_point);

        vertexB = centre_a_and_b  + 2.0*mCellRearrangementRatio*mCellRearrangementThreshold*vector_a_to_b/norm_2(vector_a_to_b);
        ChastePoint<SPACE_DIM> vertex_B_point(vertexB);
        SetNode(indexB, vertex_B_point);

        intersection = centre_a_and_b;
    }

    // Reset distances
    vector_a_to_b = this->GetVectorFromAtoB(vertexA, vertexB);
    c_vector<double,2> edge_ab_unit_vector = vector_a_to_b/norm_2(vector_a_to_b);

    // Reset the intersection away from vertices A and B to allow enough room for new nodes
    /**
     * If the intersection is within mCellRearrangementRatio^2*mCellRearrangementThreshold of vertexA or vertexB move it
     * mCellRearrangementRatio^2*mCellRearrangementThreshold away.
     *
     * Note: this distance so that there is always enough room for new nodes (if necessary).
     * \todo currently this assumes a worst case scenario of 3 nodes between A and B; could be less movement for other cases
     *       (see #2401)
     */
    if (norm_2(intersection - vertexA) < 2.0*mCellRearrangementRatio*mCellRearrangementThreshold)
    {
        intersection = vertexA + 2.0*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;
    }
    if (norm_2(intersection - vertexB) < 2.0*mCellRearrangementRatio*mCellRearrangementThreshold)
    {
        intersection = vertexB - 2.0*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;
    }
    return intersection;
}

template <unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetMeshOperationTracking(const bool track)
{
    mTrackMeshOperations = track;
    if (track)
    {
        EdgeHelper<SPACE_DIM> *pEdgeHelper(&(this->mEdges));
        mOperationRecorder.SetEdgeHelper(pEdgeHelper);
    }
}

template <unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>* MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetOperationRecorder()
{
    VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>* value(&mOperationRecorder);
    return value;
}

// Explicit instantiation
template class MutableVertexMesh<1,1>;
template class MutableVertexMesh<1,2>;
template class MutableVertexMesh<1,3>;
template class MutableVertexMesh<2,2>;
template class MutableVertexMesh<2,3>;
template class MutableVertexMesh<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MutableVertexMesh)
