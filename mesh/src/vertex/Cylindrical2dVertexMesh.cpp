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

#include "Cylindrical2dVertexMesh.hpp"
#include "Cylindrical2dMesh.hpp"

Cylindrical2dVertexMesh::Cylindrical2dVertexMesh(double width,
                                                 std::vector<Node<2>*> nodes,
                                                 std::vector<VertexElement<2, 2>*> vertexElements,
                                                 double cellRearrangementThreshold,
                                                 double t2Threshold)
    : MutableVertexMesh<2,2>(nodes, vertexElements, cellRearrangementThreshold, t2Threshold),
      mWidth(width),
      mpMeshForVtk(nullptr)
{
    // Call ReMesh() to remove any deleted nodes and relabel
    ReMesh();
}

Cylindrical2dVertexMesh::Cylindrical2dVertexMesh(Cylindrical2dMesh& rMesh)
    : mWidth(rMesh.GetWidth(0)),
      mpMeshForVtk(nullptr)
{
    mpDelaunayMesh = &rMesh;

    // Reset member variables and clear mNodes, mFaces and mElements
    Clear();

    unsigned num_elements = mpDelaunayMesh->GetNumAllNodes();
    unsigned num_nodes = mpDelaunayMesh->GetNumAllElements();

    // Allocate memory for mNodes and mElements
    this->mNodes.reserve(num_nodes);

    // Create as many elements as there are nodes in the mesh
    mElements.reserve(num_elements);

    for (unsigned elem_index=0; elem_index<num_elements; elem_index++)
    {
        VertexElement<2,2>* p_element = new VertexElement<2,2>(elem_index);
        mElements.push_back(p_element);
    }

    // Populate mNodes
    GenerateVerticesFromElementCircumcentres(rMesh);

    // Loop over all generated nodes and check they're not outside [0,mWidth]
    for (unsigned i=0; i<num_nodes; i++)
    {
        double x_location = mNodes[i]->rGetLocation()[0];
        if (x_location < 0)
        {
            mNodes[i]->rGetModifiableLocation()[0] = x_location + mWidth;
        }
        else if (x_location > mWidth)
        {
            mNodes[i]->rGetModifiableLocation()[0] = x_location - mWidth;
        }
    }

    // Loop over elements of the Delaunay mesh (which are nodes/vertices of this mesh)
    for (unsigned i=0; i<num_nodes; i++)
    {
        // Loop over nodes owned by this triangular element in the Delaunay mesh
        // Add this node/vertex to each of the 3 vertex elements
        for (unsigned local_index=0; local_index<3; local_index++)
        {
            unsigned elem_index = mpDelaunayMesh->GetElement(i)->GetNodeGlobalIndex(local_index);
            unsigned num_nodes_in_elem = mElements[elem_index]->GetNumNodes();
            unsigned end_index = num_nodes_in_elem>0 ? num_nodes_in_elem-1 : 0;

            mElements[elem_index]->AddNode(this->mNodes[i], end_index);
        }
    }

    // Reorder mNodes anticlockwise
    for (unsigned elem_index=0; elem_index<mElements.size(); elem_index++)
    {
        /**
         * Create a std::vector of pairs, where each pair comprises the angle
         * between the centre of the Voronoi element and each node with that
         * node's global index in the Voronoi mesh.
         */
        std::vector<std::pair<double, unsigned> > index_angle_list;
        for (unsigned local_index=0; local_index<mElements[elem_index]->GetNumNodes(); local_index++)
        {
            c_vector<double, 2> vectorA = mpDelaunayMesh->GetNode(elem_index)->rGetLocation();
            c_vector<double, 2> vectorB = mElements[elem_index]->GetNodeLocation(local_index);
            c_vector<double, 2> centre_to_vertex = mpDelaunayMesh->GetVectorFromAtoB(vectorA, vectorB);

            double angle = atan2(centre_to_vertex(1), centre_to_vertex(0));
            unsigned global_index = mElements[elem_index]->GetNodeGlobalIndex(local_index);

            std::pair<double, unsigned> pair(angle, global_index);
            index_angle_list.push_back(pair);
        }

        // Sort the list in order of increasing angle
        sort(index_angle_list.begin(), index_angle_list.end());

        // Create a new Voronoi element and pass in the appropriate Nodes, ordered anticlockwise
        VertexElement<2,2>* p_new_element = new VertexElement<2,2>(elem_index);
        for (unsigned count = 0; count < index_angle_list.size(); count++)
        {
            unsigned local_index = count>1 ? count-1 : 0;
            p_new_element->AddNode(mNodes[index_angle_list[count].second], local_index);
        }

        // Replace the relevant member of mElements with this Voronoi element
        delete mElements[elem_index];
        mElements[elem_index] = p_new_element;
    }

    this->mMeshChangesDuringSimulation = false;
}

Cylindrical2dVertexMesh::Cylindrical2dVertexMesh()
    : mpMeshForVtk(nullptr)
{
}

Cylindrical2dVertexMesh::~Cylindrical2dVertexMesh()
{
    if (mpMeshForVtk != nullptr)
    {
         delete mpMeshForVtk;
    }
}

c_vector<double, 2> Cylindrical2dVertexMesh::GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
{
    assert(mWidth > 0.0);

    c_vector<double, 2> vector = rLocation2 - rLocation1;
    vector[0] = fmod(vector[0], mWidth);

    // If the points are more than halfway around the cylinder apart, measure the other way
    if (vector[0] > 0.5*mWidth)
    {
        vector[0] -= mWidth;
    }
    else if (vector[0] < -0.5*mWidth)
    {
        vector[0] += mWidth;
    }
    return vector;
}

void Cylindrical2dVertexMesh::SetNode(unsigned nodeIndex, ChastePoint<2> point)
{
    double x_coord = point.rGetLocation()[0];

    // Perform a periodic movement if necessary
    if (x_coord >= mWidth)
    {
        // Move point to the left
        point.SetCoordinate(0, x_coord - mWidth);
    }
    else if (x_coord < 0.0)
    {
        // Move point to the right
        point.SetCoordinate(0, x_coord + mWidth);
    }

    // Update the node's location
    MutableVertexMesh<2,2>::SetNode(nodeIndex, point);
}

double Cylindrical2dVertexMesh::GetWidth(const unsigned& rDimension) const
{
    double width = 0.0;
    assert(rDimension==0 || rDimension==1);
    if (rDimension==0)
    {
        width = mWidth;
    }
    else
    {
        width = VertexMesh<2,2>::GetWidth(rDimension);
    }
    return width;
}

unsigned Cylindrical2dVertexMesh::AddNode(Node<2>* pNewNode)
{
    unsigned node_index = MutableVertexMesh<2,2>::AddNode(pNewNode);

    // If necessary move it to be back on the cylinder
    ChastePoint<2> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point);

    return node_index;
}

void Cylindrical2dVertexMesh::Scale(const double xScale, const double yScale, const double zScale)
{
    assert(zScale == 1.0);

    AbstractMesh<2, 2>::Scale(xScale, yScale);

    // Also rescale the width of the mesh (this effectively scales the domain)
    mWidth *= xScale;
}

VertexMesh<2, 2>* Cylindrical2dVertexMesh::GetMeshForVtk()
{
    unsigned num_nodes = GetNumNodes();

    std::vector<Node<2>*> temp_nodes(2*num_nodes);
    std::vector<VertexElement<2, 2>*> elements;

    // Create four copies of each node
    for (unsigned index=0; index<num_nodes; index++)
    {
        c_vector<double, 2> location;
        location = GetNode(index)->rGetLocation();

        // Node copy at original location
        Node<2>* p_node = new Node<2>(index, false, location[0], location[1]);
        temp_nodes[index] = p_node;

        // Node copy shifted right
        p_node = new Node<2>(num_nodes + index, false, location[0] + mWidth, location[1]);
        temp_nodes[num_nodes + index] = p_node;
    }

    // Iterate over elements
    for (VertexMesh<2,2>::VertexElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        unsigned num_nodes_in_elem = elem_iter->GetNumNodes();

        std::vector<Node<2>*> elem_nodes;

        // Compute whether the element straddles either periodic boundary
        bool element_straddles_left_right_boundary = false;

        const c_vector<double, 2>& r_this_node_location = elem_iter->GetNode(0)->rGetLocation();
        for (unsigned local_index=0; local_index<num_nodes_in_elem; local_index++)
        {
            const c_vector<double, 2>& r_next_node_location = elem_iter->GetNode((local_index+1)%num_nodes_in_elem)->rGetLocation();
            c_vector<double, 2> vector;
            vector = r_next_node_location - r_this_node_location;

            if (fabs(vector[0]) > 0.5*mWidth)
            {
                element_straddles_left_right_boundary = true;
            }
        }

        // Use the above information when duplicating the element
        for (unsigned local_index=0; local_index<num_nodes_in_elem; local_index++)
        {
            unsigned this_node_index = elem_iter->GetNodeGlobalIndex(local_index);

            // If the element straddles the left/right periodic boundary...
            if (element_straddles_left_right_boundary)
            {
                // ...and this node is located to the left of the centre of the mesh...
                bool node_is_right_of_centre = (elem_iter->GetNode(local_index)->rGetLocation()[0] - 0.5*mWidth > 0);
                if (!node_is_right_of_centre)
                {
                    // ...then choose the equivalent node to the right
                    this_node_index += num_nodes;
                }
            }

            elem_nodes.push_back(temp_nodes[this_node_index]);
        }

        VertexElement<2,2>* p_element = new VertexElement<2,2>(elem_index, elem_nodes);
        elements.push_back(p_element);
    }

    // Now delete any nodes from the mesh for VTK that are not contained in any elements
    std::vector<Node<2>*> nodes;
    unsigned count = 0;
    for (unsigned index=0; index<temp_nodes.size(); index++)
    {
        unsigned num_elems_containing_this_node = temp_nodes[index]->rGetContainingElementIndices().size();

        if (num_elems_containing_this_node == 0)
        {
            // Avoid memory leak
            delete temp_nodes[index];
        }
        else
        {
            temp_nodes[index]->SetIndex(count);
            nodes.push_back(temp_nodes[index]);
            count++;
        }
    }

    if (mpMeshForVtk != nullptr)
    {
        delete mpMeshForVtk;
    }

    mpMeshForVtk = new VertexMesh<2,2>(nodes, elements);
    return mpMeshForVtk;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Cylindrical2dVertexMesh)
