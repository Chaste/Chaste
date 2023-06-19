/*

Copyright (c) 2005-2023, University of Oxford.
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

#include "Toroidal2dVertexMesh.hpp"
#include "Toroidal2dMesh.hpp"

Toroidal2dVertexMesh::Toroidal2dVertexMesh(double width,
                                           double height,
                                           std::vector<Node<2>*> nodes,
                                           std::vector<VertexElement<2, 2>*> vertexElements,
                                           double cellRearrangementThreshold,
                                           double t2Threshold)
    : MutableVertexMesh<2,2>(nodes, vertexElements, cellRearrangementThreshold, t2Threshold),
      mWidth(width),
      mHeight(height),
      mpMeshForVtk(nullptr)
{
    // Call ReMesh() to remove any deleted nodes and relabel
    ReMesh();
}

Toroidal2dVertexMesh::Toroidal2dVertexMesh(Toroidal2dMesh& rMesh, bool isBounded)
    : mWidth(rMesh.GetWidth(0)),
      mHeight(rMesh.GetWidth(1)),
      mpMeshForVtk(nullptr)
{
    mpDelaunayMesh = &rMesh;

    // Reset member variables and clear mNodes, mFaces and mElements
    Clear();

    if (!isBounded)
    {
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

        // Loop over all generated nodes and check they're not outside [0,mWidth]x[0,mHeight]
        for (unsigned i=0; i<num_nodes; i++)
        {
            CheckNodeLocation(mNodes[i]);
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
    }
    else // Is Bounded
    {
        // First create an extended mesh to include points extended from the boundary
        std::vector<Node<2> *> nodes;
        for (typename TetrahedralMesh<2,2>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
            node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
            ++node_iter)
        {
            nodes.push_back(new Node<2>(node_iter->GetIndex(), node_iter->rGetLocation(),node_iter->IsBoundaryNode()));
        }

        // Add new nodes
        unsigned new_node_index = mpDelaunayMesh->GetNumNodes();
        for (TetrahedralMesh<2,2>::ElementIterator elem_iter = mpDelaunayMesh->GetElementIteratorBegin();
            elem_iter != mpDelaunayMesh->GetElementIteratorEnd();
            ++elem_iter)
        {
            bool bad_element = false;
            double edge_threshold = 1.5; //TODO Make settable variable!

            for (unsigned j=0; j<3; j++)
            {
                Node<2>* p_node_a = mpDelaunayMesh->GetNode(elem_iter->GetNodeGlobalIndex(j));
                Node<2>* p_node_b = mpDelaunayMesh->GetNode(elem_iter->GetNodeGlobalIndex((j+1)%3));
                if (norm_2(mpDelaunayMesh->GetVectorFromAtoB(p_node_a->rGetLocation(), p_node_b->rGetLocation()))>edge_threshold)
                {
                    bad_element = true;
                    break;
                }
            }

            if (bad_element)
            {
                for (unsigned j=0; j<3; j++)
                {
                    Node<2>* p_node_a = mpDelaunayMesh->GetNode(elem_iter->GetNodeGlobalIndex(j));
                    Node<2>* p_node_b = mpDelaunayMesh->GetNode(elem_iter->GetNodeGlobalIndex((j+1)%3));

                    c_vector<double,2> edge = mpDelaunayMesh->GetVectorFromAtoB(p_node_a->rGetLocation(), p_node_b->rGetLocation());
                    double edge_length = norm_2(edge);
                    
                    // The short edges in these elements are the boundaries of the void
                    if (edge_length<edge_threshold)
                    {
                        // Short Edge so add new node
                        c_vector<double,2> normal_vector;

                        // Outward Normal
                        normal_vector[0]= edge[1];
                        normal_vector[1]= -edge[0];

                        double dij = norm_2(normal_vector);
                        assert(dij>1e-5); //Sanity check
                        normal_vector /= dij;

                        double bound_offset = 0.5; //TODO Make settable variable!
                        c_vector<double,2> new_node_location = -bound_offset*normal_vector + p_node_a->rGetLocation() + 0.5*edge;

                        nodes.push_back(new Node<2>(new_node_index, new_node_location));
                        new_node_index++;

                        // Now add extra end nodes if appropriate.
                        unsigned num_sections = 1; // I.e. only at ends.
                        for (unsigned section=0; section <= num_sections; section++)
                        {
                            double ratio = (double)section/(double)num_sections;
                            c_vector<double,2> new_node_location = -normal_vector + ratio*p_node_a->rGetLocation() + (1-ratio)*p_node_b->rGetLocation();

                            //Check if near other nodes (could be inefficient)
                            bool node_clear = true;
                            double node_clearance = 0.3; // Make settable variable??

                            for (unsigned i=0; i<nodes.size(); i++)
                            {
                                double distance = norm_2(mpDelaunayMesh->GetVectorFromAtoB(nodes[i]->rGetLocation(), new_node_location));
                                if (distance < node_clearance)
                                {
                                    node_clear = false;
                                    break;
                                }
                            }

                            if (node_clear)
                            {
                                nodes.push_back(new Node<2>(new_node_index, new_node_location));
                                new_node_index++;
                            }
                        }
                    }
                }
            }
        }

        // Loop over all nodes and check they're not outside [0,mWidth]x[0,mHeight]
        for (unsigned i=0; i<nodes.size(); i++)
        {
            CheckNodeLocation(nodes[i]);
        }

        Toroidal2dMesh extended_mesh(mpDelaunayMesh->GetWidth(0),mpDelaunayMesh->GetWidth(1), nodes);

        unsigned num_elements = mpDelaunayMesh->GetNumAllNodes();
        unsigned num_nodes = extended_mesh.GetNumAllElements();

        // Allocate memory for mNodes and mElements
        this->mNodes.reserve(num_nodes);

        // Create as many elements as there are nodes in the mesh
        mElements.reserve(num_elements);
        for (unsigned elem_index = 0; elem_index < num_elements; elem_index++)
        {
            VertexElement<2, 2>* p_element = new VertexElement<2, 2>(elem_index);
            mElements.push_back(p_element);
        }

        // Populate mNodes
        GenerateVerticesFromElementCircumcentres(extended_mesh);

        // Loop over all generated nodes and check they're not outside [0,mWidth]x[0,mHeight]
        for (unsigned i=0; i<num_nodes; i++)
        {
            CheckNodeLocation(mNodes[i]);
        }


        // Loop over elements of the Delaunay mesh (which are nodes/vertices of this mesh)
        for (unsigned i = 0; i < num_nodes; i++)
        {
            // Loop over nodes owned by this triangular element in the Delaunay mesh
            // Add this node/vertex to each of the 3 vertex elements
            for (unsigned local_index = 0; local_index < 3; local_index++)
            {
                unsigned elem_index = extended_mesh.GetElement(i)->GetNodeGlobalIndex(local_index);

                if (elem_index < num_elements)
                {
                    unsigned num_nodes_in_elem = mElements[elem_index]->GetNumNodes();
                    unsigned end_index = num_nodes_in_elem > 0 ? num_nodes_in_elem - 1 : 0;

                    mElements[elem_index]->AddNode(this->mNodes[i], end_index);
                }
            }
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

Toroidal2dVertexMesh::Toroidal2dVertexMesh()
    : mpMeshForVtk(nullptr)
{
}

Toroidal2dVertexMesh::~Toroidal2dVertexMesh()
{
    if (mpMeshForVtk != NULL)
    {
         delete mpMeshForVtk;
    }
}

c_vector<double, 2> Toroidal2dVertexMesh::GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
{
    assert(mWidth > 0.0);
    assert(mHeight > 0.0);

    c_vector<double, 2> vector = rLocation2 - rLocation1;
    vector[0] = fmod(vector[0], mWidth);
    vector[1] = fmod(vector[1], mHeight);

    // If the points are more than halfway across the domain, measure the other way
    if (vector[0] > 0.5*mWidth)
    {
        vector[0] -= mWidth;
    }
    else if (vector[0] < -0.5*mWidth)
    {
        vector[0] += mWidth;
    }

    // If the points are more than halfway up the domain, measure the other way
    if (vector[1] > 0.5*mHeight)
    {
        vector[1] -= mHeight;
    }
    else if (vector[1] < -0.5*mHeight)
    {
        vector[1] += mHeight;
    }
    return vector;
}

void Toroidal2dVertexMesh::SetNode(unsigned nodeIndex, ChastePoint<2> point)
{
    double x_coord = point.rGetLocation()[0];
    double y_coord = point.rGetLocation()[1];

    // Perform a periodic movement if necessary
    if (x_coord >= mWidth)
    {
        // Move point left
        point.SetCoordinate(0, x_coord - mWidth);
    }
    else if (x_coord < 0.0)
    {
        // Move point right
        point.SetCoordinate(0, x_coord + mWidth);
    }
    if (y_coord >= mHeight)
    {
        // Move point down
        point.SetCoordinate(1, y_coord - mHeight);
    }
    else if (y_coord < 0.0)
    {
        // Move point up
        point.SetCoordinate(1, y_coord + mHeight);
    }

    // Update the node's location
    MutableVertexMesh<2,2>::SetNode(nodeIndex, point);
}

double Toroidal2dVertexMesh::GetWidth(const unsigned& rDimension) const
{
    assert(rDimension==0 || rDimension==1);

    double width = mWidth;
    if (rDimension == 1)
    {
        width = mHeight;
    }

    return width;
}

void Toroidal2dVertexMesh::SetHeight(double height)
{
    assert(height > 0);
    mHeight = height;
}

void Toroidal2dVertexMesh::SetWidth(double width)
{
    assert(width > 0);
    mWidth = width;
}

unsigned Toroidal2dVertexMesh::AddNode(Node<2>* pNewNode)
{
    // If necessary move it to be back onto the torus
    CheckNodeLocation(pNewNode);

    unsigned node_index = MutableVertexMesh<2,2>::AddNode(pNewNode);

    return node_index;
}

void Toroidal2dVertexMesh::CheckNodeLocation(Node<2>* pNode)
{
    double x_location = pNode->rGetLocation()[0];
    if (x_location < 0)
    {
        pNode->rGetModifiableLocation()[0] = x_location + mWidth;
    }
    else if (x_location > mWidth)
    {
        pNode->rGetModifiableLocation()[0] = x_location - mWidth;
    }
    double y_location = pNode->rGetLocation()[1];
    if (y_location < 0)
    {
        pNode->rGetModifiableLocation()[1] = y_location + mHeight;
    }
    else if (y_location > mHeight)
    {
        pNode->rGetModifiableLocation()[1] = y_location - mHeight;
    }
}

VertexMesh<2, 2>* Toroidal2dVertexMesh::GetMeshForVtk()
{
    unsigned num_nodes = GetNumNodes();

    std::vector<Node<2>*> temp_nodes;
    std::vector<VertexElement<2, 2>*> elements;

    if(!mpDelaunayMesh) // No Delaunay mesh so less copies as all too top right
    {
        temp_nodes.resize(4 * num_nodes);

        // Create four copies of each node.
        for (unsigned index = 0; index < num_nodes; index++)
        {
            c_vector<double, 2> location;
            location = GetNode(index)->rGetLocation();

            // Node copy at original location
            Node<2>* p_node = new Node<2>(index, false, location[0], location[1]);
            temp_nodes[index] = p_node;

            // Node copy shifted right
            p_node = new Node<2>(num_nodes + index, false, location[0] + mWidth, location[1]);
            temp_nodes[num_nodes + index] = p_node;

            // Node copy shifted up
            p_node = new Node<2>(2*num_nodes + index, false, location[0], location[1] + mHeight);
            temp_nodes[2*num_nodes + index] = p_node;

            // Node copy shifted right and up
            p_node = new Node<2>(3*num_nodes + index, false, location[0] + mWidth, location[1] + mHeight);
            temp_nodes[3*num_nodes + index] = p_node;
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
            bool element_straddles_top_bottom_boundary = false;

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
                if (fabs(vector[1]) > 0.5*mHeight)
                {
                    element_straddles_top_bottom_boundary = true;
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

                // If the element straddles the top/bottom periodic boundary...
                if (element_straddles_top_bottom_boundary)
                {
                    // ...and this node is located below the centre of the mesh...
                    bool node_is_above_centre = (elem_iter->GetNode(local_index)->rGetLocation()[1] - 0.5*mHeight > 0);
                    if (!node_is_above_centre)
                    {
                        // ...then choose the equivalent node above
                        this_node_index += 2*num_nodes;
                    }
                }

                elem_nodes.push_back(temp_nodes[this_node_index]);
            }

            VertexElement<2,2>* p_element = new VertexElement<2,2>(elem_index, elem_nodes);
            elements.push_back(p_element);
        }
    }
    else // Has Delaunay mesh so match elements to centres
    {
        temp_nodes.resize(9*num_nodes);

        // Create nine copies of each node.
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

            // Node copy shifted up and right
            p_node = new Node<2>(2*num_nodes + index, false, location[0] + mWidth, location[1] + mHeight);
            temp_nodes[2*num_nodes + index] = p_node;

            // Node copy shifted up
            p_node = new Node<2>(3*num_nodes + index, false, location[0], location[1] + mHeight);
            temp_nodes[3*num_nodes + index] = p_node;

            // Node copy shifted up and left
            p_node = new Node<2>(4*num_nodes + index, false, location[0] - mWidth, location[1] + mHeight);
            temp_nodes[4*num_nodes + index] = p_node;

            // Node copy shifted left
            p_node = new Node<2>(5*num_nodes + index, false, location[0] - mWidth, location[1]);
            temp_nodes[5*num_nodes + index] = p_node;

            // Node copy shifted left and down
            p_node = new Node<2>(6*num_nodes + index, false, location[0] - mWidth, location[1] - mHeight);
            temp_nodes[6*num_nodes + index] = p_node;

            // Node copy shifted down
            p_node = new Node<2>(7*num_nodes + index, false, location[0], location[1] - mHeight);
            temp_nodes[7*num_nodes + index] = p_node;

            // Node copy shifted down and right
            p_node = new Node<2>(8*num_nodes + index, false, location[0] + mWidth, location[1] - mHeight);
            temp_nodes[8*num_nodes + index] = p_node;
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
            bool element_straddles_top_bottom_boundary = false;

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
                if (fabs(vector[1]) > 0.5*mHeight)
                {
                    element_straddles_top_bottom_boundary = true;
                }
            }
            /* If this is a voronoi tesselation make sure the elements contain
            * the original Delaunay node
            */
            bool element_centre_on_right = true;
            bool element_centre_on_top = true;

            unsigned delaunay_index = this->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);
            double element_centre_x_location = this->mpDelaunayMesh->GetNode(delaunay_index)->rGetLocation()[0];
            double element_centre_y_location = this->mpDelaunayMesh->GetNode(delaunay_index)->rGetLocation()[1];

            if (element_centre_x_location < 0.5*mWidth)
            {
                element_centre_on_right = false;
            }
            if (element_centre_y_location < 0.5*mHeight)
            {
                element_centre_on_top = false;
            }

            // Use the above information when duplicating the element
            for (unsigned local_index=0; local_index<num_nodes_in_elem; local_index++)
            {
                unsigned this_node_index = elem_iter->GetNodeGlobalIndex(local_index);

                // If the element straddles the left/right periodic boundary...
                if (element_straddles_left_right_boundary && !element_straddles_top_bottom_boundary)
                {
                    // ...and this node is located to the left of the centre of the mesh...
                    bool node_is_right_of_centre = (elem_iter->GetNode(local_index)->rGetLocation()[0] - 0.5*mWidth > 0);
                    if (!node_is_right_of_centre && element_centre_on_right)
                    {
                        // ...then choose the equivalent node to the right
                        this_node_index += num_nodes;
                    }
                    else if (node_is_right_of_centre && !element_centre_on_right)
                    {
                        // ...then choose the equivalent node to the left
                        this_node_index += 5*num_nodes;
                    }
                }
                else if (!element_straddles_left_right_boundary && element_straddles_top_bottom_boundary)
                {
                    // ...and this node is located to the bottom of the centre of the mesh...
                    bool node_is_top_of_centre = (elem_iter->GetNode(local_index)->rGetLocation()[1] - 0.5*mHeight > 0);
                    if (!node_is_top_of_centre && element_centre_on_top)
                    {
                        // ...then choose the equivalent node to the top
                        this_node_index += 3*num_nodes;
                    }
                    else if (node_is_top_of_centre && !element_centre_on_top)
                    {
                        // ...then choose the equivalent node to the bottom
                        this_node_index += 7*num_nodes;
                    }
                }
                else if (element_straddles_left_right_boundary && element_straddles_top_bottom_boundary)
                {
                    bool node_is_right_of_centre = (elem_iter->GetNode(local_index)->rGetLocation()[0] - 0.5*mWidth > 0);
                    bool node_is_top_of_centre = (elem_iter->GetNode(local_index)->rGetLocation()[1] - 0.5*mHeight > 0);

                    if (!node_is_top_of_centre && element_centre_on_top)
                    {
                        if (!node_is_right_of_centre && element_centre_on_right)
                        {
                            this_node_index += 2*num_nodes;
                        }
                        else if (node_is_right_of_centre && !element_centre_on_right)
                        {
                            this_node_index += 4*num_nodes;
                        }
                        else
                        {
                            this_node_index += 3*num_nodes;
                        }
                    }
                    else if (node_is_top_of_centre && !element_centre_on_top)
                    {
                        if (!node_is_right_of_centre && element_centre_on_right)
                        {
                            this_node_index += 8*num_nodes;
                        }
                        else if (node_is_right_of_centre && !element_centre_on_right)
                        {
                            this_node_index += 6*num_nodes;
                        }
                        else
                        {
                            this_node_index += 7*num_nodes;
                        }
                    }
                    else
                    {
                        if (!node_is_right_of_centre && element_centre_on_right)
                        {
                            this_node_index += num_nodes;
                        }
                        else if (node_is_right_of_centre && !element_centre_on_right)
                        {
                            this_node_index += 5*num_nodes;
                        }
                    }
                }

                // If the element straddles the top/bottom periodic boundary...
                // if (element_straddles_top_bottom_boundary)
                // {
                //     // ...and this node is located below the centre of the mesh...
                //     bool node_is_above_centre = (elem_iter->GetNode(local_index)->rGetLocation()[1] - 0.5*mHeight > 0);
                //     if (!node_is_above_centre)
                //     {
                //         // ...then choose the equivalent node above
                //         this_node_index += 2*num_nodes;
                //     }
                // }

                elem_nodes.push_back(temp_nodes[this_node_index]);
            }

            VertexElement<2,2>* p_element = new VertexElement<2,2>(elem_index, elem_nodes);
            elements.push_back(p_element);
        }
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

void Toroidal2dVertexMesh::ConstructFromMeshReader(AbstractMeshReader<2,2>& rMeshReader, double width, double height)
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
    for (unsigned node_idx = 0 ; node_idx < num_nodes ; node_idx++)
    {
        node_data = rMeshReader.GetNextNode();
        node_data.pop_back();
        this->mNodes.push_back(new Node<2>(node_idx, node_data, false));
    }

    rMeshReader.Reset();

    // Reserve memory for elements
    mElements.reserve(rMeshReader.GetNumElements());

    // Add elements
    for (unsigned elem_idx = 0 ; elem_idx < num_elements ; elem_idx++)
    {
        // Get the data for this element
        ElementData element_data = rMeshReader.GetNextElementData();

        // Get the nodes owned by this element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j=0; j<num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        // Use nodes and index to construct this element
        VertexElement<2,2>* p_element = new VertexElement<2,2>(elem_idx, nodes);
        mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = (unsigned) element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }

    /*
     * Set width and height from function arguments, and validate by checking area is correct
     */
    this->mWidth = width;
    this->mHeight = height;

    double total_surface_area = 0.0;
    for (unsigned elem_idx = 0 ; elem_idx < num_elements ; elem_idx++)
    {
        total_surface_area += this->GetVolumeOfElement(elem_idx);
    }

    if (fabs(mWidth * mHeight - total_surface_area) > 1e-6)
    {
        EXCEPTION("Mesh width and height do not match sheet surface area.");
    }

    // Set default parameter values
    this->mCellRearrangementRatio = 1.5;
    this->mCellRearrangementThreshold = 0.01;
    this->mT2Threshold = 0.001;
    this->mMeshChangesDuringSimulation = true;
    this->mpMeshForVtk = nullptr;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Toroidal2dVertexMesh)
