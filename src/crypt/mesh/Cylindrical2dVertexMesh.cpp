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
#include "Cylindrical2dVertexMesh.hpp"

Cylindrical2dVertexMesh::Cylindrical2dVertexMesh(double width,
                                                 std::vector<Node<2>*> nodes,
                                                 std::vector<VertexElement<2, 2>*> vertexElements,
                                                 double cellRearrangementThreshold,
                                                 double edgeDivisionThreshold,
                                                 double t2Threshold)
    : MutableVertexMesh<2,2>(nodes, vertexElements, cellRearrangementThreshold, edgeDivisionThreshold, t2Threshold),
      mWidth(width)
{
    // ReMesh to remove any deleted nodes and relabel
    ReMesh();
}

Cylindrical2dVertexMesh::~Cylindrical2dVertexMesh()
{
}

c_vector<double, 2> Cylindrical2dVertexMesh::GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
{
    assert(mWidth > 0.0);

    c_vector<double, 2> vector = rLocation2 - rLocation1;
    vector[0] = fmod(vector[0], mWidth);

    // If the points are more than halfway around the cylinder apart, measure the other way
    if (vector[0] > mWidth/2.0)
    {
        vector[0] -= mWidth;
    }
    else if (vector[0] < -mWidth/2.0)
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

double Cylindrical2dVertexMesh::GetAreaOfElement(unsigned index)
{
    VertexElement<2, 2>* p_element = GetElement(index);

    c_vector<double, 2> first_node = p_element->GetNodeLocation(0);
    c_vector<double, 2> current_node;
    c_vector<double, 2> anticlockwise_node;
    c_vector<double, 2> transformed_current_node;
    c_vector<double, 2> transformed_anticlockwise_node;

    unsigned num_nodes_in_element = p_element->GetNumNodes();

    double element_area = 0;

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        // Find locations of current node and anticlockwise node
        current_node = p_element->GetNodeLocation(local_index);
        anticlockwise_node = p_element->GetNodeLocation((local_index+1)%num_nodes_in_element);

        /*
         * In order to calculate the area we map the origin to (x[0],y[0])
         * then use GetVectorFromAtoB() to get node cooordiantes
         */

        transformed_current_node = GetVectorFromAtoB(first_node, current_node);
        transformed_anticlockwise_node = GetVectorFromAtoB(first_node, anticlockwise_node);

        element_area += 0.5*(transformed_current_node[0]*transformed_anticlockwise_node[1]
                           - transformed_anticlockwise_node[0]*transformed_current_node[1]);
    }

    return element_area;
}

c_vector<double, 2> Cylindrical2dVertexMesh::GetCentroidOfElement(unsigned index)
{
    VertexElement<2, 2>* p_element = GetElement(index);

    c_vector<double, 2> centroid;
    c_vector<double, 2> transformed_centroid = zero_vector<double>(2);
    c_vector<double, 2> first_node = p_element->GetNodeLocation(0);
    c_vector<double, 2> current_node_location;
    c_vector<double, 2> next_node_location;
    c_vector<double, 2> transformed_current_node;
    c_vector<double, 2> transformed_anticlockwise_node;

    double temp_centroid_x = 0;
    double temp_centroid_y = 0;

    unsigned num_nodes_in_element = p_element->GetNumNodes();

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        // Find locations of current node and anticlockwise node
        current_node_location = p_element->GetNodeLocation(local_index);
        next_node_location = p_element->GetNodeLocation((local_index+1)%num_nodes_in_element);

        /*
         * In order to calculate the centroid we map the origin to (x[0],y[0])
         * then use  GetVectorFromAtoB() to get node cooordiantes
         */

        transformed_current_node = GetVectorFromAtoB(first_node, current_node_location);
        transformed_anticlockwise_node = GetVectorFromAtoB(first_node, next_node_location);

        temp_centroid_x += (transformed_current_node[0]+transformed_anticlockwise_node[0])*(transformed_current_node[0]*transformed_anticlockwise_node[1]-transformed_current_node[1]*transformed_anticlockwise_node[0]);
        temp_centroid_y += (transformed_current_node[1]+transformed_anticlockwise_node[1])*(transformed_current_node[0]*transformed_anticlockwise_node[1]-transformed_current_node[1]*transformed_anticlockwise_node[0]);
    }

    double vertex_area = GetAreaOfElement(index);
    double centroid_coefficient = 1.0/(6.0*vertex_area);

    transformed_centroid(0) = centroid_coefficient*temp_centroid_x;
    transformed_centroid(1) = centroid_coefficient*temp_centroid_y;

    centroid = transformed_centroid + first_node;

    return centroid;
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Cylindrical2dVertexMesh)
