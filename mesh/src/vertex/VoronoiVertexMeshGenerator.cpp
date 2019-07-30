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

#include "VoronoiVertexMeshGenerator.hpp"

#if BOOST_VERSION >= 105200

using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;
typedef boost::polygon::point_data<int> boost_point;

VoronoiVertexMeshGenerator::VoronoiVertexMeshGenerator(unsigned numElementsX,
                                                       unsigned numElementsY,
                                                       unsigned numRelaxationSteps,
                                                       double elementTargetArea)
        : mpMesh(nullptr),
          mpTorMesh(nullptr),
          mNumElementsX(numElementsX),
          mNumElementsY(numElementsY),
          mNumRelaxationSteps(numRelaxationSteps),
          mElementTargetArea(elementTargetArea),
          mTol(1e-7),
          mMaxExpectedNumSidesPerPolygon(17)
{
    this->ValidateInputAndSetMembers();
    this->GenerateVoronoiMesh();
}

VoronoiVertexMeshGenerator::~VoronoiVertexMeshGenerator()
{
    if (mpMesh)
    {
        delete mpMesh;
    }
    if (mpTorMesh)
    {
        delete mpTorMesh;
    }
}

void VoronoiVertexMeshGenerator::GenerateVoronoiMesh()
{
    /**
     * The initial points will be randomly distributed in the box [0.0, mMultiplierInX] x [0.0, mMultiplierInY], which
     * is a subset of the box [0.0, 1.0] x [0.0, 1.0].  At least one of mMultiplierInX and mMultiplierInY will equal
     * 1.0, and this will depend whether more elements were requested in the x direction or y direction.
     */

    // Get initial seed locations
    std::vector<c_vector<double, 2> > seed_locations = this->GetInitialPointLocations();
    this->ValidateSeedLocations(seed_locations);

    // We now create the initial Voronoi tessellation. This method updates mpMesh.
    this->CreateVoronoiTessellation(seed_locations);

    /**
     * Next, we perform the relaxation steps. The points used as seeds in the new Voronoi tessellation are the
     * centroids of the elements which are currently in the mesh.
     */
    for (unsigned relaxation = 0 ; relaxation < mNumRelaxationSteps ; relaxation++)
    {
        seed_locations.clear();
        seed_locations = this->GetElementCentroidsFromMesh();

        this->CreateVoronoiTessellation(seed_locations);
    }

    // We need to modify the node locations to achieve the correct target average element area
    double scale_factor = double(mMaxNumElems) * sqrt(mElementTargetArea);
    for (unsigned node_idx = 0 ; node_idx < mpMesh->GetNumNodes() ; node_idx++)
    {
        c_vector<double, 2>& node_location = mpMesh->GetNode(node_idx)->rGetModifiableLocation();
        node_location *= scale_factor;
    }
}

MutableVertexMesh<2,2>* VoronoiVertexMeshGenerator::GetMesh()
{
    return mpMesh;
}

MutableVertexMesh<2,2>* VoronoiVertexMeshGenerator::GetMeshAfterReMesh()
{
    mpMesh->ReMesh();
    return mpMesh;
}

Toroidal2dVertexMesh* VoronoiVertexMeshGenerator::GetToroidalMesh()
{
    /*
     * METHOD OUTLINE:
     *
     * 1. Copy nodes and elements from mpMesh into a new MutableVertexMesh, so no data is shared between mpMesh and the
     *    nodes and elements we will be working with.  There are no available copy constructors, so this is done from
     *    scratch.
     *
     * 2. Identify which nodes on the boundary are congruent to each other.  This is done by recursively using the
     *    helper function CheckForCongruentNodes().
     *
     * 3. Create mpTorMesh by copying the subset of remaining nodes and all elements in the new MutableVertexMesh.
     *    Some elements have been modified to replace nodes by congruent partner nodes and some nodes have been deleted.
     */

    // The width and height of the mesh for periodicity purposes
    double width  = mNumElementsX * sqrt(mElementTargetArea);
    double height = mNumElementsY * sqrt(mElementTargetArea);

    // We need to construct new nodes and elements so we don't have mpTorMesh sharing data with mpMesh
    std::vector<Node<2>*> new_nodes(mpMesh->GetNumNodes());
    std::vector<VertexElement<2,2>*> new_elems(mpMesh->GetNumElements());

    // Copy nodes
    for (unsigned node_counter = 0 ; node_counter < mpMesh->GetNumNodes() ; node_counter++)
    {
        Node<2>* p_node_to_copy = mpMesh->GetNode(node_counter);

        // Get all the information about the node we are copying
        unsigned copy_index = p_node_to_copy->GetIndex();
        bool copy_is_boundary = p_node_to_copy->IsBoundaryNode();
        const c_vector<double, 2>& copy_location = p_node_to_copy->rGetLocation();

        // There should not be any 'gaps' in node numbering, but we will assert just to make sure
        assert(copy_index < mpMesh->GetNumNodes());

        // Create a new node and place it in index order. Every node in a periodic mesh is non-boundary.
        new_nodes[copy_index] = new Node<2>(copy_index, copy_location, copy_is_boundary);
    }

    // Copy elements
    for (unsigned elem_counter = 0; elem_counter < mpMesh->GetNumElements(); elem_counter++)
    {
        VertexElement<2,2>* p_elem_to_copy = mpMesh->GetElement(elem_counter);

        // Get the information relating to the element we are copying
        unsigned copy_index     = p_elem_to_copy->GetIndex();
        unsigned copy_num_nodes = p_elem_to_copy->GetNumNodes();

        // There should not be any 'gaps' in element numbering, but we will assert just to make sure
        assert(copy_index < mpMesh->GetNumElements());

        // The vertex element is created from a vector of nodes
        std::vector<Node<2>*> nodes_this_elem;

        // Loop through the nodes in p_elem_to_copy and add the corresponding nodes that we have already copied
        for (unsigned node_local_idx = 0 ; node_local_idx < copy_num_nodes ; node_local_idx++)
        {
            Node<2>* p_local_node = p_elem_to_copy->GetNode(node_local_idx);
            unsigned local_node_global_idx = p_local_node->GetIndex();
            nodes_this_elem.push_back(new_nodes[local_node_global_idx]);
        }

        // Create a new node and place it in index order
        new_elems[copy_index] = new VertexElement<2,2>(copy_index, nodes_this_elem);
    }

    // We can now create the mesh with new_elements and the subset of new_nodes
    MutableVertexMesh<2,2>* p_temp_mesh = new MutableVertexMesh<2,2>(new_nodes, new_elems);

    /*
     * Recursively associate congruent nodes.
     *
     * For each node currently on the boundary, we loop through all other boundary nodes and check against all eight
     * possible congruent locations.  If any one coincides, we replace the identified node with its congruent partner,
     * delete and remove the identified node from the mesh, and start again from scratch.
     *
     * If a node has no congruent partners, we simply mark it as non-boundary, as it is then a node that needs to appear
     * in the final toroidal mesh.
     *
     * This recursion is guaranteed to terminate as the number of boundary nodes decreases by precisely one each time
     * CheckForCongruentNodes() is called, and CheckForCongruentNodes() returns false if there are no boundary nodes.
     */
    bool re_check = true;

    while (re_check)
    {
        re_check = this->CheckForCongruentNodes(p_temp_mesh, width, height);
    }

    // We now copy the nodes and elements into a new toroidal mesh, and delete p_temp_mesh
    new_nodes.clear();
    new_nodes.resize(p_temp_mesh->GetNumNodes());

    new_elems.clear();
    new_elems.resize(p_temp_mesh->GetNumElements());

    // Copy nodes
    for (unsigned node_counter = 0 ; node_counter < p_temp_mesh->GetNumNodes() ; node_counter++)
    {
        Node<2>* p_node_to_copy = p_temp_mesh->GetNode(node_counter);

        // Get all the information about the node we are copying
        unsigned copy_index = p_node_to_copy->GetIndex();
        const c_vector<double, 2>& copy_location = p_node_to_copy->rGetLocation();

        // No nodes should be boundary nodes
        assert(!p_node_to_copy->IsBoundaryNode());

        // There should not be any 'gaps' in node numbering, but we will assert just to make sure
        assert(copy_index < p_temp_mesh->GetNumNodes());

        // Create a new node and place it in index order. Every node in a periodic mesh is non-boundary.
        new_nodes[copy_index] = new Node<2>(copy_index, copy_location, false);
    }

    // Copy elements
    for (unsigned elem_counter = 0; elem_counter < p_temp_mesh->GetNumElements(); elem_counter++)
    {
        VertexElement<2,2>* p_elem_to_copy = p_temp_mesh->GetElement(elem_counter);

        // Get the information relating to the element we are copying
        unsigned copy_index     = p_elem_to_copy->GetIndex();
        unsigned copy_num_nodes = p_elem_to_copy->GetNumNodes();

        // There should not be any 'gaps' in element numbering, but we will assert just to make sure
        assert(copy_index < p_temp_mesh->GetNumElements());

        // The vertex element is created from a vector of nodes
        std::vector<Node<2>*> nodes_this_elem;

        // Loop through the nodes in p_elem_to_copy and add the corresponding nodes that we have already copied
        for (unsigned node_local_idx = 0 ; node_local_idx < copy_num_nodes ; node_local_idx++)
        {
            Node<2>* p_local_node = p_elem_to_copy->GetNode(node_local_idx);

            unsigned local_node_global_idx = p_local_node->GetIndex();

            nodes_this_elem.push_back(new_nodes[local_node_global_idx]);
        }

        // Create a new node and place it in index order
        new_elems[copy_index] = new VertexElement<2,2>(copy_index, nodes_this_elem);
    }

    delete p_temp_mesh;

    /*
     * We now create the mesh with new_elements and new_nodes.  We immediately call ReMesh() to tidy up any short edges.
     *
     * We then reposition all nodes to be within the box [0, width]x[0, height] to make mesh look nicer when visualised.
     */
    mpTorMesh = new Toroidal2dVertexMesh(width, height, new_nodes, new_elems);
    mpTorMesh->ReMesh();

    c_vector<double, 2> min_x_y;
    min_x_y[0] = DBL_MAX;
    min_x_y[1] = DBL_MAX;

    // First loop is to calculate the correct offset, min_x_y
    for (unsigned node_idx = 0 ; node_idx < mpTorMesh->GetNumNodes() ; node_idx++)
    {
        const c_vector<double, 2>& r_this_node_location = mpTorMesh->GetNode(node_idx)->rGetLocation();

        if (r_this_node_location[0] < min_x_y[0])
        {
            min_x_y[0] = r_this_node_location[0];
        }
        if (r_this_node_location[1] < min_x_y[1])
        {
            min_x_y[1] = r_this_node_location[1];
        }
    }

    // Second loop applies the offset, min_x_y, to each node in the mesh
    for (unsigned node_idx = 0 ; node_idx < mpTorMesh->GetNumNodes() ; node_idx++)
    {
        mpTorMesh->GetNode(node_idx)->rGetModifiableLocation() -= min_x_y;
    }

    // Third loop is to reposition any nodes that are now outside the bounding rectangle
    for (unsigned node_idx = 0 ; node_idx < mpTorMesh->GetNumNodes() ; node_idx++)
    {
        if (mpTorMesh->GetNode(node_idx)->rGetLocation()[0] >= width)
        {
            mpTorMesh->GetNode(node_idx)->rGetModifiableLocation()[0] -= width;
        }
        if (mpTorMesh->GetNode(node_idx)->rGetLocation()[1] >= height)
        {
            mpTorMesh->GetNode(node_idx)->rGetModifiableLocation()[1] -= height;
        }
    }

    return mpTorMesh;
}

bool VoronoiVertexMeshGenerator::CheckForCongruentNodes(MutableVertexMesh<2,2>* pMesh, double width, double height)
{
    // First find all the current boundary nodes in pMesh
    std::vector<Node<2>*> boundary_nodes;
    for (unsigned node_idx = 0 ; node_idx < pMesh->GetNumNodes() ; node_idx++)
    {
        Node<2>* p_node = pMesh->GetNode(node_idx);

        if (p_node->IsBoundaryNode())
        {
            boundary_nodes.push_back(p_node);
        }
    }

    // If there are no boundary nodes, return false (we're done checking congruencies)
    if (boundary_nodes.empty())
    {
        return false;
    }

    // Otherwise, calculate the eight possible congruent locations for the current node
    Node<2>* p_node_a = *(boundary_nodes.begin());
    c_vector<double, 2> node_a_pos = p_node_a->rGetLocation();
    std::vector<c_vector<double,2> > congruent_locations(8, node_a_pos);

    congruent_locations[0][0] += width;

    congruent_locations[1][1] += height;

    congruent_locations[2][0] -= width;

    congruent_locations[3][1] -= height;

    congruent_locations[4][0] += width;
    congruent_locations[4][1] += height;

    congruent_locations[5][0] -= width;
    congruent_locations[5][1] += height;

    congruent_locations[6][0] -= width;
    congruent_locations[6][1] -= height;

    congruent_locations[7][0] += width;
    congruent_locations[7][1] -= height;

    // Loop over all other boundary nodes
    for (unsigned node_b_counter = 0; node_b_counter < boundary_nodes.size(); node_b_counter++)
    {
        // Get the index of the current boundary node, and the corresponding node from the mesh
        unsigned node_b_idx = boundary_nodes[node_b_counter]->GetIndex();
        Node<2>* p_mesh_node_b = pMesh->GetNode(node_b_idx);

        if (p_node_a == p_mesh_node_b)
        {
            continue;
        }

        c_vector<double, 2> node_b_pos = p_mesh_node_b->rGetLocation();

        // Loop over the eight possible congruent locations and check for coincidence of location
        for (unsigned pos = 0 ; pos < congruent_locations.size() ; pos++)
        {
            if (norm_inf(node_b_pos - congruent_locations[pos]) < mTol)
            {
                // Once we find a congruent location, we replace that node in all containing elements
                std::set<unsigned> containing_elems = p_mesh_node_b->rGetContainingElementIndices();

                assert(containing_elems.size() > 0);

                for (std::set<unsigned>::iterator it = containing_elems.begin() ; it != containing_elems.end() ; ++it)
                {
                    VertexElement<2,2>* p_this_elem = pMesh->GetElement(*it);
                    unsigned local_idx = p_this_elem->GetNodeLocalIndex(p_mesh_node_b->GetIndex());

                    assert(local_idx < UINT_MAX);

                    p_this_elem->AddNode(p_node_a, local_idx);
                    p_this_elem->DeleteNode(local_idx);
                }

                // Delete the node_b and remove it from the mesh, and return true to start checking again from scratch
                pMesh->DeleteNodePriorToReMesh(p_mesh_node_b->GetIndex());
                pMesh->RemoveDeletedNodes();
                return true;
            }
        }
    }

    /*
     * If we have checked p_node_a against all node_b candidates and have not found a congruent location, p_node_a can
     * be re-tagged as a non-boundary node, and so will not get checked again next time this function is called.
     */
    p_node_a->SetAsBoundaryNode(false);

    return true;
}

std::vector<c_vector<double, 2> > VoronoiVertexMeshGenerator::GetInitialPointLocations()
{
    // Create a vector which contains mTotalNumElements spaces
    std::vector<c_vector<double, 2> > seed_points(mTotalNumElements);

    RandomNumberGenerator* p_rand_gen = RandomNumberGenerator::Instance();

    // Create the correct number of suitably scaled random numbers
    for (unsigned point_idx = 0 ; point_idx < mTotalNumElements ; point_idx++)
    {
        seed_points[point_idx][0] = p_rand_gen->ranf() * mMultiplierInX;
        seed_points[point_idx][1] = p_rand_gen->ranf() * mMultiplierInY;
    }

    return seed_points;
}

std::vector<c_vector<double, 2> > VoronoiVertexMeshGenerator::GetElementCentroidsFromMesh()
{
    assert(mpMesh->GetNumElements() == mNumElementsX * mNumElementsY);

    std::vector<c_vector<double, 2> > element_centroids;

    // Loop over all elements in the mesh
    for (unsigned elem_idx = 0; elem_idx < mpMesh->GetNumElements(); elem_idx++)
    {
        // Get the current centroid of the element
        c_vector<double, 2> this_centroid = mpMesh->GetCentroidOfElement(elem_idx);

        /**
         * We cannot be certain the centroid is in the correct location: while the initial seed points are all located
         * in the correct bounding rectangle, after each Voronoi tessellation step, the centroids might move out of
         * this box.  We can correct for this by exploiting the periodicity which we have because of the 3x3 tiling.
         */

        // Account for possible wrap-around in the x-direction
        if (this_centroid[0] < 0.0)
        {
            this_centroid[0] += mMultiplierInX;
        }
        else if (this_centroid[0] > mMultiplierInX)
        {
            this_centroid[0] -= mMultiplierInX;
        }

        // Account for possible wrap-around in the y-direction
        if (this_centroid[1] < 0.0)
        {
            this_centroid[1] += mMultiplierInY;
        }
        else if (this_centroid[1] > mMultiplierInY)
        {
            this_centroid[1] -= mMultiplierInY;
        }

        element_centroids.push_back(this_centroid);
    }

    return element_centroids;
}

void VoronoiVertexMeshGenerator::CreateVoronoiTessellation(std::vector<c_vector<double, 2> >& rSeedLocations)
{
    // Clear the mesh nodes and elements, as they will be replaced in this method
    std::vector<Node<2>*> nodes;
    std::vector<VertexElement<2, 2>* > elements;

    /**
     * We need to take the vector of locations and snap each element to an integer gridpoint (the boost Voronoi library
     * requires integer points).
     *
     * We then form a 3x3 tessellation of this set of points to bypass any problems associated with boundary conditions.
     */

    // Datatype is boost::polygon::point_data<int>
    std::vector<boost_point> points;

    // Take the locations and map them to integers in a subset of [0, INT_MAX/2] x [0, INT_MAX/2]
    for (unsigned point_idx = 0 ; point_idx < rSeedLocations.size() ; point_idx++)
    {
        // Calculate the correct integer gridpoints
        int point_x = int( floor( (rSeedLocations[point_idx][0] * mSamplingMultiplier) + 0.5) );
        int point_y = int( floor( (rSeedLocations[point_idx][1] * mSamplingMultiplier) + 0.5) );

        points.push_back( boost_point(point_x, point_y) );
    }

    // Define offset vectors for the 3 by 3 tessellation
    std::vector<boost_point> offsets;
    offsets.push_back( boost_point(-mMaxIntX,  mMaxIntY) );
    offsets.push_back( boost_point(        0,  mMaxIntY) );
    offsets.push_back( boost_point( mMaxIntX,  mMaxIntY) );
    offsets.push_back( boost_point( mMaxIntX,         0) );
    offsets.push_back( boost_point( mMaxIntX, -mMaxIntY) );
    offsets.push_back( boost_point(        0, -mMaxIntY) );
    offsets.push_back( boost_point(-mMaxIntX, -mMaxIntY) );
    offsets.push_back( boost_point(-mMaxIntX,         0) );

    // Add all the points in the tessellation to the vector of boost points so in total there are 9 copies of each
    // seed location, suitably tiled.
    for (unsigned rep = 0; rep < offsets.size(); rep++)
    {
        boost_point offset = offsets[rep];
        for (unsigned point_idx = 0; point_idx < rSeedLocations.size(); point_idx++)
        {
            boost_point new_point = boost_point(points[point_idx].x() + offset.x(), points[point_idx].y() + offset.y());
            points.push_back(new_point);
        }
    }

    /**
     * Now we have the 3x3 tessellation, we can create the Voronoi diagram.
     *
     * Note: The boost Voronoi implementation refers to cells.  In Chaste language, we should think of 'cells' in the
     * remainder of this function as 'elements'.
     *
     * We then loop over the cells in the Voronoi diagram corresponding only to the points from the centre of the
     * 3x3 tessellation.  These cells will correspond to the elements in the MutableVertexMesh.  From these cells,
     * we find the nodes which will correspond to nodes in the MutableVertexMesh.
     *
     * We then loop over the cells again, this time only those corresponding to points not in the centre of the 3x3
     * tessellation.  This allows us to tag boundary nodes, by finding those nodes in outer voronoi cells that coincide
     * in location with the nodes we have already identified.
     */

    // Construct the Voronoi tessellation of these 9 x mTotalNumElements points
    voronoi_diagram<double> vd;
    construct_voronoi(points.begin(), points.end(), &vd);

    // Loop over the cells in the voronoi diagram to find nodes for our mesh
    for (voronoi_diagram<double>::const_cell_iterator it = vd.cells().begin();
         it != vd.cells().end();
         ++it)
    {
        // Get a reference to the current cell
        const voronoi_diagram<double>::cell_type& cell = *it;

        // The cells we care about are exactly those whose source_index is less than the size of the locations vector
        // (i.e. those cells in the central portion of the 3x3 tessellation of source points)
        if (cell.source_index() < rSeedLocations.size())
        {
            // We create a vector of nodes, which will be used to create a MutableElement
            std::vector<Node<2>*> nodes_this_elem;

            // Loop over the edges of the current cell
            const voronoi_diagram<double>::edge_type *edge = cell.incident_edge();

            do
            {
                if (edge->is_primary())
                {
                    // Calculate the location corresponding to the location of vertex0 of the current edge
                    double x_location = (edge->vertex0()->x()) / mSamplingMultiplier;
                    double y_location = (edge->vertex0()->y()) / mSamplingMultiplier;

                    // Create a node at this location.  Default to non-boundary node; this will be updated later
                    Node<2>* p_this_node = new Node<2>(nodes.size(), false, x_location, y_location);

                    /**
                     * Check whether this node has already been created - all non-boundary nodes will be found twice
                     *
                     * When the Voronoi Tessellation has been created, there is no simple way to iterate over the Voronoi vertices
                     * without double-counting some.  We therefore loop through the vector of nodes and checks for equality (of position)
                     * with the new node (p_this_node).  If an 'equal' node is found, the following finds the index of that node.
                     */
                    unsigned existing_node_idx = UINT_MAX;

                    for (unsigned node_idx = 0; node_idx < nodes.size(); node_idx++)
                    {
                        // Grab the existing node location
                        const c_vector<double, 2>& r_existing_node_location = nodes[node_idx]->rGetLocation();

                        // Equality here is determined entirely on coincidence of position
                        if (fabs(r_existing_node_location[0] - p_this_node->rGetLocation()[0]) < DBL_EPSILON)
                        {
                            if (fabs(r_existing_node_location[1] - p_this_node->rGetLocation()[1]) < DBL_EPSILON)
                            {
                                // If the nodes match, return the existing node index
                                existing_node_idx = node_idx;
                            }
                        }
                    }

                    if (existing_node_idx < UINT_MAX)
                    {
                        // The node was already in nodes vector, and its index is the variable 'existing_node'
                        nodes_this_elem.push_back(nodes[existing_node_idx]);
                        delete p_this_node;
                    }
                    else
                    {
                        // The node does not yet exist - we add it to the nodes vector
                        nodes.push_back(p_this_node);
                        nodes_this_elem.push_back(p_this_node);
                    }

                    // Move to the next edge
                    edge = edge->next();
                }

            } while (edge != cell.incident_edge());

            // Add a new VertexElement to the mElements vector
            elements.push_back(new VertexElement<2,2>(elements.size(), nodes_this_elem));
        }
    }

    // Loop over the cells in the voronoi diagram to identify boundary nodes
    for (voronoi_diagram<double>::const_cell_iterator it = vd.cells().begin();
         it != vd.cells().end();
         ++it)
    {
        // Get a reference to the current cell
        const voronoi_diagram<double>::cell_type& cell = *it;

        // The cells we care about are exactly those whose source_index is greater than the size of the locations vector
        // (i.e. those cells in the outside eight portions of the 3x3 tessellation of source points)
        if (cell.source_index() >= rSeedLocations.size())
        {
            // Loop over the edges of the current cell
            const voronoi_diagram<double>::edge_type *edge = cell.incident_edge();

            do
            {
                /*
                 * Break out of the do-while if there is an infinite edge; we needn't care about cells at the very edge.
                 */
                if (edge->is_infinite())
                {
                    break;
                }

                if (edge->is_primary())
                {
                    c_vector<double, 2> vertex_location;

                    // Get the location of vertex0 of the current edge
                    vertex_location[0] = (edge->vertex0()->x()) / mSamplingMultiplier;
                    vertex_location[1] = (edge->vertex0()->y()) / mSamplingMultiplier;

                    /*
                     * Check whether this location coincides with one of our nodes; if it does, it must be a boundary
                     * node
                     */
                    for (unsigned node_idx = 0 ; node_idx < nodes.size() ; node_idx++)
                    {
                        // Grab the existing node location
                        const c_vector<double, 2>& r_existing_node_location = nodes[node_idx]->rGetLocation();

                        // Equality here is determined entirely on coincidence of position
                        if (fabs(r_existing_node_location[0] - vertex_location[0]) < DBL_EPSILON)
                        {
                            if (fabs(r_existing_node_location[1] - vertex_location[1]) < DBL_EPSILON)
                            {
                                // If the locations match, tag the node as being on the boundary
                                nodes[node_idx]->SetAsBoundaryNode(true);
                            }
                        }
                    }
                    // Move to the next edge
                    edge = edge->next();
                }

            } while (edge != cell.incident_edge());
        }
    }

    // Create a new mesh with the current vector of nodes and elements
    if (mpMesh)
    {
        delete mpMesh;
    }
    mpMesh = new MutableVertexMesh<2,2>(nodes, elements);
}

void VoronoiVertexMeshGenerator::ValidateInputAndSetMembers()
{
    // Validate inputs
    if ((mNumElementsX < 2) || (mNumElementsY < 2))
    {
        EXCEPTION("Need at least 2 by 2 cells");
    }

    if (mElementTargetArea <= 0.0)
    {
        EXCEPTION("Specified target area must be strictly positive");
    }

    // The total number of elements requested
    mTotalNumElements = mNumElementsX * mNumElementsY;

    // Max of the numbers of elements across and up
    mMaxNumElems = std::max(mNumElementsX, mNumElementsY);

    // The multipliers necessary in x and y to ensure all seed points are in [0,1]x[0,1]
    mMultiplierInX = double(mNumElementsX) / double(mMaxNumElems);
    mMultiplierInY = double(mNumElementsY) / double(mMaxNumElems);

    // Floating point representation of INT_MAX/2
    mSamplingMultiplier = 0.5 * double(INT_MAX);

    // The max integer that a seed point could be mapped to when discretising
    mMaxIntX = int( floor( (mMultiplierInX * mSamplingMultiplier) + 0.5 ) );
    mMaxIntY = int( floor( (mMultiplierInY * mSamplingMultiplier) + 0.5 ) );
}

void VoronoiVertexMeshGenerator::ValidateSeedLocations(std::vector<c_vector<double, 2> >& rSeedLocations)
{
    unsigned num_seeds = rSeedLocations.size();

    /*
     * Seeds at least 1.0 / mSamplingMultiplier apart from each other
     * are acceptable; the 1.5 term below could just as well be any
     * number that is strictly greater than 1.0.
     */
    double safe_distance = 1.5 / mSamplingMultiplier;

    // If we find a seed that needs to move position, we will move it and start checking again from the beginning
    bool recheck = true;

    while (recheck)
    {
        // If no seeds needs moving (as is overwhelmingly likely), we do not need to check again
        recheck = false;

        // We check each seed against every other, without double counting
        for (unsigned seed_idx_one = 0; seed_idx_one < num_seeds; seed_idx_one++)
        {
            for (unsigned seed_idx_two = seed_idx_one + 1; seed_idx_two < num_seeds; seed_idx_two++)
            {
                // We get the distance between the two seeds currently being considered
                c_vector<double, 2> one_to_two = rSeedLocations[seed_idx_two] - rSeedLocations[seed_idx_one];
                double distance = norm_2(one_to_two);

                if (distance < safe_distance)
                {
                    // We will now need to re-check
                    recheck = true;

                    /*
                     * If the distance is non-zero, we will move the second seed away along the line joining the two
                     * seeds until it's at the safe distance. If the distance is somehow zero, we will just move the
                     * x-location of the second seed by the safe distance.
                     *
                     * In all cases, we also need to ensure the new location is within the boundary box.
                     */
                    if (distance > DBL_EPSILON)
                    {
                        double multiplier = safe_distance / distance;
                        rSeedLocations[seed_idx_two] += (one_to_two * multiplier);

                        rSeedLocations[seed_idx_two][0] = fmod(rSeedLocations[seed_idx_two][0] + mMultiplierInX, mMultiplierInX);
                        rSeedLocations[seed_idx_two][1] = fmod(rSeedLocations[seed_idx_two][1] + mMultiplierInX, mMultiplierInX);
                    }
                    else
                    {
                        rSeedLocations[seed_idx_two][0] += safe_distance;
                        rSeedLocations[seed_idx_two][0] = fmod(rSeedLocations[seed_idx_two][0], mMultiplierInX);
                    }
                }
            }

            // We will re-check now rather than finishing the outer loop first
            if (recheck)
            {
                break;
            }
        }
    }
}

std::vector<double> VoronoiVertexMeshGenerator::GetPolygonDistribution()
{
    assert(mpMesh != nullptr);

    // Number of elements in the mesh
    unsigned num_elems = mpMesh->GetNumElements();

    // Store the number of each class of polygons
    std::vector<unsigned> num_polygons(mMaxExpectedNumSidesPerPolygon - 2, 0);

    // Container to return the polygon distribution
    std::vector<double> polygon_dist;

    // Loop over elements in the mesh to get the number of each class of polygon
    for (unsigned elem_idx = 0; elem_idx < num_elems; elem_idx++)
    {
        unsigned num_nodes_this_elem = mpMesh->GetElement(elem_idx)->GetNumNodes();

        // All polygons are assumed to have 3, 4, 5, ..., mMaxExpectedNumSidesPerPolygon sides
        assert(num_nodes_this_elem > 2);
        assert(num_nodes_this_elem <= mMaxExpectedNumSidesPerPolygon);

        // Increment correct place in counter - triangles in place 0, squares in 1, etc
        num_polygons[num_nodes_this_elem - 3]++;
    }

    // Loop over the vector of polygon numbers and calculate the distribution vector to return
    unsigned elems_accounted_for = 0;
    for (unsigned polygon = 0 ; polygon < num_polygons.size() ; polygon++)
    {
        elems_accounted_for += num_polygons[polygon];

        polygon_dist.push_back(static_cast<double>(num_polygons[polygon]) / static_cast<double>(num_elems));

        // Only fill the vector of polygon distributions to the point where there are none of higher class in the mesh
        if (elems_accounted_for == num_elems)
        {
            break;
        }
    }

    return polygon_dist;
}

double VoronoiVertexMeshGenerator::GetAreaCoefficientOfVariation()
{
    assert(mpMesh != nullptr);

    // Number of elements in the mesh, and check there are at least two
    unsigned num_elems = mpMesh->GetNumElements();
    assert(num_elems > 1);

    double var = 0.0;

    // Loop over elements in the mesh to get the contributions to the variance
    for (unsigned elem_idx = 0 ; elem_idx < num_elems ; elem_idx++)
    {
        double deviation = mpMesh->GetVolumeOfElement(elem_idx) - mElementTargetArea;
        var += deviation * deviation;
    }

    var /= static_cast<double>(num_elems - 1);

    return sqrt(var) / mElementTargetArea;
}

void VoronoiVertexMeshGenerator::RefreshSeedsAndRegenerateMesh()
{
    this->GenerateVoronoiMesh();
}

void VoronoiVertexMeshGenerator::SetMaxExpectedNumSidesPerPolygon(unsigned maxExpectedNumSidesPerPolygon)
{
    mMaxExpectedNumSidesPerPolygon = maxExpectedNumSidesPerPolygon;
}

unsigned VoronoiVertexMeshGenerator::GetMaxExpectedNumSidesPerPolygon()
{
    return mMaxExpectedNumSidesPerPolygon;
}

#endif // BOOST_VERSION >= 105200
#if BOOST_VERSION < 105200

/**
 * This is a fake class to suppress coverage warnings. To get the real class
 * you must build with Boost version 1.52 or above.
 */
VoronoiVertexMeshGenerator::VoronoiVertexMeshGenerator()
{
    EXCEPTION("This is a dummy class. Build with Boost version 1.52 or above for functionality.");
}

#endif // BOOST_VERSION < 105200
