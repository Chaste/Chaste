/*

Copyright (c) 2005-2015, University of Oxford.
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
        : mpMesh(NULL),
          mpTorMesh(NULL),
          mNumElementsX(numElementsX),
          mNumElementsY(numElementsY),
          mNumRelaxationSteps(numRelaxationSteps),
          mElementTargetArea(elementTargetArea),
          mTol(1e-7)
{

    this->ValidateInputAndSetMembers();

    /**
     * The initial points will be randomly distributed in the box [0.0, mMultiplierInX] x [0.0, mMultiplierInY], which
     * is a subset of the box [0.0, 1.0] x [0.0, 1.0].  At least one of mMultiplierInX and mMultiplierInY will equal
     * 1.0, and this will depend whether more elements were requested in the x direction or y direction.
     */

    // Get initial seed locations
    std::vector<c_vector<double, 2> > seed_locations = this->GetInitialPointLocations();
    this->ValidateSeedLocations(seed_locations);

    // We now create the initial Voronoi tessellation.  This method updates mpMesh.
    this->CreateVoronoiTessellation(seed_locations);

    /**
     * Next, we perform the relaxation steps.  The points used as seeds in the new Voronoi tessellation are the
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

    // We reposition every node so that its x and y coordinates are >= 0.0. This helps tagging boundary nodes
    this->RepositionNodes();

    // Finally, we tag the boundary nodes to ensure the mesh properties are correct
    this->TagBoundaryNodes();

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
     * 1. Copy nodes and elements from mpMesh so no data is shared between the MutableVertexMesh and the
     *    ToroidalVertexMesh. There are no available copy constructors, so this is done from scratch.
     *
     * 2. Identify which nodes on the boundary are congruent to each other. Nodes on the left and bottom edges of the
     *    boundary, WLOG, can replace those on the top and right edges. Those we are discarding must be removed from
     *    all elements they are contained in, and replaced by their congruent partner.
     *
     * 3. Create mpTorMesh with the subset of the node we copied from mpMesh and all of the elements we copied from
     *    mpMesh, some of which have been modified to replace nodes by their congruent partners.
     */

    // This method should only ever be called after mpMesh has been created
    assert(mpMesh != NULL);

    // The width and height of the mesh for periodicity purposes
    double width = mNumElementsX * sqrt(mElementTargetArea);
    double hight = mNumElementsY * sqrt(mElementTargetArea);

    // We need to construct new nodes and elements so we don't have mpTorMesh sharing data with mpMesh
    std::vector<Node<2>*> new_nodes(mpMesh->GetNumNodes());
    std::vector<VertexElement<2,2>*> new_elems(mpMesh->GetNumElements());

    // Only boundary nodes need be identified with their congruent partners, so we will keep track of them separately
    std::vector<Node<2>*> boundary_nodes;

    // Copy nodes
    for (unsigned node_counter = 0 ; node_counter < mpMesh->GetNumNodes() ; node_counter++)
    {
        Node<2>* p_node_to_copy = mpMesh->GetNode(node_counter);

        // Get all the information about the node we are copying
        unsigned            copy_index       = p_node_to_copy->GetIndex();
        c_vector<double, 2> copy_location    = p_node_to_copy->rGetLocation();
        bool                copy_is_boundary = p_node_to_copy->IsBoundaryNode();

        // There should not be any 'gaps' in node numbering, but we will assert just to make sure
        assert(copy_index < mpMesh->GetNumNodes());

        // Create a new node and place it in index order. Every node in a periodic mesh is non-boundary.
        new_nodes[copy_index] = new Node<2>(copy_index, copy_location, false);

        // If the original node was boundary, we will keep a separate copy for use later
        if (copy_is_boundary)
        {
            boundary_nodes.push_back(new_nodes[copy_index]);
        }
    }

    // Copy elements
    for (unsigned elem_counter = 0 ; elem_counter < mpMesh->GetNumElements() ; elem_counter++)
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

    /*
     * We now need to identify congruent boundary nodes. This process in inverse to the process used to identify
     * boundary nodes. Instead of considering those boundary nodes with locations >= (width, hight), we identify those
     * with locations <= (width, hight).  Each one of these nodes will be congruent to at least one other boundary node,
     * and will replace each of those that it's congruent to in any elements containing those congruent nodes.
     *
     * The positions must be checked up to some tolerance, and we use mTol, the same value as when identifying
     * boundary nodes.
     */

    // We need to keep a track of which nodes will make it into the final mesh
    std::vector<bool> nodes_to_keep(new_nodes.size(), true);

    // Loop through boundary nodes to decide if we need to check for congruence
    for (unsigned node_a_idx = 0 ; node_a_idx < boundary_nodes.size() ; node_a_idx++)
    {
        Node<2>* p_node_a = boundary_nodes[node_a_idx];
        c_vector<double, 2> node_a_location = p_node_a->rGetLocation();

        // We consider only boundary nodes inside [0, width] x [0, hight]
        if ( (node_a_location[0] > width - mTol) || (node_a_location[1] > hight - mTol) )
        {
            continue;
        }

        // There are three possible congruent locations for each candidate boundary node
        c_vector<double, 2> congruent_location_1 = node_a_location;
        congruent_location_1[0] += width;

        c_vector<double, 2> congruent_location_2 = node_a_location;
        congruent_location_2[1] += hight;

        c_vector<double, 2> congruent_location_3 = node_a_location;
        congruent_location_3[0] += width;
        congruent_location_3[1] += hight;

        // Loop over all other boundary nodes and keep track of any in the same position as a congruent location
        std::vector<Node<2>*> congruent_nodes;

        for (unsigned node_b_idx = 0 ; node_b_idx < boundary_nodes.size() ; node_b_idx++)
        {
            if (node_a_idx == node_b_idx)
            {
                continue;
            }

            Node<2>* p_node_b = boundary_nodes[node_b_idx];
            c_vector<double, 2> node_b_location = p_node_b->rGetLocation();

            if (norm_2(congruent_location_3 - node_b_location) < mTol)
            {
                congruent_nodes.push_back(p_node_b);
            }
            else if (norm_2(congruent_location_2 - node_b_location) < mTol)
            {
                congruent_nodes.push_back(p_node_b);
            }
            else if (norm_2(congruent_location_1 - node_b_location) < mTol)
            {
                congruent_nodes.push_back(p_node_b);
            }
        }

        // Each of these node_a should be congruent to at least one node_b
        assert(congruent_nodes.size() > 0);

        // We now replace each node_b we found by node_a in all elements containing node_b
        for (unsigned idx = 0 ; idx < congruent_nodes.size() ; idx++)
        {
            Node<2>* p_congruent_node = congruent_nodes[idx];
            unsigned congruent_idx = p_congruent_node->GetIndex();

            // Tag this node as one not to include in the final mesh
            nodes_to_keep[congruent_idx] = false;

            // Identify which elements contain the current congruent node
            std::vector<VertexElement<2,2>*> containing_elems;

            for (unsigned elem_idx = 0 ; elem_idx < new_elems.size() ; elem_idx++)
            {
                VertexElement<2,2>* p_this_elem = new_elems[elem_idx];

                // GetNodeLocalIndex() returns UINT_MAX if the test-node is not in the element
                unsigned congruent_node_local_idx = p_this_elem->GetNodeLocalIndex(congruent_idx);

                if (congruent_node_local_idx < UINT_MAX)
                {
                    containing_elems.push_back(p_this_elem);
                }
            }

            // Loop over containing elements and replace node_b with node_a
            for (unsigned containing_elem_idx = 0 ; containing_elem_idx < containing_elems.size() ; containing_elem_idx++)
            {
                VertexElement<2,2>* p_this_elem = containing_elems[containing_elem_idx];
                unsigned local_idx = p_this_elem->GetNodeLocalIndex(congruent_idx);

                assert(local_idx < UINT_MAX);

                p_this_elem->AddNode(p_node_a, local_idx);
                p_this_elem->DeleteNode(local_idx);
            }
        }
    }

    /*
     * We recreate a vector of nodes that only includes those that have not been deleted above.
     *
     * There may also be interior nodes still outside the box [0, width] x [0, hight], so we move those to their
     * congruent location.
     */
    std::vector<Node<2>* > nodes_for_mesh;

    for (unsigned node_counter = 0 ; node_counter < new_nodes.size() ; node_counter++)
    {
        Node<2>* p_this_node = new_nodes[node_counter];

        if ( nodes_to_keep[node_counter] )
        {
            if (p_this_node->rGetLocation()[0] > width)
            {
                p_this_node->rGetModifiableLocation()[0] -= width;
            }

            if (p_this_node->rGetLocation()[1] > hight)
            {
                p_this_node->rGetModifiableLocation()[1] -= hight;
            }

            nodes_for_mesh.push_back(p_this_node);
            nodes_for_mesh.back()->SetIndex(nodes_for_mesh.size() - 1);
        }
        else
        {
            delete p_this_node;
        }
    }

    // We can now create the mesh with new_elements and the subset of new_nodes
    mpTorMesh = new Toroidal2dVertexMesh(width, hight, nodes_for_mesh, new_elems);

    return mpTorMesh;
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
    for (unsigned elem_idx = 0 ; elem_idx < mpMesh->GetNumElements() ; elem_idx++)
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
    for (unsigned rep = 0 ; rep < offsets.size() ; rep++)
    {
        boost_point offset = offsets[rep];
        for (unsigned point_idx = 0 ; point_idx < rSeedLocations.size() ; point_idx++)
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
     */

    // Construct the Voronoi tessellation of these 9 x mTotalNumElements points
    voronoi_diagram<double> vd;
    construct_voronoi(points.begin(), points.end(), &vd);

    // Loop over the cells in the voronoi diagram
    for (voronoi_diagram<double>::const_cell_iterator it = vd.cells().begin();
         it != vd.cells().end();
         ++it)
    {
        // Get a reference to the current cell
        const voronoi_diagram<double>::cell_type& cell = *it;

        // The cells we care about are exactly those whose source_index is less than the size of the locations vector
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

                    for (unsigned node_idx = 0 ; node_idx < nodes.size() ; node_idx++)
                    {
                        // Grab the existing node location
                        c_vector<double, 2> existing_node_location = nodes[node_idx]->rGetLocation();

                        // Equality here is determined entirely on coincidence of position
                        if ( fabs(existing_node_location[0] - p_this_node->rGetLocation()[0]) < DBL_EPSILON )
                        {
                            if ( fabs(existing_node_location[1] - p_this_node->rGetLocation()[1]) < DBL_EPSILON )
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

    // Create a new mesh with the current vector of nodes and elements
    if (mpMesh)
    {
        delete mpMesh;
    }
    mpMesh = new MutableVertexMesh<2,2>(nodes, elements);
}

void VoronoiVertexMeshGenerator::TagBoundaryNodes()
{
    /*
     * The nodes have already been repositioned to have non-negative x and y coords.  This means that each boundary
     * node is now congruent to at least one boundary node that has x-coord or y-coord greater than the mesh width
     * or height.
     */
    double width = mNumElementsX * sqrt(mElementTargetArea);
    double hight = mNumElementsY * sqrt(mElementTargetArea);

    for (unsigned node_a_idx = 0 ; node_a_idx < mpMesh->GetNumNodes() ; node_a_idx++)
    {
        Node<2>* p_node_a = mpMesh->GetNode(node_a_idx);
        c_vector<double, 2> node_a_location = p_node_a->rGetLocation();

        if ( (node_a_location[0] < width - mTol) && (node_a_location[1] < hight - mTol) )
        {
            // Most nodes will be within the width and height
            continue;
        }

        // There are three possible congruent locations for each candidate boundary node
        c_vector<double, 2> congruent_location_1 = node_a_location;
        congruent_location_1[0] -= width;

        c_vector<double, 2> congruent_location_2 = node_a_location;
        congruent_location_2[1] -= hight;

        c_vector<double, 2> congruent_location_3 = node_a_location;
        congruent_location_3[0] -= width;
        congruent_location_3[1] -= hight;

        /*
         * Loop over all other nodes and check if any have the same location up to mTol
         */
        std::vector<Node<2>*> congruent_nodes;

        for (unsigned node_b_idx = 0 ; node_b_idx < mpMesh->GetNumNodes() ; node_b_idx++)
        {
            if (node_a_idx == node_b_idx)
            {
                continue;
            }

            Node<2>* p_node_b = mpMesh->GetNode(node_b_idx);
            c_vector<double, 2> node_b_location = p_node_b->rGetLocation();

            if (norm_2(congruent_location_1 - node_b_location) < mTol)
            {
                congruent_nodes.push_back(p_node_b);
            }
            else if (norm_2(congruent_location_2 - node_b_location) < mTol)
            {
                congruent_nodes.push_back(p_node_b);
            }
            else if (norm_2(congruent_location_3 - node_b_location) < mTol)
            {
                congruent_nodes.push_back(p_node_b);
            }
        }

        // If there is at least one congruent location, we must nodes as being on the boundary
        if (congruent_nodes.size() > 0)
        {
            p_node_a->SetAsBoundaryNode(true);

            for (unsigned idx = 0 ; idx < congruent_nodes.size() ; idx++)
            {
                congruent_nodes[idx]->SetAsBoundaryNode(true);
            }
        }
    }
}

void VoronoiVertexMeshGenerator::ValidateInputAndSetMembers()
{
    // Validate inputs
    if ( (mNumElementsX < 2) || (mNumElementsY < 2) )
    {
        EXCEPTION("Need at least 2 by 2 cells");
    }

    if ( mElementTargetArea <= 0.0 )
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

    // Seeds at least 1.0 / mSamplingMultiplier are acceptable, but we use 1.5 to be absolutely safe
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

void VoronoiVertexMeshGenerator::RepositionNodes()
{
    assert(mpMesh != NULL);
    assert(mpMesh->GetNumNodes() > 0);

    c_vector<double, 2> min_x_y;
    min_x_y[0] = DBL_MAX;
    min_x_y[1] = DBL_MAX;

    // First loop is to calculate the correct offset, min_x_y
    for (unsigned node_idx = 0 ; node_idx < mpMesh->GetNumNodes() ; node_idx++)
    {
        c_vector<double, 2> this_node_location = mpMesh->GetNode(node_idx)->rGetLocation();

        if(this_node_location[0] < min_x_y[0])
        {
            min_x_y[0] = this_node_location[0];
        }
        if(this_node_location[1] < min_x_y[1])
        {
            min_x_y[1] = this_node_location[1];
        }
    }

    // Second loop applies the offset, min_x_y, to each node in the mesh
    for (unsigned node_idx = 0 ; node_idx < mpMesh->GetNumNodes() ; node_idx++)
    {
        mpMesh->GetNode(node_idx)->rGetModifiableLocation() -= min_x_y;
    }
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
