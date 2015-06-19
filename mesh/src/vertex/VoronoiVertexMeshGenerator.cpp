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
          mNumElementsX(numElementsX),
          mNumElementsY(numElementsY),
          mNumRelaxationSteps(numRelaxationSteps),
          mElementTargetArea(elementTargetArea)
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

    // Finally, we tag the boundary nodes to ensure the mesh properties are correct
    this->TagBoundaryNodes();

}

VoronoiVertexMeshGenerator::~VoronoiVertexMeshGenerator()
{
    if (mpMesh)
    {
        delete mpMesh;
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

            do // other instances of do - while in Chaste code
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
    for (unsigned node_idx = 0 ; node_idx < mpMesh->GetNumNodes() ; node_idx++)
    {
        Node<2>* p_node = mpMesh->GetNode(node_idx);

        unsigned num_containing_elements = p_node->GetNumContainingElements();

        if (num_containing_elements < 3)
        {
            p_node->SetAsBoundaryNode(true);
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

#endif // BOOST_VERSION >= 105200