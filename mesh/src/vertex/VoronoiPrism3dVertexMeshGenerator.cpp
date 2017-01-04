/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "VoronoiPrism3dVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshCustomFunctions.hpp"
#include <map>

#include "Debug.hpp"

#if BOOST_VERSION >= 105200

using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;
typedef boost::polygon::point_data<int> boost_point;

VoronoiPrism3dVertexMeshGenerator::VoronoiPrism3dVertexMeshGenerator(unsigned numElementsX,
                                                                     unsigned numElementsY,
                                                                     double elementHeightZ,
                                                                     unsigned numRelaxationSteps,
                                                                     double elementTargetApicalArea)
        : mpMesh(NULL),
          mNumElementsX(numElementsX),
          mNumElementsY(numElementsY),
          mElementHeightZ(elementHeightZ),
          mNumRelaxationSteps(numRelaxationSteps),
          mElementTargetApicalArea(elementTargetApicalArea),
          mTol(1e-7),
          mMaxExpectedNumSidesPerPolygon(17)
{
    this->ValidateInputAndSetMembers();
    this->GenerateVoronoiMesh();
}

VoronoiPrism3dVertexMeshGenerator::~VoronoiPrism3dVertexMeshGenerator()
{
    if (mpMesh)
    {
        delete mpMesh;
    }
}

void VoronoiPrism3dVertexMeshGenerator::GenerateVoronoiMesh()
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
    double scale_factor = double(mMaxNumElems) * sqrt(mElementTargetApicalArea);
    for (unsigned node_idx = 0 ; node_idx < mpMesh->GetNumNodes() ; node_idx++)
    {
        c_vector<double, 3>& node_location = mpMesh->GetNode(node_idx)->rGetModifiableLocation();
        node_location[0] *= scale_factor;
        node_location[1] *= scale_factor;
    }
}

MutableVertexMesh<3,3>* VoronoiPrism3dVertexMeshGenerator::GetMesh()
{
    return mpMesh;
}

MutableVertexMesh<3,3>* VoronoiPrism3dVertexMeshGenerator::GetMeshAfterReMesh()
{
    mpMesh->ReMesh();
    return mpMesh;
}

std::vector<c_vector<double, 2> > VoronoiPrism3dVertexMeshGenerator::GetInitialPointLocations()
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

std::vector<c_vector<double, 2> > VoronoiPrism3dVertexMeshGenerator::GetElementCentroidsFromMesh()
{
    // Return a vector of 2D centroids as Voronoi tessellation is implemented on a 2D surface.
    // for fully 3D purpose please change accordingly
    assert(mpMesh->GetNumElements() == mNumElementsX * mNumElementsY);

    std::vector<c_vector<double, 2> > element_centroids;

    // Loop over all elements in the mesh
    for (unsigned elem_idx = 0; elem_idx < mpMesh->GetNumElements(); elem_idx++)
    {
        // Get the current centroid of the element
        c_vector<double, 3> this_centroid_3d = mpMesh->GetCentroidOfElement(elem_idx);
        // As for Voronoi tessellation algorithm from boost, 2D seed is required.
        c_vector<double, 2> this_centroid;
        this_centroid[0] = this_centroid_3d[0];
        this_centroid[1] = this_centroid_3d[1];

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

void VoronoiPrism3dVertexMeshGenerator::CreateVoronoiTessellation(std::vector<c_vector<double, 2> >& rSeedLocations)
{
    // Clear the mesh nodes, faces and elements, as they will be replaced in this method
    std::vector<Node<3>*> lower_nodes;
    std::vector<Node<3>*> upper_nodes;
    std::vector<VertexElement<2, 3>*> faces;
    std::vector<VertexElement<3, 3>*> elements;

    // These two containers are used to recycle created lateral faces
    // (map is chosen over vector as the number of nodes is not known a priori)
    std::vector<bool> does_node_already_exist;
    std::map<unsigned, std::vector<unsigned> > node_to_lateral_face_indices;

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
            // We create two vectors of nodes, which will be used to create a MutableElement
            // two vectors are used so that we can save the time to figure out the indices for the pair of upper&lower nodes
            std::vector<Node<3>*> lower_nodes_this_elem;
            std::vector<Node<3>*> upper_nodes_this_elem;

            // Loop over the edges of the current cell
            const voronoi_diagram<double>::edge_type *edge = cell.incident_edge();

            do
            {
                if (edge->is_primary())
                {
                    // Calculate the location corresponding to the location of vertex0 of the current edge
                    double x_location = (edge->vertex0()->x()) / mSamplingMultiplier;
                    double y_location = (edge->vertex0()->y()) / mSamplingMultiplier;

                    // Create a node at this location. Default to non-boundary node; this will be updated later
                    Node<3>* p_this_lower_node = new Node<3>(lower_nodes.size(), false, x_location, y_location, 0);
                    // Attribute is added so that it can be identified in simulation (as basal and apical has different forces)
                    // 1.1 instead of 1.0 as it will be casted into unsigned for simpler comparison.
                    p_this_lower_node->AddNodeAttribute(1.1);

                    /**
                     * Check whether this node has already been created - all non-boundary nodes will be found twice
                     *
                     * When the Voronoi tessellation has been created, there is no simple way to iterate over the Voronoi vertices
                     * without double-counting some.  We therefore loop through the vector of nodes and checks for equality (of position)
                     * with the new node (p_this_node).  If an 'equal' node is found, the following finds the index of that node.
                     */
                    unsigned existing_node_idx = UINT_MAX;

                    for (unsigned node_idx = 0; node_idx < lower_nodes.size(); node_idx++)
                    {
                        // Grab the existing node location
                        const c_vector<double, 3>& r_existing_node_location = lower_nodes[node_idx]->rGetLocation();

                        // Equality here is determined entirely on coincidence of position
                        if ( fabs(r_existing_node_location[0] - p_this_lower_node->rGetLocation()[0]) < DBL_EPSILON )
                        {
                            if ( fabs(r_existing_node_location[1] - p_this_lower_node->rGetLocation()[1]) < DBL_EPSILON )
                            {
                                // If the nodes match, return the existing node index
                                existing_node_idx = node_idx;
                            }
                        }
                    }

                    if (existing_node_idx < UINT_MAX)
                    {
                        // The node was already in nodes vector, and its index is the variable 'existing_node'
                        lower_nodes_this_elem.push_back(lower_nodes[existing_node_idx]);
                        upper_nodes_this_elem.push_back(upper_nodes[existing_node_idx]);
                        does_node_already_exist[existing_node_idx] = true;
                        delete p_this_lower_node;
                    }
                    else
                    {
                        // The node does not yet exist - we add it to the nodes vector
                        lower_nodes.push_back(p_this_lower_node);
                        lower_nodes_this_elem.push_back(p_this_lower_node);
                        does_node_already_exist.push_back(false);

                        // Creating the corresponding upper node
                        ///\todo check this - apparently, 2 nodes having the same index do not cause any problem yet #2850
                        Node<3>* p_this_upper_node = new Node<3>(upper_nodes.size(), false, x_location, y_location, mElementHeightZ);
                        // Attribute is added so that it can be identified in simulation (as basal and apical has different forces)
                        // 2.1 instead of 2.0 as it will be casted into unsigned for simpler comparison.
                        p_this_upper_node->AddNodeAttribute(2.1);
                        upper_nodes.push_back(p_this_upper_node);
                        upper_nodes_this_elem.push_back(p_this_upper_node);
                    }

                    // Move to the next edge
                    edge = edge->next();
                }

            } while (edge != cell.incident_edge());

            /**
             * Create faces and element.
             * According to boost documentation, with `edge = edge->next()` we should obtain the nodes in counter-clockwise order.
             * Upper and lower faces will be created as usual.
             * As for lateral faces, to avoid constantly searching through the entire faces vector and also creating the
             * exact same lateral face twice, a map which has the index of the (lower) node as key and a vector of its lateral face
             * indices as value is used.
            **/

            // Initializing vectors which are required for the generation of the VertexElement<3, 3>
            unsigned numLowerNodesThisElement = lower_nodes_this_elem.size();
            std::vector<VertexElement<2, 3>*> faces_this_elem;
            std::vector<bool> faces_orientation;

            // Creating the lower face
            VertexElement<2, 3>* p_lower_face = new VertexElement<2, 3>(faces.size(), lower_nodes_this_elem);
            // Attribute is added so that it can be identified in simulation (as basal, apical and lateral faces have different contributions)
            // 1.1 instead of 1.0 as it will be casted into unsigned for simpler comparison.
            p_lower_face->AddElementAttribute(1.1);
            faces.push_back(p_lower_face);
            faces_this_elem.push_back(p_lower_face);
            faces_orientation.push_back(true);

            // Creating the upper face
            VertexElement<2,3>* p_upper_face = new VertexElement<2,3>(faces.size(), upper_nodes_this_elem);
            // Attribute is added so that it can be identified in simulation (as basal, apical and lateral faces have different contributions)
            // 2.1 instead of 2.0 as it will be casted into unsigned for simpler comparison.
            p_upper_face->AddElementAttribute(2.1);
            faces.push_back(p_upper_face);
            faces_this_elem.push_back(p_upper_face);
            faces_orientation.push_back(false);

            // Creating all the lateral faces in CCW
            for (unsigned localNodeIndex=0; localNodeIndex<numLowerNodesThisElement; ++localNodeIndex )
            {
                unsigned node1Index = lower_nodes_this_elem[localNodeIndex]->GetIndex();
                unsigned node2Index = lower_nodes_this_elem[(localNodeIndex+1) % numLowerNodesThisElement]->GetIndex();

                // The values of the maps are called here because they will be used both existing and creating branch.
                // They are called by reference as they will be modified if they enter creating branch.
                std::vector<unsigned>& r_face1_indices = node_to_lateral_face_indices[node1Index];
                std::vector<unsigned>& r_face2_indices = node_to_lateral_face_indices[node2Index];

                // If both nodes already exist, the lateral face MIGHT have been created.
                unsigned existing_face_index = UINT_MAX;

                if (does_node_already_exist[node1Index] && does_node_already_exist[node2Index])
                {
                    // Now need to search for the same lateral face index in both vector
                    // not a too complicated and resource intensive (as r_faces_index vectors have length of at most 3 or 4
                    // therefore not using existing function in <algorithm>
                    for (unsigned i1 = 0; (i1 < r_face1_indices.size() && existing_face_index==UINT_MAX); ++i1)
                    {
                        for (unsigned i2 = 0; (i2 < r_face2_indices.size() && existing_face_index==UINT_MAX); ++i2)
                        {
                            if (r_face1_indices[i1] == r_face2_indices[i2])
                            {
                                existing_face_index = r_face1_indices[i1];
                                break;
                            }
                        }
                    }

                    if (existing_face_index != UINT_MAX) // meaning it's found
                    {
                        faces_this_elem.push_back(faces[existing_face_index]);
                        // Face orientation is false as it was created by another element. CCW for another will be CW when
                        // viewing from the other side as rotation is pseudovectorial
                        faces_orientation.push_back(true);
                    }

                }
                if (existing_face_index == UINT_MAX)
                {
                    // Create new lateral rectangular face
                    std::vector<Node<3>*> nodes_of_lateral_face;
                    nodes_of_lateral_face.push_back(lower_nodes[node1Index]);
                    nodes_of_lateral_face.push_back(lower_nodes[node2Index]);
                    nodes_of_lateral_face.push_back(upper_nodes[node2Index]);
                    nodes_of_lateral_face.push_back(upper_nodes[node1Index]);

                    unsigned newFaceIndex = faces.size();
                    VertexElement<2, 3>* p_lateral_face = new VertexElement<2, 3>(newFaceIndex, nodes_of_lateral_face);
                    // Attribute is added so that it can be identified in simulation (as basal, apical and lateral faces have different contributions)
                    // 3.1 instead of 3.0 as it will be casted into unsigned for simpler comparison.
                    p_lateral_face->AddElementAttribute(3.1);
                    faces.push_back(p_lateral_face);
                    faces_this_elem.push_back(p_lateral_face);
                    faces_orientation.push_back(false);

                    // Update node_to_lateral_face_indices
                    r_face1_indices.push_back(newFaceIndex);
                    r_face2_indices.push_back(newFaceIndex);
                }
            }

            lower_nodes_this_elem.insert(lower_nodes_this_elem.end(), upper_nodes_this_elem.begin(), upper_nodes_this_elem.end());

            // Use the contructor with std::vector<Node> so that the output will be cleaner
            VertexElement<3, 3>* p_elem = new VertexElement<3, 3>(elements.size(), faces_this_elem, faces_orientation, lower_nodes_this_elem);
//            SetElementAsMonolayer(p_elem);
MARK
            elements.push_back( p_elem );
MARK
        }
    }

    // Loop over the cells in the Voronoi diagram to identify boundary nodes
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
                    for (unsigned node_idx = 0 ; node_idx < lower_nodes.size() ; node_idx++)
                    {
                        // Grab the existing node location
                        const c_vector<double, 3>& r_existing_node_location = lower_nodes[node_idx]->rGetLocation();

                        // Equality here is determined entirely on coincidence of position
                        if ( fabs(r_existing_node_location[0] - vertex_location[0]) < DBL_EPSILON )
                        {
                            if ( fabs(r_existing_node_location[1] - vertex_location[1]) < DBL_EPSILON )
                            {
                                // If the locations match, tag the node as being on the boundary
                                lower_nodes[node_idx]->SetAsBoundaryNode(true);
                                upper_nodes[node_idx]->SetAsBoundaryNode(true);
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
    // Combine the upper and lower nodes by adding upper_nodes into lower_nodes.
    // The index of upper nodes need to be modified.
    // unsigned lowerNodeLength = lower_nodes.size();
    unsigned upperNodeLength = upper_nodes.size();

    for (unsigned upperRunningIndex=0; upperRunningIndex<upperNodeLength; ++upperRunningIndex)
    {
        upper_nodes[upperRunningIndex]->SetIndex(lower_nodes.size());
        // Alternatively if the upper nodes start with 0 as well
        // upper_nodes[upperRunningIndex]->rGetModifiableIndex() += loweNodeLength;
        lower_nodes.push_back(upper_nodes[upperRunningIndex]);
    }

    mpMesh = new MutableVertexMesh<3,3>(lower_nodes, elements);
}

void VoronoiPrism3dVertexMeshGenerator::ValidateInputAndSetMembers()
{
    // Validate inputs
    if ( (mNumElementsX < 2) || (mNumElementsY < 2) )
    {
        EXCEPTION("Need at least 2 by 2 cells");
    }

    if ( mElementHeightZ <= 0.0 )
    {
        EXCEPTION("Specified element height must be strictly positive");
    }

    if ( mElementTargetApicalArea <= 0.0 )
    {
        EXCEPTION("Specified target apical area must be strictly positive");
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

void VoronoiPrism3dVertexMeshGenerator::ValidateSeedLocations(std::vector<c_vector<double, 2> >& rSeedLocations)
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

std::vector<double> VoronoiPrism3dVertexMeshGenerator::GetPolygonDistribution()
{
    assert(mpMesh != NULL);

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
        // and since there should be an upper node for every lower node, it should be even number
        assert(num_nodes_this_elem%2 == 0);     ///\todo check if the pairs is really lower and upper after simulation #2850
        assert(num_nodes_this_elem > 2);
        assert(num_nodes_this_elem/2 <= mMaxExpectedNumSidesPerPolygon);

        // Increment correct place in counter - triangles in place 0, squares in 1, etc
        num_polygons[num_nodes_this_elem/2 - 3]++;
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

double VoronoiPrism3dVertexMeshGenerator::GetApicalAreaCoefficientOfVariation()
{
    assert(mpMesh != NULL);

    // Number of elements in the mesh, and check there are at least two
    unsigned num_elems = mpMesh->GetNumElements();
    assert(num_elems > 1);

    double var = 0.0;

    // Loop over elements in the mesh to get the contributions to the variance
    for (unsigned elem_idx = 0 ; elem_idx < num_elems ; elem_idx++)
    {
        double deviation = mpMesh->GetVolumeOfElement(elem_idx)/mElementHeightZ - mElementTargetApicalArea;
        var += deviation * deviation;
    }

    var /= static_cast<double>(num_elems - 1);

    return sqrt(var) / mElementTargetApicalArea;
}

void VoronoiPrism3dVertexMeshGenerator::RefreshSeedsAndRegenerateMesh()
{
    this->GenerateVoronoiMesh();
}

void VoronoiPrism3dVertexMeshGenerator::SetMaxExpectedNumSidesPerPolygon(unsigned maxExpectedNumSidesPerPolygon)
{
    mMaxExpectedNumSidesPerPolygon = maxExpectedNumSidesPerPolygon;
}

unsigned VoronoiPrism3dVertexMeshGenerator::GetMaxExpectedNumSidesPerPolygon()
{
    return mMaxExpectedNumSidesPerPolygon;
}

#endif // BOOST_VERSION >= 105200
#if BOOST_VERSION < 105200

/**
 * This is a fake class to suppress coverage warnings. To get the real class
 * you must build with Boost version 1.52 or above.
 */
VoronoiPrism3dVertexMeshGenerator::VoronoiPrism3dVertexMeshGenerator()
{
    EXCEPTION("This is a dummy class. Build with Boost version 1.52 or above for functionality.");
}

#endif // BOOST_VERSION < 105200
