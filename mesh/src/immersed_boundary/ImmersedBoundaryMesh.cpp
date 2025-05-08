/*

Copyright (c) 2005-2025, University of Oxford.
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

#include "ImmersedBoundaryMesh.hpp"

#include <algorithm>
#include <cmath>

#include "ImmersedBoundaryEnumerations.hpp"
#include "MeshUtilityFunctions.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"
#include "Warnings.hpp"


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                                                   std::vector<ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>*> elements,
                                                                   std::vector<ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>*> laminas,
                                                                   unsigned numGridPtsX,
                                                                   unsigned numGridPtsY)
        : mNumGridPtsX(numGridPtsX),
          mNumGridPtsY(numGridPtsY),
          mCharacteristicNodeSpacing(DOUBLE_UNSET),
          mElementDivisionSpacing(DOUBLE_UNSET),
          mNeighbourDist(0.1),
          mSummaryOfNodeLocations(0.0),
          mCellRearrangementThreshold(0.05)
{
    // Clear mNodes and mElements
    Clear();

    switch (SPACE_DIM)
    {
        case 2:
            m2dVelocityGrids.resize(extents[2][mNumGridPtsX][mNumGridPtsY]);
            break;

        case 3:
            EXCEPTION("Not implemented yet in 3D");
            break;

        //LCOV_EXCL_START
        default:
            NEVER_REACHED;
            break;
        //LCOV_EXCL_STOP
    }

    // Populate mNodes, mElements, and mLaminas
    for (unsigned node_it = 0; node_it < nodes.size(); node_it++)
    {
        Node<SPACE_DIM>* p_temp_node = nodes[node_it];
        this->mNodes.push_back(p_temp_node);
    }
    for (unsigned elem_it = 0; elem_it < elements.size(); elem_it++)
    {
        ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_temp_element = elements[elem_it];
        mElements.push_back(p_temp_element);
    }
    for (unsigned lam_it = 0; lam_it < laminas.size(); lam_it++)
    {
        ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>* p_temp_lamina = laminas[lam_it];
        mLaminas.push_back(p_temp_lamina);
    }

    // Register elements with nodes
    for (unsigned elem_it = 0; elem_it < mElements.size(); elem_it++)
    {
        ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_element = mElements[elem_it];

        unsigned element_index = p_element->GetIndex();
        unsigned num_nodes_in_element = p_element->GetNumNodes();

        for (unsigned node_idx = 0; node_idx < num_nodes_in_element; node_idx++)
        {
            p_element->GetNode(node_idx)->AddElement(element_index);
        }
    }

    // Register laminas with nodes
    //\todo is there a way we can register laminas with nodes?

    // Set characteristic node spacing to the average distance between nodes in elements
    double total_perimeter = 0.0;
    unsigned total_nodes = 0;
    for (unsigned elem_it = 0; elem_it < mElements.size(); elem_it++)
    {
        total_perimeter += this->GetSurfaceAreaOfElement(elem_it);
        total_nodes += mElements[elem_it]->GetNumNodes();
    }
    mCharacteristicNodeSpacing = total_perimeter / double(total_nodes);

    // Position fluid sources at the centroid of each element, and set strength to zero
    for (unsigned elem_it = 0; elem_it < elements.size(); elem_it++)
    {
        unsigned elem_idx = mElements[elem_it]->GetIndex();

        // Create a new fluid source at the correct location for each element
        unsigned source_idx = static_cast<unsigned>(mElementFluidSources.size());
        c_vector<double, SPACE_DIM> source_location = this->GetCentroidOfElement(elem_idx);
        mElementFluidSources.push_back(std::make_shared<FluidSource<SPACE_DIM>>(source_idx, source_location));

        // Set source parameters
        mElementFluidSources.back()->SetAssociatedElementIndex(elem_idx);
        mElementFluidSources.back()->SetStrength(0.0);

        // Associate source with element
        mElements[elem_it]->SetFluidSource(mElementFluidSources.back());
    }

    //Set up a number of sources to balance any active sources associated with elements
    double balancing_source_spacing = 2.0 * mCharacteristicNodeSpacing;

    // We start at the characteristic spacing in from the left-hand end, and place a source every 2 spacings
    double current_location = mCharacteristicNodeSpacing;

    while (current_location < 1.0)
    {
        // Create a new fluid source at the current x-location and zero y-location
        unsigned source_idx = static_cast<unsigned>(mBalancingFluidSources.size());
        mBalancingFluidSources.push_back(std::make_shared<FluidSource<SPACE_DIM>>(source_idx, current_location));

        mBalancingFluidSources.back()->SetStrength(0.0);

        // Increment the current location
        current_location += balancing_source_spacing;
    }

    // Calculate a default neighbour dist, as half the root of the average element volume
    constexpr double power = 1.0 / SPACE_DIM;
    const double total_volume_of_elems = std::accumulate(mElements.begin(), mElements.end(), 0.0,
                                                         [this](double d, ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* a)
                                                         {
                                                             return d + this->GetVolumeOfElement(a->GetIndex());
                                                         });
    mNeighbourDist = 0.5 * std::pow(total_volume_of_elems / mElements.size(), power);

    this->mMeshChangesDuringSimulation = true;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetElongationShapeFactorOfElement(
    [[maybe_unused]] unsigned index) // [[maybe_unused]] due to unused-but-set-parameter warning in GCC 7,8,9
{
    if constexpr (SPACE_DIM == 2)
    {
        c_vector<double, 3> moments = CalculateMomentsOfElement(index);

        double discriminant = sqrt((moments(0) - moments(1)) * (moments(0) - moments(1)) + 4.0 * moments(2) * moments(2));

        // Note that as the matrix of second moments of area is symmetric, both its eigenvalues are real
        double largest_eigenvalue = (moments(0) + moments(1) + discriminant) * 0.5;
        double smallest_eigenvalue = (moments(0) + moments(1) - discriminant) * 0.5;

        double elongation_shape_factor = sqrt(largest_eigenvalue / smallest_eigenvalue);
        return elongation_shape_factor;
    }
    else
    {
        NEVER_REACHED;
    }
}

bool CustomComparisonForVectorX(c_vector<double, 2> vecA, c_vector<double, 2> vecB)
{
    return vecA[0] < vecB[0];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetTortuosityOfMesh()
{
    if constexpr (SPACE_DIM == 2)
    {
        double total_length = 0.0;

        // Get the current elements
        std::vector<c_vector<double, 2> > centroids(mElements.size());
        for (unsigned elem_it = 0; elem_it < mElements.size(); elem_it++)
        {
            centroids[elem_it] = this->GetCentroidOfElement(mElements[elem_it]->GetIndex());
        }

        // Sort centroids by X
        std::sort(centroids.begin(), centroids.end(), CustomComparisonForVectorX);

        // Calculate piecewise linear length connecting centroids
        for (unsigned cent_it = 1; cent_it < centroids.size(); cent_it++)
        {
            total_length += norm_2(centroids[cent_it - 1] - centroids[cent_it]);
        }

        double straight_line_length = norm_2(centroids[0] - centroids[centroids.size() - 1]);

        // Avoid division by zero
        if (straight_line_length < 0.00001)
        {
            return 0;
        }

        return total_length / straight_line_length;
    }
    else
    {
        NEVER_REACHED;
    }
}

bool CustomComparisonForSkewnessMeasure(std::pair<unsigned, c_vector<double, 2> > pairA, std::pair<unsigned, c_vector<double, 2> > pairB)
{
    return pairA.second[0] < pairB.second[0];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetSkewnessOfElementMassDistributionAboutAxis(
    [[maybe_unused]] unsigned elemIndex, // [[maybe_unused]] due to unused-but-set-parameter warning in GCC 7,8,9
    c_vector<double, SPACE_DIM> axis)
{
    /*
     * Method outline:
     *
     * Given an arbitrary axis and a closed polygon, we calculate the skewness of the mass distribution of the polygon
     * perpendicular to the axis.  This is used as a measure of asymmetry.
     *
     * To simplify calculating the mass distribution, we translate the centroid of the element to the origin and rotate
     * about the centroid so the axis is vertical; then we sort all the nodes in ascending order of their x-coordinate.
     *
     * For each node in order, we need to know the length of the intersection through the node of the vertical line with
     * the polygon.  Once calculated, we have a piecewise-linear PDF for the mass distribution, which can be normalised
     * by the surface area of the polygon.
     *
     * By integrating the pdf directly, we can calculate exactly the necessary moments of the distribution needed for
     * the skewness.
     *
     * This method does not work for concave shapes, and emits a warning but continues anyway. Results will be close
     * to correct for mildly concave shapes or may be correct depending on the alignment of the element.
     */

    // This method only works in 2D
    if constexpr (SPACE_DIM == 2 && ELEMENT_DIM == 2)
    {
        // Get relevant info about the element
        ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_elem = this->GetElement(elemIndex);
        unsigned num_nodes = p_elem->GetNumNodes();
        double area_of_elem = this->GetVolumeOfElement(elemIndex);
        c_vector<double, SPACE_DIM> centroid = this->GetCentroidOfElement(elemIndex);

        // Get the unit axis and trig terms for rotation
        c_vector<double, SPACE_DIM> unit_axis = axis / norm_2(axis);
        double sin_theta = unit_axis[0];
        double cos_theta = unit_axis[1];

        // We need the (rotated) node locations in two orders - original and ordered left-to-right.
        // For the latter we need to keep track of index, so we store that as part of a pair.
        std::vector<c_vector<double, SPACE_DIM> > node_locations_original_order;
        std::vector<std::pair<unsigned, c_vector<double, SPACE_DIM> > > ordered_locations;

        // Get the node locations of the current element relative to its centroid, and rotate them
        for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
        {
            const c_vector<double, SPACE_DIM>& node_location = p_elem->GetNode(node_idx)->rGetLocation();

            c_vector<double, SPACE_DIM> displacement = this->GetVectorFromAtoB(centroid, node_location);

            c_vector<double, SPACE_DIM> rotated_location;
            rotated_location[0] = cos_theta * displacement[0] - sin_theta * displacement[1];
            rotated_location[1] = sin_theta * displacement[0] + cos_theta * displacement[1];

            node_locations_original_order.push_back(rotated_location);
        }

        // Fill up a vector of identical points, and sort it so nodes are ordered in ascending x value
        for (unsigned i = 0; i < node_locations_original_order.size(); i++)
        {
            ordered_locations.push_back(std::pair<unsigned, c_vector<double, SPACE_DIM> >(i, node_locations_original_order[i]));
        }

        std::sort(ordered_locations.begin(), ordered_locations.end(), CustomComparisonForSkewnessMeasure);

        /*
        * For each node, we must find every place where the axis (now rotated to be vertical) intersects the polygon:
        *
        *       |
        *     __|______
        *    /  |      \
        *   /   |       \
        *  /____|___    |
        *       |  |    |
        *  _____|__|    |
        *  \    |       |
        *   \   |      /
        *    \__|_____/
        *       |
        *       |
        *       ^
        * For instance, the number of times the vertical intersects the polygon above is 4 and, for each node, we need to
        * find all such intersections.  We can do this by checking where the dot product of the the vector a with the unit
        * x direction changes sign as we iterate over the original node locations, where a is the vector from the current
        * node to the test node.
        */

        // For each node, we keep track of all the y-locations where the vertical through the node meets the polygon
        std::vector<std::vector<double> > knots(num_nodes);

        // Iterate over ordered locations from left to right
        for (unsigned location = 0; location < num_nodes; location++)
        {
            // Get the two parts of the pair
            unsigned this_idx = ordered_locations[location].first;
            // this_location is the location of the node relative to centroid
            c_vector<double, SPACE_DIM> this_location = ordered_locations[location].second;

            // The y-coordinate of the current location is always a knot
            // because we are passing the vertical through this node
            knots[location].push_back(this_location[1]);

            // To calculate all the intersection points, we need to iterate over every other location and see, sequentially,
            // if the x-coordinate of location i+1 and i+2 crosses the x-coordinate of the current location.
            // i.e. check whether each other edge crosses the vertical through this_location/current node
            unsigned next_idx = (this_idx + 1) % num_nodes;
            c_vector<double, SPACE_DIM> to_previous = node_locations_original_order[next_idx] - this_location;
            c_vector<double, SPACE_DIM> to_next;

            for (unsigned node_idx = this_idx + 2; node_idx < this_idx + num_nodes; node_idx++)
            {
                unsigned idx = node_idx % num_nodes;

                to_next = node_locations_original_order[idx] - this_location;

                // If the segment between to_previous and to_next intersects the vertical through this_location, the clause
                // in the if statement below will be triggered
                if (to_previous[0] * to_next[0] <= 0.0)
                {
                    // Find how far between to_previous and to_next the point of intersection is
                    double interp = 0.5;
                    if (to_previous[0] - to_next[0] != 0.0)
                    {
                        interp = to_previous[0] / (to_previous[0] - to_next[0]);
                    }

                    assert(interp >= 0.0 && interp <= 1.0);

                    // Record the y-value of the intersection point
                    double new_intersection = this_location[1] + to_previous[1] + interp * (to_next[1] - to_previous[1]);
                    knots[location].push_back(new_intersection);
                }

                to_previous = to_next;
            }

            if (knots[location].size() > 2)
            {
                WARN_ONCE_ONLY("Axis intersects polygon more than 2 times (concavity) - check element is fairly convex.");
            }
        }

        // For ease, construct a vector of the x-locations of all the nodes, in order
        std::vector<double> ordered_x(num_nodes);
        for (unsigned location = 0; location < num_nodes; location++)
        {
            ordered_x[location] = ordered_locations[location].second[0];
        }

        // Calculate the mass contributions at each x-location - this is the length of the intersection of the vertical
        // through each location
        std::vector<double> mass_contributions(num_nodes);
        for (unsigned i = 0; i < num_nodes; i++)
        {
            std::sort(knots[i].begin(), knots[i].end());

            switch (knots[i].size())
            {
                case 1:
                    mass_contributions[i] = 0.0;
                    break;

                case 2:
                    mass_contributions[i] = knots[i][1] - knots[i][0];
                    break;

                default:
                    mass_contributions[i] += knots[i][knots[i].size()-1] - knots[i][0];
            }

            // Normalise, so that these lengths define a pdf
            mass_contributions[i] /= area_of_elem;
        }

        // Calculate moments. Because we just have a bunch of linear segments, we can integrate the pdf exactly
        double e_x0 = 0.0;
        double e_x1 = 0.0;
        double e_x2 = 0.0;
        double e_x3 = 0.0;

        for (unsigned i = 1; i < num_nodes; i++)
        {
            double x0 = ordered_x[i - 1];
            double x1 = ordered_x[i];

            double fx0 = mass_contributions[i - 1];
            double fx1 = mass_contributions[i];

            // We need squared, cubed, ..., order 5 for each x
            double x0_2 = x0 * x0;
            double x0_3 = x0_2 * x0;
            double x0_4 = x0_3 * x0;
            double x0_5 = x0_4 * x0;

            double x1_2 = x1 * x1;
            double x1_3 = x1_2 * x1;
            double x1_4 = x1_3 * x1;
            double x1_5 = x1_4 * x1;

            if (x1 - x0 > 0)
            {
                // Calculate y = mx + c for this section of the pdf
                double m = (fx1 - fx0) / (x1 - x0);
                double c = fx0 - m * x0;

                e_x0 += m * (x1_2 - x0_2) / 2.0 + c * (x1 - x0);
                e_x1 += m * (x1_3 - x0_3) / 3.0 + c * (x1_2 - x0_2) / 2.0;
                e_x2 += m * (x1_4 - x0_4) / 4.0 + c * (x1_3 - x0_3) / 3.0;
                e_x3 += m * (x1_5 - x0_5) / 5.0 + c * (x1_4 - x0_4) / 4.0;
            }
        }

        // Check that we have correctly defined a pdf
        if (fabs(e_x0 - 1.0) < 1e-6)
        {
            WARN_ONCE_ONLY("Mass distribution of element calculated incorrectly due to element concavity. Skewness may not be correct!");
        }

        // Calculate the standard deviation, and return the skewness
        double sd = sqrt(e_x2 - e_x1 * e_x1);
        return (e_x3 - 3.0 * e_x1 * sd * sd - e_x1 * e_x1 * e_x1) / (sd * sd * sd);
    }
    else
    {
        NEVER_REACHED;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ChasteCuboid<SPACE_DIM> ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::CalculateBoundingBoxOfElement(unsigned index)
{
    ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_elem = this->GetElement(index);

    // Get the location of node zero as a reference point
    c_vector<double, SPACE_DIM> ref_point = p_elem->GetNode(0)->rGetLocation();

    // Vector to represent the n-dimensional 'bottom left'-most node
    c_vector<double, SPACE_DIM> bottom_left = zero_vector<double>(SPACE_DIM);

    // Vector to represent the n-dimensional 'top right'-most node
    c_vector<double, SPACE_DIM> top_right = zero_vector<double>(SPACE_DIM);

    // Loop over all nodes in the element and update bottom_left and top_right, relative to node zero to account for periodicity
    for (unsigned node_idx = 0; node_idx < p_elem->GetNumNodes(); node_idx++)
    {
        c_vector<double, SPACE_DIM> vec_to_node = this->GetVectorFromAtoB(ref_point, p_elem->GetNode(node_idx)->rGetLocation());

        for (unsigned dim = 0; dim < SPACE_DIM; dim++)
        {
            if (vec_to_node[dim] < bottom_left[dim])
            {
                bottom_left[dim] = vec_to_node[dim];
            }
            else if (vec_to_node[dim] > top_right[dim])
            {
                top_right[dim] = vec_to_node[dim];
            }
        }
    }

    // Create Chaste points, rescaled by the location of node zero
    ChastePoint<SPACE_DIM> min(bottom_left + ref_point);
    ChastePoint<SPACE_DIM> max(top_right + ref_point);

    return ChasteCuboid<SPACE_DIM>(min, max);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryMesh()
{
    this->mMeshChangesDuringSimulation = false;
    Clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::~ImmersedBoundaryMesh()
{
    Clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodeIndices(unsigned nodeIndex)
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
        unsigned previous_local_index = (local_index + num_nodes - 1) % num_nodes;
        unsigned next_local_index = (local_index + 1) % num_nodes;

        // Add the global indices of these two nodes to the set of neighbouring node indices
        neighbouring_node_indices.insert(GetElement(*elem_iter)->GetNodeGlobalIndex(previous_local_index));
        neighbouring_node_indices.insert(GetElement(*elem_iter)->GetNodeGlobalIndex(next_local_index));
    }

    return neighbouring_node_indices;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size());
    return index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size());
    return index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    return index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    // Delete elements
    for (unsigned i = 0; i < mElements.size(); i++)
    {
        delete mElements[i];
    }
    mElements.clear();

    // Delete laminas
    for (auto lamina : mLaminas)
    {
        delete(lamina);
    }
    mLaminas.clear();

    // Delete nodes
    for (unsigned i = 0; i < this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }
    this->mNodes.clear();

    mBalancingFluidSources.clear();
    mElementFluidSources.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetCharacteristicNodeSpacing() const
{
    return mCharacteristicNodeSpacing;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetSpacingRatio() const
{
    return mCharacteristicNodeSpacing / (1.0 / double(mNumGridPtsX));
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<Node<SPACE_DIM>*>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGetNodes() const
{
    return this->mNodes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNumGridPtsX() const
{
    return mNumGridPtsX;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNumGridPtsY() const
{
    return mNumGridPtsY;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetNumGridPtsX(unsigned meshPointsX)
{
    mNumGridPtsX = meshPointsX;
    m2dVelocityGrids.resize(extents[2][mNumGridPtsX][mNumGridPtsY]);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetNumGridPtsY(unsigned meshPointsY)
{
    mNumGridPtsY = meshPointsY;
    m2dVelocityGrids.resize(extents[2][mNumGridPtsX][mNumGridPtsY]);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetNumGridPtsXAndY(unsigned numGridPts)
{
    mNumGridPtsX = numGridPts;
    mNumGridPtsY = numGridPts;
    m2dVelocityGrids.resize(extents[2][mNumGridPtsX][mNumGridPtsY]);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetCharacteristicNodeSpacing(double nodeSpacing)
{
    mCharacteristicNodeSpacing = nodeSpacing;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::shared_ptr<FluidSource<SPACE_DIM>>>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGetElementFluidSources()
{
    mElementFluidSources.clear();
    for (const auto& element : mElements) {
        if (element->GetFluidSource() != nullptr) {
            mElementFluidSources.push_back(element->GetFluidSource());
        }
    }
    return mElementFluidSources;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::shared_ptr<FluidSource<SPACE_DIM>>>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGetBalancingFluidSources()
{
    return mBalancingFluidSources;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const multi_array<double, 3>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGet2dVelocityGrids() const
{
    return m2dVelocityGrids;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
multi_array<double, 3>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGetModifiable2dVelocityGrids()
{
    return m2dVelocityGrids;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(
    Node<SPACE_DIM>* pNewNode)
{
    pNewNode->SetIndex(this->mNodes.size());
    this->mNodes.push_back(pNewNode);
    return this->mNodes.size() - 1;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetCellRearrangementThreshold(double cellRearrangementThreshold)
{
    mCellRearrangementThreshold = cellRearrangementThreshold;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetCellRearrangementThreshold()
{
    return mCellRearrangementThreshold;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocation1, const c_vector<double, SPACE_DIM>& rLocation2)
{
    // This code currently assumes the grid is precisely [0,1)x[0,1)
    c_vector<double, SPACE_DIM> vector = rLocation2 - rLocation1;

    /*
     * Handle the periodic condition here: if the points are more
     * than 0.5 apart in any direction, choose -(1.0-dist).
     */
    for (unsigned dim = 0; dim < SPACE_DIM; dim++)
    {
        if (fabs(vector[dim]) > 0.5)
        {
            vector[dim] = copysign(fabs(vector[dim]) - 1.0, -vector[dim]);
        }
    }

    return vector;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point)
{
    this->mNodes[nodeIndex]->SetPoint(point);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements() const
{
    return mElements.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLaminas() const
{
    return mLaminas.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>* ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetLamina(unsigned index) const
{
    assert(index < mLaminas.size());
    return mLaminas[index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbourDist() const
{
    return mNeighbourDist;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetNeighbourDist(double neighbourDist)
{
    mNeighbourDist = neighbourDist;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetCentroidOfElement(
    [[maybe_unused]] unsigned index) // [[maybe_unused]] due to unused-but-set-parameter warning in GCC 7,8,9
{
    // Only implemented in 2D
    if constexpr (SPACE_DIM == 2)
    {
        ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

        unsigned num_nodes = p_element->GetNumNodes();
        c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);

        double centroid_x = 0;
        double centroid_y = 0;

        // Note that we cannot use GetVolumeOfElement() below as it returns the absolute, rather than signed, area
        double element_signed_area = 0.0;

        // Map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
        c_vector<double, SPACE_DIM> first_node_location = p_element->GetNodeLocation(0);
        c_vector<double, SPACE_DIM> pos_1 = zero_vector<double>(SPACE_DIM);

        // Loop over vertices
        for (unsigned local_index = 0; local_index < num_nodes; local_index++)
        {
            c_vector<double, SPACE_DIM> next_node_location = p_element->GetNodeLocation((local_index + 1) % num_nodes);
            c_vector<double, SPACE_DIM> pos_2 = GetVectorFromAtoB(first_node_location, next_node_location);

            double this_x = pos_1[0];
            double this_y = pos_1[1];
            double next_x = pos_2[0];
            double next_y = pos_2[1];

            double signed_area_term = this_x * next_y - this_y * next_x;

            centroid_x += (this_x + next_x) * signed_area_term;
            centroid_y += (this_y + next_y) * signed_area_term;
            element_signed_area += 0.5 * signed_area_term;

            pos_1 = pos_2;
        }

        assert(element_signed_area != 0.0);

        // Finally, map back and employ GetVectorFromAtoB() to allow for periodicity
        centroid = first_node_location;
        centroid[0] += centroid_x / (6.0 * element_signed_area);
        centroid[1] += centroid_y / (6.0 * element_signed_area);

        centroid[0] = centroid[0] < 0 ? centroid[0] + 1.0 : fmod(centroid[0], 1.0);
        centroid[1] = centroid[1] < 0 ? centroid[1] + 1.0 : fmod(centroid[1], 1.0);

        return centroid;
    }
    else
    {
        NEVER_REACHED;
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void ImmersedBoundaryMesh<1, 1>::ConstructFromMeshReader(AbstractMeshReader<1, 1>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("ImmersedBoundaryMesh not yet supported for the specified dimensions");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void ImmersedBoundaryMesh<1, 2>::ConstructFromMeshReader(AbstractMeshReader<1, 2>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("ImmersedBoundaryMesh not yet supported for the specified dimensions");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void ImmersedBoundaryMesh<1, 3>::ConstructFromMeshReader(AbstractMeshReader<1, 3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("ImmersedBoundaryMesh not yet supported for the specified dimensions");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void ImmersedBoundaryMesh<2, 3>::ConstructFromMeshReader(AbstractMeshReader<2, 3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("ImmersedBoundaryMesh not yet supported for the specified dimensions");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void ImmersedBoundaryMesh<2, 2>::ConstructFromMeshReader(AbstractMeshReader<2, 2>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    ImmersedBoundaryMeshReader<2, 2>& rIBMeshReader = dynamic_cast<ImmersedBoundaryMeshReader<2, 2>&>(rMeshReader);

    assert(!rIBMeshReader.HasNodePermutation());

    // Store numbers of nodes and elements
    unsigned num_nodes = rIBMeshReader.GetNumNodes();
    unsigned num_elements = rIBMeshReader.GetNumElements();
    unsigned num_laminas = rIBMeshReader.GetNumLaminas();
    this->mCharacteristicNodeSpacing = rIBMeshReader.GetCharacteristicNodeSpacing();

    // Add nodes
    rIBMeshReader.Reset();
    mNodes.reserve(num_nodes);
    std::vector<double> node_data;
    for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
    {
        node_data = rIBMeshReader.GetNextNode();
        unsigned is_boundary_node = (bool)node_data[2];
        node_data.pop_back();
        this->mNodes.push_back(new Node<2>(node_idx, node_data, is_boundary_node));
    }

    // Add laminas
    rIBMeshReader.Reset();
    mLaminas.reserve(num_laminas);
    for (unsigned lam_idx = 0; lam_idx < num_laminas; lam_idx++)
    {
        // Get the data for this element
        ImmersedBoundaryElementData lamina_data = rIBMeshReader.GetNextImmersedBoundaryLaminaData();

        // Get the nodes owned by this element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes_in_lamina = lamina_data.NodeIndices.size();
        for (unsigned node_idx = 0; node_idx < num_nodes_in_lamina; node_idx++)
        {
            assert(lamina_data.NodeIndices[node_idx] < this->mNodes.size());
            nodes.push_back(this->mNodes[lamina_data.NodeIndices[node_idx]]);
        }

        // Use nodes and index to construct this element
        ImmersedBoundaryElement<1, 2>* p_lamina = new ImmersedBoundaryElement<1, 2>(lam_idx, nodes);
        mLaminas.push_back(p_lamina);

        if (rIBMeshReader.GetNumLaminaAttributes() > 0)
        {
            assert(rIBMeshReader.GetNumLaminaAttributes() == 1);
            unsigned attribute_value = lamina_data.AttributeValue;
            p_lamina->SetAttribute(attribute_value);
        }
    }

    // Add elements
    rIBMeshReader.Reset();
    mElements.reserve(num_elements);
    for (unsigned elem_idx = 0; elem_idx < num_elements; elem_idx++)
    {
        // Get the data for this element
        ImmersedBoundaryElementData element_data = rIBMeshReader.GetNextImmersedBoundaryElementData();

        // Get the nodes owned by this element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned node_idx = 0; node_idx < num_nodes_in_element; node_idx++)
        {
            assert(element_data.NodeIndices[node_idx] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[node_idx]]);
        }

        // Use nodes and index to construct this element
        ImmersedBoundaryElement<2, 2>* p_element = new ImmersedBoundaryElement<2, 2>(elem_idx, nodes);
        mElements.push_back(p_element);

        if (rIBMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rIBMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }

    // Get grid dimensions from grid file and set up grids accordingly
    this->mNumGridPtsX = rIBMeshReader.GetNumGridPtsX();
    this->mNumGridPtsY = rIBMeshReader.GetNumGridPtsY();
    m2dVelocityGrids.resize(extents[2][mNumGridPtsX][mNumGridPtsY]);

    // Construct the velocity grids from mesh reader
    for (unsigned dim = 0; dim < 2; dim++)
    {
        for (unsigned grid_row = 0; grid_row < mNumGridPtsY; grid_row++)
        {
            std::vector<double> next_row = rIBMeshReader.GetNextGridRow();
            assert(next_row.size() == mNumGridPtsX);

            for (unsigned i = 0; i < mNumGridPtsX; i++)
            {
                m2dVelocityGrids[dim][i][grid_row] = next_row[i];
            }
        }
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void ImmersedBoundaryMesh<3, 3>::ConstructFromMeshReader(AbstractMeshReader<3, 3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("ImmersedBoundaryMesh not yet supported for the specified dimensions");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetVolumeOfElement(unsigned index)
{
    if constexpr (SPACE_DIM == 2)
    {
        // Get pointer to this element
        ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

        double element_volume = 0.0;

        // We map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
        c_vector<double, SPACE_DIM> first_node_location = p_element->GetNodeLocation(0);
        c_vector<double, SPACE_DIM> pos_1 = zero_vector<double>(SPACE_DIM);

        unsigned num_nodes = p_element->GetNumNodes();
        for (unsigned local_index = 0; local_index < num_nodes; local_index++)
        {
            c_vector<double, SPACE_DIM> next_node_location = p_element->GetNodeLocation((local_index + 1) % num_nodes);
            c_vector<double, SPACE_DIM> pos_2 = GetVectorFromAtoB(first_node_location, next_node_location);

            double this_x = pos_1[0];
            double this_y = pos_1[1];
            double next_x = pos_2[0];
            double next_y = pos_2[1];

            element_volume += 0.5 * (this_x * next_y - next_x * this_y);

            pos_1 = pos_2;
        }

        // We take the absolute value just in case the nodes were really oriented clockwise
        return fabs(element_volume);
    }
    else
    {
        NEVER_REACHED;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetSurfaceAreaOfElement(unsigned index)
{
    if constexpr (SPACE_DIM == 2)
    {
        // Get pointer to this element
        ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

        double surface_area = 0.0;
        unsigned num_nodes = p_element->GetNumNodes();
        unsigned this_node_index = p_element->GetNodeGlobalIndex(0);
        for (unsigned local_index = 0; local_index < num_nodes; local_index++)
        {
            unsigned next_node_index = p_element->GetNodeGlobalIndex((local_index + 1) % num_nodes);
            surface_area += this->GetDistanceBetweenNodes(this_node_index, next_node_index);
            this_node_index = next_node_index;
        }

        return surface_area;
    }
    else
    {
        NEVER_REACHED;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetVoronoiSurfaceAreaOfElement(const unsigned elemIdx)
{
    if constexpr (SPACE_DIM == 2)
    {
        UpdateNodeLocationsVoronoiDiagramIfOutOfDate();
        ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(elemIdx);
        double surface_area = 0.0;

        // Sum contributions from each node in the element
        for (unsigned local_idx = 0; local_idx < p_element->GetNumNodes(); local_idx++)
        {
            const unsigned global_node_idx = p_element->GetNodeGlobalIndex(local_idx);

            // Get the voronoi cell that corresponds to this node
            const unsigned voronoi_cell_id = mVoronoiCellIdsIndexedByNodeIndex[global_node_idx];
            const auto& voronoi_cell = mNodeLocationsVoronoiDiagram.cells()[voronoi_cell_id];

            // Iterate over the edges of this voronoi cell
            auto p_edge = voronoi_cell.incident_edge();
            do
            {
                // The global node index corresponding to a voronoi cell is encoded in its 'color' variable
                const unsigned twin_node_idx = p_edge->twin()->cell()->color();
                const unsigned twin_elem_idx = *this->GetNode(twin_node_idx)->ContainingElementsBegin();

                // Check if the nodes are in different elements
                if (elemIdx != twin_elem_idx)
                {
                    surface_area += this->CalculateLengthOfVoronoiEdge(*p_edge);
                }

                p_edge = p_edge->next();
            } while (p_edge != voronoi_cell.incident_edge());
        }

        return surface_area;
    }
    else
    {
        NEVER_REACHED;
    }
}

// Excluded from coverage as it is impossible to construct a 1 dimensional element currently
//LCOV_EXCL_START
/**
 * Template specialisation for impossible to reach code path
 * @param elemIdx the element index
 * @returns 0.0
 */
template <>
double ImmersedBoundaryMesh<1, 1>::GetVoronoiSurfaceAreaOfElement(const unsigned elemIdx)
{
    return 0.0;
}
//LCOV_EXCL_STOP

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetAverageNodeSpacingOfElement(unsigned index, bool recalculate)
{
    if (recalculate || (this->GetElement(index)->GetAverageNodeSpacing() == DOUBLE_UNSET))
    {
        double average_node_spacing = this->GetSurfaceAreaOfElement(index) / this->GetElement(index)->GetNumNodes();
        this->GetElement(index)->SetAverageNodeSpacing(average_node_spacing);

        return average_node_spacing;
    }
    else
    {
        return this->GetElement(index)->GetAverageNodeSpacing();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetAverageNodeSpacingOfLamina(unsigned index, bool recalculate)
{
    if (recalculate || (this->GetLamina(index)->GetAverageNodeSpacing() == DOUBLE_UNSET))
    {
        ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>* p_lam = this->GetLamina(index);

        // Explicitly calculate the average node spacing
        double average_node_spacing = 0.0;
        for (unsigned node_it = 1; node_it < p_lam->GetNumNodes(); node_it++)
        {
            average_node_spacing += this->GetDistanceBetweenNodes(p_lam->GetNodeGlobalIndex(node_it),
                                                                  p_lam->GetNodeGlobalIndex(node_it - 1));
        }

        average_node_spacing /= (p_lam->GetNumNodes() - 1);

        // Set it for quick retrieval next time
        p_lam->SetAverageNodeSpacing(average_node_spacing);
        return average_node_spacing;
    }
    else
    {
        return this->GetLamina(index)->GetAverageNodeSpacing();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetElementDivisionSpacing()
{
    return mElementDivisionSpacing;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetElementDivisionSpacing(double elementDivisionSpacing)
{
    mElementDivisionSpacing = elementDivisionSpacing;
}

//////////////////////////////////////////////////////////////////////
//                        2D-specific methods                       //
//////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMomentsOfElement(unsigned index)
{
    if constexpr (SPACE_DIM == 2)
    {
        // Define helper variables
        ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);
        unsigned num_nodes = p_element->GetNumNodes();
        c_vector<double, 3> moments = zero_vector<double>(3);

        // Since we compute I_xx, I_yy and I_xy about the centroid, we must shift each vertex accordingly
        c_vector<double, SPACE_DIM> centroid = GetCentroidOfElement(index);

        c_vector<double, SPACE_DIM> this_node_location = p_element->GetNodeLocation(0);
        c_vector<double, SPACE_DIM> pos_1 = this->GetVectorFromAtoB(centroid, this_node_location);

        for (unsigned local_index = 0; local_index < num_nodes; local_index++)
        {
            unsigned next_index = (local_index + 1) % num_nodes;
            c_vector<double, SPACE_DIM> next_node_location = p_element->GetNodeLocation(next_index);
            c_vector<double, SPACE_DIM> pos_2 = this->GetVectorFromAtoB(centroid, next_node_location);

            double signed_area_term = pos_1(0) * pos_2(1) - pos_2(0) * pos_1(1);
            // Ixx
            moments(0) += (pos_1(1) * pos_1(1) + pos_1(1) * pos_2(1) + pos_2(1) * pos_2(1)) * signed_area_term;

            // Iyy
            moments(1) += (pos_1(0) * pos_1(0) + pos_1(0) * pos_2(0) + pos_2(0) * pos_2(0)) * signed_area_term;

            // Ixy
            moments(2) += (pos_1(0) * pos_2(1) + 2 * pos_1(0) * pos_1(1) + 2 * pos_2(0) * pos_2(1) + pos_2(0) * pos_1(1)) * signed_area_term;

            pos_1 = pos_2;
        }

        moments(0) /= 12;
        moments(1) /= 12;
        moments(2) /= 24;

        /*
        * If the nodes owned by the element were supplied in a clockwise rather
        * than anticlockwise manner, or if this arose as a result of enforcing
        * periodicity, then our computed quantities will be the wrong sign, so
        * we need to fix this.
        */
        if (moments(0) < 0.0)
        {
            moments(0) = -moments(0);
            moments(1) = -moments(1);
            moments(2) = -moments(2);
        }
        return moments;
    }
    else
    {
        NEVER_REACHED;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetShortAxisOfElement(unsigned index)
{
    if constexpr (SPACE_DIM == 2)
    {
        c_vector<double, SPACE_DIM> short_axis = zero_vector<double>(SPACE_DIM);

        // Calculate the moments of the element about its centroid (recall that I_xx and I_yy must be non-negative)
        c_vector<double, 3> moments = CalculateMomentsOfElement(index);

        // Normalise the moments vector to remove problem of a very small discriminant (see #2874)
        moments /= norm_2(moments);

        // If the principal moments are equal...
        double discriminant = (moments(0) - moments(1)) * (moments(0) - moments(1)) + 4.0 * moments(2) * moments(2);
        if (fabs(discriminant) < DBL_EPSILON)
        {
            // ...then every axis through the centroid is a principal axis, so return a random unit vector
            short_axis(0) = RandomNumberGenerator::Instance()->ranf();
            short_axis(1) = sqrt(1.0 - short_axis(0) * short_axis(0));
        }
        else
        {
            // If the product of inertia is zero, then the coordinate axes are the principal axes
            if (fabs(moments(2)) < DBL_EPSILON)
            {
                if (moments(0) < moments(1))
                {
                    short_axis(0) = 0.0;
                    short_axis(1) = 1.0;
                }
                else
                {
                    short_axis(0) = 1.0;
                    short_axis(1) = 0.0;
                }
            }
            else
            {
                // Otherwise we find the eigenvector of the inertia matrix corresponding to the largest eigenvalue
                double lambda = 0.5 * (moments(0) + moments(1) + sqrt(discriminant));

                short_axis(0) = 1.0;
                short_axis(1) = (moments(0) - lambda) / moments(2);

                // Normalise the short axis before returning it
                short_axis /= norm_2(short_axis);
            }
        }

        return short_axis;
    }
    else
    {
        NEVER_REACHED;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::DivideElementAlongGivenAxis(ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                                   c_vector<double, SPACE_DIM> axisOfDivision,
                                                                                   bool placeOriginalElementBelow)
{
    if constexpr (SPACE_DIM == 2 && ELEMENT_DIM == 2)
    {
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
        for (unsigned i = 0; i < num_nodes; i++)
        {
            bool is_next_node_on_left = (inner_prod(this->GetVectorFromAtoB(pElement->GetNodeLocation((i + 1) % num_nodes), centroid), perp_axis) >= 0);
            if (is_current_node_on_left != is_next_node_on_left)
            {
                intersecting_nodes.push_back(i);
            }
            is_current_node_on_left = is_next_node_on_left;
        }

        // If the axis of division does not cross two edges then we cannot proceed
        if (intersecting_nodes.size() != 2)
        {
            EXCEPTION("Cannot proceed with element division: the given axis of division does not cross two edges of the element"); // LCOV_EXCL_LINE
        }

        // Now call DivideElement() to divide the element using the nodes found above
        //unsigned new_element_index = 0;
        unsigned new_element_index = DivideElement(pElement,
                                                intersecting_nodes[0],
                                                intersecting_nodes[1],
                                                centroid,
                                                axisOfDivision);

        return new_element_index;
    }
    else
    {
        NEVER_REACHED;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::DivideElementAlongShortAxis(ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                                   bool placeOriginalElementBelow)
{
    if constexpr (SPACE_DIM == 2 && ELEMENT_DIM == 2)
    {
        c_vector<double, SPACE_DIM> short_axis = this->GetShortAxisOfElement(pElement->GetIndex());
        unsigned new_element_index = DivideElementAlongGivenAxis(pElement, short_axis, placeOriginalElementBelow);
        return new_element_index;
    }
    else
    {
        NEVER_REACHED;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::DivideElement(ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                     unsigned nodeAIndex,
                                                                     unsigned nodeBIndex,
                                                                     c_vector<double, SPACE_DIM> centroid,
                                                                     c_vector<double, SPACE_DIM> axisOfDivision)
{
    if constexpr (SPACE_DIM == 2 && ELEMENT_DIM == 2)
    {
        if (mElementDivisionSpacing == DOUBLE_UNSET)
        {
            EXCEPTION("The value of mElementDivisionSpacing has not been set.");
        }

        /*
        * Method outline:
        *
        *   Each element needs to end up with the same number of nodes as the original element, and those nodes will be
        *   equally spaced around the outline of each of the two daughter elements.
        *
        *   The two elements need to be divided by a distance of mElementDivisionSpacing, where the distance is measured
        *   perpendicular to the axis of division.
        *
        *   To achieve this, we find four 'corner' locations, each of which has a perpendicular distance from the axis of
        *   half the required spacing, and are found by using the locations from the existing element as a stencil.
        */

        double half_spacing = 0.5 * mElementDivisionSpacing;

        // Get unit vectors in the direction of the division axis, and the perpendicular
        c_vector<double, SPACE_DIM> unit_axis = axisOfDivision / norm_2(axisOfDivision);
        c_vector<double, SPACE_DIM> unit_perp;
        unit_perp[0] = -unit_axis[1];
        unit_perp[1] = unit_axis[0];

        unsigned num_nodes = pElement->GetNumNodes();

        /*
        * We first identify the start and end indices of the nodes which will form the location stencil for each daughter
        * cell.  Our starting point is the node indices already identified.
        *
        * In order to ensure the resulting gap between the elements is the correct size, we remove as many nodes as
        * necessary until the perpendicular distance between the centroid and the node is at least half the required
        * spacing.
        *
        * Finally, we move the relevant node to be exactly half the required spacing.
        */
        unsigned start_a = (nodeAIndex + 1) % num_nodes;
        unsigned end_a = nodeBIndex;

        unsigned start_b = (nodeBIndex + 1) % num_nodes;
        unsigned end_b = nodeAIndex;

        // Find correct start_a
        bool no_node_satisfied_condition_1 = true;
        for (unsigned i = start_a; i != end_a;)
        {
            c_vector<double, SPACE_DIM> centroid_to_i = this->GetVectorFromAtoB(centroid, pElement->GetNode(i)->rGetLocation());
            double perpendicular_dist = inner_prod(centroid_to_i, unit_perp);

            if (fabs(perpendicular_dist) >= half_spacing)
            {
                no_node_satisfied_condition_1 = false;
                start_a = i;

                // Calculate position so it's exactly 0.5 * elem_spacing perpendicular distance from the centroid
                c_vector<double, SPACE_DIM> new_location = pElement->GetNode(i)->rGetLocation();
                new_location -= unit_perp * copysign(fabs(perpendicular_dist) - half_spacing, perpendicular_dist);

                pElement->GetNode(i)->SetPoint(ChastePoint<SPACE_DIM>(new_location));
                break;
            }

            // Go to the next node
            i = (i + 1) % num_nodes;
        }

        // Find correct end_a
        bool no_node_satisfied_condition_2 = true;
        for (unsigned i = end_a; i != start_a;)
        {
            c_vector<double, SPACE_DIM> centroid_to_i = this->GetVectorFromAtoB(centroid, pElement->GetNode(i)->rGetLocation());
            double perpendicular_dist = inner_prod(centroid_to_i, unit_perp);

            if (fabs(perpendicular_dist) >= half_spacing)
            {
                no_node_satisfied_condition_2 = false;
                end_a = i;

                // Calculate position so it's exactly 0.5 * elem_spacing perpendicular distance from the centroid
                c_vector<double, SPACE_DIM> new_location = pElement->GetNode(i)->rGetLocation();
                new_location -= unit_perp * copysign(fabs(perpendicular_dist) - half_spacing, perpendicular_dist);

                pElement->GetNode(i)->SetPoint(ChastePoint<SPACE_DIM>(new_location));
                break;
            }

            // Go to the previous node
            i = (i + num_nodes - 1) % num_nodes; // LCOV_EXCL_LINE
        }

        // Find correct start_b
        bool no_node_satisfied_condition_3 = true;
        for (unsigned i = start_b; i != end_b;)
        {
            c_vector<double, SPACE_DIM> centroid_to_i = this->GetVectorFromAtoB(centroid, pElement->GetNode(i)->rGetLocation());
            double perpendicular_dist = inner_prod(centroid_to_i, unit_perp);

            if (fabs(perpendicular_dist) >= half_spacing)
            {
                no_node_satisfied_condition_3 = false;
                start_b = i;

                // Calculate position so it's exactly 0.5 * elem_spacing perpendicular distance from the centroid
                c_vector<double, SPACE_DIM> new_location = pElement->GetNode(i)->rGetLocation();
                new_location -= unit_perp * copysign(fabs(perpendicular_dist) - half_spacing, perpendicular_dist);

                pElement->GetNode(i)->SetPoint(ChastePoint<SPACE_DIM>(new_location));
                break;
            }

            // Go to the next node
            i = (i + 1) % num_nodes;
        }

        // Find correct end_b
        bool no_node_satisfied_condition_4 = true;
        for (unsigned i = end_b; i != start_b;)
        {
            c_vector<double, SPACE_DIM> centroid_to_i = this->GetVectorFromAtoB(centroid, pElement->GetNode(i)->rGetLocation());
            double perpendicular_dist = inner_prod(centroid_to_i, unit_perp);

            if (fabs(perpendicular_dist) >= half_spacing)
            {
                no_node_satisfied_condition_4 = false;
                end_b = i;

                // Calculate position so it's exactly 0.5 * elem_spacing perpendicular distance from the centroid
                c_vector<double, SPACE_DIM> new_location = pElement->GetNode(i)->rGetLocation();
                new_location -= unit_perp * copysign(fabs(perpendicular_dist) - half_spacing, perpendicular_dist);

                pElement->GetNode(i)->SetPoint(ChastePoint<SPACE_DIM>(new_location));
                break;
            }

            // Go to the previous node
            i = (i + num_nodes - 1) % num_nodes;
        }

        if (no_node_satisfied_condition_1 || no_node_satisfied_condition_2 || no_node_satisfied_condition_3 || no_node_satisfied_condition_4)
        {
            EXCEPTION("Could not space elements far enough apart during cell division.  Cannot currently handle this case");
        }

        /*
        * Create location stencils for each of the daughter cells
        */
        std::vector<c_vector<double, SPACE_DIM> > daughter_a_location_stencil;
        for (unsigned node_idx = start_a; node_idx != (end_a + 1) % num_nodes;)
        {
            daughter_a_location_stencil.push_back(c_vector<double, SPACE_DIM>(pElement->GetNode(node_idx)->rGetLocation()));

            // Go to next node
            node_idx = (node_idx + 1) % num_nodes;
        }

        std::vector<c_vector<double, SPACE_DIM> > daughter_b_location_stencil;
        for (unsigned node_idx = start_b; node_idx != (end_b + 1) % num_nodes;)
        {
            daughter_b_location_stencil.push_back(c_vector<double, SPACE_DIM>(pElement->GetNode(node_idx)->rGetLocation()));

            // Go to next node
            node_idx = (node_idx + 1) % num_nodes;
        }

        assert(daughter_a_location_stencil.size() > 1);
        assert(daughter_b_location_stencil.size() > 1);

        // To help calculating cumulative distances, add the first location on to the end
        daughter_a_location_stencil.push_back(daughter_a_location_stencil[0]);
        daughter_b_location_stencil.push_back(daughter_b_location_stencil[0]);

        // Calculate the cumulative distances around the stencils
        std::vector<double> cumulative_distances_a;
        std::vector<double> cumulative_distances_b;
        cumulative_distances_a.push_back(0.0);
        cumulative_distances_b.push_back(0.0);
        for (unsigned loc_idx = 1; loc_idx < daughter_a_location_stencil.size(); loc_idx++)
        {
            cumulative_distances_a.push_back(cumulative_distances_a.back() + norm_2(this->GetVectorFromAtoB(daughter_a_location_stencil[loc_idx - 1], daughter_a_location_stencil[loc_idx])));
        }
        for (unsigned loc_idx = 1; loc_idx < daughter_b_location_stencil.size(); loc_idx++)
        {
            cumulative_distances_b.push_back(cumulative_distances_b.back() + norm_2(this->GetVectorFromAtoB(daughter_b_location_stencil[loc_idx - 1], daughter_b_location_stencil[loc_idx])));
        }

        // Find the target node spacing for each of the daughter elements
        double target_spacing_a = cumulative_distances_a.back() / (double)num_nodes;
        double target_spacing_b = cumulative_distances_b.back() / (double)num_nodes;

        // Move the existing nodes into position to become daughter-A nodes
        unsigned last_idx_used = 0;
        for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
        {
            double location_along_arc = (double)node_idx * target_spacing_a;

            while (location_along_arc > cumulative_distances_a[last_idx_used + 1])
            {
                last_idx_used++;
            }

            // Interpolant is the extra distance past the last index used divided by the length of the next line segment
            double interpolant = (location_along_arc - cumulative_distances_a[last_idx_used]) / (cumulative_distances_a[last_idx_used + 1] - cumulative_distances_a[last_idx_used]);

            c_vector<double, SPACE_DIM> this_to_next = this->GetVectorFromAtoB(daughter_a_location_stencil[last_idx_used],
                                                                            daughter_a_location_stencil[last_idx_used + 1]);

            c_vector<double, SPACE_DIM> new_location_a = daughter_a_location_stencil[last_idx_used] + interpolant * this_to_next;

            pElement->GetNode(node_idx)->SetPoint(ChastePoint<SPACE_DIM>(new_location_a));
        }

        // Create new nodes at positions around the daughter-B stencil
        last_idx_used = 0;
        std::vector<Node<SPACE_DIM>*> new_nodes_vec;
        for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
        {
            double location_along_arc = (double)node_idx * target_spacing_b;

            while (location_along_arc > cumulative_distances_b[last_idx_used + 1])
            {
                last_idx_used++;
            }

            // Interpolant is the extra distance past the last index used divided by the length of the next line segment
            double interpolant = (location_along_arc - cumulative_distances_b[last_idx_used]) / (cumulative_distances_b[last_idx_used + 1] - cumulative_distances_b[last_idx_used]);

            c_vector<double, SPACE_DIM> this_to_next = this->GetVectorFromAtoB(daughter_b_location_stencil[last_idx_used],
                                                                            daughter_b_location_stencil[last_idx_used + 1]);

            c_vector<double, SPACE_DIM> new_location_b = daughter_b_location_stencil[last_idx_used] + interpolant * this_to_next;

            unsigned new_node_idx = this->mNodes.size();
            this->mNodes.push_back(new Node<SPACE_DIM>(new_node_idx, new_location_b, true));
            new_nodes_vec.push_back(this->mNodes.back());
        }

        // Copy node attributes
        for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
        {
            new_nodes_vec[node_idx]->SetRegion(pElement->GetNode(node_idx)->GetRegion());

            for (unsigned node_attribute = 0; node_attribute < pElement->GetNode(node_idx)->GetNumNodeAttributes(); node_attribute++)
            {
                new_nodes_vec[node_idx]->AddNodeAttribute(pElement->GetNode(node_idx)->rGetNodeAttributes()[node_attribute]);
            }
        }

        // Create the new element
        unsigned new_elem_idx = this->mElements.size();
        this->mElements.push_back(new ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>(new_elem_idx, new_nodes_vec));
        this->mElements.back()->RegisterWithNodes();

        // Copy any element attributes
        for (unsigned elem_attribute = 0; elem_attribute < pElement->GetNumElementAttributes(); elem_attribute++)
        {
            this->mElements.back()->AddElementAttribute(pElement->rGetElementAttributes()[elem_attribute]);
        }

        // Add the necessary corners to keep consistency with the other daughter element
        for (unsigned corner = 0; corner < pElement->rGetCornerNodes().size(); corner++)
        {
            this->mElements.back()->rGetCornerNodes().push_back(pElement->rGetCornerNodes()[corner]);
        }

        // Update fluid source location for the existing element
        pElement->GetFluidSource()->rGetModifiableLocation() = this->GetCentroidOfElement(pElement->GetIndex());

        // Add a fluid source for the new element
        c_vector<double, SPACE_DIM> new_centroid = this->GetCentroidOfElement(new_elem_idx);
        mElementFluidSources.push_back(std::make_shared<FluidSource<SPACE_DIM>>(new_elem_idx, new_centroid));

        // Set source parameters
        mElementFluidSources.back()->SetAssociatedElementIndex(new_elem_idx);
        mElementFluidSources.back()->SetStrength(0.0);

        // Associate source with element
        mElements[new_elem_idx]->SetFluidSource(mElementFluidSources.back());

        return new_elem_idx;
    }
    else
    {
        NEVER_REACHED;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh(bool randomOrder)
{
    // Iterate over elements and remesh each one
    for (auto p_elem : mElements)
    {
        ReMeshElement(p_elem, randomOrder);
    }

    // Iterate over laminas and remesh each one
    for (auto p_lam : mLaminas)
    {
        ReMeshLamina(p_lam, randomOrder);
    }

    // Reposition fluid sources to the centroid of cells
    for (auto& p_source : mElementFluidSources)
    {
        p_source->rGetModifiableLocation() = this->GetCentroidOfElement(p_source->GetAssociatedElementIndex());
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ReMeshElement(ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                 bool randomOrder)
{
    if constexpr (SPACE_DIM == 2)
    {
        const unsigned num_nodes = pElement->GetNumNodes();

        // Start at a random location in the vector of element nodes if requested
        const unsigned start_idx = randomOrder ? RandomNumberGenerator::Instance()->randMod(num_nodes) : 0u;

        // Straighten out locations to a contiguous polygon rather than wrapped around the domain due to periodicity
        std::vector<c_vector<double, SPACE_DIM>> locations_straightened;
        locations_straightened.reserve(num_nodes);

        locations_straightened.emplace_back(pElement->GetNodeLocation(start_idx));

        for (unsigned node_idx = 1; node_idx < num_nodes; ++node_idx)
        {
            const unsigned prev_idx = AdvanceMod(start_idx, node_idx - 1, num_nodes);
            const unsigned this_idx = AdvanceMod(start_idx, node_idx, num_nodes);

            const c_vector<double, SPACE_DIM>& r_last_location_added = locations_straightened.back();

            const c_vector<double, SPACE_DIM>& r_prev_location = pElement->GetNodeLocation(prev_idx);
            const c_vector<double, SPACE_DIM>& r_this_location = pElement->GetNodeLocation(this_idx);

            locations_straightened.emplace_back(r_last_location_added +
                                                this->GetVectorFromAtoB(r_prev_location, r_this_location));
        }

        assert(locations_straightened.size() == num_nodes);

        const bool closed_path = true;
        const bool permute_order = false;
        const std::size_t num_pts_to_place = num_nodes;

        std::vector<c_vector<double, SPACE_DIM>> evenly_spaced_locations = EvenlySpaceAlongPath(
                locations_straightened,
                closed_path,
                permute_order,
                num_pts_to_place
        );

        assert(evenly_spaced_locations.size() == num_nodes);

        // Conform all locations to geometry
        for (c_vector<double, SPACE_DIM>& r_loc : evenly_spaced_locations)
        {
            ConformToGeometry(r_loc);
        }

        // Update the node locations
        for (unsigned node_idx = 0; node_idx < num_nodes; ++node_idx)
        {
            const unsigned this_idx = AdvanceMod(node_idx, start_idx, num_nodes);
            pElement->GetNode(this_idx)->rGetModifiableLocation() = evenly_spaced_locations[node_idx];
        }
    }
    else
    {
        NEVER_REACHED;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ReMeshLamina(ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>* pLamina,
                                                                bool randomOrder)
{
    if constexpr (SPACE_DIM == 2)
    {
        const unsigned num_nodes = pLamina->GetNumNodes();

        // Start at a random location in the vector of element nodes
        const unsigned start_idx = randomOrder ? RandomNumberGenerator::Instance()->randMod(num_nodes) : 0u;

        // Straighten out locations to a contiguous line rather than wrapped around the domain due to periodicity
        std::vector<c_vector<double, SPACE_DIM>> locations_straightened;
        locations_straightened.reserve(1 + num_nodes);

        locations_straightened.emplace_back(pLamina->GetNodeLocation(start_idx));

        // We go to 1 + num_nodes so that we add on a node in a location congruent to the first location added
        for (unsigned node_idx = 1; node_idx < 1 + num_nodes; ++node_idx)
        {
            const unsigned prev_idx = AdvanceMod(start_idx, node_idx - 1, num_nodes);
            const unsigned this_idx = AdvanceMod(start_idx, node_idx, num_nodes);

            const c_vector<double, SPACE_DIM>& r_last_location_added = locations_straightened.back();

            const c_vector<double, SPACE_DIM>& r_prev_location = pLamina->GetNodeLocation(prev_idx);
            const c_vector<double, SPACE_DIM>& r_this_location = pLamina->GetNodeLocation(this_idx);

            locations_straightened.emplace_back(r_last_location_added +
                                                this->GetVectorFromAtoB(r_prev_location, r_this_location));
        }

        assert(locations_straightened.size() == 1 + num_nodes);

        const bool closed_path = false;
        const bool permute_order = false;
        const std::size_t num_pts_to_place = 1 + num_nodes;

        std::vector<c_vector<double, SPACE_DIM>> evenly_spaced_locations = EvenlySpaceAlongPath(
                locations_straightened,
                closed_path,
                permute_order,
                num_pts_to_place
        );

        assert(evenly_spaced_locations.size() == 1 + num_nodes);

        // Conform all locations to geometry
        for (c_vector<double, SPACE_DIM>& r_loc : evenly_spaced_locations)
        {
            ConformToGeometry(r_loc);
        }

        // Update the node locations, ignoring the very last location that was added to make the path spacing even
        for (unsigned node_idx = 0; node_idx < num_nodes; ++node_idx)
        {
            const unsigned this_idx = AdvanceMod(start_idx, node_idx, num_nodes);
            pLamina->GetNode(this_idx)->rGetModifiableLocation() = evenly_spaced_locations[node_idx];
        }
    }
    else
    {
        NEVER_REACHED;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ConformToGeometry(c_vector<double, SPACE_DIM>& rLocation)
{
    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        assert(rLocation[dim] >= -2.0);  // It's expected that the location is near the domain.  This method is
        assert(rLocation[dim] < 3.0);    // inefficient if the location is far away.

        while (rLocation[dim] < 0.0)
        {
            rLocation[dim] += 1.0;
        }

        while (rLocation[dim] >= 1.0)
        {
            rLocation[dim] -= 1.0;
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::NodesInDifferentElementOrLamina(Node<SPACE_DIM>* pNodeA,
                                                                                   Node<SPACE_DIM>* pNodeB)
{
    // If neither are lamina nodes, we can just check equality of first containing element indices
    if (pNodeA->GetRegion() != LAMINA_REGION && pNodeB->GetRegion() != LAMINA_REGION)
    {
        return *(pNodeA->ContainingElementsBegin()) != *(pNodeB->ContainingElementsBegin());
    }
    else // either one or both nodes is in lamina; assume that two laminas never interact
    {
        return true;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringElementIndices(unsigned elemIdx)
{
    if (!mLaminas.empty())
    {
        EXCEPTION("This method does not yet work in the presence of laminas");
    }

    /*
     * This method uses geometric information from a voronoi diagram of every node.  The voronoi cell of an element
     * (the union of voronoi cells of nodes in the element) provides a precise boundary between elements on which the
     * concept of 'neighbourhood' is well-defined.
     *
     * We first update the voronoi diagram, if it is out of date.
     */
    UpdateNodeLocationsVoronoiDiagramIfOutOfDate();

    std::set<unsigned> neighbouring_element_indices;

    for (unsigned node_local_idx = 0; node_local_idx < this->GetElement(elemIdx)->GetNumNodes(); ++node_local_idx)
    {
        const unsigned node_global_idx = this->GetElement(elemIdx)->GetNodeGlobalIndex(node_local_idx);
        const unsigned voronoi_cell_id = mVoronoiCellIdsIndexedByNodeIndex[node_global_idx];

        /*
         * Iterate over the edges of the voronoi cell corresponding to the current node.  Each primary edge has a twin
         * in a voronoi cell corresponding to a different node.  The element containing that node (which may be the
         * current element under consideration) is added to the set of neighbours precisely if the distance between
         * nodes is less than mNeighbourDist.
         */
        const auto voronoi_cell = mNodeLocationsVoronoiDiagram.cells()[voronoi_cell_id];
        auto p_edge = voronoi_cell.incident_edge();

        do
        {
            // The global node index corresponding to a voronoi cell cell is encoded in its 'color' variable
            const unsigned twin_node_idx = p_edge->twin()->cell()->color();
            const unsigned twin_elem_idx = *this->GetNode(twin_node_idx)->ContainingElementsBegin();

            // Only bother to check the node distances if the nodes are in different elements
            if (twin_elem_idx != elemIdx)
            {
                if (this->GetDistanceBetweenNodes(node_global_idx, twin_node_idx) < mNeighbourDist)
                {
                    neighbouring_element_indices.insert(twin_elem_idx);
                }
            }
            p_edge = p_edge->next();
        } while (p_edge != voronoi_cell.incident_edge());
    }

    return neighbouring_element_indices;
} // LCOV_EXCL_LINE

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::array<unsigned, 13> ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetPolygonDistribution()
{
    if constexpr (SPACE_DIM == 2)
    {
        std::array<unsigned, 13> polygon_dist = {{0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u}};

        for (auto elem_it = this->GetElementIteratorBegin(); elem_it != this->GetElementIteratorEnd(); ++elem_it)
        {
            if (!elem_it->IsElementOnBoundary())
            {
                // Accumulate all 12+ sided shapes
                unsigned num_neighbours = std::min<unsigned>(12u, GetNeighbouringElementIndices(elem_it->GetIndex()).size());
                polygon_dist[num_neighbours]++;
            }
        }

        return polygon_dist;
    }
    else
    {
        NEVER_REACHED;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::CalculateLengthOfVoronoiEdge(const boost::polygon::voronoi_diagram<double>::edge_type& rEdge)
{
    assert(rEdge.is_finite());

    const double d_x = rEdge.vertex1()->x() - rEdge.vertex0()->x();
    const double d_y = rEdge.vertex1()->y() - rEdge.vertex0()->y();

    return ScaleDistanceDownFromVoronoi(std::sqrt(d_x * d_x + d_y * d_y));
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetMaxNodeIndex() const
{
    if (this->mNodes.begin() == this->mNodes.end())
    {
        return UINT_MAX;
    }
    else
    {
        const auto max_idx_it = std::max_element(this->mNodes.begin(), this->mNodes.end(),
                                                 [](Node<SPACE_DIM>* const a, Node<SPACE_DIM>* const b)
                                                 {
                                                     return a->GetIndex() < b->GetIndex();
                                                 });

        return (*max_idx_it)->GetIndex();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetMaxElementIndex() const
{
    if (this->mElements.begin() == this->mElements.end())
    {
        return UINT_MAX;
    }
    else
    {
        using IbElem = ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>;

        const auto max_idx_it = std::max_element(this->mElements.begin(), this->mElements.end(),
                                                 [](IbElem* const a, IbElem* const b)
                                                 {
                                                     return a->GetIndex() < b->GetIndex();
                                                 });

        return (*max_idx_it)->GetIndex();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetMaxLaminaIndex() const
{
    if (this->mLaminas.begin() == this->mLaminas.end())
    {
        return UINT_MAX;
    }
    else
    {
        using IbLam = ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>;

        const auto max_idx_it = std::max_element(this->mLaminas.begin(), this->mLaminas.end(),
                                                 [](IbLam* const a, IbLam* const b)
                                                 {
                                                     return a->GetIndex() < b->GetIndex();
                                                 });

        return (*max_idx_it)->GetIndex();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::TagBoundaryElements()
{
    if (!mLaminas.empty())
    {
        EXCEPTION("This method does not yet work in the presence of laminas");
    }

    // Need a large double, but since we also perform calculations with it, we can't really use DBL_MAX
    const double LARGE_DOUBLE = 1e6;

    /*
     * This method uses geometric information from a voronoi diagram of every node.
     *
     * Nodes with infinite voronoi cells must be in a boundary element.
     * In addition, some nodes in boundary elements will have unusually long shared-edges, and the voronoi perimeter
     * of boundary elements will exceed the true perimeter by an unusually large amount.
     * We use a threshold on the median values of these identifiers to robustly determine boundary elements.
     */
    UpdateNodeLocationsVoronoiDiagramIfOutOfDate();

    // Get the maximum element index in the mesh
    const unsigned max_elem_idx = GetMaxElementIndex();

    // Keep track of all finite shared edge lengths (between two elements), as well as (for each element) the maximum
    // shared edge length.
    std::vector<double> max_shared_lengths(1 + max_elem_idx, 0.0);
    std::vector<double> voronoi_perimeter(1 + max_elem_idx, 0.0);

    // Loop over nodes and calculate values to fill the containers defined above
    for (const auto& p_node : this->mNodes)
    {
        const unsigned this_node_idx = p_node->GetIndex();
        const unsigned this_elem_idx = *p_node->ContainingElementsBegin();

        // Get the voronoi cell that corresponds to this node
        const unsigned voronoi_cell_id = mVoronoiCellIdsIndexedByNodeIndex[this_node_idx];
        const auto& voronoi_cell = mNodeLocationsVoronoiDiagram.cells()[voronoi_cell_id];

        // Iterate over the edges of this cell.  Note that due to the halo region none of these edges should be infinite.
        auto p_edge = voronoi_cell.incident_edge();

        // Loop over the voronoi edges associated with this node's voronoi cell
        do
        {
            // If the edge is infinite, the element containing it must be on the boundary
            if (p_edge->is_infinite())
            {
                max_shared_lengths[this_elem_idx] = LARGE_DOUBLE;
                voronoi_perimeter[this_elem_idx] += LARGE_DOUBLE;
            }
            else
            {
                // The global node index corresponding to a voronoi cell cell is encoded in its 'color' variable
                const unsigned twin_node_idx = p_edge->twin()->cell()->color();
                const unsigned twin_elem_idx = *this->GetNode(twin_node_idx)->ContainingElementsBegin();

                // If the edge is between two elements, it counts
                if (this_elem_idx != twin_elem_idx)
                {
                    const double edge_length = CalculateLengthOfVoronoiEdge(*p_edge);
                    max_shared_lengths[this_elem_idx] = std::max(max_shared_lengths[this_elem_idx], edge_length);
                    voronoi_perimeter[this_elem_idx] += edge_length;
                }
            }
            p_edge = p_edge->next();
        } while (p_edge != voronoi_cell.incident_edge());
    }

    // Calculate what proportion bigger the voronoi perimeter is than the actual one
    std::vector<double> perimeter_multiples(voronoi_perimeter.size());
    for (const auto& p_elem : this->mElements)
    {
        const unsigned idx = p_elem->GetIndex();
        perimeter_multiples[idx] = voronoi_perimeter[idx] / this->GetSurfaceAreaOfElement(idx);
    }

    // Calculate the medians of each vector.  Copy the vector accounting for possible non-contiguous element indices.
    std::vector<double> copy_max_shared_lengths;
    std::vector<double> copy_perimeter_multiples;

    for (const auto& p_elem : this->mElements)
    {
        const unsigned idx = p_elem->GetIndex();
        copy_max_shared_lengths.emplace_back(max_shared_lengths[idx]);
        copy_perimeter_multiples.emplace_back(perimeter_multiples[idx]);
    }

    const std::size_t half_way = copy_max_shared_lengths.size() / 2;
    assert(half_way == copy_perimeter_multiples.size() / 2);

    std::nth_element(copy_max_shared_lengths.begin(), copy_max_shared_lengths.begin() + half_way, copy_max_shared_lengths.end());
    std::nth_element(copy_perimeter_multiples.begin(), copy_perimeter_multiples.begin() + half_way, copy_perimeter_multiples.end());

    double median_max_edge_length = copy_max_shared_lengths[half_way];
    double median_perimeter_multiple = copy_perimeter_multiples[half_way];

    // In the event that most elements have infinite edges, for instance if there are a small number of elements,
    // a reduced median can be used
    if (median_max_edge_length == LARGE_DOUBLE || median_perimeter_multiple >= LARGE_DOUBLE)
    {
        median_max_edge_length = 0.5 * LARGE_DOUBLE;
        median_perimeter_multiple = 0.5 * LARGE_DOUBLE;
    }

    // Finally, tag the boundary elements
    for (const auto& p_elem : this->mElements)
    {
        const unsigned idx = p_elem->GetIndex();

        const bool large_shared_edge = max_shared_lengths[idx] > 1.1 * median_max_edge_length;
        const bool large_perimeter = perimeter_multiples[idx] > 1.1 * median_perimeter_multiple;

        p_elem->SetIsBoundaryElement(large_shared_edge && large_perimeter);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::UpdateNodeLocationsVoronoiDiagramIfOutOfDate()
{
    if constexpr (SPACE_DIM == 2)
    {
        using boost_point = boost::polygon::point_data<int>;

        /*
        * The voronoi diagram needs updating only if the node locations have changed (i.e. not more than once per
        * timestep).  We test a chosen summary of node locations against the cached mSummaryOfNodeLocations to determine
        * this.
        */
        double new_location_summary = this->mNodes.front()->rGetLocation()[0] +
                                    this->mNodes.front()->rGetLocation()[1] +
                                    this->mNodes.back()->rGetLocation()[0] +
                                    this->mNodes.back()->rGetLocation()[1];

        bool voronoi_needs_updating = std::fabs(mSummaryOfNodeLocations - new_location_summary) > DBL_EPSILON;

        /*
        * Method outline:
        * - Add all node locations as boost points, correctly scaled
        * - Add extra locations for any nodes within distance mVoronoiHalo of the domain edge for periodicity
        * - Calculate the voronoi diagram
        * - Populate mVoronoiCellIdsIndexedByNodeIndex for efficient mapping between node index and corresponding voronoi
        *   cell
        * - Store the corresponding node index in each voronoi cell's "color" variable for efficient reverse-lookup
        * - Tag boundary elements
        */
        if (voronoi_needs_updating)
        {
            mSummaryOfNodeLocations = new_location_summary;

            // Helper c_vectors for halo region
            const c_vector<double, SPACE_DIM> halo_up = Create_c_vector(0.0, 1.0);
            const c_vector<double, SPACE_DIM> halo_down = Create_c_vector(0.0, -1.0);
            const c_vector<double, SPACE_DIM> halo_left = Create_c_vector(-1.0, 0.0);
            const c_vector<double, SPACE_DIM> halo_right = Create_c_vector(1.0, 0.0);

            std::vector<std::pair<unsigned, c_vector<double, SPACE_DIM>>> halo_ids_and_locations;
            std::vector<unsigned> node_ids_in_source_idx_order;

            // We need to translate node locations into boost points, which are scaled to an integer grid
            std::vector<boost_point> points;

            // First add a point for each node, and calculate any additional locations needed within the halo
            for (const auto& p_node : this->mNodes)
            {
                const double x_pos = p_node->rGetLocation()[0];
                const double y_pos = p_node->rGetLocation()[1];

                // Put the integer voronoi coordinate in for this node's location
                points.emplace_back(boost_point(ScaleUpToVoronoiCoordinate(x_pos), ScaleUpToVoronoiCoordinate(y_pos)));
                node_ids_in_source_idx_order.emplace_back(p_node->GetIndex());

                // Now, for this location, decide whether any copies of this location are needed in the halo region
                const bool needed_up = y_pos < mVoronoiHalo;
                const bool needed_down = y_pos > 1.0 - mVoronoiHalo;
                const bool needed_left = x_pos > 1.0 - mVoronoiHalo;
                const bool needed_right = x_pos < mVoronoiHalo;

                if (needed_up)
                {
                    halo_ids_and_locations.emplace_back(std::make_pair(p_node->GetIndex(),
                                                                    p_node->rGetLocation() + halo_up));
                }
                if (needed_down)
                {
                    halo_ids_and_locations.emplace_back(std::make_pair(p_node->GetIndex(),
                                                                    p_node->rGetLocation() + halo_down));
                }
                if (needed_left)
                {
                    halo_ids_and_locations.emplace_back(std::make_pair(p_node->GetIndex(),
                                                                    p_node->rGetLocation() + halo_left));
                }
                if (needed_right)
                {
                    halo_ids_and_locations.emplace_back(std::make_pair(p_node->GetIndex(),
                                                                    p_node->rGetLocation() + halo_right));
                }
                if (needed_up && needed_left)
                {
                    halo_ids_and_locations.emplace_back(std::make_pair(p_node->GetIndex(),
                                                                    p_node->rGetLocation() + halo_up + halo_left));
                }
                if (needed_up && needed_right)
                {
                    halo_ids_and_locations.emplace_back(std::make_pair(p_node->GetIndex(),
                                                                    p_node->rGetLocation() + halo_up + halo_right));
                }
                if (needed_down && needed_left)
                {
                    halo_ids_and_locations.emplace_back(std::make_pair(p_node->GetIndex(),
                                                                    p_node->rGetLocation() + halo_down + halo_left));
                }
                if (needed_down && needed_right)
                {
                    halo_ids_and_locations.emplace_back(std::make_pair(p_node->GetIndex(),
                                                                    p_node->rGetLocation() + halo_down + halo_right));
                }
            }

            // Next add the additional points
            for (const auto& pair : halo_ids_and_locations)
            {
                const unsigned node_idx = pair.first;
                c_vector<double, SPACE_DIM> location = pair.second;

                const int x_coord = ScaleUpToVoronoiCoordinate(location[0]);
                const int y_coord = ScaleUpToVoronoiCoordinate(location[1]);

                points.emplace_back(boost_point(x_coord, y_coord));
                node_ids_in_source_idx_order.emplace_back(node_idx);
            }

            // Construct the voronoi diagram.  This is the costly part of this method.
            mNodeLocationsVoronoiDiagram.clear();
            construct_voronoi(std::begin(points), std::end(points), &mNodeLocationsVoronoiDiagram);

            // We need an efficient map from node global index to the voronoi cell representing it, and we can't assume
            // that the nodes are ordered sequentially.  We first identify the largest node index.
            const unsigned max_node_idx = GetMaxNodeIndex();

            mVoronoiCellIdsIndexedByNodeIndex.resize(1 + max_node_idx);

            for (unsigned vor_cell_id = 0; vor_cell_id < mNodeLocationsVoronoiDiagram.cells().size(); ++vor_cell_id)
            {
                // Source index is incrementally given to each input point, which is in order of nodes in this->mNodes.
                // We need to be able to identify each vor_cell_id by the global node index

                auto& r_this_cell = mNodeLocationsVoronoiDiagram.cells()[vor_cell_id];
                const auto source_idx = r_this_cell.source_index();
                const unsigned node_idx = node_ids_in_source_idx_order[source_idx];

                r_this_cell.color(node_idx);

                if (source_idx < this->mNodes.size())
                {
                    mVoronoiCellIdsIndexedByNodeIndex[node_idx] = vor_cell_id;
                }
            }

            // Finally, if the diagram was out-of-date, we will need to re-tag boundary elements.
            this->TagBoundaryElements();
        }
    }
    else
    {
        NEVER_REACHED;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
int ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ScaleUpToVoronoiCoordinate(const double location) const
{
    assert(location >= -mVoronoiHalo);
    assert(location <= 1.0 + mVoronoiHalo);

    constexpr auto DBL_INT_MIN = static_cast<double>(INT_MIN);
    constexpr auto DBL_INT_RANGE = static_cast<double>(INT_MAX) - static_cast<double>(INT_MIN);

    // To handle periodicity, we add (again) any nodes within a halo region of the unit square
    constexpr double scale_factor = DBL_INT_RANGE / (1.0 + 2.0 * mVoronoiHalo);

    // By construction there is no narrowing, so this static_cast might not be necessary
    return static_cast<int>(std::lround(DBL_INT_MIN + (mVoronoiHalo + location) * scale_factor));
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ScaleDistanceDownFromVoronoi(const double distance) const
{
    constexpr auto DBL_INT_RANGE = static_cast<double>(INT_MAX) - static_cast<double>(INT_MIN);
    constexpr double scale_factor = (1.0 + 2.0 * mVoronoiHalo) / DBL_INT_RANGE;

    return scale_factor * distance;
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
const boost::polygon::voronoi_diagram<double>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGetNodeLocationsVoronoiDiagram(bool update)
{
    if (update)
    {
        this->UpdateNodeLocationsVoronoiDiagramIfOutOfDate();
    }

    return mNodeLocationsVoronoiDiagram;
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
const std::vector<unsigned int>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetVoronoiCellIdsIndexedByNodeIndex() const
{
    return mVoronoiCellIdsIndexedByNodeIndex;
}

// Explicit instantiation
template class ImmersedBoundaryMesh<1, 1>;
template class ImmersedBoundaryMesh<1, 2>;
template class ImmersedBoundaryMesh<1, 3>;
template class ImmersedBoundaryMesh<2, 2>;
template class ImmersedBoundaryMesh<2, 3>;
template class ImmersedBoundaryMesh<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ImmersedBoundaryMesh)
