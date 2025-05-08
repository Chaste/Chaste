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

#include "SuperellipseGenerator.hpp"

SuperellipseGenerator::SuperellipseGenerator(unsigned numPoints,
                                             double ellipseExponent,
                                             double width,
                                             double height,
                                             double botLeftX,
                                             double botLeftY)
        : mTargetNodeSpacing(DOUBLE_UNSET),
          mHeightOfTopSurface(0.3 * height)
{
    // Validate input
    assert(numPoints > 1);
    assert(ellipseExponent > 0.0);
    assert(width > 0.0);
    assert(height > 0.0);

    // Set up member variables
    mPoints.resize(numPoints);

    /*
     * Run a high-density pass around the parametric curve first
     */
    unsigned dense_pts = 50000;

    // Vectors to store information from the loop
    std::vector<double> cumulative_arc_length(dense_pts);
    std::vector<c_vector<double, 2> > dense_locations(dense_pts);

    // Helper variables for the loop
    double dense_spacing = (2.0 * M_PI) / (double)dense_pts;
    double cumulative_dist = 0.0;
    double temp_sin, temp_cos;

    // Fill in first location by hand
    cumulative_arc_length[0] = 0.0;
    dense_locations[0][0] = width * 0.5;
    dense_locations[0][1] = 0.0;

    // Fill in all other locations
    for (unsigned point = 1; point < dense_pts; point++)
    {
        // Update cumulative distance
        cumulative_dist += dense_spacing;

        // Get the current sin and cos values
        temp_cos = cos(cumulative_dist);
        temp_sin = sin(cumulative_dist);

        // x and y are powers of cos and sin, with the sign of sin and cos
        dense_locations[point][0] = copysign(pow(fabs(temp_cos), ellipseExponent), temp_cos) * width * 0.5;
        dense_locations[point][1] = copysign(pow(fabs(temp_sin), ellipseExponent), temp_sin) * height * 0.5;

        // Check that the new point created lies within [-width/2, width/2] x [-height/2, height/2]
        assert(fabs(dense_locations[point][0]) < width * 0.5 + 1e-10);
        assert(fabs(dense_locations[point][1]) < height * 0.5 + 1e-10);

        // Fill in arc length property
        cumulative_arc_length[point] = cumulative_arc_length[point - 1] + norm_2(dense_locations[point] - dense_locations[point - 1]);
    }

    double total_arc_length = cumulative_arc_length[dense_pts - 1] + norm_2(dense_locations[dense_pts - 1] - dense_locations[0]);

    // Since our perimeter is periodic, we get back to the start
    cumulative_arc_length.push_back(total_arc_length);
    dense_locations.push_back(dense_locations[0]);

    /*
     * Decide on the best place to put the actual boundary points
     */

    // Helper variables for loop
    unsigned dense_it = 0;

    mTargetNodeSpacing = total_arc_length / double(numPoints);

    double target_arclength = 0.0;
    double interpolant;

    // Fill in first point by hand
    mPoints[0] = dense_locations[0];

    // Fill in all other locations
    for (unsigned point = 1; point < numPoints; point++)
    {
        target_arclength = double(point) * mTargetNodeSpacing;

        while (cumulative_arc_length[dense_it] < target_arclength && dense_it < dense_pts)
        {
            dense_it ++;
        }

        // Check we're one position either side of the target arclength
        assert(cumulative_arc_length[dense_it - 1] < target_arclength);
        assert(cumulative_arc_length[dense_it] >= target_arclength);

        // Find the proportion between the two known locations we are
        interpolant = (target_arclength - cumulative_arc_length[dense_it - 1]) / (cumulative_arc_length[dense_it] - cumulative_arc_length[dense_it - 1]);

        // Calculate the location of the new point
        mPoints[point] = (1.0 - interpolant) * dense_locations[dense_it - 1] + interpolant * dense_locations[dense_it];
    }

    /*
     * If the exponent is at least 1, we return the default value which is the y-coord with largest value.
     *
     * If the exponent is less than 1, the top surface is defined as the height at which maximum curvature is attained.
     * We loop through the points looking at the change in x and change in y from the previous point.  The first
     * difference in which the absolute change in x is greater than the absolute change in y will be the point with
     * maximal curvature, for a superellipse.
     */
    if (ellipseExponent < 1.0)
    {
        double delta_x = fabs(mPoints[0][0] - mPoints[1][0]);
        double delta_y = fabs(mPoints[0][1] - mPoints[1][1]);

        // The following should always hold, for an exponent less than 1.0
        if (delta_x > delta_y)
        {
            NEVER_REACHED;
        }

        for (unsigned i = 1; i < numPoints; i++)
        {
            delta_x = fabs(mPoints[i][0] - mPoints[i - 1][0]);
            delta_y = fabs(mPoints[i][1] - mPoints[i - 1][1]);

            if (delta_x > delta_y)
            {
                mHeightOfTopSurface = 0.5 * (mPoints[i][1] + mPoints[i - 1][1]);

                break;
            }

            // We should meet this condition before being half way through the list
            if (i == unsigned(numPoints / 2.0))
            {
                NEVER_REACHED;
            }
        }
    }

    /*
     * Rescale all points to match parameters
     */

    // Move perimeter from [-width/2, width/2] x [-height/2, height/2] to correct location
    c_vector<double, 2> offset;
    offset[0] = width * 0.5 + botLeftX;
    offset[1] = height * 0.5 + botLeftY;

    // Reposition the height of the top surface
    mHeightOfTopSurface += offset[1];

    for (unsigned point = 0; point < numPoints; point++)
    {
        // Reposition all points
        mPoints[point] += offset;

        // Check we're in the right place
        assert((mPoints[point][0] > botLeftX - 1e-10) && (mPoints[point][0] < botLeftX + width + 1e-10));
        assert((mPoints[point][1] > botLeftY - 1e-10) && (mPoints[point][1] < botLeftY + height + 1e-10));
    }
}

SuperellipseGenerator::~SuperellipseGenerator()
{
    mPoints.clear();
}

double SuperellipseGenerator::GetTargetNodeSpacing()
{
    return mTargetNodeSpacing;
}

double SuperellipseGenerator::GetHeightOfTopSurface()
{
    return mHeightOfTopSurface;
}

const std::vector<c_vector<double, 2> > SuperellipseGenerator::GetPointsAsVectors() const
{
    return mPoints;
}

const std::vector<ChastePoint<2> > SuperellipseGenerator::GetPointsAsChastePoints() const
{
    std::vector<ChastePoint<2> > chaste_points;

    for (unsigned point = 0; point < mPoints.size(); point++)
    {
        chaste_points.push_back(ChastePoint<2>(mPoints[point]));
    }

    return chaste_points;
}
