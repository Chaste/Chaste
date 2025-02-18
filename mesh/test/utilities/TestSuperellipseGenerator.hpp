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

#ifndef TESTSUPERELLIPSEGENERATOR_HPP_
#define TESTSUPERELLIPSEGENERATOR_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>

#include "SuperellipseGenerator.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestSuperellipseGenerator : public CxxTest::TestSuite
{
public:

    void TestCreateReferenceSuperellipse()
    {
        // Generate an ellipse (superellipse with exponent 1.0) that should sit inside [0, 0.6]x[0, 0.8]
        unsigned num_points = 100;
        double   exponent = 1.0;
        double   width = 0.6;
        double   height = 0.8;
        double   bot_left_x = 0.0;
        double   bot_left_y = 0.0;

        SuperellipseGenerator gen(num_points, exponent, width, height, bot_left_x, bot_left_y);
        std::vector<c_vector<double, 2> > points = gen.GetPointsAsVectors();

        double min_x = DBL_MAX;
        double min_y = DBL_MAX;
        double max_x = -DBL_MAX;
        double max_y = -DBL_MAX;

        // Test all points are within the unit square, and find the four extreme positions
        for (unsigned idx = 0; idx < points.size(); idx++)
        {
            double x_pos = points[idx][0];
            double y_pos = points[idx][1];

            TS_ASSERT(x_pos >= 0.0);
            TS_ASSERT(x_pos <= 0.6);
            TS_ASSERT(y_pos >= 0.0);
            TS_ASSERT(y_pos <= 0.8);

            if (x_pos > max_x)
            {
                max_x = x_pos;
            }
            else if (x_pos < min_x)
            {
                min_x = x_pos;
            }

            if (y_pos > max_y)
            {
                max_y = y_pos;
            }
            else if (y_pos < min_y)
            {
                min_y = y_pos;
            }
        }

        // Test location of extreme points
        TS_ASSERT_DELTA(min_x, 0.0, 1e-6);
        TS_ASSERT_DELTA(max_x, 0.6, 1e-6);
        TS_ASSERT_DELTA(min_y, 0.0, 1e-6);
        TS_ASSERT_DELTA(max_y, 0.8, 1e-6);

        // Test points are roughly uniformly spaced
        double previous_dist = norm_2(points[0] - points[1]);
        for (unsigned idx = 1; idx < points.size(); idx++)
        {
            double this_dist = norm_2(points[idx] - points[(idx + 1) % points.size()]);

            TS_ASSERT_DELTA(previous_dist, this_dist, 1e-6);

            previous_dist = this_dist;
        }
    }

    void TestGetTargetNodeSpacing()
    {
        /*
         * The target node spacing is the ellipse arclength divided by the number of points parameterising the ellipse.
         * In general, it is not easy to calculate the arclength analytically, so an approximation is used.  We can test
         * the method by trying special cases (circle and square).
         */

        unsigned num_points = 123;
        double   width = 1.0;
        double   height = 1.0;
        double   bot_left_x = 0.0;
        double   bot_left_y = 0.0;

        double circle_exponent = 1.0;
        double square_exponent = 0.0001;

        // The circle arclength should be approximately PI
        SuperellipseGenerator gen_circle(num_points, circle_exponent, width, height, bot_left_x, bot_left_y);
        TS_ASSERT_DELTA(gen_circle.GetTargetNodeSpacing(), M_PI / double(num_points), 1e-6);

        // The square arclength should be approximately 4.0
        SuperellipseGenerator gen_square(num_points, square_exponent, width, height, bot_left_x, bot_left_y);
        TS_ASSERT_DELTA(gen_square.GetTargetNodeSpacing(), 4.0 / double(num_points), 1e-6);
    }

    void TestGetHeightOfTopSurface()
    {
        {
            // Generate an ellipse (superellipse with exponent 1.0) that should sit inside [0, 0.6]x[0, 0.8]
            unsigned num_points = 100;
            double   exponent = 1.0;
            double   width = 0.4;
            double   height = 0.8;
            double   bot_left_x = 0.0;
            double   bot_left_y = 0.0;

            SuperellipseGenerator gen(num_points, exponent, width, height, bot_left_x, bot_left_y);

            TS_ASSERT_DELTA(gen.GetHeightOfTopSurface(), 0.8 * height, 1e-6);
        }

        //\todo implement this test
        {
            // case 2: exponent < 1
        }

    }

    void TestGetPointsAsChastePoints()
    {
        unsigned num_points = 123;
        double   exponent = 4.56;
        double   width = 0.789;
        double   height = 0.123;
        double   bot_left_x = 4.56;
        double   bot_left_y = 7.89;

        SuperellipseGenerator gen(num_points, exponent, width, height, bot_left_x, bot_left_y);

        // Get vectors
        std::vector<c_vector<double, 2> > vector_points = gen.GetPointsAsVectors();

        // Get Chaste points
        std::vector<ChastePoint<2> > chaste_points = gen.GetPointsAsChastePoints();

        // Check vectors are the same length
        TS_ASSERT_EQUALS(vector_points.size(), chaste_points.size());

        // Check locations coincide
        for (unsigned idx = 1; idx < vector_points.size(); idx++)
        {
            TS_ASSERT_DELTA(vector_points[idx][0], chaste_points[idx][0], 1e-6);
            TS_ASSERT_DELTA(vector_points[idx][1], chaste_points[idx][1], 1e-6);
        }
    }
};

#endif /*TESTSUPERELLIPSEGENERATOR_HPP_*/