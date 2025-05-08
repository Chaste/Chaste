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

#ifndef TESTMESHUTILITYFUNCTIONS_HPP_
#define TESTMESHUTILITYFUNCTIONS_HPP_

#include <cxxtest/TestSuite.h>

#include "MeshUtilityFunctions.hpp"
#include "UblasCustomFunctions.hpp"

#include <algorithm>

#include "PetscSetupAndFinalize.hpp"


class TestMeshUtilityFunctions : public CxxTest::TestSuite
{
private:
    const bool OPEN_PATH = false;
    const bool CLOSED_PATH = true;

    const bool KEEP_ORDER = false;
    const bool PERMUTE_ORDER = true;

public:

    void TestEvenlySpaceAlongPath1d()
    {
        // Simple open path with number of points
        {
            std::vector<c_vector<double, 1>> path = {
                    Create_c_vector(0.0),
                    Create_c_vector(1.0)
            };

            const unsigned num_points = 3u;

            std::vector<c_vector<double, 1>> output = EvenlySpaceAlongPath(path, OPEN_PATH, KEEP_ORDER, num_points);

            TS_ASSERT_EQUALS(output.size(), num_points);
            TS_ASSERT_DELTA(norm_inf(output[0] - Create_c_vector(0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[1] - Create_c_vector(0.5)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[2] - Create_c_vector(1.0)), 0.0, 1e-6);
        }

        // Multi-point open path with number of points
        {
            std::vector<c_vector<double, 1>> path = {
                    Create_c_vector(0.0),
                    Create_c_vector(0.1),
                    Create_c_vector(0.2),
                    Create_c_vector(0.6),
                    Create_c_vector(1.0)
            };

            const unsigned num_points = 3u;

            std::vector<c_vector<double, 1>> output = EvenlySpaceAlongPath(path, OPEN_PATH, KEEP_ORDER, num_points);

            TS_ASSERT_EQUALS(output.size(), num_points);
            TS_ASSERT_DELTA(norm_inf(output[0] - Create_c_vector(0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[1] - Create_c_vector(0.5)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[2] - Create_c_vector(1.0)), 0.0, 1e-6);
        }

        // Simple open path with target spacing
        {
            std::vector<c_vector<double, 1>> path = {
                    Create_c_vector(0.0),
                    Create_c_vector(1.0)
            };

            const double target_spacing = 0.5;

            std::vector<c_vector<double, 1>> output = EvenlySpaceAlongPath(path, OPEN_PATH, KEEP_ORDER, 0u, target_spacing);

            TS_ASSERT_EQUALS(output.size(), 3u);
            TS_ASSERT_DELTA(norm_inf(output[0] - Create_c_vector(0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[1] - Create_c_vector(0.5)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[2] - Create_c_vector(1.0)), 0.0, 1e-6);
        }

        // Simple closed path with number of points
        {
            std::vector<c_vector<double, 1>> path = {
                    Create_c_vector(0.0),
                    Create_c_vector(1.0)
            };

            const unsigned num_points = 3u;

            std::vector<c_vector<double, 1>> output = EvenlySpaceAlongPath(path, CLOSED_PATH, KEEP_ORDER, num_points);

            TS_ASSERT_EQUALS(output.size(), 3u);
            TS_ASSERT_DELTA(norm_inf(output[0] - Create_c_vector(0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[1] - Create_c_vector(2.0/3.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[2] - Create_c_vector(2.0/3.0)), 0.0, 1e-6);
        }

        // Simple closed path with target spacing
        {
            std::vector<c_vector<double, 1>> path = {
                    Create_c_vector(0.0),
                    Create_c_vector(1.0),
                    Create_c_vector(0.8),
                    Create_c_vector(0.79)
            };

            const double target_spacing = 0.2857; // spacing still 0.25 as 7 output points will have spacing too big

            std::vector<c_vector<double, 1>> output = EvenlySpaceAlongPath(path, CLOSED_PATH, KEEP_ORDER, 0u, target_spacing);

            TS_ASSERT_EQUALS(output.size(), 8u);
            TS_ASSERT_DELTA(norm_inf(output[0] - Create_c_vector(0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[1] - Create_c_vector(0.25)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[2] - Create_c_vector(0.5)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[3] - Create_c_vector(0.75)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[4] - Create_c_vector(1.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[5] - Create_c_vector(0.75)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[6] - Create_c_vector(0.5)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[7] - Create_c_vector(0.25)), 0.0, 1e-6);
        }
    }

    void TestEvenlySpaceAlongPath2d()
    {
        // Open square with 5 points
        {
            std::vector<c_vector<double, 2>> path = {
                    Create_c_vector(0.0, 0.0),
                    Create_c_vector(1.0, 0.0),
                    Create_c_vector(1.0, 1.0),
                    Create_c_vector(0.0, 1.0),
                    Create_c_vector(0.0, 0.0)
            };

            const unsigned num_points = 5;

            std::vector<c_vector<double, 2>> output = EvenlySpaceAlongPath(path, OPEN_PATH, KEEP_ORDER, num_points);

            TS_ASSERT_EQUALS(output.size(), num_points);
            TS_ASSERT_DELTA(norm_inf(output[0] - Create_c_vector(0.0, 0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[1] - Create_c_vector(1.0, 0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[2] - Create_c_vector(1.0, 1.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[3] - Create_c_vector(0.0, 1.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[4] - Create_c_vector(0.0, 0.0)), 0.0, 1e-6);
        }

        // Closed square with 4 points
        {
            std::vector<c_vector<double, 2>> path = {
                    Create_c_vector(0.0, 0.0),
                    Create_c_vector(1.0, 0.0),
                    Create_c_vector(1.0, 1.0),
                    Create_c_vector(0.0, 1.0)
            };

            const unsigned num_points = 4;

            std::vector<c_vector<double, 2>> output = EvenlySpaceAlongPath(path, CLOSED_PATH, KEEP_ORDER, num_points);

            TS_ASSERT_EQUALS(output.size(), num_points);
            TS_ASSERT_DELTA(norm_inf(output[0] - Create_c_vector(0.0, 0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[1] - Create_c_vector(1.0, 0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[2] - Create_c_vector(1.0, 1.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[3] - Create_c_vector(0.0, 1.0)), 0.0, 1e-6);
        }
    }

    void TestEvenlySpaceAlongPath3d()
    {
        // Open cube sides with 5 points
        {
            std::vector<c_vector<double, 3>> path = {
                    Create_c_vector(0.0, 0.0, 0.0),
                    Create_c_vector(0.0, 1.0, 0.0),
                    Create_c_vector(0.0, 1.0, 1.0),
                    Create_c_vector(1.0, 1.0, 1.0)
            };

            const unsigned num_points = 5;

            std::vector<c_vector<double, 3>> output = EvenlySpaceAlongPath(path, OPEN_PATH, KEEP_ORDER, num_points);

            TS_ASSERT_EQUALS(output.size(), num_points);
            TS_ASSERT_DELTA(norm_inf(output[0] - Create_c_vector(0.0, 0.0, 0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[1] - Create_c_vector(0.0, 0.75, 0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[2] - Create_c_vector(0.0, 1.0, 0.5)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[3] - Create_c_vector(0.25, 1.0, 1.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[4] - Create_c_vector(1.0, 1.0, 1.0)), 0.0, 1e-6);
        }

        // Closed curve with permutation
        {
            std::vector<c_vector<double, 3>> path = {
                    Create_c_vector(0.0, 0.0, 0.0),
                    Create_c_vector(0.0, 1.0, 0.0)
            };

            const unsigned num_points = 4;

            std::vector<c_vector<double, 3>> output = EvenlySpaceAlongPath(path, CLOSED_PATH, PERMUTE_ORDER, num_points);

            TS_ASSERT_EQUALS(output.size(), num_points);

            // One output should have y=0: find this one and re-rotate.  The output should now be in the original order.
            auto pivot_it = std::find_if(output.begin(), output.end(), [](const c_vector<double, 3>& a){return std::fabs(a[1]) < 1e-6;});
            std::rotate(output.begin(), pivot_it, output.end());

            TS_ASSERT_DELTA(norm_inf(output[0] - Create_c_vector(0.0, 0.0, 0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[1] - Create_c_vector(0.0, 0.5, 0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[2] - Create_c_vector(0.0, 1.0, 0.0)), 0.0, 1e-6);
            TS_ASSERT_DELTA(norm_inf(output[3] - Create_c_vector(0.0, 0.5, 0.0)), 0.0, 1e-6);
        }
    }
};

#endif /*TESTMESHUTILITYFUNCTIONS_HPP_*/
