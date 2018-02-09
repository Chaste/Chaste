/*

Copyright (c) 2005-2018, University of Oxford.
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

#ifndef TESTSAMPLESFROMCACHEDRANDOMFIELD_HPP_
#define TESTSAMPLESFROMCACHEDRANDOMFIELD_HPP_

#include <cxxtest/TestSuite.h>

#include <array>
#include <string>

#include "OutputFileHandler.hpp"
#include "UniformGridRandomFieldGenerator.hpp"

// These tests do not run in parallel
#include "FakePetscSetup.hpp"

class TestSamplesFromCachedRandomField : public CxxTest::TestSuite
{
public:

    /**
     * Construct a small field so that there is always a cached field in $CHASTE_TEST_OUTPUT directory
     */
    void TestConstructSmallField()
    {
        // Create a random field over a 2D rectangle
        const std::array<double, 2> lower_corner = {{0.0, 1.0}};
        const std::array<double, 2> upper_corner = {{2.0, 5.0}};

        // 10 * 12 = 120 grid points (requires calculating the eigenvalues of a 120 * 120 matrix)
        const std::array<unsigned, 2> num_grid_pts = {{10, 12}};

        // No periodicity
        const std::array<bool, 2> periodicity = {{false, false}};

        // Calculate most of the eigenvalues
        const double trace_proportion = 0.8;

        // Fairly large lengthscale compared to the domain size
        const double length_scale = 0.5;

        // Generate and cache the random field
        UniformGridRandomFieldGenerator<2> gen(lower_corner, upper_corner, num_grid_pts, periodicity, trace_proportion, length_scale);
        gen.SaveToCache();
    }

    /**
     * Sample from your cached random field to test whether the numbers are well distributed.
     *
     * Run the python file PROJECT/apps/TestHistogramOfRandomFieldSamples.py after changing the relevant file names
     * in this test and in the python file.
     *
     * If the resulting histogram does not fit well within the blue dotted line, you probably need to calculate more
     * eigenvalues when you generate the random field.
     */
    void TestOuputSamplesToFile()
    {
        // Change this vaiable to your cached random field (relative to $CHASTE_TEST_OUTPUT)
        const std::string file_path = "CachedRandomFields/xy_0.000_1.000_2.000_5.000_10_12_0_0_0.800_0.500.rfg";

        // Change this variable to the number of instances of the random field you want to sample
        const unsigned num_fields_to_sample = 100u;

        // Generate random field from cache and sample from it
        UniformGridRandomFieldGenerator<2> gen(file_path);

        OutputFileHandler results_handler("", false);
        out_stream results_file = results_handler.OpenOutputFile(file_path + ".check");

        // Get a number of random field instances and dump the values to file
        for (unsigned i = 0; i < num_fields_to_sample; ++i)
        {
            const std::vector<double> grf = gen.SampleRandomField();
            for (const double& var : grf)
            {
                (*results_file) << var << '\n';
            }
        }

        // Tidy up
        results_file->close();
    }
};

#endif /*TESTSAMPLESFROMCACHEDRANDOMFIELD_HPP_*/