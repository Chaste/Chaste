/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTGENERATEANDCACHERANDOMFIELD_HPP_
#define TESTGENERATEANDCACHERANDOMFIELD_HPP_

#include <cxxtest/TestSuite.h>

#include <array>

#include "UniformGridRandomFieldGenerator.hpp"

// These tests do not run in parallel
#include "FakePetscSetup.hpp"

class TestGenerateAndCacheRandomField : public CxxTest::TestSuite
{
public:

    /**
     * Change the variables here as required.  If the field is big this will take a long time and may produce a large
     * output file.
     *
     * For instance, a 2D grid with 128x128 using 1000 eigenvalues will generate a file ~ 130MB in size.
     *
     * You should only ever run this test in anger in Release mode, or it will take forever to generate the field.
     */
    void TestGenerateAndCache()
    {
        // Define the domain over which the field is calculated
        const std::array<double, 2> lower_corner = {{0.0, 1.0}};
        const std::array<double, 2> upper_corner = {{2.0, 5.0}};

        // Define the number of grid points required in each dimension
        const std::array<unsigned, 2> num_grid_pts = {{10, 12}};

        // Define whether there is periodicity in any dimension
        const std::array<bool, 2> periodicity = {{false, false}};

        // Define the number of eigenvalues to calculate
        const unsigned num_evals = 100u;

        // Define the correlation length for the noise in the random field
        const double length_scale = 0.5;

        // Generate and cache the random field
        UniformGridRandomFieldGenerator<2> gen(lower_corner, upper_corner, num_grid_pts, periodicity, num_evals, length_scale);
        gen.SaveToCache();
    }
};

#endif /*TESTGENERATEANDCACHERANDOMFIELD_HPP_*/