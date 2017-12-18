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

#ifndef TESTUNIFORMGRIDRANDOMFIELDGENERATOR_HPP_
#define TESTUNIFORMGRIDRANDOMFIELDGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"
#include "UniformGridRandomFieldGenerator.hpp"

#include "Debug.hpp"
// These tests do not run in parallel
#include "FakePetscSetup.hpp"

class TestUniformGridRandomFieldGenerator : public CxxTest::TestSuite
{
public:

    void TestConstructor()
    {
        const unsigned n = 12u;

        UniformGridRandomFieldGenerator<2> gen({{0.0, 0.0}}, {{1.0, 1.0}}, {{n, n}}, {{true, true}}, 0.8, 0.03);
        gen.SaveToCache();

        OutputFileHandler results_handler(".", false);
        out_stream results_file = results_handler.OpenOutputFile("a_hist.csv");

        auto p_rng = RandomNumberGenerator::Instance();

        for (unsigned i = 0; i < 500; ++i)
        {
            auto grf = gen.SampleRandomField();

            // Generate some random points
            for (unsigned j = 0; j < 100; ++j)
            {
                c_vector<double, 2> vec = Create_c_vector(p_rng->ranf(), p_rng->ranf());
                (*results_file) << gen.Interpolate(grf, vec) << '\n';
            }
        }

        // Tidy up
        results_file->close();
    }
};

#endif /*TESTUNIFORMGRIDRANDOMFIELDGENERATOR_HPP_*/