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

#ifndef TESTFFTWMETHODS_HPP_
#define TESTFFTWMETHODS_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>

// Chaste includes
#include "FileFinder.hpp"
#include "RandomNumberGenerator.hpp"
#include "Warnings.hpp"

// Needed for immersed boundary simulations
#include <boost/multi_array.hpp>
#include <complex>
#include <fftw3.h>

// This test is never run in parallel
#include "FakePetscSetup.hpp"

typedef boost::multi_array<std::complex<double>, 2> complex_array_2d;
typedef boost::multi_array<std::complex<double>, 3> complex_array_3d;

typedef boost::multi_array<double, 2> real_array_2d;
typedef boost::multi_array<double, 3> real_array_3d;

class TestFftwMethods : public CxxTest::TestSuite
{
public:
    void TestLoadWisdom()
    {
        std::string wisdom_filename = "fftw.wisdom";
        FileFinder file_finder(wisdom_filename, RelativeTo::ChasteTestOutput);

        // Verify the file can be found in its default location
        std::string wisdom_path = file_finder.GetAbsolutePath();

        if (file_finder.IsFile())
        {
            // 1 means it's read correctly, 0 indicates a failure
            int wisdom_flag = fftw_import_wisdom_from_filename(wisdom_path.c_str());
            TS_ASSERT_EQUALS(wisdom_flag, 1);
        }
        else
        {
            WARNING("fftw wisdom file not found; simulation setup may take a long time.")
        }
    }

    void Test2DR2C2R()
    {
        /*
         * We create a 2D vector of doubles and carry out a real-to-complex transform, followed by a complex-to-real
         * transform, and test that we get back to the same result.
         */

        unsigned size = 1024;
        unsigned reduced = 1 + (size / 2);

        double norm = double(size * size);

        // Create the 2D arrays we need, each of size 1024x1024
        real_array_2d original_data(boost::extents[size][size]);
        real_array_2d real_input(boost::extents[size][size]);
        real_array_2d real_output(boost::extents[size][size]);
        complex_array_2d complex_data(boost::extents[size][reduced]);

        double* p_real_input = real_input.data();
        double* p_real_output = real_output.data();
        fftw_complex* p_complex_data = reinterpret_cast<fftw_complex*>(complex_data.data());

        fftw_plan plan_f;
        plan_f = fftw_plan_dft_r2c_2d(size, size, p_real_input, p_complex_data, FFTW_ESTIMATE);

        fftw_plan plan_b;
        plan_b = fftw_plan_dft_c2r_2d(size, size, p_complex_data, p_real_output, FFTW_ESTIMATE);

        for (unsigned i = 0; i < size; i++)
        {
            for (unsigned j = 0; j < size; j++)
            {
                original_data[i][j] = RandomNumberGenerator::Instance()->ranf();
                real_input[i][j] = original_data[i][j];
            }
        }

        fftw_execute(plan_f);
        fftw_destroy_plan(plan_f);

        fftw_execute(plan_b);
        fftw_destroy_plan(plan_b);

        for (unsigned i = 0; i < size; i++)
        {
            for (unsigned j = 0; j < size; j++)
            {
                TS_ASSERT_DELTA(original_data[i][j], real_output[i][j] / norm, 1e-15);
            }
        }
    }
};

#endif /*TESTFFTWMETHODS_HPP_*/