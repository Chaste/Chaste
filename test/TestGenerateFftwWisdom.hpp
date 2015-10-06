/*

Copyright (c) 2005-2014, University of Oxford.
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

// Needed for test framework
#include <cxxtest/TestSuite.h>

// Needed for Immersed Boundary simulations
#include <complex>
#include <fftw3.h>
#include <boost/multi_array.hpp>

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestGenerateFftwWisdom : public CxxTest::TestSuite
{
public:

    void TestGenerateWisdom() throw(Exception)
    {
        /*
         * This test generates an fftw wisdom file telling fftw how to efficiently compute fourier transforms of a
         * given size.  We generate wisdom for:
         *    * 2d forward and backward complex-to-complex transforms (16x16 --> 4096x4096)
         *    * 3d forward and backward complex-to-complex transforms (16x16x16 --> 256x256x256)
         *
         * This test takes a LONG time to run if there is currently no wisdom (around 4 hours).
         */

        std::string filename = "./projects/ImmersedBoundary/src/fftw.wisdom";
        int wisdom_flag = fftw_import_wisdom_from_filename(filename.c_str());

        // 1 means it's read correctly, 0 indicates a failure
        TS_ASSERT_EQUALS(wisdom_flag, 1);

        // Create a 3D array that is 64 x 64 x 64
        typedef boost::multi_array<std::complex<double>, 2> complex_array_2d;
        typedef boost::multi_array<std::complex<double>, 3> complex_array_3d;

        // Create 2D wisdom
        for (unsigned i = 16 ; i < 5000 ; i*=2)
        {
            complex_array_2d input(boost::extents[i][i]);
            complex_array_2d output(boost::extents[i][i]);

            fftw_complex* fftw_input = reinterpret_cast<fftw_complex*>(input.data());
            fftw_complex* fftw_output = reinterpret_cast<fftw_complex*>(output.data());

            fftw_plan plan_f;
            plan_f = fftw_plan_dft_2d(i, i, fftw_input, fftw_output, FFTW_FORWARD, FFTW_EXHAUSTIVE);

            fftw_plan plan_b;
            plan_b = fftw_plan_dft_2d(i, i, fftw_input, fftw_output, FFTW_BACKWARD, FFTW_EXHAUSTIVE);
        }

        // Create 3D wisdom
        for (unsigned i = 16 ; i < 257 ; i*=2)
        {
            complex_array_3d input(boost::extents[i][i][i]);
            complex_array_3d output(boost::extents[i][i][i]);

            fftw_complex* fftw_input = reinterpret_cast<fftw_complex*>(input.data());
            fftw_complex* fftw_output = reinterpret_cast<fftw_complex*>(output.data());

            fftw_plan plan_f;
            plan_f = fftw_plan_dft_3d(i, i, i, fftw_input, fftw_output, FFTW_FORWARD, FFTW_EXHAUSTIVE);

            fftw_plan plan_b;
            plan_b = fftw_plan_dft_3d(i, i, i, fftw_input, fftw_output, FFTW_BACKWARD, FFTW_EXHAUSTIVE);
        }

        fftw_export_wisdom_to_filename(filename.c_str());
    }
};