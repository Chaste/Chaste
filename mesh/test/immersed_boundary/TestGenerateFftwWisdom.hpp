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

#ifndef TESTGENERATEFFTWWISDOM_HPP_
#define TESTGENERATEFFTWWISDOM_HPP_

// Needed for test framework
#include <cxxtest/cxxtest/TestSuite.h>

// Needed for immersed boundary simulations
#include <complex>
#include <fftw3.h>
#include <boost/multi_array.hpp>

#include "ImmersedBoundaryArray.hpp"
#include "RandomNumberGenerator.hpp"
#include "FileFinder.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestGenerateFftwWisdom : public CxxTest::TestSuite
{
private:
    /** Name of wisdom file when using a single thread */
    std::string mWisdomFilename;

    /** The max array size N for N by N arrays. Must be a multiple of 2 */
    static const unsigned mMaxArraySize = 4096;

    /** Flag to determine length of time to be taken testing potentially faster algorithms */
    unsigned mFftwFlag;

public:

    void TestSetupWisdomFileOneThread()
    {
        std::string single_thread_file_name = "fftw.wisdom";

        // Set up the file finder and get the absolute path
        FileFinder file_finder(single_thread_file_name, RelativeTo::ChasteTestOutput);
        mWisdomFilename = file_finder.GetAbsolutePath();

        // If it doesn't exists, create it with blank wisdom file
        if (!file_finder.IsFile())
        {
            void fftw_forget_wisdom();
            fftw_export_wisdom_to_filename(mWisdomFilename.c_str());
        }

        // The file should now definitely exist
        TS_ASSERT(file_finder.IsFile());

        mFftwFlag = FFTW_PATIENT;    // Faster, slightly less thorough
        //mFftwFlag = FFTW_EXHAUSTIVE; // Much slower, but guaranteed to produce fastest transforms
    }

    void TestGenerateManyR2CWisdomOneThread()
    {
        /*
         * This test generates an fftw wisdom file telling fftw how to efficiently compute fourier transforms of a
         * given size, using two threads.  We generate wisdom for two and three DFTs of data contiguous in memory:
         *
         *    * 2d forward R2C and backward C2R transforms (16x16 --> mMaxArraySize x mMaxArraySize)
         */

        // We first forget all wisdom and re-load, as threaded wisdom doesn't play well with un-threaded
        void fftw_forget_wisdom();
        int wisdom_flag = fftw_import_wisdom_from_filename(mWisdomFilename.c_str());

        // 1 means it's read correctly, 0 indicates a failure
        TS_ASSERT_EQUALS(wisdom_flag, 1);

        // Create 3D arrays that will represent two 2D arrays
        typedef boost::multi_array<std::complex<double>, 3> complex_array_2d;
        typedef boost::multi_array<double, 3> real_array_2d;

        // Create 2D wisdom with 2 transforms
        for (unsigned i = 16; i <= mMaxArraySize; i*=2)
        {
            real_array_2d input(boost::extents[2][i][i]);
            complex_array_2d output(boost::extents[2][i][(i/2) + 1]);
            real_array_2d check(boost::extents[2][i][i]);

            double* fftw_input = &input[0][0][0];
            fftw_complex* fftw_output = reinterpret_cast<fftw_complex*>(&output[0][0][0]);
            double* fftw_check = &check[0][0][0];

            // Cast i as int (stop compiler warnings)
            int j = int(i);

            // Plan variables
            int rank = 2;                   // Number of dimensions for each array
            int real_dims[] = {j, j};       // Dimensions of each real array
            int comp_dims[] = {j, 1 + j/2}; // Dimensions of each complex array
            int how_many = 2;               // Number of transforms
            int real_sep = j * j;           // How many doubles between start of first array and start of second
            int comp_sep = j * (1 + j/2);   // How many fftw_complex between start of first array and start of second
            int real_stride = 1;            // Each real array is contiguous in memory
            int comp_stride = 1;            // Each complex array is contiguous in memory
            int* real_nembed = real_dims;
            int* comp_nembed = comp_dims;

            fftw_plan plan_f;
            plan_f = fftw_plan_many_dft_r2c(rank, real_dims, how_many,
                                            fftw_input,  real_nembed, real_stride, real_sep,
                                            fftw_output, comp_nembed, comp_stride, comp_sep,
                                            mFftwFlag);

            fftw_plan plan_b;
            plan_b = fftw_plan_many_dft_c2r(rank, real_dims, how_many,
                                            fftw_output, comp_nembed, comp_stride, comp_sep,
                                            fftw_check,  real_nembed, real_stride, real_sep,
                                            mFftwFlag);

            // We now verify that the forward followed by inverse transform produces the correct result
            for (unsigned dim = 0; dim < 2; dim++)
            {
                for (unsigned x = 0; x < i; x++)
                {
                    for (unsigned y = 0; y < i; y++)
                    {
                        input[dim][x][y] = RandomNumberGenerator::Instance()->ranf();
                    }
                }
            }

            fftw_execute(plan_f);
            fftw_execute(plan_b);

            for (unsigned dim = 0; dim < 2; dim++)
            {
                for (unsigned x = 0; x < i; x++)
                {
                    for (unsigned y = 0; y < i; y++)
                    {
                        TS_ASSERT_DELTA(input[dim][x][y], check[dim][x][y]/(i*i), 1e-10);
                    }
                }
            }

            // Export each step of the loop so if the test is stopped progress is kept
            fftw_export_wisdom_to_filename(mWisdomFilename.c_str());
            std::cout << "Exported 1-thread wisdom for two contiguous arrays of size " << i << " by " << i << std::endl;
        }

        // Create 2D wisdom with 3 transforms
        for (unsigned i = 16; i <= mMaxArraySize; i*=2)
        {
            real_array_2d input(boost::extents[3][i][i]);
            complex_array_2d output(boost::extents[3][i][(i/2) + 1]);
            real_array_2d check(boost::extents[3][i][i]);

            double* fftw_input = &input[0][0][0];
            fftw_complex* fftw_output = reinterpret_cast<fftw_complex*>(&output[0][0][0]);
            double* fftw_check = &check[0][0][0];

            // Cast i as int (stop compiler warnings)
            int j = int(i);

            // Plan variables
            int rank = 2;                   // Number of dimensions for each array
            int real_dims[] = {j, j};       // Dimensions of each real array
            int comp_dims[] = {j, 1 + j/2}; // Dimensions of each complex array
            int how_many = 3;               // Number of transforms
            int real_sep = j * j;           // How many doubles between start of first array and start of second
            int comp_sep = j * (1 + j/2);   // How many fftw_complex between start of first array and start of second
            int real_stride = 1;            // Each real array is contiguous in memory
            int comp_stride = 1;            // Each complex array is contiguous in memory
            int* real_nembed = real_dims;
            int* comp_nembed = comp_dims;

            fftw_plan plan_f;
            plan_f = fftw_plan_many_dft_r2c(rank, real_dims, how_many,
                                            fftw_input,  real_nembed, real_stride, real_sep,
                                            fftw_output, comp_nembed, comp_stride, comp_sep,
                                            mFftwFlag);

            fftw_plan plan_b;
            plan_b = fftw_plan_many_dft_c2r(rank, real_dims, how_many,
                                            fftw_output, comp_nembed, comp_stride, comp_sep,
                                            fftw_check,  real_nembed, real_stride, real_sep,
                                            mFftwFlag);

            // We now verify that the forward followed by inverse transform produces the correct result
            for (unsigned dim = 0; dim < 3; dim++)
            {
                for (unsigned x = 0; x < i; x++)
                {
                    for (unsigned y = 0; y < i; y++)
                    {
                        input[dim][x][y] = RandomNumberGenerator::Instance()->ranf();
                    }
                }
            }

            fftw_execute(plan_f);
            fftw_execute(plan_b);

            for (unsigned dim = 0; dim < 3; dim++)
            {
                for (unsigned x = 0; x < i; x++)
                {
                    for (unsigned y = 0; y < i; y++)
                    {
                        TS_ASSERT_DELTA(input[dim][x][y], check[dim][x][y]/(i*i), 1e-10);
                    }
                }
            }

            // Export each step of the loop so if the test is stopped progress is kept
            fftw_export_wisdom_to_filename(mWisdomFilename.c_str());
            std::cout << "Exported 1-thread wisdom for three contiguous arrays of size " << i << " by " << i << std::endl;
        }
    }
};

#endif /*TESTGENERATEFFTWWISDOM_HPP_*/