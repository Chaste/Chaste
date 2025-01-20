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


#ifndef TESTUNIFORMGRIDRANDOMFIELDGENERATOR_HPP_
#define TESTUNIFORMGRIDRANDOMFIELDGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <type_traits>

#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"
#include "UniformGridRandomFieldGenerator.hpp"

// These tests do not run in parallel
//#include "FakePetscSetup.hpp"

class TestUniformGridRandomFieldGenerator : public CxxTest::TestSuite
{
private:

    template <typename T, unsigned DIM>
    void TestEqualityOfArrays(std::array<T, DIM> a, std::array<T, DIM> b)
    {
        if (std::is_integral<T>::value)
        {
            for (unsigned dim = 0; dim < DIM; ++dim)
            {
                TS_ASSERT_EQUALS(a[dim], b[dim]);
            }
        }
        else if (std::is_floating_point<T>::value)
        {
            for (unsigned dim = 0; dim < DIM; ++dim)
            {
                TS_ASSERT_DELTA(a[dim], b[dim], 1e-6);
            }
        }
        else
        {
            NEVER_REACHED;
        }
    }

public:

    void TestMainConstructor()
    {
        const unsigned n = 12u;
        UniformGridRandomFieldGenerator<2> gen({{0.0, 0.0}}, {{1.0, 1.0}}, {{n, n}}, {{true, true}}, 0.03);

        OutputFileHandler results_handler("", true);
        out_stream results_file = results_handler.OpenOutputFile("a_hist.csv");

        auto p_rng = RandomNumberGenerator::Instance();

        auto grf = gen.SampleRandomField();
        for (unsigned i = 0; i < 500; ++i)
        {

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

    void TestSampleRandomField()
    {
        // Correct sample numbers
        {
            // 1D
            {
                const unsigned n = 12u;
                UniformGridRandomFieldGenerator<1> gen({{0.0}}, {{1.0}}, {{n}}, {{true}}, 0.03);
                auto samples = gen.SampleRandomField();
                TS_ASSERT_EQUALS(samples.size(), n);
            }

            // 2D
            {
                const unsigned n = 12u;
                UniformGridRandomFieldGenerator<2> gen({{0.0, 0.0}}, {{1.0, 1.0}}, {{n, n}}, {{true, true}}, 0.03);
                auto samples = gen.SampleRandomField();
                TS_ASSERT_EQUALS(samples.size(), n * n);
            }

            // 3D
            {
                const unsigned n = 12u;
                UniformGridRandomFieldGenerator<3> gen({{0.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0}}, {{n, n, n}}, {{true, true, true}}, 0.03);
                auto samples = gen.SampleRandomField();
                TS_ASSERT_EQUALS(samples.size(), n * n * n);
            }
        }

        // Length scale
        {
            const unsigned n = 4;
            UniformGridRandomFieldGenerator<1> gen1({{0.0}}, {{1.0}}, {{n}}, {{true}}, 1.00);
            gen1.SetRandomSeed(123);
            auto samples1 = gen1.SampleRandomFieldAtTime(0.0);
            UniformGridRandomFieldGenerator<1> gen2({{0.0}}, {{1.0}}, {{n}}, {{true}}, 2.00);
            gen2.SetRandomSeed(123);
            auto samples2 = gen2.SampleRandomFieldAtTime(0.0);

            TS_ASSERT_EQUALS(samples1[0], samples2[0]);
            TS_ASSERT_EQUALS(samples1[2], samples2[1]);
        }
    }

    void TestGetLinearIndex()
    {
        // 1D
        {
            UniformGridRandomFieldGenerator<1> gen({{0.0}}, {{1.0}}, {{8}}, {{true}}, 0.3);

            TS_ASSERT_EQUALS(gen.GetLinearIndex({0l}), 0l);
            TS_ASSERT_EQUALS(gen.GetLinearIndex({15l}), 15l);
            TS_ASSERT_EQUALS(gen.GetLinearIndex({-50l}), -50l);
        }

        // 2D
        {
            UniformGridRandomFieldGenerator<2> gen({{0.0, 0.0}}, {{1.0, 1.0}}, {{4, 4}}, {{true, true}}, 0.3);

            TS_ASSERT_EQUALS(gen.GetLinearIndex({0l, 0l}), 0l);
            TS_ASSERT_EQUALS(gen.GetLinearIndex({15l, 0l}), 15l);
            TS_ASSERT_EQUALS(gen.GetLinearIndex({-50l, 0l}), -50l);
            TS_ASSERT_EQUALS(gen.GetLinearIndex({0l, 1l}), 4l);
            TS_ASSERT_EQUALS(gen.GetLinearIndex({3l, 2l}), 11l);
        }

        // 3D
        {
            UniformGridRandomFieldGenerator<3> gen({{0.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0}}, {{4, 4, 4}}, {{true, true, true}}, 0.3);

            TS_ASSERT_EQUALS(gen.GetLinearIndex({0l, 0l, 0l}), 0l);
            TS_ASSERT_EQUALS(gen.GetLinearIndex({15l, 0l, 0l}), 15l);
            TS_ASSERT_EQUALS(gen.GetLinearIndex({-50l, 0l, 0l}), -50l);
            TS_ASSERT_EQUALS(gen.GetLinearIndex({0l, 1l, 0l}), 4l);
            TS_ASSERT_EQUALS(gen.GetLinearIndex({3l, 2l, 0l}), 11l);
            TS_ASSERT_EQUALS(gen.GetLinearIndex({0l, 1l, 2l}), 36l);
            TS_ASSERT_EQUALS(gen.GetLinearIndex({3l, 2l, 3l}), 59l);
        }
    }

    void TestGetPositionUsingGridIndex()
    {
        // 1D
        {
            UniformGridRandomFieldGenerator<1> gen({{0.0}}, {{1.0}}, {{8}}, {{true}}, 0.3);

            TestEqualityOfArrays<double, 1>(gen.GetPositionUsingGridIndex({0l}), {{0.0}});
            TestEqualityOfArrays<double, 1>(gen.GetPositionUsingGridIndex({3l}), {{3.0 / 8.0}});
            TestEqualityOfArrays<double, 1>(gen.GetPositionUsingGridIndex({6l}), {{6.0 / 8.0}});
        }

        // 2D
        {
            UniformGridRandomFieldGenerator<2> gen({{0.0, 0.0}}, {{1.0, 1.0}}, {{4, 4}}, {{true, true}}, 0.3);

            TestEqualityOfArrays<double, 2>(gen.GetPositionUsingGridIndex({0l, 0l}), {{0.0, 0.0}});
            TestEqualityOfArrays<double, 2>(gen.GetPositionUsingGridIndex({1l, 2l}), {{0.25, 0.5}});
            TestEqualityOfArrays<double, 2>(gen.GetPositionUsingGridIndex({3l, 3l}), {{0.75, 0.75}});
        }

        // 3D
        {
            UniformGridRandomFieldGenerator<3> gen({{0.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0}}, {{4, 4, 4}}, {{true, true, true}}, 0.3);

            TestEqualityOfArrays<double, 3>(gen.GetPositionUsingGridIndex({0l, 0l, 0l}), {{0.0, 0.0, 0.0}});
            TestEqualityOfArrays<double, 3>(gen.GetPositionUsingGridIndex({1l, 2l, 3l}), {{0.25, 0.5, 0.75}});
            TestEqualityOfArrays<double, 3>(gen.GetPositionUsingGridIndex({3l, 1l, 0l}), {{0.75, 0.25, 0.0}});
        }
    }

    void TestInterpolate()
    {
        { // lower left < 0

            UniformGridRandomFieldGenerator<1> gen({{0.0}}, {{1.0}}, {{8}}, {{true}}, 0.3);

            // If the field is constant everywhere then that constant should always be the interpolated value
            {
                std::vector<double> constant_field(8, 1.23);

                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(-0.1)), 1.23, 1e-6);
            }
        }
        // 1D
        {
            UniformGridRandomFieldGenerator<1> gen({{0.0}}, {{1.0}}, {{8}}, {{true}}, 0.3);

            // If the field is constant everywhere then that constant should always be the interpolated value
            {
                std::vector<double> constant_field(8, 1.23);

                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.0)), 1.23, 1e-6);
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.543)), 1.23, 1e-6);
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.897)), 1.23, 1e-6);
            }

            // If the field is not constant everywhere, interpolation is linear
            {
                std::vector<double> field(8, 1.0);
                field[4] = 2.0;

                TS_ASSERT_DELTA(gen.Interpolate(field, Create_c_vector(3.0 / 8.0)), 1.0, 1e-6); // below
                TS_ASSERT_DELTA(gen.Interpolate(field, Create_c_vector(4.0 / 8.0)), 2.0, 1e-6); // exactly at the point we changed
                TS_ASSERT_DELTA(gen.Interpolate(field, Create_c_vector(5.0 / 8.0)), 1.0, 1e-6); // above

                TS_ASSERT_DELTA(gen.Interpolate(field, Create_c_vector(3.25 / 8.0)), 1.25, 1e-6); // quarter
                TS_ASSERT_DELTA(gen.Interpolate(field, Create_c_vector(3.5 / 8.0)), 1.5, 1e-6); // half
                TS_ASSERT_DELTA(gen.Interpolate(field, Create_c_vector(4.25 / 8.0)), 1.75, 1e-6); // one + quarter
            }
        }

        // 2D
        {
            UniformGridRandomFieldGenerator<2> gen({{0.0, 0.0}}, {{1.0, 1.0}}, {{4, 4}}, {{true, true}}, 0.3);

            // If the field is constant everywhere then that constant should always be the interpolated value
            {
                std::vector<double> constant_field(16, 2.34);

                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.0, 0.0)), 2.34, 1e-6);
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.543, 0.143)), 2.34, 1e-6);
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.897, 0.453)), 2.34, 1e-6);
            }

            // If the field is not constant everywhere, interpolation is bilinear
            {
                std::vector<double> constant_field(16, 1.0);

                // The four corners of this square, and any point in the square, should have value 2.0
                constant_field[5] = 2.0;
                constant_field[6] = 2.0;
                constant_field[9] = 2.0;
                constant_field[10] = 2.0;

                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.25, 0.25)), 2.0, 1e-6); // corner
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.25, 0.5)), 2.0, 1e-6); // corner
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.5, 0.25)), 2.0, 1e-6); // corner
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.5, 0.5)), 2.0, 1e-6); // corner
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.3, 0.4)), 2.0, 1e-6); // inside

                // Interpolation should be linear along the edges of the square, with a quadratic surface inside the square
                constant_field[5] = 2.0;
                constant_field[6] = 3.0;
                constant_field[9] = 6.0;
                constant_field[10] = 9.0;

                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.3, 0.25)), 2.2, 1e-6); // bottom side
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.25, 0.4)), 4.4, 1e-6); // left side
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.5, 0.35)), 5.4, 1e-6); // right side
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.375, 0.5)), 7.5, 1e-6); // top side


                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.375, 0.375)), 5.0, 1e-6); // inside
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.4, 0.3)), 3.64, 1e-6); // inside
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.3, 0.4)), 4.84, 1e-6); // inside
            }
        }

        // 3D
        {
            UniformGridRandomFieldGenerator<3> gen({{0.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0}}, {{4, 4, 4}}, {{true, true, true}}, 0.3);

            // If the field is constant everywhere then that constant should always be the interpolated value
            {
                std::vector<double> constant_field(64, 3.45);

                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.0, 0.0, 0.0)), 3.45, 1e-6);
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.543, 0.143, 0.0)), 3.45, 1e-6);
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.897, 0.453, 0.876)), 3.45, 1e-6);
            }

            // If the field is not constant everywhere, interpolation is a simple nearest-neighbour
            {
                std::vector<double> constant_field(64, 3.45);
                constant_field[21] = 4.56;

                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.0, 0.0, 0.0)), 3.45, 1e-6); // Lower corner of supported region
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.25, 0.25, 0.25)), 4.56, 1e-6); // Center of supported region
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.5, 0.5, 0.5)), 3.45, 1e-6); // top corner of supported region
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.125, 0.25, 0.25)), 4.005, 1e-6); // mid 4.56 region
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.25, 0.25, 0.125)), 4.005, 1e-6); // mid 4.56 region
                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(0.25, 0.375, 0.25)), 4.005, 1e-6); // mid 4.56 region

                TS_ASSERT_DELTA(gen.Interpolate(constant_field, Create_c_vector(1.0, 0.75, 1.0)), 3.45, 1e-6); // outside 4.56 region (above)
            }
        }
    }
};

#endif /*TESTUNIFORMGRIDRANDOMFIELDGENERATOR_HPP_*/