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


#ifndef TESTOFFLATTICERANDOMFIELDGENERATOR_HPP_
#define TESTOFFLATTICERANDOMFIELDGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>

#include "Node.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"
#include "OffLatticeRandomFieldGenerator.hpp"

// These tests do not run in parallel
#include <sys/stat.h>

#include "FakePetscSetup.hpp"

class TestOffLatticeRandomFieldGenerator : public CxxTest::TestSuite
{
private:
    std::tuple<double, double> CalculateMeanAndVariance(const std::vector<double>& vec)
    {
        const double f_size = static_cast<double>(vec.size());

        // Calculate mean (using std::reduce for a modern C++ approach)
        const double mean = std::reduce(vec.begin(), vec.end(), 0.0) / f_size;

        // Calculate variance using std::transform_reduce
        const double variance = std::transform_reduce(
                                    vec.begin(), vec.end(), // Range of the first sequence
                                    vec.begin(), // Start of the second sequence, reused vec as dummy to match transform_reduce signature
                                    0.0, // Initial value for the reduction
                                    std::plus<>(), // Reduction operation (summing up)
                                    [mean](double val, double /*unused*/)
                                    { return (val - mean) * (val - mean); } // Transformation
                                    )
            / f_size;

        return std::make_tuple(mean, variance);
    }

    double CalculateMeanSquaredError(const std::vector<double>& vecA, const std::vector<double>& vecB)
    {
        return std::transform_reduce(
                   vecA.begin(), vecA.end(), vecB.begin(), 0.0,
                   std::plus<>(), // Reduction operation (summing up)
                   [](double a, double b)
                   { return (a - b) * (a - b); } // Transformation operation (squared difference)
                   )
            / static_cast<double>(vecA.size());
    }

public:

    void TestStatsUtilFunctions()
    {
        const std::vector<double> a = {1.2, 2.3, 4.0};
        const std::vector<double> b = {3.4, 5.6, 8.7};

        auto [mean_a, var_a] = CalculateMeanAndVariance(a);
        auto [mean_b, var_b] = CalculateMeanAndVariance(b);

        const double mse_aa = CalculateMeanSquaredError(a, a);
        const double mse_bb = CalculateMeanSquaredError(b, b);
        const double mse_ab = CalculateMeanSquaredError(a, b);
        const double mse_ba = CalculateMeanSquaredError(b, a);

        TS_ASSERT_DELTA(mean_a, 2.5, 1e-15);
        TS_ASSERT_DELTA(mean_b, 5.9, 1e-15);

        TS_ASSERT_DELTA(var_a, 3.98 / 3.0, 1e-15);
        TS_ASSERT_DELTA(var_b, 14.18 / 3.0, 1e-15);

        TS_ASSERT_DELTA(mse_aa, 0.0, 1e-15);
        TS_ASSERT_DELTA(mse_bb, 0.0, 1e-15);
        TS_ASSERT_DELTA(mse_ab, mse_ba, 1e-15);
        TS_ASSERT_DELTA(mse_ab, ((3.4 - 1.2) * (3.4 - 1.2) + (5.6 - 2.3) * (5.6 - 2.3) + (8.7 - 4.0) * (8.7 - 4.0)) / 3.0, 1e-15);
    }

    void TestConstructor()
    {
        const std::array<double, 2> lower_corner {{0.0, 0.0}};
        const std::array<double, 2> upper_corner {{10.0, 10.0}};
        const std::array<bool, 2> periodicity {{false, false}};
        const double lengthscale = 0.1;

        OffLatticeRandomFieldGenerator<2> gen(
                lower_corner,
                upper_corner,
                periodicity,
                lengthscale
        );
    }

    void TestSetRandomSeed()
    {
        const std::array<double, 1> lower_corner {{0.0}};
        const std::array<double, 1> upper_corner {{10.0}};
        const std::array<bool, 1> periodicity {{false}};
        const double lengthscale = 0.1;

        OffLatticeRandomFieldGenerator<1> gen(
                lower_corner,
                upper_corner,
                periodicity,
                lengthscale
        );

        auto p_gen = RandomNumberGenerator::Instance();
        std::vector<Node<1>*> nodes(10);
        for (unsigned node_idx = 0; node_idx < 10; ++node_idx)
        {
            nodes[node_idx] = new Node<1>(node_idx, Create_c_vector(10.0 * p_gen->ranf()));
        }

        gen.SetRandomSeed(10);
        auto res1 = gen.SampleRandomFieldAtTime(nodes, 0.0);
        gen.SetRandomSeed(11);
        auto res2 = gen.SampleRandomFieldAtTime(nodes, 0.0);

        TS_ASSERT_DIFFERS(res1[0], res2[0]);

        for (auto& p_node : nodes)
        {
            delete p_node;
        }
    }

    void TestSampleFromRandomField()
    {
        { // Without specifying time
            auto p_gen = RandomNumberGenerator::Instance();

            const unsigned n = 1'000;

            const std::array<double, 1> lower_corner {{0.0}};
            const std::array<double, 1> upper_corner {{10.0}};
            const std::array<bool, 1> periodicity {{false}};
            const double lengthscale = 1.0;

            OffLatticeRandomFieldGenerator<1> gen(
                    lower_corner,
                    upper_corner,
                    periodicity,
                    lengthscale
            );
            gen.SetRandomSeed(10);

            // Generate some nodes
            std::vector<Node<1>*> nodes(n);
            for (unsigned node_idx = 0; node_idx < n; ++node_idx)
            {
                // Nodes spaced out substantially, so that the random field valeus are very uncorellated
                nodes[node_idx] = new Node<1>(node_idx, Create_c_vector(1e4 * p_gen->ranf()));
            }

            auto random_field_a = gen.SampleRandomField(nodes);
            auto random_field_b = gen.SampleRandomField(nodes);
            auto [mean_a, var_a] = CalculateMeanAndVariance(random_field_a);
            auto [mean_b, var_b] = CalculateMeanAndVariance(random_field_b);
            const double mse = CalculateMeanSquaredError(random_field_a, random_field_b);

            // The means should be close to zero
            TS_ASSERT_DELTA(mean_a, 0.0, 0.1);
            TS_ASSERT_DELTA(mean_b, 0.0, 0.1);

            // The variances should be close to one another
            TS_ASSERT_DELTA(var_a, var_b, 0.05);

            // The MSE should be big
            TS_ASSERT(mse > 0.1);

            for (auto& p_node : nodes)
            {
                delete p_node;
            }
        }
        { // 1D
            auto p_gen = RandomNumberGenerator::Instance();

            const unsigned n = 1'000;

            const std::array<double, 1> lower_corner {{0.0}};
            const std::array<double, 1> upper_corner {{10.0}};
            const std::array<bool, 1> periodicity {{false}};
            const double lengthscale = 2.0;

            OffLatticeRandomFieldGenerator<1> gen(
                    lower_corner,
                    upper_corner,
                    periodicity,
                    lengthscale
            );

            // Generate some nodes
            std::vector<Node<1>*> nodes(n);
            for (unsigned node_idx = 0; node_idx < n; ++node_idx)
            {
                nodes[node_idx] = new Node<1>(node_idx, Create_c_vector(10.0 * p_gen->ranf()));
            }

            // These three fields should be quite closely correlated
            auto random_field_a = gen.SampleRandomFieldAtTime(nodes, 0.0);
            auto random_field_b = gen.SampleRandomFieldAtTime(nodes, 0.01);
            auto random_field_c = gen.SampleRandomFieldAtTime(nodes, 0.02);

            auto [mean_a, var_a] = CalculateMeanAndVariance(random_field_a);
            const double mse_ab = CalculateMeanSquaredError(random_field_a, random_field_b);
            const double mse_ac = CalculateMeanSquaredError(random_field_a, random_field_c);

            // The mean should be close to zero
            TS_ASSERT_DELTA(mean_a, 0.0, 0.1);

            // The MSE between a and b should be less than between a and c, as a and c are further apart in time
            TS_ASSERT(mse_ac > mse_ab);

            for (auto& p_node : nodes)
            {
                delete p_node;
            }
        }
        { // 2D
            auto p_gen = RandomNumberGenerator::Instance();

            const unsigned n = 1'000;

            const std::array<double, 2> lower_corner {{0.0, 0.0}};
            const std::array<double, 2> upper_corner {{10.0, 10.0}};
            const std::array<bool, 2> periodicity {{false, false}};
            const double lengthscale = 2.0;

            OffLatticeRandomFieldGenerator<2> gen(
                    lower_corner,
                    upper_corner,
                    periodicity,
                    lengthscale
            );

            // Generate some nodes very closely spaced; the field should be nearly constant.
            std::vector<Node<2>*> nodes(n);
            for (unsigned node_idx = 0; node_idx < n; ++node_idx)
            {
                nodes[node_idx] = new Node<2>(node_idx, Create_c_vector(1e-8 * p_gen->ranf(), 1e-8 * p_gen->ranf()));
            }

            auto random_field = gen.SampleRandomFieldAtTime(nodes, 0.0);
            auto [mean, var] = CalculateMeanAndVariance(random_field);

            // Mean should be close to any value, and variance should be near zero
            TS_ASSERT_DELTA(mean, random_field.at(0), 1e-6);
            TS_ASSERT_DELTA(var, 0.0, 1e-12);

            for (auto& p_node : nodes)
            {
                delete p_node;
            }
        }
        { // 3D
            auto p_gen = RandomNumberGenerator::Instance();

            const unsigned n = 1'000;

            const std::array<double, 3> lower_corner {{0.0, 0.0, 0.0}};
            const std::array<double, 3> upper_corner {{10.0, 10.0, 10.0}};
            const std::array<bool, 3> periodicity {{false, false, false}};

            // The expected average spacing between n points in the unit cube
            const double lengthscale = std::cbrt(1.0 / n);

            OffLatticeRandomFieldGenerator<3> gen(
                    lower_corner,
                    upper_corner,
                    periodicity,
                    lengthscale
            );

            // Generate some nodes
            std::vector<Node<3>*> nodes_close(n);
            for (unsigned node_idx = 0; node_idx < n; ++node_idx)
            {
                const double x = lengthscale * p_gen->ranf();
                const double y = lengthscale * p_gen->ranf();
                const double z = lengthscale * p_gen->ranf();
                nodes_close[node_idx] = new Node<3>(node_idx, Create_c_vector(x, y, z));
            }
            std::vector<Node<3>*> nodes_far(n);
            for (unsigned node_idx = 0; node_idx < n; ++node_idx)
            {
                nodes_far[node_idx] = new Node<3>(node_idx, 2.0 * nodes_close[node_idx]->rGetLocation());
            }

            auto random_field_a = gen.SampleRandomFieldAtTime(nodes_close, 0.0);
            auto random_field_b = gen.SampleRandomFieldAtTime(nodes_far, 0.0);

            auto [mean_a, var_a] = CalculateMeanAndVariance(random_field_a);
            auto [mean_b, var_b] = CalculateMeanAndVariance(random_field_b);

            // The nearer nodes should generate a field with lower variance
            TS_ASSERT(var_b > var_a);

            for (auto& p_node : nodes_close)
            {
                delete p_node;
            }
            for (auto& p_node : nodes_far)
            {
                delete p_node;
            }
        }
    }

    void TestSampleFromRandomFieldWithLengthscaleZero()
    {
        auto p_gen = RandomNumberGenerator::Instance();

        const unsigned n = 100;

        const std::array<double, 2> lower_corner {{0.0, 0.0}};
        const std::array<double, 2> upper_corner {{10.0, 10.0}};
        const std::array<bool, 2> periodicity {{false, false}};
        const double lengthscale = 0.0;

        OffLatticeRandomFieldGenerator<2> gen(
                lower_corner,
                upper_corner,
                periodicity,
                lengthscale
        );

        // Generate some nodes
        std::vector<Node<2>*> nodes(n);
        for (unsigned node_idx = 0; node_idx < n; ++node_idx)
        {
            nodes[node_idx] = new Node<2>(node_idx, Create_c_vector(10.0 * p_gen->ranf(), 10.0 * p_gen->ranf()));
        }

        const std::vector<double> grf = gen.SampleRandomFieldAtTime(nodes, 0.0);

        for (auto& p_node : nodes)
        {
            delete p_node;
        }

        for (const auto& val : grf)
        {
            TS_ASSERT_EQUALS(val, 0.0);
        }
    }
};

#endif /*TESTOFFLATTICERANDOMFIELDGENERATOR_HPP_*/