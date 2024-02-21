/*

Copyright (c) 2005-2024, University of Oxford.
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

#include <array>
#include <vector>
#include <algorithm>

#include "Node.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"
#include "OffLatticeRandomFieldGenerator.hpp"

// These tests do not run in parallel
#include "FakePetscSetup.hpp"

class TestOffLatticeRandomFieldGenerator : public CxxTest::TestSuite
{
public:

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

            const unsigned n = 100;

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
            gen.SetRandomSeed(10);

            // Generate some nodes
            std::vector<Node<1>*> nodes(n);
            for (unsigned node_idx = 0; node_idx < n; ++node_idx)
            {
                nodes[node_idx] = new Node<1>(node_idx, Create_c_vector(10.0 * p_gen->ranf()));
            }

            auto random_field = gen.SampleRandomFieldAtTime(nodes, 0.5);

            std::transform(random_field.begin(), random_field.end(), random_field.begin(), [] (const double& v) { return std::abs(v); });
            auto sum = std::accumulate(random_field.begin(), random_field.end(), 0.0);

            for (auto& p_node : nodes)
            {
                delete p_node;
            }
            TS_ASSERT_DELTA(sum, 28.0820, 0.1)
        }
        { // 1D
            auto p_gen = RandomNumberGenerator::Instance();

            const unsigned n = 100;

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

            auto random_field = gen.SampleRandomFieldAtTime(nodes, 0.0);

            std::transform(random_field.begin(), random_field.end(), random_field.begin(), [] (const double& v) { return std::abs(v); });
            auto sum = std::accumulate(random_field.begin(), random_field.end(), 0.0);
            for (auto& p_node : nodes)
            {
                delete p_node;
            }
            TS_ASSERT(sum > 0.0)
            TS_ASSERT_DELTA(sum, 30.6619, 0.1)
        }
        { // 2D
            auto p_gen = RandomNumberGenerator::Instance();

            const unsigned n = 100;

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

            // Generate some nodes
            std::vector<Node<2>*> nodes(n);
            for (unsigned node_idx = 0; node_idx < n; ++node_idx)
            {
                nodes[node_idx] = new Node<2>(node_idx, Create_c_vector(10.0 * p_gen->ranf(), 10.0 * p_gen->ranf()));
            }

            auto random_field = gen.SampleRandomFieldAtTime(nodes, 0.0);

            std::transform(random_field.begin(), random_field.end(), random_field.begin(), [] (const double& v) { return std::abs(v); });
            auto sum = std::accumulate(random_field.begin(), random_field.end(), 0.0);
            for (auto& p_node : nodes)
            {
                delete p_node;
            }
            TS_ASSERT(sum > 0.0)
            TS_ASSERT_DELTA(sum, 24.5821, 0.6)
        }
        { // 3D
            auto p_gen = RandomNumberGenerator::Instance();

            const unsigned n = 100;

            const std::array<double, 3> lower_corner {{0.0, 0.0, 0.0}};
            const std::array<double, 3> upper_corner {{10.0, 10.0, 10.0}};
            const std::array<bool, 3> periodicity {{false, false, false}};
            const double lengthscale = 2.0;

            OffLatticeRandomFieldGenerator<3> gen(
                    lower_corner,
                    upper_corner,
                    periodicity,
                    lengthscale
            );

            // Generate some nodes
            std::vector<Node<3>*> nodes(n);
            for (unsigned node_idx = 0; node_idx < n; ++node_idx)
            {
                nodes[node_idx] = new Node<3>(node_idx, Create_c_vector(10.0 * p_gen->ranf(), 10.0 * p_gen->ranf(), 10.0 * p_gen->ranf()));
            }

            auto random_field = gen.SampleRandomFieldAtTime(nodes, 0.0);

            std::transform(random_field.begin(), random_field.end(), random_field.begin(), [] (const double& v) { return std::abs(v); });
            auto sum = std::accumulate(random_field.begin(), random_field.end(), 0.0);
            for (auto& p_node : nodes)
            {
                delete p_node;
            }
            TS_ASSERT(sum > 0.0)
            TS_ASSERT_DELTA(sum, 16.1120, 0.1)
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