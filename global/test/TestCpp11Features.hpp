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

#ifndef TESTCPP11FEATURES_HPP_
#define TESTCPP11FEATURES_HPP_

#include <cxxtest/TestSuite.h>

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "FakePetscSetup.hpp"

/**
 * Verify that certain C++11 features will compile correctly.
 *
 * There is no real 'testing' in this class; just verification that the compiler can cope with various new language
 * features.
 */
class TestCpp11Features : public CxxTest::TestSuite
{
public:

    void TestAutoKeyword()
    {
        auto i = 0;
        ++i; // Stop the compiler complaining about unused variable
    }

    void TestSmartPointers()
    {
        std::shared_ptr<double> p_shared;
        std::unique_ptr<double> p_unique;

        auto p_made_shared = std::make_shared<double>();
    }

    void TestStdArray()
    {
        std::array<unsigned, 3> my_array;
        ++my_array[0]; // Stop the compiler complaining about unused variable
    }

    void TestStdInitializerList()
    {
        std::vector<unsigned> my_vec = {1, 2, 3, 4, 5};
    }

    void TestRangeFor()
    {
        // Traversing a vector
        {
            std::vector<unsigned> my_vec = {0, 1, 2, 3, 4, 5};

            // Access by const reference
            for (const unsigned &i : my_vec)
            {
                std::cout << i;
            }

            // Access by value.  The type of i is unsigned
            for (auto i : my_vec)
            {
                std::cout << i;
            }

            // Access by reference.  The type of i is unsigned&
            for (auto &&i : my_vec)
            {
                std::cout << i;
            }

            // Initializer may be generated directly
            for (int i : {0, 1, 2, 3, 4, 5})
            {
                std::cout << i;
            }
        }

        // Traversing a map
        {
            std::map<std::string, double> my_map = {{"muck", 1.23}, {"scoop", 2.34}, {"roley", 3.45}};

            for (const auto& i : my_map)
            {
                std::cout << i.first << " " << i.second << std::endl;
            }
        }
    }

    void TestLambdas()
    {
        // Sort a vector of arrays by their x-location
        const unsigned DIM = 3;
        std::vector<std::array<double, DIM>> my_vec(20);

        // Standard mersenne_twister_engine seeded with 0
        std::mt19937 gen(0);
        std::uniform_real_distribution<double> dis(0.0, 1.0);
        auto rand_gen = std::bind(dis, std::ref(gen));

        // Use algorithm std::generate to fill the arrays with random numbers
        for (auto&& array : my_vec)
        {
            std::generate(array.begin(), array.end(), rand_gen);
        }

        // Use algorithm std::sort, with a lambda, to sort smallest-to-largest x-val
        auto sort_x = [](std::array<double, DIM> a, std::array<double, DIM> b) -> bool
        {
            return a[0] < b[0];
        };

        std::sort(my_vec.begin(), my_vec.end(), sort_x);
    }

    void TestTuple()
    {
        auto my_tuple = std::make_tuple(1.23, 5u, "a_string");
        ++std::get<1>(my_tuple); // Stop the compiler complaining about unused variable
    }

    void TestNullptrType()
    {
        int* i = nullptr;
        ++i;
    }
};

#endif /*TESTCPP11FEATURES_HPP_*/
