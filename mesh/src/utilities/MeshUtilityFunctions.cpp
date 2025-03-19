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

#include "MeshUtilityFunctions.hpp"

#include "RandomNumberGenerator.hpp"

#include <algorithm>
#include <cassert>

template <std::size_t DIM>
std::vector<c_vector<double, DIM>> EvenlySpaceAlongPath(
        std::vector<c_vector<double, DIM>>& path,
        const bool closedPath,
        const bool permuteOrder,
        std::size_t numPointsToPlace,
        const double targetSpacing
)
{
    assert(!path.empty());
    assert(numPointsToPlace > 0 || targetSpacing < DBL_MAX);

    if (closedPath)
    {
        if (permuteOrder)
        {
            // Rotate the closed path around a random pivot so we start at a random location along it
            const unsigned random_pivot = RandomNumberGenerator::Instance()->randMod(path.size());
            std::rotate(path.begin(), path.begin() + random_pivot, path.end());
        }

        path.emplace_back(path.front());
    }

    // Calculate a vector of edge vectors, and a vector of edge lengths
    std::vector<c_vector<double, DIM>> edge_unit_vectors;
    edge_unit_vectors.reserve(path.size() - 1);

    std::vector<double> edge_length_partial_sums;
    edge_length_partial_sums.reserve(path.size() - 1);

    for (unsigned i = 1; i < path.size(); ++i)
    {
        const c_vector<double, DIM> un_normalised_vec = path[i] - path[i - 1];
        const double edge_length = norm_2(un_normalised_vec);

        edge_unit_vectors.emplace_back(un_normalised_vec / edge_length);

        if (edge_length_partial_sums.empty())
        {
            edge_length_partial_sums.emplace_back(edge_length);
        }
        else
        {
            edge_length_partial_sums.emplace_back(edge_length_partial_sums.back() + edge_length);
        }
    }

    // Calculate the spacing, allowing the targetSpacing to override the number of points to place.
    // If the path is open, we want to hit the start and end, so subtract one from the numPointsToPlace.
    if (targetSpacing < DBL_MAX)
    {
        numPointsToPlace = static_cast<unsigned>(!closedPath + std::ceil(edge_length_partial_sums.back() / targetSpacing));
    }
    const double spacing = edge_length_partial_sums.back() / (numPointsToPlace - static_cast<unsigned>(!closedPath));

    // The vector of new locations to add to
    std::vector<c_vector<double, DIM>> new_locations;
    new_locations.reserve(numPointsToPlace);

    double remainder = 0.0;
    for (unsigned edge = 0; edge < edge_length_partial_sums.size(); ++edge)
    {
        unsigned new_locs_this_edge = 0u;
        while(spacing * new_locations.size() < edge_length_partial_sums[edge])
        {
            new_locations.emplace_back(path[edge] + (new_locs_this_edge * spacing - remainder) * edge_unit_vectors[edge]);
            new_locs_this_edge++;
        }
        remainder = edge_length_partial_sums[edge] - spacing * new_locations.size();
    }

    if (!closedPath)
    {
        new_locations.emplace_back(path.back());
    }

    return new_locations;
}

// Explicit instantiation
template std::vector<c_vector<double, 1>> EvenlySpaceAlongPath(std::vector<c_vector<double, 1>>&, bool, bool, std::size_t, double);
template std::vector<c_vector<double, 2>> EvenlySpaceAlongPath(std::vector<c_vector<double, 2>>&, bool, bool, std::size_t, double);
template std::vector<c_vector<double, 3>> EvenlySpaceAlongPath(std::vector<c_vector<double, 3>>&, bool, bool, std::size_t, double);
