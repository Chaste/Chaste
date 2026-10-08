/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef SEMLATTICEGEOMETRY_HPP_
#define SEMLATTICEGEOMETRY_HPP_

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "Exception.hpp"
#include "SemEnumerations.hpp"
#include "UblasVectorInclude.hpp"

/**
 * Free functions shared by the SEM mesh generators for laying points out on a
 * regular lattice. Both the subcellular nodes within a single element and the
 * elements within a multi-element mesh are positioned this way, so the lattice
 * arithmetic lives here rather than being duplicated in each generator.
 */

/**
 * Generate the positions of a block of lattice points.
 *
 * Points are generated with the first dimension varying fastest, so the flat
 * index of the point at grid coordinate (i, j, k) is
 * i + numPerDim[0]*(j + numPerDim[1]*k) for every lattice type. Close-packed
 * lattices displace points from their grid coordinates but do not reorder them,
 * so a caller may continue to infer a point's grid coordinate from its index.
 *
 * The point at grid coordinate zero is always placed at the origin, and all
 * positions have non-negative components.
 *
 * For SemLatticeType::SEM_LATTICE_CLOSE_PACKED the returned lattice is a true close packing
 * only when every entry of rSpacing is equal; unequal spacings stretch the
 * lattice along the corresponding axes.
 *
 * @param rNumPerDim number of points in each dimension, each of which must be at least one
 * @param rSpacing distance between adjacent points in each dimension, each of which must be positive
 * @param lattice the lattice on which to place the points
 * @return the position of each lattice point
 */
template <unsigned DIM>
std::vector<c_vector<double, DIM>> GenerateSemLatticePositions(
    const std::array<unsigned, DIM>& rNumPerDim,
    const std::array<double, DIM>& rSpacing,
    SemLatticeType::Value lattice)
{
    unsigned num_points = 1u;
    for (unsigned dim = 0; dim < DIM; ++dim)
    {
        if (rNumPerDim[dim] == 0u)
        {
            EXCEPTION("GenerateSemLatticePositions: each entry of numPerDim must be >= 1");
        }
        if (rSpacing[dim] <= 0.0)
        {
            EXCEPTION("GenerateSemLatticePositions: each entry of spacing must be positive");
        }
        num_points *= rNumPerDim[dim];
    }

    // Ratios giving unit nearest-neighbour distance on a close-packed lattice: successive
    // rows are a half spacing apart in x and sqrt(3)/2 apart in y, and successive layers
    // are a third of a row apart in y and sqrt(2/3) apart in z.
    const double row_offset_x = 0.5;
    const double row_spacing_y = 0.5 * std::sqrt(3.0);
    const double layer_offset_y = 1.0 / 3.0;
    const double layer_spacing_z = std::sqrt(2.0 / 3.0);

    const bool close_packed = (lattice == SemLatticeType::SEM_LATTICE_CLOSE_PACKED);

    std::vector<c_vector<double, DIM>> positions;
    positions.reserve(num_points);

    if constexpr (DIM == 1)
    {
        // A line of points admits only one packing, so the lattice type is immaterial
        for (unsigned i = 0; i < rNumPerDim[0]; ++i)
        {
            c_vector<double, 1> v;
            v[0] = static_cast<double>(i) * rSpacing[0];
            positions.push_back(v);
        }
    }
    else if constexpr (DIM == 2)
    {
        for (unsigned j = 0; j < rNumPerDim[1]; ++j)
        {
            const double x_offset = close_packed ? row_offset_x * static_cast<double>(j % 2u) : 0.0;
            const double y_scaling = close_packed ? row_spacing_y : 1.0;

            for (unsigned i = 0; i < rNumPerDim[0]; ++i)
            {
                c_vector<double, 2> v;
                v[0] = (static_cast<double>(i) + x_offset) * rSpacing[0];
                v[1] = static_cast<double>(j) * y_scaling * rSpacing[1];
                positions.push_back(v);
            }
        }
    }
    else if constexpr (DIM == 3)
    {
        for (unsigned k = 0; k < rNumPerDim[2]; ++k)
        {
            const double y_offset = close_packed ? layer_offset_y * static_cast<double>(k % 2u) : 0.0;
            const double y_scaling = close_packed ? row_spacing_y : 1.0;
            const double z_scaling = close_packed ? layer_spacing_z : 1.0;

            for (unsigned j = 0; j < rNumPerDim[1]; ++j)
            {
                // Alternate layers stagger in x as well as y, giving ABAB stacking
                const double x_offset = close_packed ? row_offset_x * static_cast<double>((j + k) % 2u) : 0.0;

                for (unsigned i = 0; i < rNumPerDim[0]; ++i)
                {
                    c_vector<double, 3> v;
                    v[0] = (static_cast<double>(i) + x_offset) * rSpacing[0];
                    v[1] = (static_cast<double>(j) + y_offset) * y_scaling * rSpacing[1];
                    v[2] = static_cast<double>(k) * z_scaling * rSpacing[2];
                    positions.push_back(v);
                }
            }
        }
    }
    else
    {
        NEVER_REACHED;
    }

    return positions;
}

/**
 * Calculate the extent of a set of positions along each axis, i.e. the side
 * lengths of the smallest axis-aligned bounding box containing them.
 *
 * @param rPositions the positions to measure, which must not be empty
 * @return the extent along each axis
 */
template <unsigned DIM>
c_vector<double, DIM> GetSemLatticeExtent(const std::vector<c_vector<double, DIM>>& rPositions)
{
    if (rPositions.empty())
    {
        EXCEPTION("GetSemLatticeExtent: positions must not be empty");
    }

    c_vector<double, DIM> min_position = rPositions.front();
    c_vector<double, DIM> max_position = rPositions.front();

    for (const auto& r_position : rPositions)
    {
        for (unsigned dim = 0; dim < DIM; ++dim)
        {
            min_position[dim] = std::min(min_position[dim], r_position[dim]);
            max_position[dim] = std::max(max_position[dim], r_position[dim]);
        }
    }

    return max_position - min_position;
}

/**
 * Calculate the spacing at which identical blocks of the given extent may be placed on a
 * close-packed lattice while keeping at least the given clearance between them.
 *
 * A close-packed lattice staggers its rows and layers, so a block placed on it approaches its
 * neighbours along directions that mix the axes; spacing each axis by its own extent, as is
 * sufficient for a cubic lattice, would let the blocks interlock. A single spacing is therefore
 * used in every direction, chosen as the largest that any neighbour direction demands.
 *
 * For each neighbouring lattice site the width of the block along that direction is measured, and
 * the spacing must be large enough that the two blocks are separated by at least the clearance.
 * In 1D, where a close-packed lattice is just an evenly spaced line, this reduces to the extent
 * plus the clearance.
 *
 * @param rExtent the extent of the block along each axis
 * @param clearance the minimum separation to leave between neighbouring blocks
 * @return the lattice spacing to use in every dimension
 */
template <unsigned DIM>
double GetSemClosePackedSpacing(const c_vector<double, DIM>& rExtent, double clearance)
{
    if (clearance <= 0.0)
    {
        EXCEPTION("GetSemClosePackedSpacing: clearance must be positive");
    }

    // A close-packed lattice of unit spacing, whose central site's neighbours give the directions
    // along which adjacent blocks approach one another. Two shells of neighbours are enough to
    // cover every direction that can be limiting, however elongated the block.
    const unsigned num_per_dim_1d = 5u;
    const unsigned centre_index_1d = 2u;

    std::array<unsigned, DIM> num_per_dim;
    std::array<double, DIM> unit_spacing;
    num_per_dim.fill(num_per_dim_1d);
    unit_spacing.fill(1.0);

    const std::vector<c_vector<double, DIM>> lattice
        = GenerateSemLatticePositions<DIM>(num_per_dim, unit_spacing, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);

    unsigned centre_index = 0u;
    unsigned stride = 1u;
    for (unsigned dim = 0; dim < DIM; ++dim)
    {
        centre_index += centre_index_1d * stride;
        stride *= num_per_dim_1d;
    }

    double spacing = 0.0;
    for (const auto& r_site : lattice)
    {
        const c_vector<double, DIM> offset = r_site - lattice[centre_index];
        const double offset_norm = norm_2(offset);

        if (offset_norm < 1e-12)
        {
            continue;
        }

        // Width of the block along this direction, i.e. the support width of its bounding box
        double block_width = 0.0;
        for (unsigned dim = 0; dim < DIM; ++dim)
        {
            block_width += rExtent[dim] * fabs(offset[dim]) / offset_norm;
        }

        spacing = std::max(spacing, (block_width + clearance) / offset_norm);
    }

    return spacing;
}

#endif // SEMLATTICEGEOMETRY_HPP_
