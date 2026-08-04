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

#ifndef TESTSEMLATTICEGEOMETRY_HPP_
#define TESTSEMLATTICEGEOMETRY_HPP_

#include <cxxtest/TestSuite.h>

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <vector>

#include "SemEnumerations.hpp"
#include "SemLatticeGeometry.hpp"
#include "UblasCustomFunctions.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"


class TestSemLatticeGeometry : public CxxTest::TestSuite
{
private:

    /** Half the height of an equilateral triangle of unit side, i.e. the row spacing of a close packing */
    static constexpr double mRowSpacing = 0.86602540378443865;

    /** The spacing between layers of a three-dimensional close packing, sqrt(2/3) */
    static constexpr double mLayerSpacing = 0.81649658092772603;

    /**
     * Calculate the distance between two axis-aligned boxes of identical extent whose centres are
     * separated by the given offset, which is zero if the boxes overlap.
     *
     * @param rCentreOffset the offset between the two box centres
     * @param rExtent the extent of each box along each axis
     * @return the distance between the boxes
     */
    template <unsigned DIM>
    double GetBoxSeparation(const c_vector<double, DIM>& rCentreOffset, const c_vector<double, DIM>& rExtent)
    {
        double sum_of_squares = 0.0;
        for (unsigned dim = 0; dim < DIM; ++dim)
        {
            // Each box contributes half its extent, so the two together span a full extent
            const double gap = fabs(rCentreOffset[dim]) - rExtent[dim];
            if (gap > 0.0)
            {
                sum_of_squares += gap * gap;
            }
        }
        return sqrt(sum_of_squares);
    }

public:

    void TestCubicLatticePositions()
    {
        // In 1D a lattice is just evenly spaced points
        std::vector<c_vector<double, 1>> positions_1d
            = GenerateSemLatticePositions<1>({4u}, {0.5}, SemLatticeType::SEM_LATTICE_CUBIC);

        TS_ASSERT_EQUALS(positions_1d.size(), 4u);
        for (unsigned i = 0; i < 4u; ++i)
        {
            TS_ASSERT_DELTA(positions_1d[i][0], 0.5 * static_cast<double>(i), 1e-12);
        }

        // In 2D the first dimension varies fastest, so index = i + numPerDim[0]*j
        std::vector<c_vector<double, 2>> positions_2d
            = GenerateSemLatticePositions<2>({3u, 2u}, {1.0, 1.0}, SemLatticeType::SEM_LATTICE_CUBIC);

        TS_ASSERT_EQUALS(positions_2d.size(), 6u);
        for (unsigned j = 0; j < 2u; ++j)
        {
            for (unsigned i = 0; i < 3u; ++i)
            {
                const c_vector<double, 2> expected = Create_c_vector(static_cast<double>(i), static_cast<double>(j));
                TS_ASSERT_DELTA(norm_2(positions_2d[i + 3u * j] - expected), 0.0, 1e-12);
            }
        }

        // In 3D index = i + numPerDim[0]*(j + numPerDim[1]*k)
        std::vector<c_vector<double, 3>> positions_3d
            = GenerateSemLatticePositions<3>({2u, 3u, 4u}, {1.0, 1.0, 1.0}, SemLatticeType::SEM_LATTICE_CUBIC);

        TS_ASSERT_EQUALS(positions_3d.size(), 24u);
        for (unsigned k = 0; k < 4u; ++k)
        {
            for (unsigned j = 0; j < 3u; ++j)
            {
                for (unsigned i = 0; i < 2u; ++i)
                {
                    const c_vector<double, 3> expected = Create_c_vector(static_cast<double>(i),
                        static_cast<double>(j), static_cast<double>(k));
                    TS_ASSERT_DELTA(norm_2(positions_3d[i + 2u * (j + 3u * k)] - expected), 0.0, 1e-12);
                }
            }
        }
    }

    void TestClosePackedLatticeMatchesCubicIn1d()
    {
        // There is only one way to space points along a line
        std::vector<c_vector<double, 1>> cubic
            = GenerateSemLatticePositions<1>({6u}, {0.25}, SemLatticeType::SEM_LATTICE_CUBIC);
        std::vector<c_vector<double, 1>> close_packed
            = GenerateSemLatticePositions<1>({6u}, {0.25}, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);

        TS_ASSERT_EQUALS(close_packed.size(), cubic.size());
        for (unsigned i = 0; i < cubic.size(); ++i)
        {
            TS_ASSERT_DELTA(close_packed[i][0], cubic[i][0], 1e-12);
        }
    }

    void TestClosePackedLatticePositionsIn2d()
    {
        std::vector<c_vector<double, 2>> positions
            = GenerateSemLatticePositions<2>({3u, 3u}, {1.0, 1.0}, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);

        TS_ASSERT_EQUALS(positions.size(), 9u);

        // Odd-numbered rows are offset half a spacing in x, and rows are sqrt(3)/2 apart in y
        for (unsigned j = 0; j < 3u; ++j)
        {
            for (unsigned i = 0; i < 3u; ++i)
            {
                const double expected_x = static_cast<double>(i) + (j % 2u == 1u ? 0.5 : 0.0);
                const double expected_y = static_cast<double>(j) * mRowSpacing;
                const c_vector<double, 2> expected = Create_c_vector(expected_x, expected_y);

                TS_ASSERT_DELTA(norm_2(positions[i + 3u * j] - expected), 0.0, 1e-12);
            }
        }

        // The point at grid coordinate zero is at the origin, and no position is negative
        TS_ASSERT_DELTA(norm_2(positions[0]), 0.0, 1e-12);
        for (const auto& r_position : positions)
        {
            TS_ASSERT_LESS_THAN_EQUALS(0.0, r_position[0]);
            TS_ASSERT_LESS_THAN_EQUALS(0.0, r_position[1]);
        }
    }

    void TestClosePackedLatticePositionsIn3d()
    {
        std::vector<c_vector<double, 3>> positions
            = GenerateSemLatticePositions<3>({3u, 3u, 3u}, {1.0, 1.0, 1.0}, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);

        TS_ASSERT_EQUALS(positions.size(), 27u);

        // Layers stagger in both x and y, giving ABAB stacking, and are sqrt(2/3) apart in z
        for (unsigned k = 0; k < 3u; ++k)
        {
            for (unsigned j = 0; j < 3u; ++j)
            {
                for (unsigned i = 0; i < 3u; ++i)
                {
                    const double expected_x = static_cast<double>(i) + ((j + k) % 2u == 1u ? 0.5 : 0.0);
                    const double expected_y = (static_cast<double>(j) + (k % 2u == 1u ? 1.0 / 3.0 : 0.0)) * mRowSpacing;
                    const double expected_z = static_cast<double>(k) * mLayerSpacing;
                    const c_vector<double, 3> expected = Create_c_vector(expected_x, expected_y, expected_z);

                    TS_ASSERT_DELTA(norm_2(positions[i + 3u * (j + 3u * k)] - expected), 0.0, 1e-12);
                }
            }
        }
    }

    void TestClosePackedLatticeHasUnitNearestNeighbourDistance()
    {
        // The defining property of a close packing: with unit spacing, no two points are closer
        // than one spacing, and a point in the interior has six neighbours in 2D and twelve in 3D
        std::vector<c_vector<double, 2>> positions_2d
            = GenerateSemLatticePositions<2>({5u, 5u}, {1.0, 1.0}, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);
        std::vector<c_vector<double, 3>> positions_3d
            = GenerateSemLatticePositions<3>({5u, 5u, 5u}, {1.0, 1.0, 1.0}, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);

        double min_separation_2d = DBL_MAX;
        unsigned num_neighbours_2d = 0u;
        const c_vector<double, 2> centre_2d = positions_2d[2u + 5u * 2u];
        for (unsigned i = 0; i < positions_2d.size(); ++i)
        {
            const double distance_to_centre = norm_2(positions_2d[i] - centre_2d);
            if (distance_to_centre > 1e-12 && fabs(distance_to_centre - 1.0) < 1e-9)
            {
                num_neighbours_2d++;
            }
            for (unsigned j = i + 1u; j < positions_2d.size(); ++j)
            {
                min_separation_2d = std::min(min_separation_2d, norm_2(positions_2d[i] - positions_2d[j]));
            }
        }

        TS_ASSERT_DELTA(min_separation_2d, 1.0, 1e-9);
        TS_ASSERT_EQUALS(num_neighbours_2d, 6u);

        double min_separation_3d = DBL_MAX;
        unsigned num_neighbours_3d = 0u;
        const c_vector<double, 3> centre_3d = positions_3d[2u + 5u * (2u + 5u * 2u)];
        for (unsigned i = 0; i < positions_3d.size(); ++i)
        {
            const double distance_to_centre = norm_2(positions_3d[i] - centre_3d);
            if (distance_to_centre > 1e-12 && fabs(distance_to_centre - 1.0) < 1e-9)
            {
                num_neighbours_3d++;
            }
            for (unsigned j = i + 1u; j < positions_3d.size(); ++j)
            {
                min_separation_3d = std::min(min_separation_3d, norm_2(positions_3d[i] - positions_3d[j]));
            }
        }

        TS_ASSERT_DELTA(min_separation_3d, 1.0, 1e-9);
        TS_ASSERT_EQUALS(num_neighbours_3d, 12u);
    }

    void TestUnequalSpacingsStretchTheLattice()
    {
        // Spacings apply per axis, so an unequal spacing stretches rather than shears the lattice
        std::vector<c_vector<double, 2>> positions
            = GenerateSemLatticePositions<2>({2u, 2u}, {2.0, 4.0}, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);

        TS_ASSERT_DELTA(norm_2(positions[0] - Create_c_vector(0.0, 0.0)), 0.0, 1e-12);
        TS_ASSERT_DELTA(norm_2(positions[1] - Create_c_vector(2.0, 0.0)), 0.0, 1e-12);
        TS_ASSERT_DELTA(norm_2(positions[2] - Create_c_vector(1.0, 4.0 * mRowSpacing)), 0.0, 1e-12);
        TS_ASSERT_DELTA(norm_2(positions[3] - Create_c_vector(3.0, 4.0 * mRowSpacing)), 0.0, 1e-12);
    }

    void TestSinglePointLattice()
    {
        // A lattice of one point sits at the origin, whichever lattice is requested
        std::vector<c_vector<double, 3>> cubic
            = GenerateSemLatticePositions<3>({1u, 1u, 1u}, {1.0, 1.0, 1.0}, SemLatticeType::SEM_LATTICE_CUBIC);
        std::vector<c_vector<double, 3>> close_packed
            = GenerateSemLatticePositions<3>({1u, 1u, 1u}, {1.0, 1.0, 1.0}, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);

        TS_ASSERT_EQUALS(cubic.size(), 1u);
        TS_ASSERT_EQUALS(close_packed.size(), 1u);
        TS_ASSERT_DELTA(norm_2(cubic[0]), 0.0, 1e-12);
        TS_ASSERT_DELTA(norm_2(close_packed[0]), 0.0, 1e-12);
    }

    void TestLatticeExtent()
    {
        // An arbitrary set of points
        std::vector<c_vector<double, 2>> positions;
        positions.push_back(Create_c_vector(1.0, -2.0));
        positions.push_back(Create_c_vector(-3.0, 4.0));
        positions.push_back(Create_c_vector(0.5, 0.5));

        const c_vector<double, 2> extent = GetSemLatticeExtent<2>(positions);
        TS_ASSERT_DELTA(extent[0], 4.0, 1e-12);
        TS_ASSERT_DELTA(extent[1], 6.0, 1e-12);

        // A single point has zero extent
        std::vector<c_vector<double, 3>> single_position;
        single_position.push_back(Create_c_vector(7.0, 8.0, 9.0));
        TS_ASSERT_DELTA(norm_2(GetSemLatticeExtent<3>(single_position)), 0.0, 1e-12);

        // A cubic lattice spans one fewer spacing than it has points in each dimension
        const c_vector<double, 3> cubic_extent = GetSemLatticeExtent<3>(
            GenerateSemLatticePositions<3>({4u, 3u, 2u}, {1.0, 2.0, 3.0}, SemLatticeType::SEM_LATTICE_CUBIC));
        TS_ASSERT_DELTA(cubic_extent[0], 3.0, 1e-12);
        TS_ASSERT_DELTA(cubic_extent[1], 4.0, 1e-12);
        TS_ASSERT_DELTA(cubic_extent[2], 3.0, 1e-12);

        // Staggering widens a close-packed lattice in x by half a spacing, and shortens it in y
        const c_vector<double, 2> packed_extent = GetSemLatticeExtent<2>(
            GenerateSemLatticePositions<2>({3u, 3u}, {1.0, 1.0}, SemLatticeType::SEM_LATTICE_CLOSE_PACKED));
        TS_ASSERT_DELTA(packed_extent[0], 2.5, 1e-12);
        TS_ASSERT_DELTA(packed_extent[1], 2.0 * mRowSpacing, 1e-12);
    }

    void TestClosePackedSpacingIn1d()
    {
        // A close-packed lattice in 1D is an evenly spaced line, so blocks need only their own
        // extent plus the clearance
        const c_vector<double, 1> extent = Create_c_vector(3.0);
        TS_ASSERT_DELTA(GetSemClosePackedSpacing<1>(extent, 0.5), 3.5, 1e-12);
    }

    void TestClosePackedSpacingIn2d()
    {
        // A unit square approached along the diagonal to a staggered neighbour presents a width of
        // (0.5 + sqrt(3)/2), which is the limiting direction, so the spacing is that plus the clearance
        const c_vector<double, 2> extent = Create_c_vector(1.0, 1.0);
        TS_ASSERT_DELTA(GetSemClosePackedSpacing<2>(extent, 1.0), 1.5 + mRowSpacing, 1e-9);

        // The calculation is homogeneous: scaling the extent and the clearance scales the spacing
        TS_ASSERT_DELTA(GetSemClosePackedSpacing<2>(2.0 * extent, 2.0),
                        2.0 * GetSemClosePackedSpacing<2>(extent, 1.0), 1e-9);
    }

    void TestClosePackedSpacingOfPointBlocks()
    {
        // Blocks of no extent are points, which need only be one clearance apart
        const c_vector<double, 3> extent = zero_vector<double>(3);
        TS_ASSERT_DELTA(GetSemClosePackedSpacing<3>(extent, 1.5), 1.5, 1e-9);
    }

    void TestClosePackedSpacingLeavesClearanceBetweenBlocks()
    {
        // The contract: blocks of the given extent, placed on a close-packed lattice at the
        // returned spacing, are never closer to one another than the requested clearance. The
        // spacing is derived from the width of the block's bounding box along each neighbour
        // direction, which is a lower bound on the separation actually achieved, so cubic and
        // elongated blocks alike end up at least the clearance apart.
        const double clearance = 0.75;
        const c_vector<double, 3> extents[2] = {Create_c_vector(1.0, 1.0, 1.0), Create_c_vector(2.0, 3.0, 1.0)};

        for (const c_vector<double, 3>& r_extent : extents)
        {
            const double spacing = GetSemClosePackedSpacing<3>(r_extent, clearance);
            std::array<double, 3> lattice_spacing;
            lattice_spacing.fill(spacing);

            const std::vector<c_vector<double, 3>> sites
                = GenerateSemLatticePositions<3>({4u, 4u, 4u}, lattice_spacing, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);

            double min_separation = DBL_MAX;
            for (unsigned i = 0; i < sites.size(); ++i)
            {
                for (unsigned j = i + 1u; j < sites.size(); ++j)
                {
                    min_separation = std::min(min_separation, GetBoxSeparation<3>(sites[i] - sites[j], r_extent));
                }
            }

            TS_ASSERT_LESS_THAN_EQUALS(clearance - 1e-9, min_separation);
        }
    }

    void TestExceptions()
    {
        std::array<double, 2> valid_spacing = {1.0, 1.0};
        std::array<unsigned, 2> valid_num_per_dim = {2u, 2u};

        // A lattice must have at least one point in each dimension
        std::array<unsigned, 2> zero_num_per_dim = {2u, 0u};
        TS_ASSERT_THROWS_THIS(GenerateSemLatticePositions<2>(zero_num_per_dim, valid_spacing, SemLatticeType::SEM_LATTICE_CUBIC),
            "GenerateSemLatticePositions: each entry of numPerDim must be >= 1");

        // Spacings must be positive
        std::array<double, 2> zero_spacing = {1.0, 0.0};
        TS_ASSERT_THROWS_THIS(GenerateSemLatticePositions<2>(valid_num_per_dim, zero_spacing, SemLatticeType::SEM_LATTICE_CLOSE_PACKED),
            "GenerateSemLatticePositions: each entry of spacing must be positive");

        std::array<double, 2> negative_spacing = {-1.0, 1.0};
        TS_ASSERT_THROWS_THIS(GenerateSemLatticePositions<2>(valid_num_per_dim, negative_spacing, SemLatticeType::SEM_LATTICE_CUBIC),
            "GenerateSemLatticePositions: each entry of spacing must be positive");

        // There is no extent to measure without any positions
        std::vector<c_vector<double, 2>> no_positions;
        TS_ASSERT_THROWS_THIS(GetSemLatticeExtent<2>(no_positions),
            "GetSemLatticeExtent: positions must not be empty");

        // A non-positive clearance would allow blocks to touch or overlap
        const c_vector<double, 2> extent = Create_c_vector(1.0, 1.0);
        TS_ASSERT_THROWS_THIS(GetSemClosePackedSpacing<2>(extent, 0.0),
            "GetSemClosePackedSpacing: clearance must be positive");

        TS_ASSERT_THROWS_THIS(GetSemClosePackedSpacing<2>(extent, -1.0),
            "GetSemClosePackedSpacing: clearance must be positive");
    }
};

#endif /*TESTSEMLATTICEGEOMETRY_HPP_*/
