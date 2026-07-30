
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

#include "SemSphericalElementMeshGenerator.hpp"

#include <algorithm>
#include <array>
#include <cmath>

#include "Exception.hpp"
#include "SemEnumerations.hpp"
#include "SemLatticeGeometry.hpp"
#include "SemMesh.hpp"


template <unsigned DIM>
SemSphericalElementMeshGenerator<DIM>::SemSphericalElementMeshGenerator(
    unsigned numNodes,
    double cellRadius,
    double boundaryThickness)
    : mpMesh{std::make_shared<SemMesh<DIM>>()},
      mNumNodes{numNodes},
      mCellRadius{cellRadius},
      mBoundaryThickness{boundaryThickness},
      mNodeSpacing{0.0}
{
    // Two nodes are the fewest that can be scaled to put an outermost node at the requested
    // radius; a single node would sit at the centre, giving an element of radius zero
    if (mNumNodes < 2u)
    {
        EXCEPTION("SemSphericalElementMeshGenerator: numNodes must be >= 2");
    }

    if (cellRadius <= 0.0)
    {
        EXCEPTION("SemSphericalElementMeshGenerator: cellRadius must be positive");
    }

    if (boundaryThickness < 0.0)
    {
        EXCEPTION("SemSphericalElementMeshGenerator: boundaryThickness must be non-negative");
    }

    std::vector<c_vector<double, DIM>> positions = GenerateNodePositions();
    GenerateMesh(positions);
}

template <unsigned DIM>
std::shared_ptr<SemMesh<DIM>> SemSphericalElementMeshGenerator<DIM>::GetMesh()
{
    return mpMesh;
}

template <unsigned DIM>
double SemSphericalElementMeshGenerator<DIM>::GetNodeSpacing() const
{
    return mNodeSpacing;
}

template <unsigned DIM>
std::vector<c_vector<double, DIM>> SemSphericalElementMeshGenerator<DIM>::GenerateNodePositions()
{
    // Positions are compared on a fixed grid rather than directly, so that lattice sites which
    // are equidistant from the centre are ordered identically whatever rounding the norm incurs
    const double quantum = 1e-9;

    // Every point of space lies within one unit of some site of a close packing at unit spacing:
    // the covering radius is 1/sqrt(2) for the lattice in 3D, 1/sqrt(3) in 2D and 1/2 in 1D. The
    // sites lying within a radius R of the centre therefore cover a ball of radius R - 1 between
    // them, so there are at least (R - 1)^DIM of them. Taking R - 1 to be the DIM'th root of the
    // requested node count bounds the radius of the ball the retained sites occupy.
    const double covering_radius = 1.0;
    const double ball_radius
        = covering_radius + std::pow(static_cast<double>(mNumNodes), 1.0 / DIM);

    // Size the block of sites so that ball is cut cleanly from it. The block spans at least
    // (n - 1) * sqrt(2/3) in each dimension, sqrt(2/3) being the closest any close-packed lattice
    // brings successive rows or layers, and the site taken as the centre is within the covering
    // radius of the centre of the block.
    const double min_lattice_spacing = std::sqrt(2.0 / 3.0);
    const unsigned num_per_dim_1d = 1u + static_cast<unsigned>(
        std::ceil(2.0 * (ball_radius + covering_radius) / min_lattice_spacing));

    std::array<unsigned, DIM> num_per_dim;
    num_per_dim.fill(num_per_dim_1d);

    std::array<double, DIM> unit_spacing;
    unit_spacing.fill(1.0);

    std::vector<c_vector<double, DIM>> sites
        = GenerateSemLatticePositions<DIM>(num_per_dim, unit_spacing, SEM_LATTICE_CLOSE_PACKED);

    // GenerateSemLatticePositions() places the block in the positive orthant with a corner at the
    // origin, so its extent is also the position of its far corner
    const c_vector<double, DIM> extent = GetSemLatticeExtent<DIM>(sites);
    const c_vector<double, DIM> block_centre = 0.5 * extent;

    // Centre the ball on a lattice site, rather than on the centre of the block, so that a node
    // sits at the centre of the cell and the packing around it is symmetric
    c_vector<double, DIM> centre = sites.front();
    double centre_offset = norm_2(sites.front() - block_centre);
    for (const auto& r_site : sites)
    {
        const double offset = norm_2(r_site - block_centre);
        if (offset < centre_offset - quantum)
        {
            centre_offset = offset;
            centre = r_site;
        }
    }

    // Order the sites by distance from that centre, breaking ties by coordinate so that the
    // retained set does not depend on how equidistant sites happen to be ordered beforehand
    auto compare = [&centre, quantum](const c_vector<double, DIM>& rSiteA,
                                      const c_vector<double, DIM>& rSiteB)
    {
        const long long radius_a = std::llround(norm_2(rSiteA - centre) / quantum);
        const long long radius_b = std::llround(norm_2(rSiteB - centre) / quantum);
        if (radius_a != radius_b)
        {
            return radius_a < radius_b;
        }

        for (unsigned dim = 0; dim < DIM; ++dim)
        {
            const long long coord_a = std::llround(rSiteA[dim] / quantum);
            const long long coord_b = std::llround(rSiteB[dim] / quantum);
            if (coord_a != coord_b)
            {
                return coord_a < coord_b;
            }
        }

        return false;
    };

    const auto num_to_keep = static_cast<std::ptrdiff_t>(mNumNodes);
    std::partial_sort(sites.begin(), sites.begin() + num_to_keep, sites.end(), compare);

    std::vector<c_vector<double, DIM>> positions(sites.begin(), sites.begin() + num_to_keep);
    for (auto& r_position : positions)
    {
        r_position -= centre;
    }

    const double max_radius = norm_2(positions.back());

    // The retained sites are the nearest sites of the infinite lattice only if every site within
    // max_radius of the centre was present in the block, which the sizing above guarantees
    double distance_to_face = extent[0] - centre[0];
    for (unsigned dim = 0; dim < DIM; ++dim)
    {
        distance_to_face = std::min(distance_to_face, centre[dim]);
        distance_to_face = std::min(distance_to_face, extent[dim] - centre[dim]);
    }

    if (max_radius > distance_to_face + quantum)
    {
        NEVER_REACHED;
    }

    // The lattice was generated at unit spacing, so scaling the outermost node onto the requested
    // radius scales the nearest-neighbour distance by the same factor
    mNodeSpacing = mCellRadius / max_radius;
    for (auto& r_position : positions)
    {
        r_position *= mNodeSpacing;
    }

    return positions;
}

template <unsigned DIM>
bool SemSphericalElementMeshGenerator<DIM>::IsBoundaryNode(const c_vector<double, DIM>& rPosition) const
{
    const double tolerance = 1e-9 * mCellRadius;
    return norm_2(rPosition) > mCellRadius - mBoundaryThickness * mNodeSpacing - tolerance;
}

template <unsigned DIM> void SemSphericalElementMeshGenerator<DIM>::GenerateMesh(
    std::vector<c_vector<double, DIM>> positions)
{
    unsigned int new_element_id = mpMesh->GetNumElements();
    auto new_element = new SemElement<DIM>(new_element_id, {});

    for (unsigned i = 0; i < positions.size(); ++i)
    {
        unsigned new_node_index = mpMesh->GetNumNodes();
        Node<DIM>* new_node = new Node<DIM>(new_node_index, positions[i]);
        new_node->SetRegion(IsBoundaryNode(positions[i]) ? SEM_BOUNDARY_REGION : SEM_INTERIOR_REGION);
        new_node->AddElement(new_element_id);

        mpMesh->AddNode(new_node);
        new_element->AddNode(new_node);
    }

    mpMesh->AddElement(new_element);
    NodeMap map(mpMesh->GetNumNodes());
    mpMesh->ReMesh(map);
}

template class SemSphericalElementMeshGenerator<1>;
template class SemSphericalElementMeshGenerator<2>;
template class SemSphericalElementMeshGenerator<3>;
