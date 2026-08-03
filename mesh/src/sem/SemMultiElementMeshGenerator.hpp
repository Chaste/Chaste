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

#ifndef SEMMULTIELEMENTMESHGENERATOR_HPP_
#define SEMMULTIELEMENTMESHGENERATOR_HPP_


#include <vector>
#include <array>
#include <memory>

#include "SemEnumerations.hpp"
#include "SemMesh.hpp"
#include "UblasVectorInclude.hpp"


/**
 * Sem mesh generator that creates multiple SEM elements, in either 1D, 2D, or 3D,
 * on a regular lattice. The subcellular nodes within an element and the elements
 * themselves are laid out independently, so either may be cubic or close packed.
 */
template <unsigned DIM> class SemMultiElementMeshGenerator
{
private:
    /** A pointer to the mesh this class creates */
    std::shared_ptr<SemMesh<DIM>> mpMesh;

    /** Number of nodes per element in each dimension 0 <= dim < DIM */
    std::array<unsigned, DIM> mNumNodesPerElem;

    /** Number of elements in each dimension 0 <= dim < DIM */
    std::array<unsigned, DIM> mNumElems;

    /** The target diameter of a single element in the x-direction */
    double mScaleFactor;

    /** The node spacing necessary to achieve the target diameter */
    double mNodeSpacing;

    /** Spacing between elements in each dimension 0 <= dim < DIM */
    std::array<double, DIM> mElemSpacing;

    /** The lattice on which the subcellular nodes within each element are placed */
    SemLatticeType::Value mNodeLattice;

    /** The lattice on which the elements are placed */
    SemLatticeType::Value mElementLattice;

    /**
     * Generate uniformly spaced node positions.
     *
     * @return vector of node coordinates (each as c_vector<double, DIM>)
     */
    std::vector<c_vector<double, DIM>> GenerateNodePositions() const;

    /**
     * Generate element offsets.
     *
     * @return vector of offsets to translate from element 0 to each other element
     */
    std::vector<c_vector<double, DIM>> GenerateElementOffsets() const;

    /**
     * Generate the mesh from the node positions
     *
     * @param positions a vector of positions for each node in element 0
     * @param offsets a vector of offsets from element 0 to each other element
     */
    void GenerateMesh(std::vector<c_vector<double, DIM>> positions, std::vector<c_vector<double, DIM>> offsets);

    /**
     * Return whether the node at the given flat grid index is on the boundary
     * of the per-element node grid. Region 0 = interior, 1 = boundary.
     *
     * @param flatIndex linear index into the per-element node array
     * @return true if the node lies on any face of the per-element grid, false otherwise
     */
    bool IsBoundaryNode(unsigned flatIndex) const;

public:
    /**
     * Constructor.
     *
     * @param numNodesPerElem number of nodes per element in each direction (x, y, z)
     * @param numElems number of elements in each direction (x, y, z)
     * @param scaleFactor the target diameter of each element in the x-direction
     * @param nodeLattice the lattice on which to place the nodes within each element (defaults to SemLatticeType::SEM_LATTICE_CUBIC)
     * @param elementLattice the lattice on which to place the elements (defaults to SemLatticeType::SEM_LATTICE_CUBIC)
     */
    SemMultiElementMeshGenerator(const std::array<unsigned, DIM>& numNodesPerElem,
        const std::array<unsigned, DIM>& numElems, double scaleFactor = 1.0,
        SemLatticeType::Value nodeLattice = SemLatticeType::SEM_LATTICE_CUBIC,
        SemLatticeType::Value elementLattice = SemLatticeType::SEM_LATTICE_CUBIC);

    /**
     * Default constructor
     */
    SemMultiElementMeshGenerator() = default;

    /**
     * Default destructor
     */
    virtual ~SemMultiElementMeshGenerator() = default;

    /**
     * @return a shared view into the mesh
     */
    virtual std::shared_ptr<SemMesh<DIM>> GetMesh();
};

#endif /*SEMMULTIELEMENTMESHGENERATOR_HPP_*/
