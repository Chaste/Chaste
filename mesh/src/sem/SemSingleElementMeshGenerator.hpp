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

#ifndef SEMSINGLEELEMENTMESHGENERATOR_HPP_
#define SEMSINGLEELEMENTMESHGENERATOR_HPP_


#include <vector>
#include <array>
#include <memory>

#include "UblasVectorInclude.hpp"

// Forward declarations
template <unsigned DIM> class SemMesh;


/**
 * Sem mesh generator that creates a single SEM element, in either 1D, 2D, or 3D
 */
template <unsigned DIM> class SemSingleElementMeshGenerator
{
private:
    /** A pointer to the mesh this class creates */
    std::shared_ptr<SemMesh<DIM>> mpMesh;

    /** Number of nodes in each dimension 0 <= dim < DIM */
    std::array<unsigned, DIM> mNumNodes;

    /** Number of nodes in total */
    unsigned mNumAllNodes;

    /** The target diameter of the element in the x-direction */
    double mScaleFactor;

    /** The node spacing necessary to achieve the target diameter */
    double mNodeSpacing;

    /**
     * Generate uniformly spaced node positions.
     *
     * @return vector of node coordinates (each as c_vector<double, DIM>)
     */
    std::vector<c_vector<double, DIM>> GenerateNodePositions() const;

    /**
     * Generate the mesh from the node positions
     */
    void GenerateMesh(std::vector<c_vector<double, DIM>> positions);

public:
    /**
     * Constructor.
     *
     * @param numNodes number of nodes in each direction (x, y, z)
     * @param scaleFactor the target diameter of the element in the x-direction
     */
    SemSingleElementMeshGenerator(const std::array<unsigned, DIM>& numNodes, double scaleFactor = 1.0);

    /**
     * Default constructor
     */
    SemSingleElementMeshGenerator() = default;

    /**
     * Default destructor
     */
    virtual ~SemSingleElementMeshGenerator() = default;

    /**
     * @return a shared view into the mesh
     */
    virtual std::shared_ptr<SemMesh<DIM>> GetMesh();
};

#endif /*SEMSINGLEELEMENTMESHGENERATOR_HPP_*/
