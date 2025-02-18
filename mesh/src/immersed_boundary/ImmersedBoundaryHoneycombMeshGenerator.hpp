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

#ifndef IMMERSEDBOUNDARYHONEYCOMBMESHGENERATOR_HPP_
#define IMMERSEDBOUNDARYHONEYCOMBMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "ImmersedBoundaryMesh.hpp"

/**
 * Creates a honeycomb of immersed boundary elements.
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
class ImmersedBoundaryHoneycombMeshGenerator
{

private:

    /** A pointer to the mesh this class creates */
    ImmersedBoundaryMesh<2,2>* mpMesh;

public:

    /**
     * Default constructor.
     *
     * @param numElementsX  the number of cells from left to right along the domain
     * @param numElementsY  the number of cells from top to bottom up the domain
     * @param numNodesPerEdge  the number of nodes per cell (defaults to 100)
     * @param proportionalGap  the proportion of space between elements
     * @param padding  the minimum padding around the edge of the generated mesh
     */
    ImmersedBoundaryHoneycombMeshGenerator(unsigned numElementsX,
                                           unsigned numElementsY,
                                           unsigned numNodesPerEdge,
                                           double proportionalGap,
                                           double padding);

    /**
     * Null constructor for derived classes to call.
     */
    ImmersedBoundaryHoneycombMeshGenerator()
    {
    }

    /**
     * Destructor.
     *
     * Deletes the mesh object and pointer.
     */
    virtual ~ImmersedBoundaryHoneycombMeshGenerator();

    /**
     * @return a 2D honeycomb mesh based on a 2D plane
     */
    ImmersedBoundaryMesh<2,2>* GetMesh();

    /**
     * Helper method for the constructor that calculates locations around a unit hexagon centred at the origin.
     *
     * @param numPtsPerSide the number of locations along each of the six sides
     * @return a vector of locations around the unit hexagon
     */
    std::vector<c_vector<double, 2> > GetUnitHexagon(unsigned numPtsPerSide);
};

#endif /*IMMERSEDBOUNDARYHONEYCOMBMESHGENERATOR_HPP_*/
