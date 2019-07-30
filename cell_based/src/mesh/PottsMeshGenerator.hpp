/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef POTTSMESHGENERATOR_HPP_
#define POTTSMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "PottsMesh.hpp"

/**
 * Generator of regular Potts meshes, used as starting points for many simulations.
 *
 * This class takes in options such as width, height,
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
template<unsigned DIM>
class PottsMeshGenerator
{
protected:

    /** A pointer to the mesh this class creates */
    PottsMesh<DIM>* mpMesh;

public:

    /**
     * Constructor.
     *
     * @param numNodesAcross  The number of columns of nodes in the mesh
     * @param numElementsAcross  The number of columns of elements in the mesh
     * @param elementWidth  Width of the elements
     * @param numNodesUp  The number of rows of nodes in the mesh (defaults to 1)
     * @param numElementsUp  The number of rows of elements in the mesh (defaults to 1)
     * @param elementHeight  Height of the elements (defaults to 1)
     * @param numNodesDeep  The number nodes deep for this mesh (defaults to 1)
     * @param numElementsDeep  The number of elements deep for this mesh (defaults to 1)
     * @param elementDepth  The number of rows of nodes in each element (defaults to 1)
     * @param startAtBottomLeft  If true then the mesh starts in the bottom left corner
     *     of the domain rather than the centre, used for simple tests (defaults to false)
     * @param isPeriodicInX  If true then the mesh is periodic in the x dimension (defaults to false)
     * @param isPeriodicInY  If true then the mesh is periodic in the y dimension (defaults to false)
     * @param isPeriodicInZ  If true then the mesh is periodic in the y dimension (defaults to false)
     */
    PottsMeshGenerator(unsigned numNodesAcross,
                       unsigned numElementsAcross,
                       unsigned elementWidth,
                       unsigned numNodesUp=1u,
                       unsigned numElementsUp=1u,
                       unsigned elementHeight=1u,
                       unsigned numNodesDeep=1u,
                       unsigned numElementsDeep=1u,
                       unsigned elementDepth=1u,
                       bool startAtBottomLeft = false,
                       bool isPeriodicInX = false,
                       bool isPeriodicInY = false,
                       bool isPeriodicInZ = false);

    /**
     * Destructor - deletes the mesh object and pointer.
     */
    virtual ~PottsMeshGenerator();

    /**
     * Helper method to calculate the Moore and Von Neumann Neighbourhoods of all nodes
     *
     * @param isPeriodicInX  If true then the mesh is periodic in the x dimension
     */
    void CaclulateNeighbouringNodeIndices(bool isPeriodicInX);

    /**
     * @return a Cuboid or rectangular Potts mesh.
     */
    virtual PottsMesh<DIM>* GetMesh();
};

#endif /*POTTSMESHGENERATOR_HPP_*/
