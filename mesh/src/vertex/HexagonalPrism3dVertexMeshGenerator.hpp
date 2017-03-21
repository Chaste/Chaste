/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef HEXAGONALPRISM3DVERTEXMESHGENERATOR_HPP_
#define HEXAGONALPRISM3DVERTEXMESHGENERATOR_HPP_

#include "MutableVertexMesh.hpp"

/**
 * Honeycomb mesh generator that creates a 3D 'prismatic honeycomb' mesh for use in vertex simulations.
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
class HexagonalPrism3dVertexMeshGenerator
{
protected:
    /** A pointer to the mesh that this class creates. */
    MutableVertexMesh<3, 3>* mpMesh;

public:
    /**
     * Constructor.
     *
     * @param numElementsInXDirection he number of rows of elements in the x direction
     * @param numElementsInYDirection the number of rows of elements in the y direction
     * @param elementApicalArea the apical area of each element
     * @param elementHeight the height of each element in the z direction
     */
    HexagonalPrism3dVertexMeshGenerator(unsigned numElementsInXDirection,
                                        unsigned numElementsInYDirection,
                                        double elementApicalArea = 1.0,
                                        double elementHeight = 1.0);

    /**
     * Null constructor for derived classes to call.
     */
    HexagonalPrism3dVertexMeshGenerator()
    {
    }

    /**
     * Destructor - deletes the mesh object and pointer.
     */
    virtual ~HexagonalPrism3dVertexMeshGenerator();

    /**
     * @return a 3D mesh whose elements are hexagonal prisms
     */
    virtual MutableVertexMesh<3, 3>* GetMesh();
};

#endif /*HEXAGONALPRISM3DVERTEXMESHGENERATOR_HPP_*/
