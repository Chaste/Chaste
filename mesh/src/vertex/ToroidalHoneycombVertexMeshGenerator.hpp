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

#ifndef TOROIDALHONEYCOMBVERTEXMESHGENERATOR_HPP_
#define TOROIDALHONEYCOMBVERTEXMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "HoneycombVertexMeshGenerator.hpp"
#include "Toroidal2dVertexMesh.hpp"

/**
 * Honeycomb mesh generator that creates a 2D "toroidal" mesh (one in which
 * periodicity is imposed on the left and right and top and bottom boundaries)
 * for use with vertex-based simulations.
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
class ToroidalHoneycombVertexMeshGenerator : HoneycombVertexMeshGenerator
{
public:

    /**
     * Constructor.
     *
     * @param numElementsAcross  The number of columns of elements in the mesh.  This MUST be an even number.
     * @param numElementsUp  The number of rows of elements in the mesh.   This MUST be an even number.
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangement (defaults to 0.01)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     */
    ToroidalHoneycombVertexMeshGenerator(unsigned numElementsAcross,
                                         unsigned numElementsUp,
                                         double cellRearrangementThreshold=0.01,
                                         double t2Threshold=0.001);
    /**
     * @return a 2D honeycomb mesh
     */
    MutableVertexMesh<2,2>* GetMesh();

    /**
     * @return a 2D honeycomb mesh with periodic left/right and top/bottom boundaries
     */
    Toroidal2dVertexMesh* GetToroidalMesh();
};

#endif /*TOROIDALHONEYCOMBVERTEXMESHGENERATOR_HPP_*/
