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

#include "HexagonalPrism3dVertexMeshGenerator.hpp"

#include "HoneycombVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshGenerator.hpp"

HexagonalPrism3dVertexMeshGenerator::HexagonalPrism3dVertexMeshGenerator(unsigned numElementsInXDirection,
    unsigned numElementsInYDirection,
    double elementApicalArea,
    double elementHeight)
{
    assert(numElementsInXDirection > 0);
    assert(numElementsInYDirection > 0);
    assert(elementApicalArea > 0.0);
    assert(elementHeight > 0.0);

    HoneycombVertexMeshGenerator generator2 (numElementsInXDirection, numElementsInYDirection, false, 1, 1, elementApicalArea);
    MutableVertexMesh<2, 2>* p_mesh2 = generator2.GetMesh();
    MonolayerVertexMeshGenerator generator("", false);
    mpMesh = generator.MakeMeshUsing2dMesh(*p_mesh2, elementHeight);
}

HexagonalPrism3dVertexMeshGenerator::~HexagonalPrism3dVertexMeshGenerator()
{
    delete mpMesh;
}

MutableVertexMesh<3,3>* HexagonalPrism3dVertexMeshGenerator::GetMesh()
{
    return mpMesh;
}