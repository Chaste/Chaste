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


#ifndef CUBOIDMESHCONSTRUCTOR_HPP_
#define CUBOIDMESHCONSTRUCTOR_HPP_

#include "TetrahedralMesh.hpp"

/**
 * Helper class for constructing cuboidal meshes.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class CuboidMeshConstructor
{
private:

    double mMeshWidth; /**< Width of the mesh.*/
    unsigned mNumElements; /**< Number of elements in the mesh. */
    unsigned mNumNodes; /**< Number of nodes in the mesh. */

public:

    /**
     * Construct the mesh.
     *
     * @param rMesh  Input a blank mesh in which to construct the result
     * @param meshRefinementNum  Index for the mesh starting at 0 (4 elements in each space dimension)
     * @param meshWidth  Width of the mesh (in cm)
     *
     * (used to return a path to the mesh stored on the disk)
     *
     */
    void Construct(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, unsigned meshRefinementNum, double meshWidth);

    /**
     * @return the width of the mesh.
     */
    double GetWidth();

    /**
     * @return the number of elements in the mesh.
     */
    unsigned GetNumElements();

    /**
     * @return the number of nodes in the mesh.
     */
    unsigned GetNumNodes();
};

#endif /*CUBOIDMESHCONSTRUCTOR_HPP_*/
