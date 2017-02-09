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

// Created by bartosz on 07/02/17.

#ifndef CHASTE_IMMERSEDBOUNDARYAVEMIGRATIONMESHGENERATOR_HPP
#define CHASTE_IMMERSEDBOUNDARYAVEMIGRATIONMESHGENERATOR_HPP

#include <cmath>
#include <vector>

#include "ImmersedBoundaryMesh.hpp"
#include "SuperellipseGenerator.hpp"
#include "UblasCustomFunctions.hpp"

/**
 * Creates a collection of immersed boundary elements to model distal tip of a mouse embryo as an array of cells
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 * NOTE: the user should change mesh parameters like fluid grid spacing, as these are not altered from defaults here.
 */

class ImmersedBoundaryAVEMigrationMeshGenerator
{
protected:
    /** A pointer to the mesh this class creates. */
    ImmersedBoundaryMesh<2,2>* mpMesh;

public:

    /**
     * Default constructor.
     * @param numCellsVertically - number of cells vertically
     * @param numCellsHorizontally - number of cells horizontally
     * @param numNodesPerCell - the numebr of nodes in a cell
     * @param ellipseExponent - how round the cells are
     * @param spacing - the spacing between the cells in the array of cells
     * @param shift - a factor by which to shift the cells in every other row
     */

    ImmersedBoundaryAVEMigrationMeshGenerator(unsigned numCellsHorizontally=5,
                                              unsigned numCellsVertically=5,
                                              unsigned numNodesPerCell=100,
                                              double ellipseExponent=0.01,
                                              double spacing=0.05,
                                              double shift=0.0);

    /**
     * Null constructor for derived classes to call.
     */
    ImmersedBoundaryAVEMigrationMeshGenerator()
    {
    }

    /**
     * Destructor.
     *
     * Deletes the mesh object and pointer.
     */
    virtual ~ImmersedBoundaryAVEMigrationMeshGenerator();

    /**
     * @return a 2D honeycomb mesh based on a 2D plane
     */
    ImmersedBoundaryMesh<2, 2>* GetMesh();

};


#endif //CHASTE_IMMERSEDBOUNDARYAVEMIGRATIONMESHGENERATOR_HPP
